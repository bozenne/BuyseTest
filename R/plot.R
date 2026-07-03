### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 10 2021 (09:34) 
## Version: 
## Last-Updated: jul  3 2026 (14:40) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot - sensitivity
#' @rdname plot-sensitivity
#' @method plot S3sensitivity
#' @export
plot.S3sensitivity <- function(x, plot = TRUE, ...){

    out <- autoplot(x, ...)
    if(plot){
        print(out)
    }
    return(invisible(list(plot = out,
                          data = out$data)))
}

## * autoplot - sensitivity
##' @title Graphical Display for Sensitivity Analysis 
##' @description Display the statistic of interest across various threshold values, possibly with confidence intervals.
##' Currently only works when varying thresholds relative to one or two variables.
##' @rdname plot-sensitivity
##' 
##' @param object,x output of the sensitivity method
##' @param plot [logical] should the graph be displayed in a graphical window
##' @param col [character vector] color used to identify the thresholds relative to a second variable.
##' @param ci [logical] should the confidence intervals be displayed?
##' @param band [logical] should the simulatenous confidence intervals be displayed?
##' @param label [character] text used before the name of the variables in the legend.
##' @param size.line [numeric] width of the line connecting the point estimates.
##' @param size.point [numeric] size of the point representing the point estimates.
##' @param size.ci [numeric] width of the lines representing the confidence intervals.
##' @param alpha [numeric] transparency for the area representing the simultaneous confidence intervals.
##' @param position relative position of the error bars for a given x value. Can for instance be \code{position_dodge(width = 5)}.
##' @param ... not used. For compatibility with the generic method.
##' 
##' @details The \code{autoplot} and \code{plot} methods are very similar. The main difference is that the former returns a ggplot2 object whereas the later automatically display the figure in a graphical window and returns an (invible) list with the plot and the data.
##'
##' @return a ggplot2 object
##' @method autoplot S3sensitivity
##' @keywords hplot
##' 
##' @export
autoplot.S3sensitivity <- function(object, col = NULL, ci = TRUE, band = TRUE, label = "Threshold for", 
                                   position = NULL, size.line = 1, size.point = 1.75, size.ci = 0.5, alpha = 0.1, ...){
    grid <- attr(object,"gridRed")
    statistic <- switch(attr(object,"statistic"),
                        "netBenefit" = "Net benefit",
                        "winRatio" = "Win ratio",
                        "favorable" = "Proportion of favorable pairs",
                        "unfavorable" = "Proportion of unfavorable pairs")
                        
    if(NCOL(grid)>2){
        stop("No graphical display available when the sensitivity analysis is performed on more than 2 thresholds\n")
    }
    nU.var <- apply(grid,2,function(x){length(unique(x))})
    name.var <- names(sort(nU.var, decreasing = TRUE))
    n.var <- length(name.var)

    name.col <- name.var
    if(n.var==1 || (!is.null(col) && all(is.na(col)))){
        
        if("XXindexXX" %in% names(object)){
            stop("No endpoint should be named \"XXindexXX\" as this name is used internally. \n")
        }
        name.col[2] <- "XXindexXX"
        object <- data.frame(XXindexXX = "1", object)
    }else{
        object[[name.var[2]]] <- factor(object[[name.var[2]]], levels = sort(unique(object[[name.var[2]]])))
    }

    ## ** display
    ## error bar in the legend
    draw_key.save <- ggplot2::GeomErrorbar$draw_key
    GeomErrorbar$draw_key  <-   function (data, params, size) { ## https://stackoverflow.com/questions/53490654/adding-the-errorbar-icon-in-legend-in-ggplot
        .pt <- get(".pt", envir = as.environment("package:ggplot2"))
        data$linetype[is.na(data$linetype)] <- 0
        out <- grid::segmentsGrob(c(0.2, 0.2, 0.5), c(0.2, 0.8, 0.2), c(0.8, 0.8, 0.5), c(0.2, 0.8, 0.8),
                                  gp = grid::gpar(col = alpha(data$colour, data$alpha), lwd = data$linewidth * .pt, lty = data$linetype, lineend = "butt"), arrow = params$arrow)
        return(out)
    }
    on.exit(GeomErrorbar$draw_key  <-   draw_key.save)
    if(length(name.var)==1){
        gg <- ggplot2::ggplot(data = object, mapping = ggplot2::aes(x = .data[[name.var[1]]], y = .data$estimate))
    }else{
        gg <- ggplot2::ggplot(data = object, mapping = ggplot2::aes(x = .data[[name.var[1]]], y = .data$estimate, group = .data[[name.var[2]]]))
    }
    if(band && "lower.band" %in% names(object) && "upper.band" %in% names(object)){
        gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower.band, ymax = .data$upper.band, fill = .data[[name.col[2]]]), alpha = alpha)
    }else{
        band <- FALSE
    }
    gg <- gg + ggplot2::geom_point(ggplot2::aes(color = .data[[name.col[2]]]), size = size.point) + ggplot2::geom_line(ggplot2::aes(color = .data[[name.col[2]]]), linewidth = size.line)
    gg <- gg + ggplot2::xlab(paste(label,name.var[1],sep=" "))
    gg <- gg + ggplot2::ylab(statistic) + ggplot2::theme(legend.position = "bottom")

    if(ci && "lower.ci" %in% names(object) && "upper.ci" %in% names(object)){
        if(!is.null(position)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$lower.ci, ymax = .data$upper.ci, color = .data[[name.col[2]]]), size = size.ci, position = position)
        }else{
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$lower.ci, ymax = .data$upper.ci, color = .data[[name.col[2]]]), size = size.ci)
        }
    }else{
        ci <- FALSE
    }
    if(n.var==1){

        if(is.null(col) || all(is.na(col))){
            col <- "black"
        }else if(length(col)!=1){
            stop("Argument \'col\' should have lenght one when the sensitivity analysis is performed on one threshold. \n")
        }
        if(ci && "lower.ci" %in% names(object) && "upper.ci" %in% names(object)){
            gg <- gg + ggplot2::scale_color_manual("CIs", values = col, labels = "")
        }else{
            gg <- gg + ggplot2::scale_color_manual("Point estimate", values = col, labels = "")
        }
        if(band){
            gg <- gg + ggplot2::scale_fill_manual("Simulatenous CIs", values = col, labels = "")
        }
        
    }else if(n.var==2){
        if(!is.null(col) && all(is.na(col))){
            Ulevel.var2 <- unique(object[[name.var[2]]])
            label_facet <- setNames(unique(paste(label,name.var[[2]]," : ",Ulevel.var2,sep=" ")), Ulevel.var2)
            gg <- gg + ggplot2::facet_grid(as.formula(paste0("~",name.var[2])), labeller = ggplot2::as_labeller(label_facet))

            if(ci){
                gg <- gg + ggplot2::scale_color_manual("CIs", values = "black", labels = "")
            }else{
                gg <- gg + ggplot2::scale_color_manual("Point estimate", values = "black", labels = "")
            }
            if(band){
                gg <- gg + ggplot2::scale_fill_manual("Simulatenous CIs", values = "black", labels = "")
            }
        }else if(is.null(col)){
            if(ci){
                gg <- gg + ggplot2::labs(color = paste0("CIs \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"))
            }else{
                gg <- gg + ggplot2::labs(color = paste0("Point estimate \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"))
            }
            if(band){
                gg <- gg + ggplot2::labs(fill = paste0("Simulatenous CIs \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"))
            }
        }else{
            if(length(col)!=nU.var[[name.var[2]]]){
                stop("Argument \'col\' should have lenght ",nU.var[[name.var[2]]],", the number of unique thresholds relative to the endpoint \"",name.var[2],"\". \n")
            }
            if(ci){
                gg <- gg + ggplot2::scale_color_manual(paste0("CIs \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"), values = col)
            }else{
                gg <- gg + ggplot2::scale_color_manual(paste0("Point estimate \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"), values = col)
            }
            if(band){
                gg <- gg + ggplot2::scale_fill_manual(paste0("Simulatenous CIs \n (",paste(c(tolower(label),name.col[2]),collapse=" "),")"), values = col)
            }
        }
    }

    return(gg)
}

## * plot S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-plot 
#' @title Graphical Display for GPC
#' @aliases plot,S4BuyseTest,ANY-method
#' @include S4-BuyseTest.R
#' 
#' @description Graphical display of the percentage of favorable, unfavorable, neutral, and uninformative pairs per endpoint.
#' 
#' @param x an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param type [character] type of plot: histogram (\code{"hist"}), pie chart (\code{"pie"}), or nested pie charts (\code{"racetrack"}).
#' @param strata [character vector] strata(s) relative to which the percentage should be displayed.
#' @param endpoint [character vector] endpoint(s) relative to which the percentage should be displayed.
#' @param label.strata [character vector] new labels for the strata levels. Should match the length of argument \code{strata}.
#' @param label.endpoint [character vector] new labels for the endpoints. Should match the length of argument \code{endpoint}.
#' @param color [character vector] colors used to display the percentages for each type of pair.
#' @param plot [logical] should the graphic be displayed in a graphical window.
#' @param ... not used, for compatibility with the generic function.
#'
#' @return an invisible list containing the data and the ggplot object used for graphical display.
#' @keywords hplot
#' 
#' @examples 
#' if(require(ggplot2)){
#' 
#' ## simulate data
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' ff1 <- treatment ~ bin(toxicity) + TTE(eventtime, status = status,
#'                                       restriction = 1, threshold = 0.5)
#' BT1 <- BuyseTest(ff1, data= df.data)
#' plot(BT1, type = "hist")
#' plot(BT1, type = "pie")
#' plot(BT1, type = "racetrack")
#'
#' ff2 <- update(ff1, ~.+cont(score))
#' BT2 <- BuyseTest(ff2, data= df.data)
#' plot(BT2, type = "hist")
#' plot(BT2, type = "pie")
#' plot(BT2, type = "racetrack")
#'
#' }
#' @export
setMethod(f = "plot",
          signature = "S4BuyseTest",
          definition = function(x, type = "hist", strata = "global", endpoint = NULL, 
                                label.strata = NULL, label.endpoint = NULL,
                                plot = TRUE, color = c("#7CAE00", "#F8766D", "#C77CFF", "#00BFC4"), ...){

              gg <- autoplot(x, type = type, strata = strata, endpoint = endpoint, 
                             label.strata = label.strata, label.endpoint = label.endpoint,
                             color = color)

              ## ** display
              if(plot){
                  print(gg)
              }

              ## ** export
              return(invisible(list(plot = gg,
                                    data = gg$data)))
              
})



##----------------------------------------------------------------------
### plot.R ends here
