### S4-BuyseTest-plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 23 2023 (10:44) 
## Version: 
## Last-Updated: jun 28 2023 (14:04) 
##           By: Brice Ozenne
##     Update #: 117
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## 

## * Documentation - plot
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
#' \dontrun{
#' require(ggplot2)
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
#' }

## * Method - plot
#' @export
setMethod(f = "plot",
          signature = "S4BuyseTest",
          definition = function(x, type = "hist", strata = "global", endpoint = NULL, 
                                label.strata = NULL, label.endpoint = NULL,
                                plot = TRUE, color = c("#7CAE00", "#F8766D", "#C77CFF", "#00BFC4"), ...){

              objectS <- model.tables(x, percentage = TRUE)
              Ustrata <- setdiff(unique(objectS$strata),"global")
              Uendpoint <- x@endpoint
              objectS$endpoint2 <- paste0(objectS$endpoint,
                                          ifelse(!is.na(objectS$restriction),paste0("_r",objectS$restriction),""),
                                          ifelse(objectS$threshold>1e-12,paste0("_t",objectS$threshold),""))

              
              ## ** normalize arguments
              ## type
              type <- match.arg(type, c("hist","pie","racetrack"))

              ## strata
              if(!is.null(strata)){
                  if(is.numeric(strata)){
                      if(any(strata %in% 1:length(Ustrata)==FALSE)){
                          stop("Incorrect argument \'strata\': when numeric should be an integer vector with values between 1 and ",length(Ustrata),".\n", sep ="")
                      }
                      strata <- Ustrata[strata]
                  }else{
                      strata <- match.arg(strata, c("global",Ustrata), several.ok = TRUE)
                  }                  
              }else{
                  strata <- unique(objectS$strata)
              }

              ## endpoint
              if(!is.null(endpoint)){
                  if(is.numeric(endpoint)){
                      if(any(endpoint %in% 1:length(Uendpoint)==FALSE)){
                          stop("Incorrect argument \'endpoint\': when numeric should be an integer vector with values between 1 and ",length(endpoint),".\n", sep ="")
                      }
                      endpoint <- names(Uendpoint)[endpoint]
                  }else if(all(!duplicated(Uendpoint)) && all(endpoint %in% Uendpoint)){
                      endpoint <- names(Uendpoint)[match(endpoint, Uendpoint)]
                  }else{
                      endpoint <- match.arg(endpoint, names(Uendpoint), several.ok = TRUE)
                  }                  
              }else{
                  endpoint <- names(Uendpoint)
              }
              objectSS <- objectS[objectS$strata %in% strata & objectS$endpoint2 %in% endpoint,c("endpoint2","strata","total","favorable","unfavorable","neutral","uninf"),drop=FALSE]
              if(!is.null(label.strata)){
                  if(length(label.strata) != length(strata)){
                      stop("Length of argument \'label.strata\' must match the length of argument \'stata\' (here ",length(strata),").\n")
                  }
                  objectSS$strata <- factor(objectSS$strata, levels = strata, labels = label.strata)
              }else{
                  objectSS$strata <- factor(objectSS$strata, levels = strata)
              }
              if(!is.null(label.endpoint)){
                  if(length(label.endpoint) != length(endpoint)){
                      stop("Length of argument \'label.endpoint\' must match the length of argument \'endpoint\' (here ",length(endpoint),").\n")
                  }
                  objectSS$endpoint <- factor(objectSS$endpoint2, levels = endpoint, labels = label.endpoint)
              }else{
                  objectSS$endpoint <- factor(objectSS$endpoint2, levels = endpoint)
              }

              ## ** reshape data
              objectSSL <- stats::reshape(objectSS,
                                                idvar = c("endpoint","strata","total"), direction = "long",
                                                varying  = c("favorable","unfavorable","neutral","uninf"),
                                                times = c("favorable","unfavorable","neutral","uninf"),
                                                timevar = "type",
                                                v.names = "percentage")
              rownames(objectSSL) <- NULL

              objectSSL$strata <- factor(objectSSL$strata, levels = levels(objectSS$strata))
              objectSSL$endpoint <- factor(objectSSL$endpoint, levels = levels(objectSS$endpoint))
              objectSSL$type <- factor(objectSSL$type, c("favorable","unfavorable","neutral","uninf"))
                  

              ## ** graphical display
              if(type == "pie"){
                  gg <- ggplot2::ggplot(objectSSL, ggplot2::aes(x="", y=.data$percentage, fill=.data$type))
                  gg <- gg + ggplot2::geom_bar(stat="identity", width = 1)
                  gg <- gg + ggplot2::coord_polar("y", start=0) + ggplot2::labs(x = "", y = "", fill = "Pair (%)")
                  gg

                  if(length(strata)>1 && length(endpoint)>1){
                      gg <- gg + ggplot2::facet_grid(strata~endpoint)
                  }else if(length(strata)>1){
                      gg <- gg + ggplot2::facet_wrap(~strata)
                  }else if(length(endpoint)>1){
                      gg <- gg + ggplot2::facet_wrap(~endpoint)
                  }
              }else if(type == "racetrack"){
                  gg <- ggplot2::ggplot(objectSSL, ggplot2::aes(x=rev(.data$endpoint), y=.data$percentage))
                  gg <- gg + ggplot2::geom_bar(stat="identity", width = 1, ggplot2::aes(fill=.data$type))
                  gg <- gg + ggplot2::coord_polar("y", start=0) + ggplot2::labs(x = "", y = "", fill = "Pair (%)")
                  ## gg <- gg + ggplot2::theme(axis.text.y=ggplot2::element_blank())
                  ## gg <- gg + ggplot2::geom_text(data = data.frame(endpoint=levels(objectSSL$endpoint),percentage=0), mapping = ggplot2::aes(label = .data$endpoint), size = size.text)
                  ## gg <- gg + scale_x_discrete(expand = expand_scale(add = c(1,1)))
                  if(length(strata)>1){
                      gg <- gg + ggplot2::facet_wrap(~strata)
                  }
              }else if(type == "hist"){

                  gg <- ggplot2::ggplot(objectSSL, ggplot2::aes(x=.data$endpoint, y=.data$percentage/100))
                  gg <- gg + ggplot2::geom_bar(ggplot2::aes(fill=.data$type), position = "stack", stat = "identity")
                  gg <- gg + ggplot2::scale_y_continuous(labels=scales::percent)
                  gg <- gg + ggplot2::labs(x = "", y = "", fill = "Pair")
                  if(length(strata)>1){
                      gg <- gg + ggplot2::facet_wrap(~strata)
                  }
              }
              gg <- gg + ggplot2::scale_fill_manual(values = color)

              ## ** display
              if(plot){
                  print(gg)
              }

              ## ** export
              return(invisible(list(plot = gg,
                                    data = gg$data)))
              
})

##----------------------------------------------------------------------
### S4-BuyseTest-plot.R ends here
