### autoplot-S4BuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 29 2023 (09:27) 
## Version: 
## Last-Updated: jun  4 2024 (10:34) 
##           By: Brice Ozenne
##     Update #: 26
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * autoplot (documentation)
#' @title Graphical Display for GPC
#' @description Graphical display of the percentage of favorable, unfavorable, neutral, and uninformative pairs per endpoint.
#' @rdname autoplot-S4BuyseTest
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param type [character] type of plot: histogram (\code{"hist"}), pie chart (\code{"pie"}), or nested pie charts (\code{"racetrack"}).
#' @param strata [character vector] strata(s) relative to which the percentage should be displayed.
#' @param endpoint [character vector] endpoint(s) relative to which the percentage should be displayed.
#' @param label.strata [character vector] new labels for the strata levels. Should match the length of argument \code{strata}.
#' @param label.endpoint [character vector] new labels for the endpoints. Should match the length of argument \code{endpoint}.
#' @param color [character vector] colors used to display the percentages for each type of pair.
#' @param ... not used, for compatibility with the generic function.
#'
#' @return a ggplot object.
#' @method autoplot S4BuyseTest
#' @keywords hplot
#' @export
autoplot.S4BuyseTest <- function(object, type = "hist", strata = "global", endpoint = NULL, 
                                 label.strata = NULL, label.endpoint = NULL,
                                 color = c("#7CAE00", "#F8766D", "#C77CFF", "#00BFC4"), ...){


              objectS <- model.tables(object, percentage = TRUE,
                                      columns = c("endpoint","threshold","restriction","strata","total","favorable","unfavorable","neutral","uninf"))
              Ustrata <- slot(object,"level.strata")
              Uendpoint <- slot(object,"endpoint")
              objectS$endpoint2 <- attr(objectS,"endpoint")
              
              ## ** normalize arguments
              ## type
              type <- match.arg(type, c("hist","pie","racetrack"))

              ## strata
              level.strata <- object@level.strata
              if(is.null(strata)){
                  if(length(level.strata)==1){
                      strata <- "global"                      
                  }else{
                      strata <- c("global", level.strata)
                  }
              }else if(identical(strata,FALSE)){
                  strata <- "global"
              }else if(identical(strata,TRUE)){
                  strata <- level.strata
              }else if(is.numeric(strata)){
                  validInteger(strata,
                               name1 = "strata",
                               valid.length = NULL,
                               min = 1,
                               max = length(level.strata),
                               refuse.NULL = TRUE,
                               refuse.duplicates = TRUE,
                               method = "autoplot[S4BuyseTest]")
                  strata <- level.strata[strata]
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",level.strata),
                                 refuse.NULL = FALSE,
                                 method = "autoplot[S4BuyseTest]")
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
                  gg <- gg + ggplot2::labs(x = "", y = "", fill = "Pair (%)")
                  if(length(strata)>1){
                      gg <- gg + ggplot2::facet_wrap(~strata)
                  }
              }
              gg <- gg + ggplot2::scale_fill_manual(values = color)

    ## ** export
    return(gg)
}

##----------------------------------------------------------------------
### autoplot-S4BuyseTest.R ends here
