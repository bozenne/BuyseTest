### S4-BuyseTest-plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 23 2023 (10:44) 
## Version: 
## Last-Updated: Jun 29 2023 (10:25) 
##           By: Brice Ozenne
##     Update #: 131
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## 

## * plot (documentation)
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

## * plot (code)
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
### S4-BuyseTest-plot.R ends here
