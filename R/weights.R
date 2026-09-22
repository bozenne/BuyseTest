### weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 22 2026 (16:01) 
## Version: 
## Last-Updated: sep 22 2026 (16:07) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * weights S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-weights
#' @title Extract Weigths from GPC
#' @aliases weights,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Extract weights used by the GPC procedure.
#' 
#' @param object a \code{S4BuyseTest} object, output of \code{\link{BuyseTest}}.
#' @param type [character] should the weights used to pool summary statistics across strata (\code{"strata"}),
#' the weights used when summing across endpoints (\code{"endpoint"}),
#' or the weights of each observation (\code{"obs"})
#' be output.
#' @param ... ignored.
#' 
#' @return A numeric vector.
#' @keywords method
#' @author Brice Ozenne
#' @exportMethod weights
setMethod(f = "weights",
          signature = "S4BuyseTest",
          definition = function(object,
                                type = "strata",
                                ...){


              ## ** check user input
              type <- match.arg(type, c("strata","endpoint","obs"))

              ## ** export
              if(type == "strata"){
                  out <- object@weightStrata
              }else if(type == "endpoint"){
                  out <- object@weightEndpoint
              }else if(type == "obs"){
                  out <- object@weightObs
              }
              
              ## ** export
              return(out)

          })


##----------------------------------------------------------------------
### weights.R ends here
