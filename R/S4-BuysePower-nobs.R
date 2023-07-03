### S4-BuysePower-nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  3 2023 (10:00) 
## Version: 
## Last-Updated: Jul  3 2023 (10:34) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - print
#' @docType methods
#' @name S4BuysePower-nobs
#' @title Sample Size for Class "S4BuysePower"
#' @aliases nobs,S4BuysePower-method
#' @include S4-BuysePower.R
#' 
#' @description Display the sample size in each treatmnet arm as well as the number of pairs.
#' 
#' @param object an \R object of class \code{S4BuysePower}, i.e., output of \code{\link{powerBuyseTest}}
#' @param ... no used, for compatibility with the generic method.
#' 
#' @return A data.frame with two colunms, one for each treatment group, and as many rows as sample sizes used for the simulation.
#' 
#' @keywords methods
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuysePower-nobs
#' @exportMethod nobs
setMethod(f = "nobs",
          signature = "S4BuysePower",
          definition = function(object, ...){

              return(object@sample.size)
          }
)


##----------------------------------------------------------------------
### S4-BuysePower-nobs.R ends here
