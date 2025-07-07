### S4-BuyseTest-nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  3 2023 (10:00) 
## Version: 
## Last-Updated: Jul  7 2025 (12:11) 
##           By: Brice Ozenne
##     Update #: 52
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - update
#' @docType methods
#' @name S4BuyseTest-update
#' @title Re-run two-group GPC
#' @aliases update,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Re-run Generalized Pairwise Comparisons (GPC) between two groups.
#' 
#' @param object an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' @param ... additional arguments to the call, or arguments with changed values.
#' 
#' @return an \R object of class \code{S4BuyseTest}
#' 
#' @keywords methods
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuyseTest-update
#' @exportMethod update
setMethod(f = "update",
          signature = "S4BuyseTest",
          definition = function(object, ...){

              ## ** retrieve args
              old.args <- object@call
              user.args <- list(...)

              common.args <- intersect(names(old.args),names(user.args))
              new.args <- setdiff(names(user.args),names(old.args))
              
              old.args[common.args] <- user.args[common.args]

              ## ** re-run GPC
              out <- do.call(BuyseTest, args = c(old.args, user.args[new.args]))

              ## ** export
              return(out)
          }
          )


##----------------------------------------------------------------------
### S4-BuyseTest-nobs.R ends here
