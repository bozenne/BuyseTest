### update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:44) 
## Version: 
## Last-Updated: jul  3 2026 (14:44) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * update S4BuyseTest
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
### update.R ends here
