### show.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:36) 
## Version: 
## Last-Updated: jul  3 2026 (14:43) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * show S4BuysePower
#' @docType methods
#' @name S4BuysePower-show
#' @title Show Method for Class "S4BuysePower"
#' @aliases show,S4BuysePower-method
#' @include S4-BuysePower.R print.R
#' 
#' @description Display the main results stored in a \code{S4BuysePower} object.
#' 
#' @param object an \R object of class \code{S4BuysePower}, i.e., output of \code{\link{powerBuyseTest}}
#' 
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing power calculation based on GPC. \cr
#'   \code{\link{S4BuysePower-summary}} for a more detailed presentation of the \code{S4BuysePower} object.
#'  
#' @return invisible \code{NULL}
#' @keywords print
#' 
#' @author Brice Ozenne
#' @exportMethod show
setMethod(f = "show",
          signature = "S4BuysePower",
          definition = function(object){

              print(object)

              return(invisible(NULL))
          }
          )

## * show S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-show
#' @title Show Method for Class "S4BuyseTest"
#' @aliases show,S4BuyseTest-method
#' @include S4-BuyseTest.R summary.R
#' 
#' @description Display the main results stored in a \code{S4BuyseTest} object.
#' 
#' @param object an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#'  
#' @keywords print
#' @author Brice Ozenne
#' @exportMethod show
setMethod(f = "show",
          signature = "S4BuyseTest",
          definition = function(object){

              ## compute summary statistics
              outSummary <- summary(object, print = FALSE, strata = "global")

              ## only keep certain columns
              type.display <- BuyseTest.options()$print.display
              vec.tfunu <- c("total","favorable","unfavorable","neutral","uninformative")
              if(any(vec.tfunu %in% type.display)){
                  type.display[type.display %in% vec.tfunu] <- paste0(type.display[type.display %in% vec.tfunu],"(%)")
              }
              if("CI" %in% type.display){
                  type.display <- c(setdiff(type.display,"CI"),grep("^CI",names(outSummary),value=TRUE))
              }
              type.display <- intersect(names(outSummary),type.display)
              
              ## display
              table.print <- outSummary[,type.display,drop=FALSE]
              if("significance" %in% names(table.print)){
                  names(table.print)[names(table.print) == "significance"] <- ""
              }
              print(table.print, row.names = FALSE)
           
              return(invisible(NULL))
          }
          )

##----------------------------------------------------------------------
### show.R ends here
