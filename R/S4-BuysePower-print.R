## * Documentation - print
#' @docType methods
#' @name S4BuysePower-print
#' @title Print Method for Class "S4BuysePower"
#' @aliases print,S4BuysePower-method
#' @include S4-BuysePower.R S4-BuysePower-summary.R
#' 
#' @description Display the main results stored in a \code{S4BuysePower} object.
#' 
#' @param x an \R object of class \code{S4BuysePower}, i.e., output of \code{\link{powerBuyseTest}}
#' @param ... additional arguments passed to the summary method.
#' 
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing power calculation based on GPC. \cr
#'   \code{\link{S4BuysePower-summary}} for a more detailed presentation of the \code{S4BuysePower} object.
#'  
#' @return invisible table
#' @keywords print
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuysePower-print
#' @exportMethod print
setMethod(f = "print",
          signature = "S4BuysePower",
          definition = function(x, ...){

              ## compute summary statistics
              outSummary <- summary(x, print = FALSE, ...)
              print(outSummary[[1]], row.names = FALSE, quote = FALSE)
              return(invisible(outSummary$table))
          }
)
