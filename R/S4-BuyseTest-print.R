## * Documentation - print
#' @docType methods
#' @name S4BuyseTest-print
#' @title Print Method for Class "S4BuyseTest"
#' @aliases print,S4BuyseTest-method
#' @include S4-BuyseTest.R S4-BuyseTest-summary.R
#' 
#' @description Display the main results stored in a \code{S4BuyseTest} object.
#' 
#' @param x an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' @param ... additional arguments passed to the summary method.
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#'  
#' @keywords print
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuyseTest-print
#' @exportMethod print
setMethod(f = "print",
          signature = "S4BuyseTest",
          definition = function(x, ...){

              ## compute summary statistics
              outSummary <- summary(x, print = FALSE, ...)
              ## remove significance column
              table.print <- outSummary[,setdiff(names(outSummary), "significance"),drop=FALSE]
              print(table.print, row.names = FALSE)
           
              return(invisible(table.print))
          }
)
