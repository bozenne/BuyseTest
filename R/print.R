### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  6 2022 (14:56) 
## Version: 
## Last-Updated: jul  3 2026 (14:41) 
##           By: Brice Ozenne
##     Update #: 61
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.auc
#' @exportMethod print
print.BuyseTestAuc <- function(x, ...){
    ##    if(attr(x,"n.fold")==0){
        print.data.frame(x, ...)
    ## }else{
    ##     label.upper <- paste0(attr(x,"contrast")[2],">",attr(x,"contrast")[1])
    ##     label.lower <- paste0(attr(x,"contrast")[1],">",attr(x,"contrast")[2])
    ##     x$direction <- sapply(x$direction, function(iD){
    ##         if(iD==">"){return(label.upper)}else if(iD=="<"){return(label.lower)}else{return(iD)}
    ##     })
    ##     print.data.frame(x[x$fold == "global",c("direction","estimate","se","lower","upper","p.value")], row.names = FALSE)
    ## }
}

## * print.BuyseTestBrier
##' @exportMethod print
print.BuyseTestBrier <- function(x, ...){
    print.data.frame(x)
}

## * print.performance
##' @method print performance
##' @export
print.performance <- function(x, ...){
    out <- summary(x)
    return(invisible(NULL))
}

## * print S4BuysePower
#' @docType methods
#' @name S4BuysePower-print
#' @title Print Method for Class "S4BuysePower"
#' @aliases print,S4BuysePower-method
#' @include S4-BuysePower.R summary.R
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

## * print S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-print
#' @title Print Method for Class "S4BuyseTest"
#' @aliases print,S4BuyseTest-method
#' @include S4-BuyseTest.R summary.R
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
##----------------------------------------------------------------------
### print.R ends here
