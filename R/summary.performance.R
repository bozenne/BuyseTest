### summary.performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  6 2022 (14:56) 
## Version: 
## Last-Updated: apr 12 2022 (14:42) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.performance
##' @title Summary Method for Performance Objects
##' @description Summary of the performance of binary classifiers
##'
##' @param object output of performance.
##' @param digits [numeric vector of length 2] number of digits used for the estimates and p-values.
##' @param print [logical] should the performance be printed in the console.
##' @param ... not used.
##'
##' @method summary performance
##' @export
summary.performance <- function(object, digits = c(3,3), print = TRUE, ...){
    df.print <- object$performance

    df.print$p.value <- base::format.pval(df.print$p.value, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value[is.na(object$performance$p.value)] <- ""
    df.print$p.value_comp <- base::format.pval(df.print$p.value_comp, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value_comp[is.na(object$performance$p.value_comp)] <- ""
    df.print <- df.print[,union(names(which(colSums(!is.na(object$performance))>0)),"estimate")]
    print(df.print, digits = digits[1])

    return(invisible(object$performance))
}

## * summary.performance
##' @method print performance
##' @export
print.performance <- function(x, ...){
    out <- summary(x)
    return(invisible(NULL))
}

##----------------------------------------------------------------------
### summary.performance.R ends here
