### S4-BuyseTest-nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  3 2023 (10:00) 
## Version: 
## Last-Updated: Jul  3 2023 (10:17) 
##           By: Brice Ozenne
##     Update #: 22
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
#' @name S4BuyseTest-nobs
#' @title Sample Size for Class "S4BuyseTest"
#' @aliases nobs,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Display the sample size in each treatmnet arm as well as the number of pairs.
#' 
#' @param object an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' @param stratified [logical] should strata-specific sample size be output?
#' Otherwise output the sample size and number of pairs over all strata.
#' @param ... no used, for compatibility with the generic method.
#' 
#' @return A vector (when argument \code{stratified} is \code{FALSE}) or a data.frame (when argument \code{stratified} is \code{TRUE}). In the latter case each line correspond to a strata.
#' 
#' @keywords methods
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuyseTest-nobs
#' @exportMethod nobs
setMethod(f = "nobs",
          signature = "S4BuyseTest",
          definition = function(object, stratified = NULL, ...){

              ## ** normalize arguments
              indexC <- attr(object@level.treatment,"indexC")
              indexT <- attr(object@level.treatment,"indexT")

              if(!is.null(stratified)){
                  level.strata <- object@level.strata
                  if(identical(stratified,TRUE)){
                      strata <- level.strata
                  }else{
                      strata <- match.arg(stratified, level.strata, several.ok = TRUE)
                  }
                  index.strata <- attr(level.strata,"index")
              }

              ## ** extract
              if(is.null(stratified)){
                  out <- c(length(indexC),
                           length(indexT),
                           sum(object@n.pairs))
                  names(out) <- c(object@level.treatment, "pairs")
              }else{
                  out <- data.frame(sapply(lapply(index.strata[match(strata, level.strata)], intersect, indexC), length),
                                    sapply(lapply(index.strata[match(strata, level.strata)], intersect, indexT), length),
                                    object@n.pairs[strata])
                  
                  names(out) <- c(object@level.treatment, "pairs")
                  rownames(out) <- strata
              }
              return(out)
          }
)


##----------------------------------------------------------------------
### S4-BuyseTest-nobs.R ends here
