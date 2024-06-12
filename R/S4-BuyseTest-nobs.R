### S4-BuyseTest-nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  3 2023 (10:00) 
## Version: 
## Last-Updated: jun  4 2024 (10:34) 
##           By: Brice Ozenne
##     Update #: 40
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
#' @param strata [character vector] the strata relative to which the number of pairs should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the total number of pairs (i.e. across all strata),
#' or \code{TRUE} to output each strata-specific number of pairs.
#' @param simplify [logical] should the result be coerced to the lowest possible dimension?
#' @param ... no used, for compatibility with the generic method.
#' 
#' @return A vector (when argument \code{strata} is \code{FALSE}) or a matrix (when argument \code{strata} is \code{TRUE}). In the latter case each line correspond to a strata.
#' 
#' @keywords methods
#' @author Brice Ozenne

## * Method - print
#' @rdname S4BuyseTest-nobs
#' @exportMethod nobs
setMethod(f = "nobs",
          signature = "S4BuyseTest",
          definition = function(object, strata = FALSE, simplify = TRUE, ...){

              ## ** normalize arguments
              indexC <- attr(object@level.treatment,"indexC")
              indexT <- attr(object@level.treatment,"indexT")

              level.strata <- object@level.strata
              index.strata <- attr(level.strata,"index")
              if(is.null(strata)){
                  if(length(level.strata)==1){
                      strata <- "global"                      
                  }else{
                      strata <- c("global", level.strata)
                  }
              }else if(identical(strata,FALSE)){
                  strata <- "global"
              }else if(identical(strata,TRUE)){
                  strata <- level.strata
              }else if(is.numeric(strata)){
                  validInteger(strata,
                               name1 = "strata",
                               valid.length = NULL,
                               min = 1,
                               max = length(level.strata),
                               refuse.NULL = TRUE,
                               refuse.duplicates = TRUE,
                               method = "autoplot[S4BuyseTest]")
                  strata <- level.strata[strata]
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",level.strata),
                                 refuse.NULL = FALSE,
                                 method = "coef[S4BuyseTest]")
              }

              ## ** extract
              Mout <- rbind(c(length(indexC), length(indexT), sum(object@n.pairs)),
                            cbind(sapply(lapply(index.strata, intersect, indexC), length),
                                  sapply(lapply(index.strata, intersect, indexT), length),
                                  object@n.pairs))
              rownames(Mout) <- c("global", level.strata)
              colnames(Mout) <- c(object@level.treatment, "pairs")
              
              
              ## ** export
              out <- Mout[strata,,drop=simplify]
              if(is.matrix(out)){
                  return(as.data.frame(out))
              }else{
                  return(out)
              }
          }
          )


##----------------------------------------------------------------------
### S4-BuyseTest-nobs.R ends here
