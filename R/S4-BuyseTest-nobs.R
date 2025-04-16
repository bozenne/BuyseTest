### S4-BuyseTest-nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  3 2023 (10:00) 
## Version: 
## Last-Updated: apr  3 2025 (13:41) 
##           By: Brice Ozenne
##     Update #: 47
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
#' @description Display the sample size in each treatment arm as well as the number of pairs.
#' 
#' @param object an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' @param strata [character vector] the strata relative to which the number of pairs should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the total number of pairs (i.e. across all strata),
#' or \code{TRUE} to output each strata-specific number of pairs.
#' @param resampling [logical] should the sample size of each bootstrap or permutation sample be output?
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
          definition = function(object, strata = FALSE, resampling = FALSE, simplify = TRUE, ...){

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
                               method = "nobs[S4BuyseTest]")
                  strata <- level.strata[strata]
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",level.strata),
                                 refuse.NULL = FALSE,
                                 method = "nobs[S4BuyseTest]")
              }

              ## *** resampling
              if(length(resampling)==1 && is.logical(resampling)){
                  if(resampling==TRUE){
                      resampling <- 1:dim(object@nResampling)[1]
                  }
              }else if(length(resampling)==1&&resampling==0){
                  resampling <- FALSE
              }else{
                  validInteger(resampling,
                               name1 = "resampling",
                               valid.length = NULL,
                               min = 1,
                               max = dim(object@nResampling)[1],
                               refuse.NULL = TRUE,
                               refuse.duplicates = TRUE,
                               method = "nobs[S4BuyseTest]")
              }
              

              ## ** extract
              if(any(resampling>0)){
                  Mout <- do.call(rbind,lapply(level.strata, function(iStrata){
                      data.frame(sample = 1:NROW(object@nResampling),
                                 strata = iStrata,
                                 object@nResampling[,,iStrata],
                                 pairs = object@nResampling[,1,iStrata]*object@nResampling[,2,iStrata])
                  }))
                  Mout <- rbind(Mout,do.call(rbind,by(Mout, Mout$sample, function(iDF){
                      data.frame(sample = iDF[1,"sample"], strata = "global", as.list(colSums(iDF[-(1:2)])))
                  })))
              }else{
                  Mout <- rbind(c(length(indexC), length(indexT), sum(object@n.pairs)),
                                cbind(sapply(lapply(index.strata, intersect, indexC), length),
                                      sapply(lapply(index.strata, intersect, indexT), length),
                                      object@n.pairs))
                  rownames(Mout) <- c("global", level.strata)
                  colnames(Mout) <- c(object@level.treatment, "pairs")
              }
              
              ## ** export
              if(any(resampling>0)){
                  out <- Mout[Mout$strata %in% strata & Mout$sample %in% resampling,,drop=FALSE]
                  if(simplify){
                      if(length(strata)==1 && length(resampling)==1){
                          out <- unlist(as.vector(out[-(1:2)]))
                      }else if(length(strata)==1){
                          out <- as.matrix(out[-(1:2)])
                          rownames(out) <- resampling
                      }else if(length(strata)==1){
                          out <- as.matrix(out[-(1:2)])
                          rownames(out) <- strata
                      }
                  }
                  return(out)
              }else{
                  out <- Mout[strata,,drop=simplify]
                  if(is.matrix(out)){
                      return(as.data.frame(out))
                  }else{
                      return(out)
                  }
              }
          }
          )


##----------------------------------------------------------------------
### S4-BuyseTest-nobs.R ends here
