### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:26) 
## Version: 
## Last-Updated: jul  3 2026 (14:44) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .vcov.wglm
#' @title Variance-covariance matrix for Logistic Regressions
#' @description Compute the variance covariance matrix (i.e. inverse of the information) for logistic regressions.
#' @noRd
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total variance-covariance?
#' @param center [logical] should the individual contribution be centered around the average?
#' 
.vcov.logit <- function(object, indiv, center){
    X <- stats::model.matrix(object)
    n.obs <- NROW(X)
    n.param <- NCOL(X)
    pi <- stats::predict(object, type = "response")
    W <- object$prior.weights
    tXWpi <- t(.colMultiply_cpp(X, scale = W*pi*(1-pi)))
    Sigma <- solve(tXWpi %*% X)
    
    if(indiv){
        if(center){
            out <- array(unlist(lapply(1:n.obs, function(i){Sigma %*% tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE] %*% Sigma - Sigma/n.obs})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }else{
            out <- array(unlist(lapply(1:n.obs, function(i){Sigma %*% tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE] %*% Sigma})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }
        return(out)
        ## apply(out, MARGIN = 1:2, FUN = sum) - tXWpi %*% X
        
    }else{
        out <- Sigma
        rownames(out) <- colnames(out)
        return(out)
    }
}

## * vcov S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-vcov
#' @title Extract Uncertainty from GPC
#' @aliases vcov,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'  
#' @description Extract uncertainty about the summary statistics (net benefit, win ratio, ...) from GPC.
#' 
#' @param object a \code{S4BuyseTest} object, output of \code{\link{BuyseTest}}.
#' @param statistic [character] the statistic summarizing the pairwise comparison relative to which the variance is to be output: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param endpoint [character] for which endpoint(s) the variance-covariance matrix should be output?
#' If \code{NULL} consider all endpoints.
#' @param strata [character vector] the strata relative to which the variance-covariance matrix should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the statistic pooled over all strata.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param ... ignored.
#'
#' @details When consider a second order H-decomposition, it is only used to estimate variance: the correlation between endpoints is evaluated using the first order H-decomposition.
#' 
#' @return A numeric matrix.
#' 
#' @keywords method
#' @author Brice Ozenne
#' @exportMethod vcov
setMethod(f = "vcov",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                strata = FALSE,
                                cumulative = TRUE,
                                ...){

              ## ** normalize arguments
              option <- BuyseTest.options()
              mycall <- match.call()
              dots <- list(...)
              if(length(dots)>0){
                  stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
              }
              method.inference <- object@method.inference
              permutation <- attr(method.inference,"permutation")
              bootstrap <- attr(method.inference,"bootstrap")
              hprojection <- attr(method.inference,"hprojection")
              
              ## ** vcov
              if(permutation){
                  stop("Cannot evaluate the variance-covariance matrix under the alternative hypothesis when using permutations. \n")
              }else if(bootstrap){
                  ## *** extract bootstrap estimates
                  object.boot <- coef(object, endpoint = endpoint, statistic = statistic, strata = strata, 
                                      cumulative = cumulative, resampling = TRUE)
                  out <- stats::var(object.boot)
              }else{
                  ## *** extract influence function
                  object.iid <- getIid(object, endpoint = endpoint, statistic = statistic, strata = strata, 
                                       cumulative = cumulative, type = "all", simplify = FALSE)
                  if(length(object.iid)>1){
                      stop("Method vcov cannot output the variance-covariance matrix relative to multiple strata at once. \n")
                  }else{
                      colnames(object.iid[[1]]) <- names(coef(object, endpoint = endpoint, statistic = statistic, strata = strata, cumulative = cumulative))
                  }

                  ## *** evaluate variance-covariance matrix
                  out <- crossprod(object.iid[[1]])
                  if(hprojection!=1){
                      if(cumulative == TRUE && (all(strata=="global") || all(strata==FALSE))){ ## rescale 
                          se2 <- confint(object, endpoint = endpoint, statistic = statistic, strata = strata, cumulative = cumulative)$se
                          out <- stats::cov2cor(out) * tcrossprod(se2)
                      }else{
                          ## Warning already output by getIid: "Inference will be performed using a first order H projection"
                      }
                  }
              }

              ## ** export
              return(out)

          })

##----------------------------------------------------------------------
### vcov.R ends here
