### rbind.performanceResample.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 21 2022 (10:05) 
## Version: 
## Last-Updated: apr 22 2022 (10:08) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rbind.performanceResample (documentation)
##' @title Combine Resampling Results For Performance Objects
##' @description Combine permutation or bootstrap samples.
##' Useful to run parallel calculations (see example below).
##'
##' @param ... performance objects.
##' @param tolerance [numeric] maximum acceptable difference between the point estimates.
##' Can be \code{NA} to skip this sanity check.
##'
##' @examples
##' if(FALSE){
##'
##' #### simulate data ####
##' set.seed(10)
##' n <- 100
##' df.train <- data.frame(Y = rbinom(n, prob = 0.5, size = 1),
##'                        X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n), X4 = rnorm(n),
##'                        X5 = rnorm(n), X6 = rnorm(n), X7 = rnorm(n), X8 = rnorm(n),
##'                        X9 = rnorm(n), X10 = rnorm(n))
##' df.train$Y <- rbinom(n, size = 1,
##'                      prob = 1/(1+exp(-df.train$X5 - df.train$X6 - df.train$X7)))
##'
##' #### fit models ####
##' e.null <- glm(Y~1, data = df.train, family = binomial(link="logit"))
##' e.logit <- glm(Y~X1+X2, data = df.train, family = binomial(link="logit"))
##' e.logit2 <- glm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = df.train,
##'                family = binomial(link="logit"))
##'
##' #### evaluate model (same seed) ####
##' fold.repetition <- 5 ## 0: internal perf (fast)
##'                      ## >0: 10 fold CV repeated (slow)
##' test <- performanceResample(list(e.logit,e.logit2), seed = 10,
##'                              fold.repetition = fold.repetition, n.resampling = 100)
##' test.1 <- performanceResample(list(e.logit,e.logit2), seed = 10,
##'                              fold.repetition = fold.repetition, n.resampling = 1:50)
##' test.2 <- performanceResample(list(e.logit,e.logit2), seed = 10,
##'                              fold.repetition = fold.repetition, n.resampling = 51:100)
##' rbind(test.1,test.2)
##' test
##'
##' ## Note: when the prediction model call RNG then test.1 and test.2 may not give test 
##' 
##' #### evaluate model (different seed) ####
##' test.3 <- performanceResample(list(e.logit,e.logit2), seed = 11,
##'                              fold.repetition = fold.repetition, n.resampling = 1:50)
##' test.4 <- performanceResample(list(e.logit,e.logit2), seed = 12,
##'                              fold.repetition = fold.repetition, n.resampling = 51:100)
##' rbind(test.3,test.4, tolerance = NA) ## does not check equality of the point estimate
##'                                      ## between test.3 and test.4
##' test
##' }

## * rbind.performanceResample (code)
##' @export
rbind.performance <- function(..., tolerance = 1e-5){

    ls.perf <- list(...)
    out <- ls.perf[[1]]
    if(length(ls.perf)==1){
        return(out)
    }
    if(any(sapply(ls.perf, function(iX){is.null(iX[["resampling"]])}))){
        stop("Performance objects should contain permutation or boostrap samples. \n",
             "(i.e. be an output of the performanceResample function) \n")
    }
    if(!is.na(tolerance)){
        test <- lapply(ls.perf[-1], function(x){
            all.equal(x$performance[,c("method","metric","model","estimate")],
                      ls.perf[[1]]$performance[,c("method","metric","model","estimate")],
                      tolerance = 1e-5)
        })
        test.logical <- sapply(test, function(iT){any(iT!=TRUE)})
        if(any(test.logical)){
            stop("Discrepancy between the point estimates: 1 vs. ",paste(which(test.logical)+1, collapse = ", ")," \n")
        }
        
    }
    ls.resampling <- lapply(ls.perf, "[[", "resampling")
    if(any(duplicated(unlist(lapply(ls.resampling, function(iL){unique(iL$sample)}))))){
        stop("Same sample name among different performance objects. \n")
    }

    out$resampling <- do.call(rbind,ls.resampling)
    out$resampling$sample <- as.numeric(factor(out$resampling$sample))
    out$resampling <- out$resampling[order(out$resampling$sample),]
    out$args$n.resampling <- max(out$resampling$sample)
    out$performance <- .performanceResample_inference(performance = out$performance[,c("method","metric","model","estimate")],
                                                      resampling = out$resampling,
                                                      type.resampling = out$args$type.resampling,
                                                      conf.level = out$args$conf.level)

    return(out)
}

##----------------------------------------------------------------------
### rbind.performance.R ends here
