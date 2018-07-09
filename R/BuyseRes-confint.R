### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: jul  9 2018 (13:00) 
##           By: Brice Ozenne
##     Update #: 141
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - confint
#' @docType methods
#' @name BuyseRes-confint
#' @title  Confidence Intervals for Model Parameters
#' @aliases confing confint,BuyseRes-method
#' @include BuyseRes-object.R BuyseRes-summary.R
#' 
#' @description Computes confidence intervals for net chance statistic or the win ratio statisitc.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netChance"} displays the net chance in favor of treatment, as described in Buyse (2010) and Peron et al. (2016)),
#' whereas \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param method.boot [character] the method used to compute the boostrap confidence intervals and p-values.
#' Can be \code{"percentile"} for computing the CI using the quantiles of the boostrap distribution or
#' \code{"gaussian"} for using a Gaussian approximation to compute the CI where the standard error is computed using the bootstrap samples.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-summary}} for a more detailed presentation of the \code{BuyseRes} object.
#' 
#' @details 
#' When using a permutation test, the uncertainty associated with the estimator is computed under the null hypothesis.
#' Thus the confidence interval may not be valid if the null hypothesis is false. \cr
#' More precisely, the quantiles of the distribution of the statistic are computed under the null hypothesis and then shifted by the punctual estimate of the statistic.
#' Therefore it is possible that the limits of the confidence interval
#' are estimated outside of the interval of definition of the statistic (e.g. outside [-1,1] for the proportion in favor of treatment).
#'
#'
#' @return A matrix containing a column for the estimated statstic (over all strata),
#' the lower bound and upper bound of the confidence intervals, and the associated p-values.
#' When using resampling methods,
#' an attribute \code{n.resampling} specified how many samples have been used to compute the confidence intervals and the p-values.
#' 
#' 
#' @keywords confint BuyseRes-method

## * Method - confint
#' @rdname BuyseRes-confint
#' @exportMethod confint
setMethod(f = "confint",
          signature = "BuyseRes",
          definition = function(object,
                                statistic = BuyseTest.options()$statistic,
                                conf.level = 0.95,
                                alternative = "two.sided",
                                method.boot = "percentile"){

              ## ** normalize and check arguments
              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netchance" = "netChance",
                                  "winratio" = "winRatio",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netChance","winRatio"),
                             valid.length = 1,
                             method = "confint[BuyseRes]")

              method.boot <- tolower(method.boot)
              validCharacter(method.boot,
                             name1 = "method.boot",
                             valid.values = c("percentile","gaussian"),
                             valid.length = 1,
                             method = "confint[BuyseRes]")

              validNumeric(conf.level,
                           name1 = "conf.level",
                           min = 0, max = 1,
                           refuse.NA = FALSE,
                           valid.length = 1,
                           method = "confint[BuyseRes]")

              validCharacter(alternative,
                             name1 = "alternative",
                             valid.values = c("two.sided","less","greater"),
                             valid.length = 1,
                             method = "confint[BuyseRes]")

              ## ** extract information
              method.inference <- object@method.inference
              if(is.na(conf.level)){
                  method.inference <- "none"
              }

              Delta <- slot(object, name = paste0("Delta.",statistic))
              Delta.permutation <- slot(object, name = paste0("DeltaResampling.",statistic))
              covariance <- object@covariance
              endpoint <- object@endpoint
              alpha <- 1-conf.level
              
              ## ** null hypothesis
              null <- switch(statistic,
                             "netChance" = 0,
                             "winRatio" = 1)


              method.confint <- switch(method.inference,
                                       "permutation" = confint_permutation,
                                       "stratified permutation" = confint_permutation,
                                       "bootstrap" = if(method.boot == "percentile"){confint_percentileBootstrap}else{confint_gaussianBootstrap},
                                       "stratified bootstrap" = if(method.boot == "percentile"){confint_percentileBootstrap}else{confint_gaussianBootstrap},
                                       "asymptotic" = confint_Ustatistic,
                                       "none" = confint_none 
                                       )

              ## ** compute the confidence intervals
              outConfint <- do.call(method.confint, args = list(Delta = Delta,
                                                                Delta.permutation = Delta.permutation,
                                                                covariance = covariance,
                                                                alternative = alternative,
                                                                null = null,
                                                                alpha = alpha,
                                                                endpoint = endpoint))

              ## ** number of permutations
              if(method.inference %in%  c("permutation","stratified permutation","bootstrap","stratified bootstrap")){
                  attr(outConfint, "n.resampling")  <- rowSums(!is.na(Delta.permutation))
              }else{
                  attr(outConfint, "n.resampling")  <- setNames(rep(as.numeric(NA), length(endpoint)), endpoint)
              }

              ## ** export              
              return(outConfint)
              
          })

## * confint_permutation (called by confint)
confint_permutation <- function(Delta, Delta.permutation,
                                null, alternative, alpha,
                                endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 4,
                       dimnames = list(endpoint, c("estimate","lower.ci","upper.ci","p.value")))

    ## ** punctual estimate
    outTable[,"estimate"] <- Delta

    ## ** computations
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified
        
        ## *** confidence interval
        qDelta_H0 <- switch(alternative,
                            "two.sided" = stats::quantile(Delta.permutation[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                            "less" = c(stats::quantile(Delta.permutation[iE,], probs = alpha,na.rm = TRUE), Inf),
                            "greater" = c(-Inf,stats::quantile(Delta.permutation[iE,], probs = 1 - alpha,na.rm = TRUE))
                            )
        outTable[iE,c("lower.ci","upper.ci")] <- Delta[iE] + (qDelta_H0 - null)

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                                         "two.sided" = mean(abs(Delta[iE] - null) < abs(Delta.permutation[iE,] - null)),
                                         "less" = mean((Delta[iE] - null) > (Delta.permutation[iE,] - null)),
                                         "greater" = mean((Delta[iE] - null) < (Delta.permutation[iE,] - null))
                                         )
    
    }
    
    ## ** export
    return(outTable)
}

## * confint_percentileBootstrap (called by confint)
confint_percentileBootstrap <- function(Delta, Delta.permutation,
                                        null, alternative, alpha,
                                        endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 4,
                       dimnames = list(endpoint, c("estimate","lower.ci","upper.ci","p.value")))

    ## ** punctual estimate
    outTable[,"estimate"] <- Delta

    ## ** computations
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = stats::quantile(Delta.permutation[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                                                        "less" = c(stats::quantile(Delta.permutation[iE,], probs = alpha,na.rm = TRUE), Inf),
                                                        "greater" = c(-Inf,stats::quantile(Delta.permutation[iE,], probs = 1 - alpha,na.rm = TRUE))
                                                        )
        ## *** p.values
        outTable[iE, "p.value"] <- boot2pvalue(Delta.permutation[iE,], null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI)
        ## quantileCI(Delta.permutation[iE,], alternative = "two.sided", p.value = 0.64, sign.estimate = 1)
        ## quantileCI(Delta.permutation[iE,], alternative = "two.sided", p.value = 0.66, sign.estimate = 1)

    }

    ## ** export
    return(outTable)
}


## * confint_gaussianBootstrap (called by confint)
confint_gaussianBootstrap <- function(Delta, Delta.permutation,
                                      null, alternative, alpha,
                                      endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 4,
                       dimnames = list(endpoint, c("estimate","lower.ci","upper.ci","p.value")))

    ## ** punctual estimate
    outTable[,"estimate"] <- Delta
    
    ## ** CI + p
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified
        
        ## *** standard error
        iSE <- stats::sd(Delta.permutation[iE,], na.rm = TRUE)
        
        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = Delta[iE] + stats::qnorm(c(alpha/2,1 - alpha/2)) * iSE,
                                                        "less" = c(Delta[iE] + stats::qnorm(alpha) * iSE, Inf),
                                                        "greater" = c(-Inf,Delta[iE] + stats::qnorm(1-alpha) * iSE)
                                                        )

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative,
                                         "two.sided" = 2*(1-stats::pnorm(abs(Delta[iE]/iSE - null))), ## 2*(1-pnorm(1.96))
                                         "less" = stats::pnorm(Delta[iE]/iSE - null), ## pnorm(1.96)
                                         "greater" = 1-stats::pnorm(Delta[iE]/iSE - null)
                                         )
    }


    ## ** export
    return(outTable)
}
                  

## * confint_Ustatistic (called by confint)
confint_Ustatistic <- function(Delta, covariance,
                               null, alternative, alpha,
                               endpoint, ...){

    warning("In development - do not trust the results \n")
    
    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 4,
                       dimnames = list(endpoint, c("estimate","lower.ci","upper.ci","p.value")))

    ## ** punctual estimate
    outTable[,"estimate"] <- Delta
    
    ## ** CI + p
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** standard error
        iSE <- sqrt(covariance[iE,"favorable"] + covariance[iE,"unfavorable"] - 2*covariance[iE,"covariance"])
        
        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = Delta[iE] + stats::qnorm(c(alpha/2,1 - alpha/2)) * iSE,
                                                        "less" = c(Delta[iE] + stats::qnorm(alpha) * iSE, Inf),
                                                        "greater" = c(-Inf,Delta[iE] + stats::qnorm(1-alpha) * iSE)
                                                        )

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative,
                                         "two.sided" = 2*(1-stats::pnorm(abs(Delta[iE]/iSE - null))), ## 2*(1-pnorm(1.96))
                                         "less" = stats::pnorm(Delta[iE]/iSE - null), ## pnorm(1.96)
                                         "greater" = 1-stats::pnorm(Delta[iE]/iSE - null)
                                         )
    }

    ## ** export
    return(outTable)
}

## * confint_none (called by confint)
confint_none <- function(Delta, endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(NA, nrow = n.endpoint, ncol = 4,
                       dimnames = list(endpoint, c("estimate","lower.ci","upper.ci","p.value")))

    ## ** punctual estimate
    outTable[,"estimate"] <- Delta

    ## ** return
    return(outTable)

    
}
##----------------------------------------------------------------------
### BuyseRes-confint.R ends here
