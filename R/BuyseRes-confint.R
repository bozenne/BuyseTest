### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: jan 15 2019 (15:47) 
##           By: Brice Ozenne
##     Update #: 244
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
#' @include BuyseRes-object.R
#' 
#' @description Computes confidence intervals for net benefit statistic or the win ratio statistic.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' whereas \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016).
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}.
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
#' More precisely, the quantiles of the distribution of the statistic are computed under the null hypothesis and then shifted by the point estimate of the statistic.
#' Therefore it is possible that the limits of the confidence interval
#' are estimated outside of the interval of definition of the statistic (e.g. outside [-1,1] for the proportion in favor of treatment).
#'
#'
#' @return A matrix containing a column for the estimated statistic (over all strata),
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
                                statistic = NULL,
                                conf.level = NULL,
                                alternative = NULL,
                                method.boot = "percentile",
                                transformation = NULL){

              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              if(is.null(transformation)){
                  transformation <- option$transformation
              }
              if(is.null(conf.level)){
                  conf.level <- option$conf.level
              }
              if(is.null(alternative)){
                  alternative <- option$alternative
              }
              
              ## ** normalize and check arguments
              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio"),
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

              validLogical(transformation,
                           name1 = "transformation",
                           valid.length = 1,
                           method = "confint[BuyseRes]")
              
              ## ** extract information
              method.inference <- object@method.inference
              if(is.na(conf.level)){
                  method.inference <- "none"
              }

              count.favorable <- slot(object, name = paste0("count.favorable"))
              count.unfavorable <- slot(object, name = paste0("count.unfavorable"))
              n.pairs <- slot(object, name = paste0("n.pairs"))
              Delta <- slot(object, name = paste0("Delta.",statistic))
              Delta.resampling <- slot(object, name = paste0("DeltaResampling.",statistic))
              covariance <- object@covariance
              endpoint <- object@endpoint
              alpha <- 1-conf.level

              if(object@method.inference %in% c("asymptotic","asymptotic-beta")){
                  if(object@method.tte == "Peron"){
                      warning("The current implementation of the asymptotic distribution is not valid for method.tte=\"Peron\" \n",
                              "Standard errors / confidence intervals / p-values should not be trusted \n")
                  }else if(object@correction.uninf > 0){
                      warning("The current implementation of the asymptotic distribution is not valid when a correction is used \n",
                              "Standard errors / confidence intervals / p-values should not be trusted \n")
                  }
              }

              ## ** null hypothesis
              null <- switch(statistic,
                             "netBenefit" = 0,
                             "winRatio" = 1)


              method.confint <- switch(method.inference,
                                       "permutation" = confint_permutation,
                                       "stratified permutation" = confint_permutation,
                                       "bootstrap" = if(method.boot == "percentile"){confint_percentileBootstrap}else{confint_gaussianBootstrap},
                                       "stratified bootstrap" = if(method.boot == "percentile"){confint_percentileBootstrap}else{confint_gaussianBootstrap},
                                       "asymptotic" = confint_Ustatistic,
                                       "asymptotic-bebu" = confint_Ustatistic,
                                       "none" = confint_none 
                                       )

              ## ** compute the confidence intervals
              outConfint <- do.call(method.confint, args = list(Delta = Delta,
                                                                Delta.resampling = Delta.resampling,
                                                                pc.favorable = count.favorable / n.pairs,
                                                                pc.unfavorable = count.unfavorable / n.pairs,
                                                                statistic = statistic,
                                                                covariance = covariance,
                                                                alternative = alternative,
                                                                null = null,
                                                                alpha = alpha,
                                                                endpoint = endpoint,
                                                                transformation = transformation,
                                                                continuity.correction = option$continuity.correction,
                                                                n.pairs = n.pairs))

              ## ** number of permutations
              if(method.inference %in%  c("permutation","stratified permutation","bootstrap","stratified bootstrap")){
                  attr(outConfint, "n.resampling")  <- rowSums(!is.na(Delta.resampling))
              }else{
                  attr(outConfint, "n.resampling")  <- setNames(rep(as.numeric(NA), length(endpoint)), endpoint)
              }

              ## ** export              
              return(outConfint)
              
          })

## * confint_permutation (called by confint)
confint_permutation <- function(Delta, Delta.resampling,
                                null, alternative, alpha,
                                endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta

    ## ** computations
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** standard error
        outTable[iE,"se"] <- stats::sd(Delta.resampling[iE,], na.rm = TRUE)
        
        ## *** confidence interval
        qDelta_H0 <- switch(alternative,
                            "two.sided" = stats::quantile(Delta.resampling[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                            "less" = c(stats::quantile(Delta.resampling[iE,], probs = alpha,na.rm = TRUE), Inf),
                            "greater" = c(-Inf,stats::quantile(Delta.resampling[iE,], probs = 1 - alpha,na.rm = TRUE))
                            )
        outTable[iE,c("lower.ci","upper.ci")] <- Delta[iE] + (qDelta_H0 - null)

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
                                         "two.sided" = mean(abs(Delta[iE] - null) < abs(Delta.resampling[iE,] - null)),
                                         "less" = mean((Delta[iE] - null) > (Delta.resampling[iE,] - null)),
                                         "greater" = mean((Delta[iE] - null) < (Delta.resampling[iE,] - null))
                                         )
    
    }
    
    ## ** export
    return(outTable)
}

## * confint_percentileBootstrap (called by confint)
confint_percentileBootstrap <- function(Delta, Delta.resampling,
                                        null, alternative, alpha,
                                        endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta

    ## ** computations
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** standard error
        outTable[iE,"se"] <- stats::sd(Delta.resampling[iE,], na.rm = TRUE)

        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = stats::quantile(Delta.resampling[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                                                        "less" = c(stats::quantile(Delta.resampling[iE,], probs = alpha,na.rm = TRUE), Inf),
                                                        "greater" = c(-Inf,stats::quantile(Delta.resampling[iE,], probs = 1 - alpha,na.rm = TRUE))
                                                        )
        ## *** p.values
        outTable[iE, "p.value"] <- boot2pvalue(Delta.resampling[iE,], null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI)
        ## quantileCI(Delta.resampling[iE,], alternative = "two.sided", p.value = 0.64, sign.estimate = 1)
        ## quantileCI(Delta.resampling[iE,], alternative = "two.sided", p.value = 0.66, sign.estimate = 1)

    }

    ## ** export
    return(outTable)
}


## * confint_gaussianBootstrap (called by confint)
confint_gaussianBootstrap <- function(Delta, Delta.resampling,
                                      null, alternative, alpha,
                                      endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta
    
    ## ** CI + p
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified
        
        ## *** standard error
        iSE <- stats::sd(Delta.resampling[iE,], na.rm = TRUE)
        outTable[iE,"se"] <- iSE
        
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
confint_Ustatistic <- function(Delta, pc.favorable, pc.unfavorable, covariance, statistic, null,
                               alternative, alpha,
                               endpoint, transformation, continuity.correction, n.pairs, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta
    
    ## ** CI + p
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** standard error
        if(statistic == "netBenefit"){
            outTable[iE,"se"] <- sqrt(covariance[iE,"favorable"] + covariance[iE,"unfavorable"] - 2 * covariance[iE,"covariance"])

            if(continuity.correction){
                Delta[iE] <- Delta[iE] - sign(Delta[iE])/(2*n.pairs)
            }

            if(transformation){ ## atanh transform (also called fisher transform)
                iSE <- outTable[iE,"se"] / (1-Delta[iE]^2)
                iDelta <-  atanh(Delta[iE])
                backtransform <- tanh
                null <- atanh(null)
            }else{ ## on the original scale
                iSE <- outTable[iE,"se"]
                iDelta <- Delta[iE] ## 0.5/(n.y * n.x)
                backtransform <- function(x){x}
            }
            
            ## on the logit scale
        }else if(statistic == "winRatio"){
            
            outTable[iE,"se"] <- sqrt(covariance[iE,"favorable"]/pc.unfavorable[iE]^2 + covariance[iE,"unfavorable"]*Delta[iE]^2/pc.unfavorable[iE]^2 - 2 * covariance[iE,"covariance"]*Delta[iE]/pc.unfavorable[iE]^2)

            if(transformation){ ## log transform
                iSE <- outTable[iE,"se"] / Delta[iE]
                iDelta <-  log(Delta[iE])
                backtransform <- exp
                null <- log(null)
            }else{ ## on the original scale
                iSE <- outTable[iE,"se"]
                iDelta <- Delta[iE]
                backtransform <- function(x){x}
            }

        }

        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- backtransform(switch(alternative,
                                                                      "two.sided" = iDelta + stats::qnorm(c(alpha/2,1 - alpha/2)) * iSE,
                                                                      "less" = c(iDelta + stats::qnorm(alpha) * iSE, Inf),
                                                                      "greater" = c(-Inf,iDelta + stats::qnorm(1-alpha) * iSE)
                                                                      ))

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative,
                                         "two.sided" = 2*(1-stats::pnorm(abs((iDelta-null)/iSE))), ## 2*(1-pnorm(1.96))
                                         "less" = stats::pnorm((iDelta-null)/iSE), ## pnorm(1.96)
                                         "greater" = 1-stats::pnorm((iDelta-null)/iSE)
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

    ## ** point estimate
    outTable[,"estimate"] <- Delta

    ## ** return
    return(outTable)

    
}
##----------------------------------------------------------------------
### BuyseRes-confint.R ends here
