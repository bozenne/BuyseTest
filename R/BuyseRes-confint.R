### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: feb 25 2019 (19:33) 
##           By: Brice Ozenne
##     Update #: 335
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
#' @param method.ci [character] the method used to compute the confidence intervals and p-values.
#' \code{"percentile"} uses the quantiles of the empirical distribution,
#' \code{"gaussian"} uses the quantiles of a Gaussian distribution,
#' and \code{"student"} uses the quantiles of a Student's t-distribution.
#' Only relevant when using bootstrap resampling.
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
                                method.ci = NULL,
                                transformation = NULL){

              option <- BuyseTest.options()
              method.inference <- object@method.inference
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

              if(is.null(method.ci)){                  
                  if(attr(method.inference,"bootstrap")){
                      if(attr(method.inference,"studentized")){
                          method.ci <- "studentized"
                      }else{
                          method.ci <- "percentile"
                      }
                  }
              }else{
                  method.ci <- tolower(method.ci)
              }
              validCharacter(method.ci,
                             name1 = "method.ci",
                             valid.values = c("percentile","gaussian","studentized"),
                             valid.length = 1,
                             refuse.NULL = FALSE,                             
                             method = "confint[BuyseRes]")

              if(!is.null(method.ci)){
                  if(!attr(method.inference,"bootstrap")){
                      warning("Argument \'method.ci\' is disregarded when not using bootstrap resampling\n")
                  }else if(method.ci == "studentized" && !attr(method.inference,"studentized")){
                      stop("Argument \'method.ci\' cannot be set to \'studentized\' unless a studentized bootstrap has been performed\n",
                           "Consider setting \'method.ci\' to \"percentile\" or \"gaussian\" \n",
                           "or setting \'method.inference\' to \"studentized bootstrap\" or \"studentized stratified bootstrap\" when calling BuyseTest. \n")
                  }
              }

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
              if(is.na(conf.level)){
                  method.inference <- "none"
              }

              endpoint <- paste0(object@endpoint,"_",object@threshold)
              Delta <- slot(object, name = paste0("Delta.",statistic))
              Delta.resampling <- slot(object, name = paste0("DeltaResampling.",statistic))
              varianceDelta <- object@covariance
              varianceDelta.resampling <- object@covarianceResampling
              alpha <- 1-conf.level

              ## safety
              if(method.inference %in% c("asymptotic","asymptotic-bebu")){
                  if(object@method.tte == "Peron"){
                      warning("The current implementation of the asymptotic distribution is not valid for method.tte=\"Peron\" \n",
                              "Standard errors / confidence intervals / p-values will not be displayed \n")
                  }else if(object@correction.uninf > 0){
                      warning("The current implementation of the asymptotic distribution is not valid when a correction is used \n",
                              "Standard errors / confidence intervals / p-values will not be displayed \n")
                  }
              }
              
              ## ** null hypothesis
              null <- switch(statistic,
                             "netBenefit" = 0,
                             "winRatio" = 1)

              if(method.inference == "none"){
                  method.confint <- confint_none
              }else if(attr(method.inference,"ustatistic")){
                  method.confint <- confint_Ustatistic
              }else if(attr(method.inference,"permutation")){
                  method.confint <- confint_permutation
              }else if(attr(method.inference,"bootstrap")){
                  method.confint <- switch(method.ci,
                                           "percentile" = confint_percentileBootstrap,
                                           "gaussian" = confint_gaussian,
                                           "studentized" = confint_student)
              } 

              ## ** compute the confidence intervals
              outConfint <- do.call(method.confint, args = list(Delta = Delta,
                                                                Delta.resampling = Delta.resampling,
                                                                statistic = statistic,
                                                                varianceDelta = varianceDelta,
                                                                varianceDelta.resampling = varianceDelta.resampling,
                                                                alternative = alternative,
                                                                null = null,
                                                                alpha = alpha,
                                                                endpoint = endpoint,
                                                                transformation = transformation))

              ## ** number of permutations
              if(method.inference %in%  c("permutation","stratified permutation","bootstrap","stratified bootstrap")){
                  attr(outConfint, "n.resampling")  <- rowSums(!is.na(Delta.resampling))
              }else{
                  attr(outConfint, "n.resampling")  <- setNames(rep(as.numeric(NA), length(endpoint)), endpoint)
              }
              attr(outConfint,"method.ci.resampling") <- method.ci

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
        outTable[iE,"se"] <- stats::sd(Delta.resampling[,iE], na.rm = TRUE)
        
        ## *** confidence interval
        qDelta_H0 <- switch(alternative,
                            "two.sided" = stats::quantile(Delta.resampling[,iE], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                            "less" = c(-Inf,stats::quantile(Delta.resampling[,iE], probs = 1 - alpha,na.rm = TRUE),
                            "greater" = c(stats::quantile(Delta.resampling[,iE], probs = alpha,na.rm = TRUE), Inf))
                            )
        outTable[iE,c("lower.ci","upper.ci")] <- Delta[iE] + (qDelta_H0 - null)

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
                                         "two.sided" = mean(abs(Delta[iE] - null) < abs(Delta.resampling[,iE] - null)),
                                         "less" = mean((Delta[iE] - null) < (Delta.resampling[,iE] - null)),
                                         "greater" = mean((Delta[iE] - null) > (Delta.resampling[,iE] - null))
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
        outTable[iE,"se"] <- stats::sd(Delta.resampling[,iE], na.rm = TRUE)

        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = stats::quantile(Delta.resampling[,iE], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                                                        "less" = c(-Inf,stats::quantile(Delta.resampling[,iE], probs = 1 - alpha,na.rm = TRUE),
                                                        "greater" = c(stats::quantile(Delta.resampling[,iE], probs = alpha,na.rm = TRUE), Inf))
                                                        )
        ## *** p.values
        outTable[iE, "p.value"] <- boot2pvalue(Delta.resampling[,iE], null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI)
        ## quantileCI(Delta.resampling[,iE], alternative = "two.sided", p.value = 0.64, sign.estimate = 1)
        ## quantileCI(Delta.resampling[,iE], alternative = "two.sided", p.value = 0.66, sign.estimate = 1)

    }

    ## ** export
    return(outTable)
}


## * confint_gaussian (called by confint)
confint_gaussian <- function(Delta, Delta.resampling,
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
        iSE <- stats::sd(Delta.resampling[,iE], na.rm = TRUE)
        outTable[iE,"se"] <- iSE
        
        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- switch(alternative,
                                                        "two.sided" = Delta[iE] + stats::qnorm(c(alpha/2,1 - alpha/2)) * iSE,
                                                        "less" = c(-Inf,Delta[iE] + stats::qnorm(1-alpha) * iSE),
                                                        "greater" = c(Delta[iE] + stats::qnorm(alpha) * iSE, Inf)
                                                        )

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative,
                                         "two.sided" = 2*(1-stats::pnorm(abs(Delta[iE]/iSE - null))), ## 2*(1-pnorm(1.96))
                                         "less" = stats::pnorm(Delta[iE]/iSE - null),
                                         "greater" = 1-stats::pnorm(Delta[iE]/iSE - null) ## pnorm(1.96)
                                         )
    }


    ## ** export
    return(outTable)
}

## * confint_student (called by confint)
confint_student <- function(Delta, varianceDelta, Delta.resampling, varianceDelta.resampling,
                            null, alternative, alpha,
                            endpoint, statistic, transformation, ...){

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
            outTable[iE,"se"] <- sqrt(varianceDelta[iE,"netBenefit"])
            
            if(transformation){ ## atanh transform (also called fisher transform)
                iSE <- outTable[iE,"se"] / (1-Delta[iE]^2)               
                iDelta <-  atanh(Delta[iE])
                iSE.resampling <- sqrt(varianceDelta.resampling[,iE,statistic]) / (1-Delta.resampling[,iE]^2)
                iDelta.resampling <- atanh(Delta.resampling[,iE]) - mean(atanh(Delta.resampling[,iE]))
                backtransform <- tanh
                null <- atanh(null)
            }else{ ## on the original scale
                iSE <- outTable[iE,"se"]
                iDelta <- Delta[iE] ## 0.5/(n.y * n.x)
                iSE.resampling <- sqrt(varianceDelta.resampling[,iE,statistic])
                iDelta.resampling <- Delta.resampling[,iE]-mean(Delta.resampling[,iE])
                backtransform <- function(x){x}
            }
                        
        }else if(statistic == "winRatio"){
            outTable[iE,"se"] <- sqrt(varianceDelta[iE,"winRatio"])
            
            if(transformation){ ## log transform
                iSE <- outTable[iE,"se"] / Delta[iE]
                iDelta <-  log(Delta[iE])
                iSE.resampling <- sqrt(varianceDelta.resampling[,iE,statistic]) / Delta.resampling[,iE]
                iDelta.resampling <- log(Delta.resampling[,iE]) - mean(log(Delta.resampling[,iE]))
                backtransform <- exp
                null <- log(null)
            }else{ ## on the original scale
                iSE <- outTable[iE,"se"]
                iDelta <- Delta[iE]
                iSE.resampling <- sqrt(varianceDelta.resampling[,iE,statistic])
                iDelta.resampling <- Delta.resampling[,iE]-mean(Delta.resampling[,iE])
                backtransform <- function(x){x}
            }

        }

        ## *** critical quantile
        ## plot(iDelta.resampling,iSE.resampling)
        qBoot <- switch(alternative,
                        "two.sided" = quantile(iDelta.resampling/iSE.resampling, probs = c(alpha/2,1 - alpha/2)),
                        "less" = quantile(iDelta.resampling/iSE.resampling, probs = alpha),
                        "greater" = quantile(iDelta.resampling/iSE.resampling, probs = 1 - alpha/2)
                        )
        ## qBoot
       

        
        ## *** confidence interval
        outTable[iE,c("lower.ci","upper.ci")] <- backtransform(switch(alternative,
                                                                      "two.sided" = Delta[iE] + qBoot * iSE,
                                                                      "less" = c(-Inf,Delta[iE] + qBoot * iSE),
                                                                      "greater" = c(Delta[iE] + qBoot * iSE, Inf)
                                                                      ))

        ## *** p.value
        quantileCI2 <- function(x, alternative, p.value, sign.estimate, ...){
            probs <- switch(alternative,
                            "two.sided" = c(p.value/2,1-p.value/2)[2-sign.estimate], ## if positive p.value/2 otherwise 1-p.value/2
                            "less" = 1-p.value,
                            "greater" = p.value)
            iQ <- quantile(x, probs = probs)
            return(Delta[iE] + iQ * iSE)
        }
        ## quantileCI2(iDelta.resampling/iSE.resampling, alternative = alternative, p.value = 0.05, sign.estimate = (iDelta>=0))
        outTable[iE, "p.value"] <- boot2pvalue(iDelta.resampling/iSE.resampling, null = null, estimate = iDelta,
                                               alternative = alternative, FUN.ci = quantileCI2)
    }


    ## ** export
    return(outTable)
}


## * confint_Ustatistic (called by confint)
confint_Ustatistic <- function(Delta, varianceDelta, statistic, null,
                               alternative, alpha,
                               endpoint, transformation, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta
    
    ## ** CI + p
    for(iE in 1:n.endpoint){
        if(is.infinite(Delta[iE]) || is.na(Delta[iE])){next} ## do not compute CI or p-value when the estimate has not been identified

        ## *** standard error
        outTable[iE,"se"] <- sqrt(varianceDelta[iE])

        if(statistic == "netBenefit"){

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
                        
        }else if(statistic == "winRatio"){

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
                                                                      "less" = c(-Inf,iDelta + stats::qnorm(1-alpha) * iSE),
                                                                      "greater" = c(iDelta + stats::qnorm(alpha) * iSE, Inf)
                                                                      ))

        ## *** p.value
        outTable[iE,"p.value"] <- switch(alternative,
                                         "two.sided" = 2*(1-stats::pnorm(abs((iDelta-null)/iSE))), ## 2*(1-pnorm(1.96))
                                         "less" = stats::pnorm((iDelta-null)/iSE),
                                         "greater" = 1-stats::pnorm((iDelta-null)/iSE) ## pnorm(1.96)
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
