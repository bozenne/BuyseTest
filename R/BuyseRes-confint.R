### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: mar  9 2019 (10:39) 
##           By: Brice Ozenne
##     Update #: 473
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
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param method.ci.boot [character] the method used to compute the confidence intervals and p-values.
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
#'
#' @return A matrix containing a column for the estimated statistic (over all strata),
#' the lower bound and upper bound of the confidence intervals, and the associated p-values.
#' When using resampling methods,
#' an attribute \code{n.resampling} specified how many samples have been used to compute the confidence intervals and the p-values.
#' When using boostrap,
#' an attribute \code{method.ci.boot} method used to compute the confidence intervals and p-values. 
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
                                method.ci.boot = NULL,
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

              if(is.null(method.ci.boot)){                  
                  if(attr(method.inference,"bootstrap")){
                      if(attr(method.inference,"studentized")){
                          method.ci.boot <- "studentized"
                      }else{
                          method.ci.boot <- "percentile"
                      }
                  }
              }else{
                  method.ci.boot <- tolower(method.ci.boot)
              }
              validCharacter(method.ci.boot,
                             name1 = "method.ci.boot",
                             valid.values = c("percentile","gaussian","studentized"),
                             valid.length = 1,
                             refuse.NULL = FALSE,                             
                             method = "confint[BuyseRes]")

              if(!is.null(method.ci.boot)){
                  if(!attr(method.inference,"bootstrap")){
                      warning("Argument \'method.ci.boot\' is disregarded when not using bootstrap resampling\n")
                  }else if(method.ci.boot == "studentized" && !attr(method.inference,"studentized")){
                      stop("Argument \'method.ci.boot\' cannot be set to \'studentized\' unless a studentized bootstrap has been performed\n",
                           "Consider setting \'method.ci.boot\' to \"percentile\" or \"gaussian\" \n",
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
                  method.inference[] <- "none" ## uses [] to not remove the attributees of method.inference
              }

              endpoint <- paste0(object@endpoint,"_",object@threshold)
              Delta <- slot(object, name = paste0("Delta.",statistic))
              Delta.resampling <- slot(object, name = paste0("DeltaResampling.",statistic))

              if(sum(dim(object@covariance))>0){
                  Delta.se <- sqrt(object@covariance[,statistic])
              }else{
                  Delta.se <- NULL
              }
              if(sum(dim(object@covarianceResampling))>0){
                  Delta.se.resampling <- sqrt(object@covarianceResampling[,,statistic])
              }else{
                  Delta.se.resampling <- NULL
              }
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

              ## ** method
              if(method.inference == "none"){
                  method.confint <- confint_none
                  transformation <- FALSE
              }else if(attr(method.inference,"ustatistic")){
                  method.confint <- confint_Ustatistic
              }else if(attr(method.inference,"permutation")){
                  method.confint <- confint_permutation
                  transformation <- (statistic=="winRatio")
              }else if(attr(method.inference,"bootstrap")){
                  method.confint <- switch(method.ci.boot,
                                           "percentile" = confint_percentileBootstrap,
                                           "gaussian" = confint_gaussian,
                                           "studentized" = confint_student)
                  if(method.ci.boot=="percentile"){
                      transformation <- FALSE
                  }
              } 

              ## ** transformation
              if(transformation){
                  trans.delta <- switch(statistic,
                                        "netBenefit" = atanh,
                                        "winRatio" = log)
                  itrans.delta <- switch(statistic,
                                         "netBenefit" = tanh,
                                         "winRatio" = exp)                  
                  trans.se.delta <- switch(statistic,
                                           "netBenefit" = function(x,se){if(!is.null(se)){se/(1-x^2)}else{se}},
                                           "winRatio" = function(x,se){if(!is.null(se)){se/x}else{se}})
                  itrans.se <- switch(statistic,
                                           "netBenefit" = function(x,se){if(!is.null(se)){se*(1-itrans.delta(x)^2)}else{se}},
                                           "winRatio" = function(x,se){if(!is.null(se)){se*itrans.delta(x)}else{se}})
              }else{
                  trans.delta <- function(x){x}
                  itrans.delta <- function(x){x}
                  trans.se.delta <- function(x,se){se}
                  itrans.se <- function(x,se){se}
              }

              ## ** compute the confidence intervals
              outConfint <- do.call(method.confint,
                                    args = list(Delta = trans.delta(Delta),
                                                Delta.resampling = trans.delta(Delta.resampling),
                                                Delta.se = trans.se.delta(Delta, se = Delta.se),
                                                Delta.se.resampling = trans.se.delta(Delta.resampling, se = Delta.se.resampling),
                                                alternative = alternative,
                                                null = trans.delta(null),
                                                alpha = alpha,
                                                endpoint = endpoint,
                                                backtransform.delta = itrans.delta,
                                                backtransform.se = itrans.se))

              ## do not output CI or p-value when the estimate has not been identified
              index.NA <- union(which(is.infinite(outConfint[,"estimate"])),which(is.na(outConfint[,"estimate"])))
              if(length(index.NA)>0){
                  outConfint[index.NA,c("se","lower.ci","upper.ci","p.value")] <- NA
              }

              ## ** number of permutations
              if(method.inference != "none" && (attr(method.inference,"permutation") || attr(method.inference,"bootstrap"))){
                  attr(outConfint, "n.resampling")  <- colSums(!is.na(Delta.resampling))
              }else{
                  attr(outConfint, "n.resampling")  <- setNames(rep(as.numeric(NA), length(endpoint)), endpoint)
              }
              attr(outConfint,"method.ci.boot") <- method.ci.boot

              
              ## ** export
              if(attr(method.inference,"permutation")){
                  attr(outConfint,"warning") <- "Confidence intervals are computed under the null hypothesis"
              }
              return(outConfint)
              
          })

## * confint_permutation (called by confint)
confint_permutation <- function(Delta, Delta.resampling,
                                null, alternative, alpha,
                                backtransform.delta, endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))
    
    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- apply(backtransform.delta(Delta.resampling), MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    ## ** confidence interval
    outTable[,"lower.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE),
                                                        "less" = -Inf,
                                                        "greater" = Delta + apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = alpha, na.rm = TRUE)
                                                        ))
    
    outTable[,"upper.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2, na.rm = TRUE),
                                                        "less" = Delta + apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha, na.rm = TRUE),
                                                        "greater" = Inf
                                                        ))

    ## ** p-value
    outTable[,"p.value"] <- sapply(1:n.endpoint, FUN = function(iE){ ## iE <- 1
        switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
               "two.sided" = mean(abs(Delta[iE] - null) <= abs(Delta.resampling[,iE] - null)),
               "less" = mean((Delta[iE] - null) >= (Delta.resampling[,iE] - null)),
               "greater" = mean((Delta[iE] - null) <= (Delta.resampling[,iE] - null))
               )
    })

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

    ## ** standard error
    outTable[,"se"] <- apply(Delta.resampling, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    ## ** confidence interval
    outTable[,"lower.ci"] <- switch(alternative,
                                    "two.sided" = apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE),
                                    "less" = -Inf,
                                    "greater" = apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = alpha, na.rm = TRUE)
                                    )
    
    outTable[,"upper.ci"] <- switch(alternative,
                                    "two.sided" = apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2, na.rm = TRUE),
                                    "less" = apply(Delta.resampling, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha, na.rm = TRUE),
                                    "greater" = Inf
                                    )

    ## ** p.values
    for(iE in 1:n.endpoint){
        outTable[iE, "p.value"] <- boot2pvalue(na.omit(Delta.resampling[,iE]), null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI)
    }
    ## quantileCI(Delta.resampling[,iE], alternative = "two.sided", p.value = 0.64, sign.estimate = 1)
    ## quantileCI(Delta.resampling[,iE], alternative = "two.sided", p.value = 0.66, sign.estimate = 1)

    

    ## ** export
    return(outTable)
}


## * confint_gaussian (called by confint)
confint_gaussian <- function(Delta, Delta.resampling,
                             null, alternative, alpha,
                             endpoint, backtransform.delta, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    Delta.se <- apply(Delta.resampling, MARGIN = 2, FUN = stats::sd, na.rm = TRUE) ## computed based on the sample
    outTable[,"se"] <- apply(backtransform.delta(Delta.resampling), MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    ## ** confidence interval
    outTable[,"lower.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + stats::qnorm(alpha/2) * Delta.se,
                                                        "less" = -Inf,
                                                        "greater" = Delta + stats::qnorm(alpha) * Delta.se
                                                        ))
    
    outTable[,"upper.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + stats::qnorm(1-alpha/2) * Delta.se,
                                                        "less" = Delta + stats::qnorm(1-alpha) * Delta.se,
                                                        "greater" = Inf
                                                        ))

    ## ** p-value
    outTable[,"p.value"] <- switch(alternative,
                                   "two.sided" = 2*(1-stats::pnorm(abs((Delta-null)/Delta.se))), 
                                   "less" = stats::pnorm((Delta-null)/Delta.se),
                                   "greater" = 1-stats::pnorm((Delta-null)/Delta.se) 
                                   )

    ## ** export
    return(outTable)
}

## * confint_student (called by confint)
confint_student <- function(Delta, Delta.se, Delta.resampling, Delta.se.resampling,
                            null, alternative, alpha,
                            endpoint, backtransform.delta, backtransform.se, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** critical quantile
    Delta.stat.resampling <- Delta.resampling/Delta.se.resampling
    Delta.statH0.resampling <- apply(Delta.stat.resampling, MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)
    Delta.qInf <- switch(alternative,
                         "two.sided" = apply(Delta.statH0.resampling, MARGIN = 2, FUN = stats::quantile, na.rm = TRUE, probs = alpha/2),
                         "less" = -Inf,
                         "greater" = apply(Delta.statH0.resampling, MARGIN = 2, FUN = stats::quantile, na.rm = TRUE, probs = alpha)
                         )
    Delta.qSup <- switch(alternative,
                         "two.sided" = apply(Delta.statH0.resampling, MARGIN = 2, FUN = stats::quantile, na.rm = TRUE, probs = 1-alpha/2),
                         "less" = apply(Delta.statH0.resampling, MARGIN = 2, FUN = stats::quantile, na.rm = TRUE, probs = 1-alpha),
                         "greater" = Inf
                         )    

    ## ** confidence interval
    outTable[,"lower.ci"] <- backtransform.delta(Delta + Delta.qInf * Delta.se)
    outTable[,"upper.ci"] <- backtransform.delta(Delta + Delta.qSup * Delta.se)

    ## ** p.value
    quantileCI2 <- function(x, alternative, p.value, sign.estimate, ...){
        probs <- switch(alternative,
                        "two.sided" = c(p.value/2,1-p.value/2)[2-sign.estimate], ## if positive p.value/2 otherwise 1-p.value/2
                        "less" = 1-p.value,
                        "greater" = p.value)
        iQ <- stats::quantile(x - mean(x), probs = probs, na.rm = TRUE)
        return(Delta[iE] + iQ * Delta.se[iE])
    }

    for(iE in 1:n.endpoint){
        outTable[iE, "p.value"] <- boot2pvalue(Delta.stat.resampling[,iE], null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI2)
    }

    ## ** export
    return(outTable)



}


## * confint_Ustatistic (called by confint)
confint_Ustatistic <- function(Delta, Delta.se, statistic, null,
                               alternative, alpha,
                               endpoint, backtransform.delta, backtransform.se, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 5,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** confidence interval
    outTable[,"lower.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + stats::qnorm(alpha/2) * Delta.se,
                                                        "less" = -Inf,
                                                        "greater" = Delta + stats::qnorm(alpha) * Delta.se
                                                        ))
    
    outTable[,"upper.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + stats::qnorm(1-alpha/2) * Delta.se,
                                                        "less" = Delta + stats::qnorm(1-alpha) * Delta.se,
                                                        "greater" = Inf
                                                        ))

    ## ** p-value
    outTable[,"p.value"] <- switch(alternative,
                                   "two.sided" = 2*(1-stats::pnorm(abs((Delta-null)/Delta.se))), 
                                   "less" = stats::pnorm((Delta-null)/Delta.se),
                                   "greater" = 1-stats::pnorm((Delta-null)/Delta.se) 
                                   )

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
