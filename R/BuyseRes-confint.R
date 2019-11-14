### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: nov 14 2019 (14:41) 
##           By: Brice Ozenne
##     Update #: 584
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
#' @aliases confint confint,BuyseRes-method
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
#' @param order.Hprojection [integer, 1-2] order of the H-decomposition used to compute the variance.
#' @param method.ci.resampling [character] the method used to compute the confidence intervals and p-values when using bootstrap or permutation.
#' \code{"percentile"} uses the quantiles of the empirical distribution,
#' \code{"gaussian"} uses the quantiles of a Gaussian distribution,
#' and \code{"student"} uses the quantiles of a Student's t-distribution (only available for bootstrap).
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
#' When using resampling methods:
#' \itemize{
#' \item an attribute \code{n.resampling} specified how many samples have been used to compute the confidence intervals and the p-values.
#' \item an attribute \code{method.ci.resampling} method used to compute the confidence intervals and p-values. 
#' }
#' 
#' @keywords confint BuyseRes-method
#' @author Brice Ozenne

## * Method - confint
#' @rdname BuyseRes-confint
#' @exportMethod confint
setMethod(f = "confint",
          signature = "BuyseRes",
          definition = function(object,
                                statistic = NULL,
                                conf.level = NULL,
                                alternative = NULL,
                                method.ci.resampling = NULL,
                                order.Hprojection = NULL,
                                transformation = NULL){

              option <- BuyseTest.options()
              D <- length(object@endpoint)
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

              if(is.null(method.ci.resampling)){                  
                  if(attr(method.inference,"bootstrap")){
                      if(attr(method.inference,"studentized")){
                          method.ci.resampling <- "studentized"
                      }else{
                          method.ci.resampling <- "percentile"
                      }
                  }else if(attr(method.inference,"permutation")){
                      method.ci.resampling <- "percentile"
                  }
              }else{
                  method.ci.resampling <- tolower(method.ci.resampling)
              }
              validCharacter(method.ci.resampling,
                             name1 = "method.ci.resampling",
                             valid.values = c("percentile","gaussian","studentized"),
                             valid.length = 1,
                             refuse.NULL = FALSE,                             
                             method = "confint[BuyseRes]")

              if(!is.null(method.ci.resampling)){
                  if(!attr(method.inference,"bootstrap") && !attr(method.inference,"permutation")){
                      warning("Argument \'method.ci.resampling\' is disregarded when not using resampling\n")
                  }else if(method.ci.resampling == "studentized" && !attr(method.inference,"studentized")){
                      stop("Argument \'method.ci.resampling\' cannot be set to \'studentized\' unless a studentized bootstrap has been performed\n",
                           "Consider setting \'method.ci.resampling\' to \"percentile\" or \"gaussian\" \n",
                           "or setting \'method.inference\' to \"studentized bootstrap\" when calling BuyseTest. \n")
                  }
              }

              if(attr(method.inference,"ustatistic") && !is.null(order.Hprojection) && order.Hprojection != attr(method.inference,"hprojection")){

                  validInteger(order.Hprojection,
                               name1 = "order.Hprojection",
                               min = 1, max = 2, valid.length = 1,
                               method = "confint[BuyseRes]")
                  
                  if(order.Hprojection > attr(method.inference,"hprojection")){
                      stop("Cannot find the second order of the H-decomposition. \n",
                           "Consider setting order.Hprojection to 2 in BuyseTest.options before calling BuyseTest. \n") 
                  }else{ ## modify covariance
                      ls.iid <- iid(object, endpoint = 1:D)
                      delta.favorable <- colSums(object@count.favorable)/object@n.pairs
                      delta.unfavorable <- colSums(object@count.unfavorable)/object@n.pairs
                      keep.names <- dimnames(object@covariance)
                          
                      object@covariance <- do.call(rbind, lapply(1:D, function(iE){ ## iE <- 1
                          iVar <- crossprod(ls.iid[[iE]])[c(1,4,2)]
                          return( c(iVar,
                                    iVar[1] + iVar[2] - 2*iVar[3],
                                    iVar[1]/delta.unfavorable^2 + iVar[2] * delta.favorable^2/delta.unfavorable^4 - 2*iVar[3] * delta.favorable/delta.unfavorable^3)
                                 )
                      }))
                      dimnames(object@covariance) <- keep.names
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
              test.model.tte <- all(unlist(lapply(object@iidNuisance,dim))==0)
              if(method.inference %in% c("u-statistic","u-statistic-bebu") && object@correction.uninf > 0){
                  warning("The current implementation of the asymptotic distribution has not been validated when using a correction. \n",
                          "Standard errors / confidence intervals / p-values may not be correct. \n",
                          "Consider using a resampling approach or checking the control of the type 1 error with powerBuyseTest. \n")
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
                  method.confint <- switch(method.ci.resampling,
                                           "percentile" = confint_permutation,
                                           "gaussian" = confint_gaussian)
                  transformation <- (statistic=="winRatio")
              }else if(attr(method.inference,"bootstrap")){
                  method.confint <- switch(method.ci.resampling,
                                           "percentile" = confint_percentileBootstrap,
                                           "gaussian" = confint_gaussian,
                                           "studentized" = confint_student)
                  if(method.ci.resampling=="percentile"){
                      transformation <- FALSE
                  }else if(transformation){
                      test <- switch(statistic,
                                     "netBenefit"= any(abs(Delta.resampling-1)<1e-12),
                                     "winRatio"= any(abs(Delta.resampling)<1e-12))
                      if(test){
                          warning("Extreme values of the statistic \n",
                                  "Consider setting the argument \'transformation\' to FALSE \n")
                      }
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
                                           "netBenefit" = function(x,se){
                                               if(is.null(se)){
                                                   out <- se
                                               }else{
                                                   out <- se/(1-x^2)
                                                   if(any(se==0)){
                                                       out[se==0] <- 0
                                                   }
                                               }
                                               return(out)
                                           },
                                           "winRatio" = function(x,se){
                                               if(is.null(se)){
                                                   out <- se
                                               }else{
                                                   out <- se/x
                                                   if(any(se==0)){
                                                       out[se==0] <- 0
                                                   }
                                               }
                                               return(out)
                                           })
                  itrans.se.delta <- switch(statistic,
                                            "netBenefit" = function(x,se){
                                                if(is.null(se)){
                                                    out <- se
                                                }else{
                                                    out <- se*(1-itrans.delta(x)^2)
                                                    if(any(se==0)){
                                                        out[se==0] <- 0
                                                    }
                                                }
                                                return(out)
                                            },
                                            "winRatio" = function(x,se){
                                                if(is.null(se)){
                                                    out <- se
                                                }else{
                                                    out <- se*itrans.delta(x)
                                                    if(any(se==0)){
                                                        out[se==0] <- 0
                                                    }
                                                }
                                                return(out)
                                            })
              }else{
                  trans.delta <- function(x){x}
                  itrans.delta <- function(x){x}
                  trans.se.delta <- function(x,se){se}
                  itrans.se.delta <- function(x,se){se}
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
                                                backtransform.se = itrans.se.delta))

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
              attr(outConfint,"method.ci.resampling") <- method.ci.resampling

              
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
    ## Delta.statH0.resampling <- sweep(Delta.resampling, MARGIN = 2, FUN = "-", STATS = Delta)/Delta.se.resampling
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
        iQ <- stats::quantile(x, probs = probs, na.rm = TRUE)
        return(Delta[iE] + iQ * Delta.se[iE])
    }

    for(iE in 1:n.endpoint){ ## iE <- 1
        ## if(sign(mean(Delta.resampling[,iE]))!=sign(Delta[iE])){
            ## warning("the estimate and the average bootstrap estimate do not have same sign \n")
        ## }
        outTable[iE, "p.value"] <- boot2pvalue(Delta.statH0.resampling[,iE], null = null, estimate = Delta[iE],
                                               alternative = alternative, FUN.ci = quantileCI2, checkSign = FALSE)
    }

    ## special case
    if(any(Delta.se==0)){
        index0 <- which(Delta.se==0)
        outTable[index0,"lower.ci"] <- outTable[index0,"estimate"]
        outTable[index0,"upper.ci"] <- outTable[index0,"estimate"]
        outTable[index0,"p.value"] <- as.numeric(outTable[index0,"estimate"]==null)
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

    ## special case with no variability
    if(any((Delta==null)*(Delta.se==0) == 1)){
        outTable[(Delta==null)*(Delta.se==0) == 1,"p.value"] <- 1
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
