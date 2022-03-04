### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
## Last-Updated: mar  4 2022 (10:16) 
##           By: Brice Ozenne
##     Update #: 861
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
#' @name S4BuyseTest-confint
#' @title  Confidence Intervals for Model Parameters
#' @aliases confint,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Computes confidence intervals for net benefit statistic or the win ratio statistic.
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' \code{"unfavorable"} displays the proportion in favor of the control.
#' Default value read from \code{BuyseTest.options()}.
#' @param endpoint [character] for which endpoint(s) the confidence intervals should be output?
#' If \code{NULL} returns the confidence intervals for all endpoints.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param null [numeric] right hand side of the null hypothesis (used for the computation of the p-value).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param order.Hprojection [integer, 1-2] order of the H-decomposition used to compute the variance.
#' @param method.ci.resampling [character] the method used to compute the confidence intervals and p-values when using bootstrap or permutation (\code{"percentile"}, \code{"gaussian"}, \code{"student"}).
#' See the details section.
#' @param cluster [numeric vector] Group of observations for which the iid assumption holds .
#'  
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#' 
#' @details
#' \bold{statistic}: when considering a single endpoint and denoting
#' \eqn{Y} the endpoint in the treatment group,
#' \eqn{X} the endpoint in the control group,
#' and \eqn{\tau} the threshold of clinical relevance,
#' the net benefit is \eqn{P[Y \ge X + \tau] - P[X \ge Y + \tau]},
#' the win ratio is \eqn{\frac{P[Y \ge X + \tau]}{P[X \ge Y + \tau]}},
#' the proportion in favor of treatment is \eqn{P[Y \ge X + \tau]},
#' the proportion in favor of control is \eqn{P[X \ge Y + \tau]}.
#' 
#' \bold{method.ci.resampling}: when using bootstrap/permutation, p-values and confidence intervals are computing as follow: \itemize{
#' \item \code{percentile} (bootstrap): compute the confidence interval using the quantiles of the bootstrap estimates.
#' Compute the p-value by finding the confidence level at which a bound of the confidence interval equals the null hypothesis.
#' 
#' \item \code{percentile} (permutation): apply the selected transformation to the estimate and permutation estimates.
#' Compute the confidence interval by (i) shfiting the estimate by the quantiles of the centered permutation estimates and (ii) back-transforming .
#' Compute the p-value as the relative frequency at which the estimate are less extreme than the permutation estimates.
#'
#' \item \code{gaussian} (bootstrap and permutation): apply the selected transformation to the estimate and bootstrap/permutation estimates.
#' Estimate the variance of the estimator using the empirical variance of the transformed boostrap/permutation estimates.
#' Compute confidence intervals and p-values under the normality assumption and back-transform the confidence intervals.
#' 
#' \item \code{student} (bootstrap): apply the selected transformation to the estimate, its standard error, the bootstrap estimates, and their standard error.
#' Compute the studentized bootstrap estimates by dividing the centered bootstrap estimates by their standard error. 
#' Compute the confidence interval based on the standard error of the estimate and the quantiles of the studentized bootstrap estimates, and back-transform.
#' Compute the p-value by finding the confidence level at which a bound of the confidence interval equals the null hypothesis.
#' 
#' \item \code{student} (permutation): apply the selected transformation to the estimate, its standard error, the permutation estimates, and their standard error.
#' Compute the studentized permutation estimates by dividing the centered permutation estimates by their standard error.
#' Compute the confidence interval based on the standard error of the estimate and the quantiles of the studentized permutation estimates, and back-transform.
#' Compute the p-value as the relative frequency at which the studentized estimate are less extreme than the permutation studentized estimates.
#'
#' }
#' 
#' \bold{WARNING}: when using a permutation test, the uncertainty associated with the estimator is computed under the null hypothesis.
#' Thus the confidence interval may not be valid if the null hypothesis is false. \cr
#'
#' @return A matrix containing a column for the estimated statistic (over all strata),
#' the lower bound and upper bound of the confidence intervals, and the associated p-values.
#' When using resampling methods:
#' \itemize{
#' \item an attribute \code{n.resampling} specified how many samples have been used to compute the confidence intervals and the p-values.
#' \item an attribute \code{method.ci.resampling} method used to compute the confidence intervals and p-values. 
#' }
#' 
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937 \cr
#'
#' @keywords confint S4BuyseTest-method
#' @author Brice Ozenne

## * Method - confint
#' @rdname S4BuyseTest-confint
#' @exportMethod confint
setMethod(f = "confint",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                cumulative = TRUE,
                                null = NULL,
                                conf.level = NULL,
                                alternative = NULL,
                                method.ci.resampling = NULL,
                                order.Hprojection = NULL,
                                transformation = NULL,
                                cluster = NULL){

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
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                             valid.length = 1,
                             method = "confint[S4BuyseTest]")

              if(attr(method.inference,"permutation") || attr(method.inference,"bootstrap")){
                  if(is.null(method.ci.resampling)){                  
                      if(attr(method.inference,"studentized")){
                          method.ci.resampling <- "studentized"
                      }else{
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
                                 method = "confint[S4BuyseTest]")

                  if(method.ci.resampling == "studentized" && !attr(method.inference,"studentized")){
                      stop("Argument \'method.ci.resampling\' cannot be set to \'studentized\' unless a studentized bootstrap/permutation has been performed.\n",
                           "Consider setting \'method.ci.resampling\' to \"percentile\" or \"gaussian\" \n",
                           "or setting \'method.inference\' to \"studentized bootstrap\" or \"studentized permutation\" when calling BuyseTest. \n")
                  }
                  if(method.ci.resampling == "studentized" && cumulative == FALSE){
                      stop("Endpoint specific confidence intervals are not available with studentized bootstrap/permutation. \n",
                           "Consider applying the BuyseTest function separately to each endpoint \n",
                           "or set \'method.inference\' to \"studentized bootstrap\" or \"studentized permutation\" when calling BuyseTest. \n")
                  }
              }else if(!is.null(method.ci.resampling)){
                  warning("Argument \'method.ci.resampling\' is disregarded when not using resampling\n")                  
              }

              if(attr(method.inference,"ustatistic")){
                  if(!is.null(order.Hprojection) && order.Hprojection != attr(method.inference,"hprojection") || !is.null(cluster)){

                      if(!is.null(order.Hprojection)){
                          validInteger(order.Hprojection,
                                       name1 = "order.Hprojection",
                                       min = 1, max = 2, valid.length = 1,
                                       method = "confint[S4BuyseTest]")
                  
                          if(order.Hprojection > attr(method.inference,"hprojection")){
                              stop("Cannot find the second order of the H-decomposition. \n",
                                   "Consider setting order.Hprojection to 2 in BuyseTest.options before calling BuyseTest. \n")
                          }
                      }
                      if(identical(order.Hprojection,2) && !is.null(cluster)){
                          warning("Inference will be performed using a first order H projection. \n")
                      }
                      ls.iid <- getIid(object, endpoint = 1:D, cluster = cluster, statistic = c("favorable","unfavorable"), cumulative = cumulative)
                      delta.favorable <- coef(object, statistic = "favorable", cumulative = cumulative)
                      delta.unfavorable <- coef(object, statistic = "unfavorable", cumulative = cumulative)
                      keep.names <- dimnames(object@covariance)

                      ls.iid <- setNames(lapply(1:D, function(iD){
                          cbind(ls.iid[[iD]],
                                "covariance" = ls.iid[[iD]][,"favorable"]*ls.iid[[iD]][,"unfavorable"],
                                "netBenefit" = ls.iid[[iD]][,"favorable"] - ls.iid[[iD]][,"unfavorable"],
                                "winRatio" = ls.iid[[iD]][,"favorable"]/delta.unfavorable[iD] - ls.iid[[iD]][,"unfavorable"]*delta.favorable[iD]/delta.unfavorable[iD]^2 
                                )
                      }), names(ls.iid))
                      object@covariance <- do.call(rbind,lapply(ls.iid,function(iIID){
                          c(colSums(iIID[,c("favorable","unfavorable")]^2),covariance = sum(iIID[,c("covariance")]),colSums(iIID[,c("netBenefit","winRatio")]^2))
                      }))
                      dimnames(object@covariance) <- keep.names
                      ## object@covariance <- do.call(rbind, lapply(1:D, function(iE){ ## iE <- 1
                      ##     iVar <- crossprod(ls.iid[[iE]])[c(1,4,2)]
                      ##     return( c(iVar,
                      ##               iVar[1] + iVar[2] - 2*iVar[3],
                      ##               iVar[1]/delta.unfavorable[iE]^2 + iVar[2] * delta.favorable[iE]^2/delta.unfavorable[iE]^4 - 2*iVar[3] * delta.favorable[iE]/delta.unfavorable[iE]^3)
                      ##            )
                      ## }))
                      Delta.iid <- lapply(ls.iid, function(iIID){iIID[,statistic,drop=FALSE]})
                  }else{
                      Delta.iid <- getIid(object, endpoint = 1:D, cluster = cluster, statistic = statistic, cumulative = cumulative)
                  }
                  
              }
                       
              validNumeric(conf.level,
                           name1 = "conf.level",
                           min = 0, max = 1,
                           refuse.NA = FALSE,
                           valid.length = 1,
                           method = "confint[S4BuyseTest]")

              validCharacter(alternative,
                             name1 = "alternative",
                             valid.values = c("two.sided","less","greater"),
                             valid.length = 1,
                             method = "confint[S4BuyseTest]")

              validLogical(transformation,
                           name1 = "transformation",
                           valid.length = 1,
                           method = "confint[S4BuyseTest]")

              ## endpoint
              if(!is.null(endpoint)){
                  valid.endpoint <- names(object@endpoint)
                  n.endpoint <- length(valid.endpoint)
                  if(is.numeric(endpoint)){
                      validInteger(endpoint,
                                   name1 = "endpoint",
                                   min = 1, max = length(valid.endpoint),
                                   valid.length = NULL,
                                   method = "iid[BuyseTest]")
                      endpoint <- valid.endpoint[endpoint]
                  }else{
                      validCharacter(endpoint,
                                     valid.length = 1:length(valid.endpoint),
                                     valid.values = valid.endpoint,
                                     refuse.NULL = FALSE)
                  }
              }

              ## ** extract information
              if(is.na(conf.level)){
                  method.inference[] <- "none" ## uses [] to not remove the attributees of method.inference
              }

              all.endpoint <- names(object@endpoint)
              Delta <- coef(object, statistic = statistic, cumulative = cumulative, stratified = FALSE)

              if(attr(method.inference,"permutation") || attr(method.inference,"bootstrap")){
                  Delta.resampling <- slot(object, name = "DeltaResampling")[,,statistic]
                  if(!is.matrix(Delta.resampling)){
                      Delta.resampling <- matrix(Delta.resampling, ncol = length(all.endpoint),
                                                 dimnames = list(NULL, all.endpoint))
                  }
                  if(attr(method.inference,"studentized")){
                      Delta.se.resampling <- sqrt(object@covarianceResampling[,,statistic])
                      if(!is.matrix(Delta.se.resampling)){
                          Delta.se.resampling <- matrix(Delta.se.resampling, ncol = length(all.endpoint),
                                                        dimnames = list(NULL, all.endpoint))
                      }
                  }else{
                      Delta.se.resampling <- NULL
                  }
                  if(cumulative==FALSE &&  length(all.endpoint)>1){
                      ## only for point estimate as studentized returns an error
                      for(iE in 2:length(all.endpoint)){
                          if(statistic %in% c("netbenefit","favorable","unfavorable")){
                              Delta.resampling[,iE] <- Delta.resampling[,iE] - Delta.resampling[,iE-1]
                          }else if(statistic == "winratio"){
                              Delta.resampling[,iE] <- (slot(object, name = "DeltaResampling")[,iE,"favorable"]-slot(object, name = "DeltaResampling")[,iE-1,"favorable"])/(slot(object, name = "DeltaResampling")[,iE,"unfavorable"]-slot(object, name = "DeltaResampling")[,iE-1,"unfavorable"])
                          }
                      }
                  }
              }else{
                  Delta.resampling <- NULL
                  Delta.se.resampling <- NULL
              }
              
              if(attr(method.inference,"ustatistic") || attr(method.inference,"studentized")){
                  if(cumulative){
                      Delta.se <- sqrt(object@covariance[,statistic])
                  }else{
                      if(attr(method.inference,"hprojection")==2){
                          warning("Inference will be performed using a first order H projection. \n")
                      }
                      Delta.se <- sqrt(unlist(lapply(Delta.iid, function(iIID){sum(iIID^2)})))
                  }
              }else{
                  Delta.se <- NULL
              }              
              alpha <- 1-conf.level

              ## safety
              test.model.tte <- all(unlist(lapply(object@iidNuisance,dim))==0)
              if(method.inference %in% c("u-statistic","u-statistic-bebu") && object@correction.uninf > 0){
                  warning("The current implementation of the asymptotic distribution is not valid when using a correction. \n",
                          "Standard errors / confidence intervals / p-values may not be correct. \n",
                          "Consider using a resampling approach or checking the control of the type 1 error with powerBuyseTest. \n")
              }
              
              ## ** null hypothesis
              if(is.null(null)){
                  null <- switch(statistic,
                                 "netBenefit" = 0,
                                 "winRatio" = 1,
                                 "favorable" = 1/2,
                                 "unfavorable" = 1/2)
              }else{
                  validNumeric(null, valid.length = 1,
                               min = if("statistic"=="netBenefit"){-1}else{0},
                               max = if("statistic"=="winRatio"){Inf}else{1})
              }

              ## ** method
              if(method.inference == "none"){
                  method.confint <- confint_none
                  transformation <- FALSE
              }else if(attr(method.inference,"ustatistic")){
                  method.confint <- confint_Ustatistic
              }else if(attr(method.inference,"permutation")){
                  method.confint <- switch(method.ci.resampling,
                                           "percentile" = confint_percentilePermutation,
                                           "gaussian" = confint_gaussian,
                                           "studentized" = confint_studentPermutation)
              }else if(attr(method.inference,"bootstrap")){
                  method.confint <- switch(method.ci.resampling,
                                           "percentile" = confint_percentileBootstrap,
                                           "gaussian" = confint_gaussian,
                                           "studentized" = confint_studentBootstrap)
                  if(method.ci.resampling=="percentile"){
                      transformation <- FALSE
                  }
              }
              
              ## ** transformation
              if(transformation){
                  if(object@hierarchical){
                      trans.weight <- 1
                  }else{
                      trans.weight <- sum(object@weight)
                  }

                  trans.delta <- switch(statistic,
                                        "netBenefit" = function(x){if(is.null(x)){x}else{atanh(x/trans.weight)}},
                                        "winRatio" = function(x){if(is.null(x)){x}else{log(x)}},
                                        "favorable" = function(x){if(is.null(x)){x}else{atanh(2*(x/trans.weight-1/2))}},
                                        "unfavorable" = function(x){if(is.null(x)){x}else{atanh(2*(x/trans.weight-1/2))}}
                                        )
                  itrans.delta <- switch(statistic,                                         
                                         "netBenefit" = function(x){if(is.null(x)){x}else{trans.weight*tanh(x)}}, 
                                         "winRatio" = function(x){if(is.null(x)){x}else{exp(x)}},
                                         "favorable" = function(x){if(is.null(x)){x}else{trans.weight*(tanh(x)/2+1/2)}},
                                         "unfavorable" = function(x){if(is.null(x)){x}else{trans.weight*(tanh(x)/2+1/2)}}
                                         )                  
                  trans.se.delta <- switch(statistic,
                                           "netBenefit" = function(x,se){
                                               if(is.null(se)){
                                                   out <- se
                                               }else{
                                                   out <- (se/trans.weight)/(1-(x/trans.weight)^2)
                                                   if(any(na.omit(se)==0)){
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
                                                   if(any(na.omit(se)==0)){
                                                       out[se==0] <- 0
                                                   }
                                               }
                                               return(out)
                                           },
                                           "favorable" = function(x,se){
                                               if(is.null(se)){
                                                   out <- se
                                               }else{
                                                   out <- 2*(se/trans.weight)/(1-(2*(x/trans.weight-1/2))^2)
                                                   if(any(na.omit(se)==0)){
                                                       out[se==0] <- 0
                                                   }
                                               }
                                               return(out)
                                           },
                                           "unfavorable" = function(x,se){
                                               if(is.null(se)){
                                                   out <- se
                                               }else{
                                                   out <- 2*(se/trans.weight)/(1-(2*(x/trans.weight-1/2))^2)
                                                   if(any(na.omit(se)==0)){
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
                                                    out <- trans.weight*se*(1-(itrans.delta(x)/trans.weight)^2)
                                                    if(any(na.omit(se)==0)){
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
                                                    if(any(na.omit(se)==0)){
                                                        out[se==0] <- 0
                                                    }
                                                }
                                                return(out)
                                            },
                                            "favorable" = function(x,se){
                                                if(is.null(se)){
                                                    out <- se
                                                }else{
                                                    out <- trans.weight*(se/2)*(1-(2*(itrans.delta(x)/trans.weight-1/2))^2)
                                                    if(any(na.omit(se)==0)){
                                                        out[se==0] <- 0
                                                    }
                                                }
                                                return(out)
                                            },
                                            "unfavorable" = function(x,se){
                                                if(is.null(se)){
                                                    out <- se
                                                }else{
                                                    out <- trans.weight*(se/2)*(1-(2*(itrans.delta(x)/trans.weight-1/2))^2)
                                                    if(any(na.omit(se)==0)){
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
                                                endpoint = all.endpoint,
                                                backtransform.delta = itrans.delta,
                                                backtransform.se = itrans.se.delta))

              ## do not output CI or p-value when the estimate has not been identified
              index.NA <- union(which(is.infinite(outConfint[,"estimate"])),which(is.na(outConfint[,"estimate"])))
              if(length(index.NA)>0){
                  outConfint[index.NA,c("se","lower.ci","upper.ci","p.value")] <- NA
              }
              if(!is.null(endpoint)){
                  outConfint <- as.data.frame(outConfint[all.endpoint %in% endpoint,,drop=FALSE])
              }else{
                  outConfint <- as.data.frame(outConfint)
              }
              
              ## ** number of permutations
              if(method.inference != "none" && (attr(method.inference,"permutation") || attr(method.inference,"bootstrap"))){
                  attr(outConfint, "n.resampling")  <- colSums(!is.na(Delta.resampling))
              }else{
                  attr(outConfint, "n.resampling")  <- stats::setNames(rep(as.numeric(NA), D), all.endpoint)
              }
              attr(outConfint,"method.ci.resampling") <- method.ci.resampling

              ## ** iid
              if(attr(method.inference,"ustatistic")){
                  if(!is.null(endpoint)){
                      attr(outConfint,"iid") <- do.call(cbind,lapply(which(all.endpoint %in% endpoint), function(iD){ ## iD <- 1
                          if(transformation){
                              return(trans.se.delta(Delta[iD], se = Delta.iid[[iD]]))
                          }else{
                              return(Delta.iid[[iD]])
                          }
                      }))
                      colnames(attr(outConfint,"iid")) <- all.endpoint[all.endpoint %in% endpoint]
                  }else{
                      attr(outConfint,"iid") <- do.call(cbind,lapply(1:D, function(iD){ ## iD <- 1
                          if(transformation){
                              return(trans.se.delta(Delta[iD], se = Delta.iid[[iD]]))
                          }else{
                              return(Delta.iid[[iD]])
                          }
                      }))
                      colnames(attr(outConfint,"iid")) <- all.endpoint
                  }
                  attr(outConfint,"transform") <- trans.delta
                  attr(outConfint,"backtransform") <- itrans.delta
              }
                  
              ## ** export
              if(attr(method.inference,"permutation")){
                  if(is.null(attr(conf.level,"warning.permutation")) || !identical(attr(conf.level,"warning.permutation"),FALSE)){
                      warning("Confidence intervals are computed under the null hypothesis and therefore may not be valid. \n")
                  }
                  attr(outConfint,"warning") <- "Confidence intervals are computed under the null hypothesis"
              }
              return(outConfint)
              
          })

## * confint_percentilePermutation (called by confint)
confint_percentilePermutation <- function(Delta, Delta.resampling,
                                          null, alternative, alpha,
                                          endpoint, backtransform.delta, ...){

    
    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))
    
    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- apply(backtransform.delta(Delta.resampling), MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    ## ** confidence interval
    Delta.resamplingH0 <- apply(Delta.resampling, MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)
    outTable[,"lower.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + apply(Delta.resamplingH0, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE),
                                                        "less" = -Inf,
                                                        "greater" = Delta + apply(Delta.resamplingH0, MARGIN = 2, FUN = stats::quantile, probs = alpha, na.rm = TRUE)
                                                        ))
    
    outTable[,"upper.ci"] <- backtransform.delta(switch(alternative,
                                                        "two.sided" = Delta + apply(Delta.resamplingH0, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2, na.rm = TRUE),
                                                        "less" = Delta + apply(Delta.resamplingH0, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha, na.rm = TRUE),
                                                        "greater" = Inf
                                                        ))

    ## ** p-value
    outTable[,"null"] <- backtransform.delta(null)
    if(BuyseTest.options()$add.1.pperm){
        outTable[,"p.value"] <- sapply(1:n.endpoint, FUN = function(iE){ ## iE <- 1
            switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
                   "two.sided" = (1+sum(abs(Delta[iE] - null) <= abs(Delta.resampling[,iE] - null), na.rm = TRUE))/(1+sum(!is.na(abs(Delta[iE] - null) <= abs(Delta.resampling[,iE] - null)))),
                   "less" = (1+sum(Delta[iE] >= Delta.resampling[,iE], na.rm = TRUE))/(1+sum(!is.na(Delta[iE] >= Delta.resampling[,iE]))),
                   "greater" = (1+sum(Delta[iE] <= Delta.resampling[,iE], na.rm = TRUE))/(1+sum(!is.na(Delta[iE] <= Delta.resampling[,iE])))
                   )
        })
    }else{
        outTable[,"p.value"] <- sapply(1:n.endpoint, FUN = function(iE){ ## iE <- 1
            switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
                   "two.sided" = mean(abs(Delta[iE] - null) <= abs(Delta.resampling[,iE] - null), na.rm = TRUE),
                   "less" = mean(Delta[iE] >= Delta.resampling[,iE], na.rm = TRUE),
                   "greater" = mean(Delta[iE] <= Delta.resampling[,iE], na.rm = TRUE)
                   )
        })
    }

    ## ** export
    return(outTable)
}

## * confint_percentileBootstrap (called by confint)
confint_percentileBootstrap <- function(Delta, Delta.resampling,
                                        null, alternative, alpha,
                                        endpoint, backtransform.delta, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

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
    if(length(null)==1 && n.endpoint>1){
        null <- rep(null, n.endpoint)
    }
    for(iE in which(!is.na(null))){
        outTable[iE, "null"] <- backtransform.delta(null[iE])
        outTable[iE, "p.value"] <- boot2pvalue(na.omit(Delta.resampling[,iE]), null = null[iE], estimate = Delta[iE],
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
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

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
    outTable[,"null"] <- backtransform.delta(null)
    outTable[,"p.value"] <- switch(alternative,
                                   "two.sided" = 2*(1-stats::pnorm(abs((Delta-null)/Delta.se))), 
                                   "less" = stats::pnorm((Delta-null)/Delta.se),
                                   "greater" = 1-stats::pnorm((Delta-null)/Delta.se) 
                                   )

    ## ** export
    return(outTable)
}

## * confint_studentPermutation (called by confint)
confint_studentPermutation <- function(Delta, Delta.se, Delta.resampling, Delta.se.resampling,
                                       null, alternative, alpha,
                                       endpoint, backtransform.delta, backtransform.se, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** critical quantile
    Delta.statH0.resampling <- apply(Delta.resampling, MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)/Delta.se.resampling

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
    outTable[,"null"] <- backtransform.delta(null)
    Delta.stat <- (Delta-null)/Delta.se
    Delta.stat.resampling <- (Delta.resampling-null)/Delta.se.resampling

    outTable[,"p.value"] <- sapply(1:n.endpoint, FUN = function(iE){ ## iE <- 1
        switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the point estimate
               "two.sided" = mean(abs(Delta.stat[iE]) <= abs(Delta.stat.resampling[,iE]), na.rm = TRUE),
               "less" = mean(Delta.stat[iE] >= Delta.stat.resampling[,iE], na.rm = TRUE),
               "greater" = mean(Delta.stat[iE] <= Delta.stat.resampling[,iE], na.rm = TRUE)
               )
    })    

    ## special case
    if(any(Delta.se==0)){
        index0 <- which(Delta.se==0)
        outTable[index0,"lower.ci"] <- outTable[index0,"estimate"]
        outTable[index0,"upper.ci"] <- outTable[index0,"estimate"]
        outTable[index0,"p.value"] <- as.numeric(NA)
    }

    
    ## ** export
    return(outTable)



}
## * confint_studentBootstrap (called by confint)
confint_studentBootstrap <- function(Delta, Delta.se, Delta.resampling, Delta.se.resampling,
                                     null, alternative, alpha,
                                     endpoint, backtransform.delta, backtransform.se, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** critical quantile
    Delta.statH0.resampling <- apply(Delta.resampling, MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)/Delta.se.resampling  ## center around the null

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

    outTable[, "null"] <- backtransform.delta(null)
    for(iE in 1:n.endpoint){ ## iE <- 1
        ## if(sign(mean(Delta.resampling[,iE]))!=sign(Delta[iE])){
        ## warning("the estimate and the average bootstrap estimate do not have same sign \n")
        ## }
        outTable[iE, "p.value"] <- boot2pvalue(Delta.statH0.resampling[,iE], null = null, estimate = Delta[iE], ## note: estimate is not used to produce the ci, just for knowing the sign
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
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

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
    outTable[,"null"] <- backtransform.delta(null)
    outTable[,"p.value"] <- switch(alternative,
                                   "two.sided" = 2*(1-stats::pnorm(abs((Delta-null)/Delta.se))), 
                                   "less" = stats::pnorm((Delta-null)/Delta.se),
                                   "greater" = 1-stats::pnorm((Delta-null)/Delta.se) 
                                   )

    ## special case with no variability
    if(any(na.omit((Delta==null)*(Delta.se==0)) == 1)){
        outTable[(Delta==null)*(Delta.se==0) == 1,"p.value"] <- 1
    }

    ## ** export
    return(outTable)
}

## * confint_none (called by confint)
confint_none <- function(Delta, endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(NA, nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))

    ## ** point estimate
    outTable[,"estimate"] <- Delta

    ## ** return
    return(outTable)

    
}

##----------------------------------------------------------------------
### S4BuyseTest-confint.R ends here
