### BuyseTest-confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 19 2018 (23:37) 
## Version: 
##           By: Brice Ozenne
##     Update #: 1278
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
#' @title Extract Confidence Interval from GPC
#' @aliases confint,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Extract confidence intervals for summary statistics (net benefit, win ratio, ...) estimated by GPC.
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param statistic [character] the statistic summarizing the pairwise comparison: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param endpoint [character] for which endpoint(s) the confidence intervals should be output?
#' If \code{NULL} returns the confidence intervals for all endpoints.
#' @param strata [character] the strata relative to which the statistic should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the statistic pooled over all strata,
#' or \code{TRUE} to output each strata-specific statistic.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param null [numeric] right hand side of the null hypothesis (used for the computation of the p-value).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the inverse hyperbolic tangent scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param order.Hprojection [integer, 1-2] order of the H-decomposition used to compute the variance.
#' @param method.ci.resampling [character] the method used to compute the confidence intervals and p-values when using bootstrap or permutation (\code{"percentile"}, \code{"gaussian"}, \code{"student"}).
#' See the details section.
#' @param cluster [numeric vector] Group of observations for which the iid assumption holds .
#' @param sep [character] character string used to separate the endpoint and the strata when naming the statistics.
#'  
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#' 
#' @details 
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
#' @keywords method
#' @author Brice Ozenne

## * Method - confint
#' @rdname S4BuyseTest-confint
#' @exportMethod confint
setMethod(f = "confint",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                strata = FALSE,
                                cumulative = TRUE,
                                null = NULL,
                                conf.level = NULL,
                                alternative = NULL,
                                method.ci.resampling = NULL,
                                order.Hprojection = NULL,
                                transformation = NULL,
                                cluster = NULL,
                                sep="."){

              option <- BuyseTest.options()
              sep <- "."

              D <- length(object@endpoint)
              method.inference <- object@method.inference
              add.halfNeutral <- object@add.halfNeutral
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              if(is.null(conf.level)){
                  if(!attr(method.inference,"permutation")){
                      conf.level <- option$conf.level
                  }else{
                      conf.level <- NA
                  }
              }
              if(is.null(alternative)){
                  alternative <- option$alternative
              }
              
              ## ** normalize and check arguments
              ## statistic
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

              ## strata
              level.strata <- object@level.strata
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
                                 method = "confint[S4BuyseTest]")
              }
              if(attr(object@weightStrata,"type")=="standardization" && attr(method.inference,"ustatistic") && any(strata %in% level.strata)){
                  ## note: error message for studentized method is triggered later (end of method.ci checks)
                  stop("Cannot output the strata-specific uncertainty based on the H-decomposition when using standardization. \n",
                       "Consider using a resampling method by specifying the argument \'method.inference\' when calling the \"BuyseTest\" function. \n")
              }
              if(attr(slot(object,"scoring.rule"), "test.match") && any(level.strata %in% strata)){
                  stop("Cannot output p-values or confidence intervals for stratified statistics with matched data. \n")
              }

              ## method.ci
              if(attr(method.inference,"permutation") || attr(method.inference,"bootstrap")){
                  if(is.null(method.ci.resampling)){
                      if(method.inference == "varexact permutation"){
                          method.ci.resampling <- "gaussian"
                      }else if(attr(method.inference,"studentized")){
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

                  if(method.ci.resampling != "gaussian" && method.inference == "varexact permutation"){
                      stop("Argument \'method.ci.resampling\' must be set to \'gaussian\' if argument \'method.inference\' has been set to \"varexact permutation\" when calling BuyseTest. \n")
                  }
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
                  if(is.null(transformation)){
                      if(method.ci.resampling=="percentile" || method.inference == "varexact permutation"){
                          transformation <- FALSE ## ensures consistency between p-values for different statistics as transformation may lead to numerical unaccuracies when comparing resampling to observed
                      }else{
                          transformation <- option$transformation
                      }
                  }else if(transformation && method.inference == "varexact permutation"){
                      transformation <- FALSE
                      message("Argument \'transformation\' has been set to FALSE. \n",
                              "Transformation is not available if argument \'method.inference\' has been set to \"varexact permutation\" when calling BuyseTest. \n")
                  }              
              }else{
                  if(is.null(transformation)){
                      transformation <- option$transformation
                  }              
                  if(!is.null(method.ci.resampling)){
                      warning("Argument \'method.ci.resampling\' is disregarded when not using resampling\n")                  
                  }
              }
              if(attr(method.inference,"studentized") && (any(strata != "global") || (cumulative!=TRUE)) ){
                  stop("Can only perform statistical inference based on studentized resampling for global cumulative effects. \n",
                       "Consider setting argument \'strata\' to FALSE and argument \'cumulative\' to TRUE. \n")
              }

              ## order.Hprojection
              if(attr(method.inference,"ustatistic")){
                  if(!is.null(order.Hprojection) && order.Hprojection != attr(method.inference,"hprojection")){

                      validInteger(order.Hprojection,
                                   name1 = "order.Hprojection",
                                   min = 1, max = 2, valid.length = 1,
                                   method = "confint[S4BuyseTest]")
                  
                      if(order.Hprojection > attr(method.inference,"hprojection")){
                          stop("Cannot find the second order of the H-decomposition. \n",
                               "Consider setting order.Hprojection to 2 in BuyseTest.options before calling BuyseTest. \n")
                      }

                      object.hprojection <- FALSE ## move from order H-decomposition of order 2 to H-decomposition of order 1

                  }else{
                      object.hprojection <- TRUE
                  }                  
              }else{
                  object.hprojection <- TRUE
              }


              ## conf.level
              validNumeric(conf.level,
                           name1 = "conf.level",
                           min = 0, max = 1,
                           refuse.NA = FALSE,
                           valid.length = 1,
                           method = "confint[S4BuyseTest]")
              alpha <- 1-conf.level

              ## alternative
              validCharacter(alternative,
                             name1 = "alternative",
                             valid.values = c("two.sided","less","greater"),
                             valid.length = 1,
                             method = "confint[S4BuyseTest]")

              ## transformation
              validLogical(transformation,
                           name1 = "transformation",
                           valid.length = 1,
                           method = "confint[S4BuyseTest]")

              ## endpoint
              valid.endpoint <- names(object@endpoint)
              if(!is.null(endpoint)){
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
              }else{
                  endpoint <- valid.endpoint
              }
              n.endpoint <- length(endpoint)
              
              ## safety
              test.model.tte <- all(unlist(lapply(object@iidNuisance,dim))==0)
              if(method.inference %in% c("u statistic","u statistic bebu") && object@correction.uninf > 0){
                  warning("The current implementation of the asymptotic distribution is not valid when using a correction. \n",
                          "Standard errors / confidence intervals / p-values may not be correct. \n",
                          "Consider using a resampling approach or checking the control of the type 1 error with powerBuyseTest. \n")
              }

              ## weight
              if(!is.null(cluster) && any(object@weightObs!=1)){
                  stop("Cannot handle clustered observations when observations are weighted. \n")
              }

              ## ** extract estimate
              all.endpoint <- names(object@endpoint)
              DeltaW <- coef(object, endpoint = endpoint, statistic = statistic, strata = strata, cumulative = cumulative, resampling = FALSE, simplify = FALSE)
              if(length(strata)==1 && all(strata=="global")){
                  Delta <- stats::setNames(DeltaW["global",], endpoint)
              }else{
                  DeltaL <- stats::reshape(data.frame(strata = strata, DeltaW), direction = "long", varying = endpoint,
                                           times = endpoint, v.names = "statistic")
                  Delta <- stats::setNames(DeltaL$statistic, paste(DeltaL$time, DeltaL$strata, sep = sep))
              }

              if(((attr(method.inference,"permutation") && method.inference != "varexact permutation")) || attr(method.inference,"bootstrap")){
                  DeltaW.resampling <- coef(object, endpoint = endpoint, statistic = statistic, strata = strata, cumulative = cumulative, resampling = TRUE, simplify = FALSE)
                  if(length(strata)==1 && all(strata=="global")){
                      Delta.resampling <- matrix(DeltaW.resampling[,"global",], ncol = length(endpoint), dimnames = list(NULL, endpoint))
                  }else{
                      Delta.resampling <- do.call(cbind,apply(DeltaW.resampling, MARGIN = 3 , FUN = base::identity, simplify = FALSE))
                      colnames(Delta.resampling) <- unlist(lapply(endpoint, paste, strata, sep = sep))
                      Delta.resampling <- Delta.resampling[,names(Delta),drop=FALSE]
                  }
              }else{
                  Delta.resampling <- NULL
              }

              ## ** extract standard error
              if(attr(method.inference,"ustatistic") || attr(method.inference,"studentized")){

                  if(object.hprojection && is.null(cluster) && (length(strata)==1 && all(strata=="global")) && cumulative == TRUE){
                      Delta.se <- sqrt(object@covariance[endpoint,statistic])
                      if(attr(method.inference,"studentized")){
                          Delta.se.resampling <- matrix(sqrt(object@covarianceResampling[,endpoint,statistic]),
                                                        ncol = n.endpoint,
                                                        dimnames = list(NULL, endpoint))
                      }else{
                          Delta.se.resampling <- NULL
                      }
                  }else{    
                      if(identical(order.Hprojection,2)){
                          warning("Inference will be performed using a first order H projection. \n")
                      }
                      ls.Delta.iid <- getIid(object, statistic = statistic, cumulative = cumulative, endpoint = endpoint, strata = strata, cluster = cluster, simplify = FALSE)
                      if(length(strata) == 1 && all(strata=="global")){
                          Delta.iid <- ls.Delta.iid$global[,names(Delta),drop=FALSE]
                      }else{
                          for(iS in 1:length(strata)){ ## iS <- 1
                              colnames(ls.Delta.iid[[iS]]) <- paste(colnames(ls.Delta.iid[[iS]]), strata[iS], sep = sep)
                          }
                          if(length(strata)==1){
                              Delta.iid <- ls.Delta.iid[[1]][,names(Delta),drop=FALSE]
                          }else{
                              Delta.iid <- do.call(cbind,ls.Delta.iid)[,names(Delta),drop=FALSE]
                          }
                      }
                      if(is.null(cluster) && any(object@weightObs!=1)){
                          Delta.iid <- .colMultiply_cpp(Delta.iid, sqrt(object@weightObs))
                      }
                      M.se <- sqrt(colSums(Delta.iid^2))
                      Delta.se <- stats::setNames(as.double(M.se), names(Delta))
                      Delta.se.resampling <- NULL
                  }
              }else{
                  if(!is.null(cluster)){
                      message("BuyseTest: argument \'cluster\' ignored when evaluating uncertainty with resampling methods. \n")
                  }
                      
                  if(method.inference == "varexact permutation"){
                      if(statistic == "winRatio"){
                          stop("BuyseTest: cannot evaluate the exact variance of the permutation distribution for the win ratio. \n",
                               "Consider using the net benefit instead (argument statistic = \"netBenefit\"). \n")
                      }
                      if(cumulative == FALSE){
                          stop("BuyseTest: cannot evaluate the exact variance of the permutation distribution for each endpoint separately. \n",
                               "Consider setting the argument \'cumulative\' to TRUE. \n")
                      }
                      if((length(strata)==1 && all(strata=="global"))){
                          Delta.se <- sqrt(object@covariance[endpoint,statistic])
                      }else{
                          Delta.se <- unlist(stats::setNames(lapply(endpoint, function(iE){ ## iE <- endpoint[1]
                              c("global" = sqrt(object@covariance[iE,statistic]), sqrt(attr(object@covariance,"strata")[[iE]][,statistic]))[strata]
                          }),endpoint))
                      }
                  }else{
                      Delta.se <- NULL
                  }
                  Delta.se.resampling <- NULL
              }

              ## ** null hypothesis
              if(attr(method.inference,"permutation") && !add.halfNeutral && is.null(null) && statistic %in% c("favorable","unfavorable")){
                  null <- NA
              }else if(is.null(null)){
                  null <- switch(statistic,
                                 "netBenefit" = 0,
                                 "winRatio" = 1,
                                 "favorable" = ifelse(add.halfNeutral,1/2,NA),
                                 "unfavorable" = ifelse(add.halfNeutral,1/2,NA))
              }else {
                  validNumeric(null, valid.length = 1,
                               refuse.NA = !attr(method.inference,"permutation"),
                               min = if(statistic=="netBenefit"){-1}else{0},
                               max = if(statistic=="winRatio"){Inf}else{1})
              }
              null <- rep(null, length(Delta))

              ## ** method
              if(method.inference == "none"){
                  method.confint <- confint_none
                  transformation <- FALSE
              }else if(method.inference == "varexact permutation"){
                  method.confint <- confint_varexactPermutation
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
                      trans.weight <- sum(object@weightEndpoint)
                  }
                  trans.name <- switch(statistic,
                                       "netBenefit" = "atanh",
                                       "winRatio" = "log",
                                       "favorable" = "atanh",
                                       "unfavorable" = "atanh"
                                       )
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
                  trans.name <- "id"
                  trans.delta <- function(x){x}
                  itrans.delta <- function(x){x}
                  trans.se.delta <- function(x,se){se}
                  itrans.se.delta <- function(x,se){se}                  
              }

              ## ** compute the confidence intervals
              if(statistic=="winRatio" && transformation==FALSE){
                  attr(null, "type")  <- "relative"
              }else{
                  attr(null, "type")  <- "absolute"
              }
              outConfint <- do.call(method.confint,
                                    args = list(Delta = trans.delta(Delta),
                                                Delta.resampling = trans.delta(Delta.resampling),
                                                Delta.se = trans.se.delta(Delta, se = Delta.se),
                                                Delta.se.resampling = trans.se.delta(Delta.resampling, se = Delta.se.resampling),
                                                alternative = alternative,
                                                null = trans.delta(null),
                                                alpha = alpha,
                                                endpoint = names(Delta),
                                                backtransform.delta = itrans.delta,
                                                backtransform.se = itrans.se.delta))

              ## do not output CI or p-value when the estimate has not been identified
              index.NA <- union(which(is.infinite(outConfint[,"estimate"])),which(is.na(outConfint[,"estimate"])))
              if(length(index.NA)>0){
                  outConfint[index.NA,c("se","lower.ci","upper.ci","p.value")] <- NA
              }
              outConfint <- as.data.frame(outConfint)

              ## ** number of permutations
              if(method.inference != "none" && ((attr(method.inference,"permutation") && (method.inference!="varexact permutation")) || attr(method.inference,"bootstrap"))){
                  attr(outConfint, "n.resampling")  <- colSums(!is.na(Delta.resampling))
              }else{
                  attr(outConfint, "n.resampling")  <- stats::setNames(rep(as.numeric(NA), D), all.endpoint)
              }
              attr(outConfint,"method.ci.resampling") <- method.ci.resampling

              ## ** transform
              if(attr(method.inference,"ustatistic")){                 
                  attr(outConfint,"nametransform") <- trans.name
                  attr(outConfint,"transform") <- trans.delta
                  attr(outConfint,"backtransform") <- itrans.delta
              }
                  
              ## ** export
              if(attr(method.inference,"permutation") && !is.na(conf.level)){
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
    if(!is.na(alpha)){
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
    }

    ## ** null
    if(any(is.na(null))){
        null[is.na(null)] <- apply(Delta.resampling,2,stats::median)[is.na(null)]
    }
    outTable[,"null"] <- backtransform.delta(null)
    
    ## ** p-value
    add.1 <- BuyseTest.options()$add.1.presample
    outTable[,"p.value"] <- sapply(1:n.endpoint, FUN = function(iE){ ## iE <- 2
    ## rounding is here to mitigate p-value mismatch between netBenefit and winRatio due to finite numeric precision

        if(alternative == "two.sided"){
            if(attr(null,"type")=="relative"){ ## win ratio without transformation
                ## H0 WR=1 so if hat(WR)=3/2 more extreme is above 3/2 or below 2/3
                test.alternative <- round(pmax(Delta.resampling[,iE]/null[iE],null[iE]/Delta.resampling[,iE]),10)/round(max(Delta[iE]/null[iE],null[iE]/Delta[iE]),10) >= 1
                ## test.alternative <- abs(log(Delta[iE]/null)) <= abs(log(Delta.resampling[,iE]/null)) ## try to avoid log-transformation
            }else if(attr(null,"type")=="absolute"){
                test.alternative <- round(abs(Delta[iE]-null[iE]),10) <= round(abs(Delta.resampling[,iE]-null[iE]),10)
            }

        }else{
            test.alternative <- switch(alternative, 
                                       "less" = round(Delta[iE],10) >= round(Delta.resampling[,iE],10),
                                       "greater" = round(Delta[iE],10) <= round(Delta.resampling[,iE],10)
                                       )
        }

        p.alternative <- (add.1 + sum(test.alternative, na.rm = TRUE)) / (add.1 + sum(!is.na(test.alternative), na.rm = TRUE))
        return(p.alternative)        
    })

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
    if(!is.na(alpha)){
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
    }

    ## ** p.values
    outTable[,"null"] <- backtransform.delta(null)

    add.1 <- BuyseTest.options()$add.1.presample
    for(iE in which(!is.na(null))){
        outTable[iE, "p.value"] <- boot2pvalue(stats::na.omit(Delta.resampling[,iE]),
                                               null = null[iE],
                                               estimate = Delta[iE],
                                               alternative = alternative,
                                               add.1 = add.1)
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
    outTable[,"se"] <- apply(backtransform.delta(Delta.resampling), MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    if(any(is.infinite(Delta.resampling))){
        txt.range <- NULL
        pc.infinite <- 100*mean(is.infinite(Delta.resampling))
        if(any(Delta.resampling[is.infinite(Delta.resampling)]>0)){
            Delta.resampling.max <- max(Delta.resampling[!is.infinite(Delta.resampling)])
            if(Delta.resampling.max<0){Delta.resampling.max <- 1}
            txt.range <- paste(txt.range, signif(Delta.resampling.max), collapse = " and ")
            Delta.resampling[is.infinite(Delta.resampling) & Delta.resampling > 0] <- 1.1*Delta.resampling.max
        }
        if(any(Delta.resampling[is.infinite(Delta.resampling)]<0)){
            Delta.resampling.min <- min(Delta.resampling[!is.infinite(Delta.resampling)])
            if(Delta.resampling.min>0){Delta.resampling.min <- -1}
            txt.range <- paste(txt.range, signif(Delta.resampling.min), collapse = " and ")
            Delta.resampling[is.infinite(Delta.resampling) & Delta.resampling < 0] <- 1.1*Delta.resampling.min
        }

        warning("Infinite statistic value after transformation in ",round(pc.infinite,2),"% of the bootstrap samples. \n",
                "Will be set to",txt.range," when evaluating CIs or p-value under Gaussian approximation. \n",
                "(110% of the most extreme, non-infinite, value or +/- 1 if not finite value of the same sign)", sep = "")
    }
    Delta.se <- apply(Delta.resampling, MARGIN = 2, FUN = stats::sd, na.rm = TRUE) ## computed based on the sample

    ## ** confidence interval
    if(!is.na(alpha)){
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
    }

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

    ## identify special case (no variability in the estimate)
    test.variability <- colSums(Delta.se.resampling!=0)+(apply(Delta.resampling,2,function(iDelta){length(unique(iDelta))})>1)+(Delta.se!=0)
    index.novar <- which(test.variability==0)
    index.var <- which(test.variability!=0)

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** critical quantile
    if(!is.na(alpha) && length(index.var)>0){
        
        Delta.statH0.resampling <- .rowCenter_cpp(Delta.resampling[,index.var,drop=FALSE],null[index.var])/Delta.se.resampling[,index.var,drop=FALSE]

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
    }

    ## ** confidence interval
    if(!is.na(alpha) && length(index.var)>0){
        outTable[index.var,"lower.ci"] <- backtransform.delta(Delta[index.var] + Delta.qInf * Delta.se[index.var])
        outTable[index.var,"upper.ci"] <- backtransform.delta(Delta[index.var] + Delta.qSup * Delta.se[index.var])
    }

    if(length(index.novar)>0){
        outTable[index.novar,"lower.ci"] <- backtransform.delta(Delta[index.novar])
        outTable[index.novar,"upper.ci"] <- backtransform.delta(Delta[index.novar])
    }
    
    ## ** null
    if(any(is.na(null))){
        null[is.na(null)] <- apply(Delta.resampling,2,stats::median)[is.na(null)]
    }
    outTable[,"null"] <- backtransform.delta(null)

    ## ** p.value
    if(length(index.var)>0){
        add.1 <- BuyseTest.options()$add.1.presample
        Delta.stat <- (Delta-null)/Delta.se
        Delta.stat.resampling <- (Delta.resampling-null)/Delta.se.resampling
        if(any(is.infinite(Delta.resampling))){
            Delta.stat.resampling[is.infinite(Delta.resampling)] <- Delta.resampling[is.infinite(Delta.resampling)]
        }
        outTable[index.var,"p.value"] <- sapply(index.var, FUN = function(iE){ ## iE <- 1

            ## rounding is here to mitigate p-value mismatch between netBenefit and winRatio due to finite numeric precision
            test.alternative <- switch(alternative,
                                       "two.sided" = round(abs(Delta.stat[iE]),10) <= round(abs(Delta.stat.resampling[,iE]),10),
                                       "less" = round(Delta.stat[iE],10) >= round(Delta.stat.resampling[,iE],10),
                                       "greater" = round(Delta.stat[iE],10) <= round(Delta.stat.resampling[,iE],10)
                                       )

            p.alternative <- (add.1 + sum(test.alternative, na.rm = TRUE)) / (add.1 + sum(!is.na(test.alternative)))
            return(p.alternative)        
        })
    }
    
    if(length(index.novar)>0){
        outTable[index.novar,"p.value"] <- as.numeric(abs(outTable[index.novar,"estimate"]-outTable[index.novar,"null"])<1e-10)
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

    ## identify special case (no variability in the estimate)
    test.variability <- colSums(Delta.se.resampling!=0)+(apply(Delta.resampling,2,function(iDelta){length(unique(iDelta))})>1)+(Delta.se!=0)
    index.novar <- which(test.variability==0)
    index.var <- which(test.variability!=0)

    ## ** point estimate
    outTable[,"estimate"] <- backtransform.delta(Delta)

    ## ** standard error
    outTable[,"se"] <- backtransform.se(Delta, se = Delta.se)

    ## ** critical quantile
    ## z-transformation: center around estimate and divide by estimated se
    if(length(index.var)>0){
        Delta.statH0.resampling <- sweep(Delta.resampling[,index.var,drop=FALSE], MARGIN = 2, FUN = "-", STATS = Delta)/Delta.se.resampling[,index.var,drop=FALSE]  
        ## Delta.statH0.resampling <- apply(Delta.resampling[,index.var,drop=FALSE], MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)/Delta.se.resampling[,index.var,drop=FALSE]  

        if(!is.na(alpha)){
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
        }

    }

    ## ** confidence interval
    ## normal case
    if(!is.na(alpha) && length(index.var)>0){
        outTable[index.var,"lower.ci"] <- backtransform.delta(Delta[index.var] + Delta.qInf * Delta.se[index.var])
        outTable[index.var,"upper.ci"] <- backtransform.delta(Delta[index.var] + Delta.qSup * Delta.se[index.var])
    }

    ## special case
    if(!is.na(alpha) && length(index.novar)>0){
        test.diff <- colSums(Delta.resampling[,index.novar,drop=FALSE] != matrix(Delta, nrow = NROW(Delta.resampling), ncol = length(Delta), byrow = TRUE))
        outTable[index.novar[test.diff],c("lower.ci","upper.ci")] <- NA
        outTable[index.novar[test.diff==0],c("lower.ci","upper.ci")] <- backtransform.delta(Delta[index.novar[test.diff==0]])
    }


    ## ** p.value
    outTable[, "null"] <- backtransform.delta(null)

    add.1 <- BuyseTest.options()$add.1.presample
    for(iE in index.var){ ## iE <- 1
        outTable[iE, "p.value"] <- boot2pvalue(stats::na.omit(Delta[iE] + Delta.se[iE] * Delta.statH0.resampling[,iE]),
                                               null = null[iE],
                                               estimate = Delta[iE], ## note: estimate is not used to produce the ci, just for knowing the sign
                                               alternative = alternative,
                                               add.1 = add.1)
    }
    
    ## special case
    if(length(index.novar)>0){
        outTable[index.novar[test.diff],c("p.value")] <- NA
        outTable[index.novar[test.diff==0],c("p.value")] <- (null==Delta) + add.1*(null!=Delta)/(NROW(Delta.resampling)+1)
    }
    
    ## ** export
    return(outTable)



}


## * confint_varexactPermutation (called by confint)
confint_varexactPermutation <- function(Delta, Delta.se, null,
                                        alternative, alpha,
                                        endpoint, ...){

    n.endpoint <- length(endpoint)
    outTable <- matrix(as.numeric(NA), nrow = n.endpoint, ncol = 6,
                       dimnames = list(endpoint, c("estimate","se","lower.ci","upper.ci","null","p.value")))
    ## Note: no transformation

    ## ** point estimate
    outTable[,"estimate"] <- Delta

    ## ** standard error
    outTable[,"se"] <- Delta.se

    ## ** confidence interval
    ## No CI because se estimated under H0 instead of H1
    
    ## ** p-value
    outTable[,"null"] <- null
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

## * confint_Ustatistic (called by confint)
confint_Ustatistic <- function(Delta, Delta.se, null,
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
    if(!is.na(alpha)){
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
    }

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
