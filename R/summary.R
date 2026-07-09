### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:27) 
## Version: 
## Last-Updated: jul  6 2026 (12:33) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.BuyseTTEM
#' @title Fitted Survival Probabilities for all Strata
#' @description Display fitted survival probabilities along with their confidence intervals for all strata.
#' 
#' @param object \code{BuyseTTEM} object that is also a \code{prodlim} object.
#' @param times [numeric vector] timepoints at which to evaluate the probabilities.
#' @param iid [logical] Should the se and iid decomposition of the fitted probabilities be output. 
#' @param ... For compatibility with the generic method.
#'
#' @return A data.frame
#' @keywords print
#' @method summary BuyseTTEM
#' @export
summary.BuyseTTEM <- function(object, times = NULL, iid = FALSE, ...){


    if(!inherits(object,"prodlim") & (iid || !is.null(times))){
        stop("Argument \'times\' and \'iid\' for the summary method of BuyseTTEM objects only developped when using prodlim. \n")
    }
    if(iid && length(object$peron$iid.cif[[1]][[1]][[1]])==0){ ## iStrata iTreat iCause
        stop("Using argument \'iid\' equal to TRUE require that the iid decomposition has already been stored in the object. \n",
             "Consider setting the argument \'iid\' to TRUE when calling BuyseTTEM. \n")
    }

    ## ** remove BuyseTTEM class
    object.original <- object
    class(object.original) <- setdiff(class(object.original), "BuyseTTEM")

    ## ** get standard summary
    if(is.null(times)){
        sobject <- summary(object.original)
    }else{
        sobject <- summary(object.original, times = times)
    }

    if(!inherits(object.original,"prodlim")){
        return(sobject)
    }

    ## ** add influence function
    treatment <- object$peron$treatment
    strata <- attr(object$peron$strata.var,"original")

    sobject <- as.data.frame(sobject)
    if(is.null(strata)){
        sobject$..strata.. <- 1
    }else{
        var2strata <- object$peron$X
        var2strata[[treatment]] <- factor(var2strata[[treatment]], levels = 0:1, object$peron$level.treatment)
        sobject <- merge(sobject, var2strata, by = c(treatment,strata))
    }

    ## find levels
    grid.treatmentstrata <- unique(sobject[c(treatment,strata)])
    levels.treatmentstrata <- interaction(grid.treatmentstrata[c(treatment,strata)], drop=TRUE)
    n.treatmentstrata <- length(levels.treatmentstrata)
    sobject$levels <- factor(interaction(sobject[c(treatment,strata)]), levels = levels(levels.treatmentstrata))

    ## find times
    ls.time <- by(sobject$time, sobject$levels, FUN = identity, simplify = FALSE)
    
    ## extract prediction and merge
    attr(sobject, "iid") <- NULL
    for(iL in 1:n.treatmentstrata){
        if(is.null(strata)){
            iPred <- predict(object, times = ls.time[[iL]], treatment = grid.treatmentstrata[iL,treatment], iid = iid)
        }else{
            iPred <- predict(object, times = ls.time[[iL]], treatment = grid.treatmentstrata[iL,treatment], strata = grid.treatmentstrata[iL,strata,drop=FALSE], iid = iid)
        }
        sobject[sobject$levels==levels.treatmentstrata[iL],"surv"] <- iPred$survival ## fix automatic NA attribution when beyond the last observation even though survival is at 0
        if(iid){
            sobject[sobject$levels==levels.treatmentstrata[iL],"se.surv"] <- iPred$survival.se
            sobject[sobject$levels==levels.treatmentstrata[iL],"lower"] <- iPred$survival + stats::qnorm(0.025)*iPred$survival.se
            sobject[sobject$levels==levels.treatmentstrata[iL],"upper"] <- iPred$survival + stats::qnorm(0.975)*iPred$survival.se
            attr(sobject, "iid") <- cbind(attr(sobject, "iid"), iPred$survival.iid)
        }
    }

    ## ** export
    if(!iid){
        sobject$se.surv <- NULL
        sobject$lower <- NULL
        sobject$upper <- NULL
    }
    sobject$levels <- NULL
    sobject$..strata.. <- NULL
    return(sobject)
}

## * summary.performance
##' @title Summary Method for Performance Objects
##' @description Summary of the performance of binary classifiers
##'
##' @param object output of performance.
##' @param digits [numeric vector of length 2] number of digits used for the estimates and p-values.
##' @param print [logical] should the performance be printed in the console.
##' @param order.model [character vector] ordering of the models.
##' @param ... not used.
##'
##' @keywords print
##' @method summary performance
##' @export
summary.performance <- function(object, order.model = NULL, digits = c(3,3), print = TRUE, ...){

    ## ** re-order models
    if(!is.null(order.model)){
        if(any(duplicated(order.model))){
            stop("Argument \'order.model\' should not contain duplicated values. \n")
        }
        Umodel <- unique(object$performance$model)
        if(is.numeric(order.model)){
            if(!identical(sort(order.model), 1:length(order.model))){
                stop("Argument \'order.model\' should contain integers from 1 to ",length(order.model)," when numeric. \n")
            }
            order.model <- Umodel[order.model]            
        }else{
            if(any(order.model %in% Umodel == FALSE)){
                stop("Unknown value \"",paste(order.model[order.model %in% Umodel == FALSE], collapse = "\" \""),"\" in argument \'order.model\'. \n")
            }
            if(any(unique(object$performance$model) %in% order.model == FALSE)){
                stop("Missing model \"",paste(Umodel[order.model %in% Umodel == FALSE], collapse = "\" \""),"\" in argument \'order.model\'. \n")
            }
        }
        if(!is.null(object$prediction$internal)){
            object$prediction$internal <- object$prediction$internal[,order.model,drop=FALSE]
        }
        if(!is.null(object$prediction$external)){
            object$prediction$external <- object$prediction$external[,order.model,drop=FALSE]
        }
        if(!is.null(object$prediction$cv)){
            object$prediction$cv <- object$prediction$cv[,order.model,,drop=FALSE]
        }
        if(!is.null(object$resampling)){
            object$resampling[,c("model") := factor(.SD$model,levels = order.model)]
            data.table::setkeyv(object$resampling, c("sample","method","metric","model"))
            object$performance <- .performanceResample_inference(performance = object$performance[order(factor(object$performance$model, levels = order.model)),
                                                                                                  c("method","metric","model","estimate")],
                                                                 resampling = object$resampling,
                                                                 type.resampling = object$args$type.resampling,
                                                                 conf.level = object$args$conf.level)
        }else{
            object$performance <- object$performance[order(factor(object$performance$method, levels = c("internal","external","cv")),
                                                           factor(object$performance$metric, levels = c("auc","brier")),
                                                           factor(object$performance$model,levels=order.model)),,drop=FALSE]
            if(length(order.model)>1){
                object$performance$p.value_comp <- NA

                for(iMethod in unique(object$performance$method)){ ## iMethod <- "internal"
                    for(iMetric in c("auc","brier")){ ## iMetric <- "auc"
                        iIndex <- which(object$performance$method==iMethod & object$performance$metric==iMetric)
                        iBeta <- object$performance[iIndex,"estimate"]
                        iIID <- object[[paste0("iid.",iMetric)]][[iMethod]][,order.model]
                        iStat <- c(NA,diff(iBeta)) / c(NA,sqrt(colSums((iIID[,1:(NCOL(iIID)-1),drop=FALSE]-iIID[,-1,drop=FALSE])^2)))
                        object$performance[iIndex,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                    }
                }
            }
        }
        object$auc <-  lapply(object$auc, function(iL){iL[order.model]})
        object$brier <-  lapply(object$brier, function(iL){iL[order.model]})
    }


    ## ** display results
    df.print <- object$performance

    df.print$p.value <- base::format.pval(df.print$p.value, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value[is.na(object$performance$p.value)] <- ""
    df.print$p.value_comp <- base::format.pval(df.print$p.value_comp, digits = digits[1], eps = 10^(-digits[2]))
    df.print$p.value_comp[is.na(object$performance$p.value_comp)] <- ""
    df.print <- df.print[,union(names(which(colSums(!is.na(object$performance))>0)),"estimate")]
    print(df.print, digits = digits[1])

    return(invisible(object$performance))
}


## * summary S4BuysePower
#' @docType methods
#' @name S4BuysePower-summary
#' @title Summary Method for Class "S4BuysePower"
#' @aliases summary,S4BuysePower-method
#' @include S4-BuysePower.R
#' 
#' @description Summarize the results from the \code{\link{powerBuyseTest}} function.
#' 
#' @param object output of \code{\link{powerBuyseTest}}
#' @param statistic [character] the statistic summarizing the pairwise comparison: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param endpoint [character vector] the endpoints to be displayed: must be the name of the endpoint followed by an underscore and then by the threshold.
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#' @param print [logical] Should the table be displayed?.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param legend [logical] should explainations about the content of each column be displayed? 
#' @param col.rep [logical] should the number of successful simulations be displayed? 
#' @param ... Not used. For compatibility with the generic method.
#'
#' @seealso 
#' \code{\link{powerBuyseTest}} for performing a simulation study for generalized pairwise comparison. \cr
#'
#' @return data.frame
#' @keywords print
#' @author Brice Ozenne
#' @exportMethod summary
setMethod(f = "summary",
          signature = "S4BuysePower",
          definition = function(object, statistic = NULL, endpoint = NULL, order.Hprojection = NULL, transformation = NULL,
                                print = TRUE, legend = TRUE, col.rep = FALSE, digit = 4, ...){

              ## ** normalize and check arguments
              option <- BuyseTest.options()
              dots <- list(...)
              if(length(dots)>0){
                  stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
              }
                  
              validLogical(print,
                           name1 = "print",
                           valid.length = 1,
                           method = "summary[S4BuysePower]")

              args <- slot(object, name = "args")
              null <- args$null
              power <- args$power
              alpha <- 1-args$conf.level
              n.rep <- args$n.rep
              restriction <- args$restriction
              method.inference <- args$method.inference
              max.sample.size <- args$max.sample.size

              dtS.res <- model.tables(object,
                                      statistic = statistic,
                                      endpoint = endpoint,
                                      order.Hprojection = order.Hprojection,
                                      transformation = transformation)
              col.value <- intersect(names(dtS.res),c("mean.estimate","sd.estimate","mean.se","rejection.rate","rep.estimate","rep.se"))
              statistic <- unique(dtS.res$statistic)
              order.Hprojection <- attr(dtS.res,"order.Hprojection")
              transformation <- attr(dtS.res,"transformation")

              nobs <- nobs(object)

              ## ** print
              ls.df.print <- stats::setNames(lapply(statistic, function(iStat){ ## iStat <- dtS.res$statistic[1]
                  iDF <- as.data.frame(dtS.res[dtS.res$statistic == iStat])
                  iDF$statistic <- NULL
                  iDF[,col.value] <- round(iDF[,col.value], digits = digit)
                  if(col.rep == FALSE){
                      iDF$rep.estimate <- NULL
                      iDF$rep.se <- NULL
                  }
                  iDF[duplicated(iDF[,c("endpoint","restriction","threshold")]),c("endpoint","restriction","threshold")] <- as.character(NA)
                  if(all(is.na(iDF$restriction))){
                      iDF$restriction <- NULL
                  }
                  if(all(is.na(iDF$threshold))){
                      iDF$threshold <- NULL
                  }
                  iDF[] <- lapply(iDF, as.character)
                  iDF[is.na(iDF)] <- ""
                  return(iDF)                      
              }), statistic)

              if(print){

                  if(!is.null(power)){
                      range.sampleC <- c(ceiling(min(attr(nobs,"sample")[,"C"])),
                                         ceiling(max(attr(nobs,"sample")[,"C"])))
                      range.sampleT <- c(ceiling(min(attr(nobs,"sample")[,"T"])),
                                         ceiling(max(attr(nobs,"sample")[,"T"])))

                      cat("        Sample size calculation with Generalized pairwise comparison\n", sep = "")
                      cat("        for a power of ",power," and type 1 error rate of ",alpha," \n\n", sep = "")
                      
                      cat(" - estimated sample size (mean [min;max]): ",nobs[,"C"]," [",range.sampleC[1],";",range.sampleC[2],"] controls\n",
                          "                                           ",nobs[,"T"]," [",range.sampleT[1],";",range.sampleT[2],"] treated\n\n",sep="")
                          
                  }else{
                      cat("        Simulation study with Generalized pairwise comparison\n", sep = "")
                      cat("        with ",n.rep[1]," samples\n\n", sep = "")
                  }
                  rm.duplicate <- c("n.T", "n.C", "rep.estimate", "rep.se", "mean.estimate", "sd.estimate")
                  
                  for(iStatistic in statistic){
                      if(all(is.na(restriction))){
                          name.statistic <- switch(iStatistic,
                                                   "netBenefit" = "net benefit",
                                                   "winRatio" = "win ratio",
                                                   "favorable" = "proportion in favor of treatment",
                                                   "unfavorable" = "proportion in favor of control"
                                                   )
                      }else{
                          name.statistic <- switch(iStatistic,
                                                   "netBenefit" = "restricted net benefit",
                                                   "winRatio" = "restricted win ratio",
                                                   "favorable" = "restricted proportion in favor of treatment",
                                                   "unfavorable" = "restricted proportion in favor of control"
                                                   )
                      }
                      cat(" - ",name.statistic," statistic (null hypothesis Delta=",null[statistic],")\n", sep = "")

                      print(ls.df.print[[iStatistic]], row.names = FALSE, quote = FALSE)
                      cat("\n")
                  }
                  
                  if(legend){
                      M <- rbind(c(" n.T",":","number of observations in the treatment group"),
                                 c(" n.C",":","number of observations in the control group"),
                                 c(" mean.estimate",":","average estimate over simulations"),
                                 c(" sd.estimate",":","standard deviation of the estimate over simulations"))
                      if(method.inference != "none"){                          
                          M <- rbind(M,
                                     c(" mean.se",":","average estimated standard error of the estimate over simulations"),
                                     c(" rejection",":","frequency of the rejection of the null hypothesis over simulations")
                                     )
                          txt.note <- paste0("(standard error: H-projection of order ",order.Hprojection,"| p-value:")
                          if(!is.null(transformation)){
                              txt.note <- paste0(txt.note," after transformation) \n", sep="")
                          }else{
                              txt.note <- paste0(txt.note," original scale) \n", sep="")
                          }
                      }else{
                          txt.note <- NULL
                      }
                      if(col.rep){
                          M <- rbind(M,
                                     c(" rep.estimate",":","number of sucessful simulations for the point estimation"),
                                     c(" rep.se",":","number of sucessful simulations for the estimation of the standard error"),
                                     )
                      }
                      
                      nchar.1 <- sapply(M[,1],nchar)
                      M[,1] <- paste0(M[,1],
                                      sapply(max(nchar.1) - nchar.1, function(iX){paste0(rep(" ",time = iX),collapse = "")}))
                      txt.legend <- apply(M, 1, function(iRow){paste(iRow[1],iRow[2]," ",iRow[3],"\n",sep = "")})
                      cat(txt.legend,sep ="")
                      cat(txt.note,sep ="")
                      cat("\n")
                  }
              }
              
              ## ** export
              return(invisible(ls.df.print))
          }
)

## * summary S4BuyseTest
#' @docType methods
#' @name S4BuyseTest-summary
#' @title Summary Method for Class "S4BuyseTest"
#' @aliases summary,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param print [logical] Should the results be displayed in the console?
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param strata [logical] should the strata-specific results be displayed or the results pooled across strata?
#' Can also be \code{NULL} to display both.
#' @param type.display [numeric or character] the results/summary statistics to be displayed.
#' Either an integer indicating refering to a type of display in \code{BuyseTest.options()}
#' or the name of the column to be output (e.g. \code{c("strata","Delta","p.value")}).
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param ... arguments to be passed to \code{\link{S4BuyseTest-confint}}
#'
#' @details
#' \bold{Content of the output} \cr
#' The "results" table in the output show the result of the GPC at each endpoint, as well as its contribution to the global statistics.
#' More precisely, the column:
#' \itemize{
#'   \item \code{endpoint} lists the endpoints, by order of priority.
#'   \item \code{threshold} lists the threshold associated to each endpoint.
#'   \item \bold{weight:} lists the weight of each priority.
#'   \item \bold{strata:} list the strata relative to which the results of the priority are displayed. If \code{"global"}, then the results are over all strata at a given priority.
#'   \item \code{total} (or \code{total(\%)}) lists the number (or percentage) of pairs to be analyzed at the current priority (or strata).
#'   \item \code{favorable} (or \code{favorable(\%)}) lists the number (or percentage) of pairs classified in favor of the treatment group at the current priority (or strata).
#'   \item \code{unfavorable} (or \code{unfavorable(\%)}) lists the number (or percentage) of pairs classified in favor of the control group at the current priority (or strata).
#'   \item \code{neutral} (or \code{neutral(\%)}) lists the number (or percentage) of pairs classified as neither in favor of the treatment group nor in favor of the control group at the current priority (or strata).
#'   \item \code{uninf} (or \code{uninf(\%)}) lists the number (or percentage) of pairs that could not be classified at the current priority (or strata) due to missing values/censoring.
#'   \item \code{delta} lists the value of the priority-specific statistic (e.g. net benefit or win ratio), i.e. computed on the pairs analyzed at the current priority only.
#'   \item \code{Delta} lists the value of the cumulative statistic (e.g. net benefit or win ratio), i.e. computed on all the pairs analyzed up to the current priority.
#'   \item \code{Delta(\%)} lists the relative statistic (i.e. statistic up to the current priority divided by the final statistic).
#'   \item \code{information(\%)} lists the information fraction (i.e. number of favorable and unfavorable pairs up to the current priority divided by the final number of favorable and unfavorable pairs).
#'   \item \code{CI} lists the confidence intervals for \code{Delta} (not adjusted for multiple comparison).
#'   \item \code{null} lists the null hypothesis (\code{Delta=null}).
#'   \item \code{p.value} p-value relative to the null hypothesis (not adjusted for multiple comparison).
#'   \item \code{resampling} number of samples used to compute the confidence intervals or p-values from permutations or bootstrap samples.
#' Only displayed if some bootstrap samples have been discarded, for example, they did not lead to sample any case or control.
#' }
#' Note: when using the Peron scoring rule or a correction for uninformative pairs, the columns \code{total}, \code{favorable}, \code{unfavorable}, \code{neutral}, and \code{uninf} are computing by summing the contribution of the pairs. This may lead to a decimal value.
#' 
#' \bold{Statistic}: when considering a single endpoint and denoting
#' \eqn{Y} the endpoint in the treatment group,
#' \eqn{X} the endpoint in the control group,
#' and \eqn{\tau} the threshold of clinical relevance,
#' the net benefit is \eqn{P[Y \ge X + \tau] - P[X \ge Y + \tau]},
#' the win ratio is \eqn{\frac{P[Y \ge X + \tau]}{P[X \ge Y + \tau]}},
#' the proportion in favor of treatment is \eqn{P[Y \ge X + \tau]},
#' the proportion in favor of control is \eqn{P[X \ge Y + \tau]}.
#'
#' \bold{Statistical inference} \cr
#' When the interest is in obtaining p-values, we recommand the use of a permutation test.
#' However, when using a permutation test confidence intervals are not displayed in the summary.
#' This is because there is no (to the best of our knowledge) straightforward way to obtain good confidence intervals with permutations. 
#' An easy way consist in using the quantiles of the permutation distribution and then shift by the point estimate of the statistic.
#' This is what is output by \code{\link{S4BuyseTest-confint}}.
#' However this approach leads to a much too high coverage when the null hypothesis is false.
#' The limits of the confidence interval can also end up being outside of the interval of definition of the statistic
#' (e.g. outside [-1,1] for the proportion in favor of treatment).
#' Therefore, for obtaining confidence intervals, we recommand the boostrap method or the u-statistic method.
#'
#' \bold{Win ratio} \cr
#' For the win ratio, the proposed implementation enables the use of thresholds and endpoints that are not time to events
#' as well as the correction proposed in Peron et al. (2016) to account for censoring. 
#' These development have not been examined by Wang et al. (2016), or in other papers (to the best of our knowledge).
#' They are only provided here by implementation convenience.
#'
#' \bold{Competing risks} \cr
#' In presence of competing risks, looking at the net benefit/win ratio computed with respect to the event of interest
#' will likely not give a full picture of the difference between the two groups.
#' For instance a treatment may decrease the risk of the event of interest (i.e. increase the net benefit for this event)
#' by increasing the risk of the competing event. If the competing event is death, this is not desirable. It is therefore advised to
#' taking into consideration the risk of the competing event, e.g. by re-running BuyseTest where cause 1 and 2 have been inverted.
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuyseTest-model.tables}} to obtain the table displayed at the end of the summary method in a \code{data.frame} format.
#'   \code{\link{S4BuyseTest-confint}} to output estimate, standard errors, confidence interval and p-values.
#'   \code{\link{S4BuyseTest-plot}} for a graphical display of the scoring of the pairs.
#'   \code{\link{BuyseMultComp}} for efficient adjustment for multiple comparisons.
#' 
#' @examples
#' library(data.table)
#' 
#' dt <- simBuyseTest(1e2, n.strata = 3)
#' 
#'  \dontrun{
#'  BT <- BuyseTest(treatment ~ TTE(eventtime, status = status) + Bin(toxicity), data=dt)
#'  }
#'  \dontshow{
#'  BT <- BuyseTest(treatment ~ TTE(eventtime, status = status) + Bin(toxicity), data=dt, n.resampling = 10, trace = 0)
#'  }
#'  summary(BT)
#'  summary(BT, percentage = FALSE)
#'  summary(BT, statistic = "winRatio")
#'
#' @keywords print
#' @author Brice Ozenne
#' @exportMethod summary
setMethod(f = "summary",
          signature = "S4BuyseTest",
          definition = function(object, print = TRUE, percentage = TRUE, statistic = NULL, conf.level = NULL, strata = NULL,
                                type.display = 1, digit = c(2,4,5), ...){

              ## ** normalize and check arguments
              mycall <- match.call()
              option <- BuyseTest.options()
              
              if(length(digit) == 1){digit <- rep(digit,3)}
              validInteger(digit,
                           name1 = "digit",
                           min = 0,
                           valid.length = 3,
                           method = "summary[S4BuyseTest]")

              if(is.numeric(type.display)){
                  validInteger(type.display,
                               name1 = "type.display",
                               min = 1,
                               max = 9, ## limitation in model.tables
                               valid.length = 1)
                  type.display.original <- paste0("summary",type.display)
                  type.display <- option$summary.display[[type.display]]
              }else{
                  type.display.original <- NA
                  validCharacter(type.display,
                                 name1 = "type.display",
                                 valid.values = c("endpoint","restriction","threshold","strata","weight","total","favorable","unfavorable","neutral","uninf","information(%)",
                                                  "delta","Delta","Delta(%)",
                                                  "p.value","CI","significance"),
                                 valid.length = NULL)
              }
              if(is.null(conf.level)){
                  conf.level <- option$conf.level
              }
              alpha <- 1-conf.level              
              method.inference <- slot(object,"method.inference")
              hierarchical <- slot(object,"hierarchical")
              scoring.rule <- slot(object,"scoring.rule")

              scoring.rule_like_peron <- scoring.rule %in% c("peron","efron","latta",names(survival::survreg.distributions))


              if(is.null(statistic)){
                  statistic <- option$statistic
              }else{
                  statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                      "netbenefit" = "netBenefit",
                                      "winratio" = "winRatio",
                                      "favorable" = "favorable",
                                      "unfavorable" = "unfavorable",
                                      statistic)
              }
              endpoint <- slot(object,"endpoint")
              n.endpoint <- length(endpoint)
              n.strata <- length(slot(object,"level.strata"))
              if(attr(scoring.rule,"test.match") && "strata" %in% names(mycall) == FALSE){
                  strata <- FALSE
              }

              ## ** build table
              if(attr(method.inference,"permutation")){
                  attr(conf.level,"warning.permutation") <- FALSE
              }
              columns.tables <- union(setdiff(type.display,"significance"),"null")
              if(method.inference == "none"){
                  columns.tables <- setdiff(columns.tables, c("CI","p.value"))
              }else if(attr(method.inference,"permutation")){
                  columns.tables <- setdiff(columns.tables, "CI")
              }else if("CI" %in% type.display){
                  columns.tables <- c(columns.tables[1:which(type.display=="CI")],
                                      columns.tables[which(type.display=="CI"):length(columns.tables)])
                  columns.tables[columns.tables=="CI"] <- c("lower.ci","upper.ci")
              }
              table.print <- model.tables(object, percentage = percentage, statistic = statistic, conf.level = conf.level, strata = strata,
                                          columns = columns.tables, ...)
              null <- table.print$null
              endpoint.restriction.threshold <- attr(table.print,"endpoint")
              transform <- attr(table.print,"transform")
              n.resampling <- attr(table.print,"n.resampling")
              method.ci.resampling <- attr(table.print,"method.ci.resampling")

              ## ** remove unnecessary columns
              if("null" %in% type.display == FALSE){
                  table.print$null <- NULL
              }
              name.print <- names(table.print)

              if(!is.na(type.display.original) && "restriction" %in% name.print && all(is.na(table.print[["restriction"]]))){
                  table.print$restriction <- NULL
                  name.print <- setdiff(name.print, "restriction")
              }
              if(!is.na(type.display.original) && "threshold" %in% name.print && all(table.print[["threshold"]]<1e-10)){
                  table.print$threshold <- NULL
                  name.print <- setdiff(name.print, "threshold")
              }
              if(!is.na(type.display.original) && "weight" %in% name.print && all(table.print[["weight"]]==1)){
                  table.print$weight <- NULL
                  name.print <- setdiff(name.print, "weight")
              }
              if(!is.na(type.display.original) && "strata" %in% name.print && n.strata==1){
                  table.print$strata <- NULL
                  name.print <- setdiff(name.print, "strata")
              }
              if(!is.na(type.display.original) && "delta" %in% name.print && n.endpoint==1 && n.strata == 1){
                  table.print$delta <- NULL
                  name.print <- setdiff(name.print, "delta")
              }

              ## CI when the estimate is not defined
              if("lower.ci" %in% name.print && "upper.ci" %in% name.print){
                  index.ci <- intersect(which(!is.na(table.print$lower.ci)),
                                        which(!is.na(table.print$upper.ci)))
                  if("Delta" %in% name.print){
                      index.ci <- intersect(index.ci,
                                            intersect(which(!is.infinite(table.print$Delta)),
                                                      which(!is.na(table.print$Delta))))                      
                  }
                  
              }else{
                  index.ci <- NULL
              }

              ## ** reformat table
              ## *** add column with stars
              if(("p.value" %in% name.print) && ("significance" %in% type.display)){
                  colStars <- sapply(table.print[,"p.value"],function(x){
                      if(is.na(x)){""}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
                  })
                  table.print[,"significance"] <- colStars
              }

              ## *** rounding
              ## counts
              col.pairs <- intersect(name.print, c("total","favorable","unfavorable","neutral","uninf"))
              if(!is.na(digit[1]) && length(col.pairs)>0){
                  table.print[,col.pairs] <- sapply(table.print[,col.pairs,drop=FALSE], round, digits = digit[1])
              }
              col.inference <- intersect(name.print, c("delta","Delta","lower.ci","upper.ci","Delta(%)","information(%)"))
              if(!is.na(digit[2]) && length(col.inference)>0){
                      table.print[,col.inference] <- sapply(table.print[,col.inference,drop=FALSE], round, digits = digit[2])
              }
              if(!is.na(digit[3]) && ("p.value" %in% name.print)){
                  table.print[!is.na(table.print$p.value),"p.value"] <- format.pval(table.print[!is.na(table.print$p.value),"p.value"], digits = digit[3])                      
              }

              ## *** set Inf to NA in summary
              ## e.g. in the case of no unfavorable pairs the win ratio is Inf
              ##      this is not a valid estimate and it is set to NA
              if("delta" %in% name.print){
                  table.print[is.infinite(table.print$delta), "delta"] <- NA
                  table.print[is.nan(table.print$delta), "delta"] <- NA
              }
              if("Delta" %in% name.print){
                  table.print[is.infinite(table.print$Delta), "Delta"] <- NA
                  table.print[is.nan(table.print$Delta), "Delta"] <- NA
              }
              
              ## *** convert NA to ""
              if("threshold" %in% name.print){
                  threshold.operator <- ifelse(attr(slot(object,"threshold"),"multiplicative"),"*","+")
                  threshold.rm <- which(table.print$threshold<=1e-12)
                  table.print$threshold <- paste0(threshold.operator,table.print$threshold)
                  table.print$threshold[threshold.rm] <- ""                  
              }
              if("restriction" %in% name.print){
                  table.print[is.na(table.print$restriction), "restriction"] <- ""
              }
              if("Delta" %in% name.print){
                  table.print[is.na(table.print$Delta), "Delta"] <- ""
              }
              if("lower.ci" %in% name.print){
                  table.print[is.na(table.print$lower), "lower.ci"] <- ""
              }
              if("upper.ci" %in% name.print){
                  table.print[is.na(table.print$upper), "upper.ci"] <- ""
              }
              if("p.value" %in% name.print){
                  table.print[is.na(table.print$p.value), "p.value"] <- ""
              }

              ## *** remove duplicated values in endpoint/threshold
              if("endpoint" %in% name.print){                  
                  table.print$endpoint[duplicated(endpoint.restriction.threshold)] <- ""
              }
              if("threshold" %in% name.print){                  
                  table.print$threshold[duplicated(endpoint.restriction.threshold)] <- ""
              }
              if("restriction" %in% name.print){                  
                  table.print$restriction[duplicated(endpoint.restriction.threshold)] <- ""
              }
              if("weight" %in% name.print){                  
                  table.print$weight[duplicated(endpoint.restriction.threshold)] <- ""
              }

              ## *** merge CI inf and CI sup column
              if("CI" %in% type.display && "lower.ci" %in% name.print && "upper.ci" %in% name.print){
                  qInf <- round(100*alpha/2, digits = digit[2])
                  qSup <- round(100*(1-alpha/2), digits = digit[2])
                  name.ci <- paste0("CI [",qInf,"% ; ",qSup,"%]")
                  table.print$upper.ci[index.ci] <- paste0("[",table.print[index.ci,"lower.ci"],
                                                              ";",table.print[index.ci,"upper.ci"],"]")

                  names(table.print)[names(table.print) == "upper.ci"] <- name.ci
                  table.print$lower.ci <- NULL
              }

              ## *** rename percentage columns
              vec.tfunu <- intersect(c("total","favorable","unfavorable","neutral","uninf"), name.print)
              if(identical(percentage,TRUE) & length(vec.tfunu)>0){
                  names(table.print)[match(vec.tfunu, names(table.print))] <- paste0(names(table.print)[names(table.print) %in% vec.tfunu],"(%)")
              }

              ## ** display
              ## *** additional text
              if(n.endpoint>1){
                  txt.endpoint <- paste0("with ",n.endpoint," ",ifelse(hierarchical, "prioritized ", ""),"endpoints", sep = "")
              }else{
                  txt.endpoint <- paste0("with 1 endpoint")
              }
              txt.strata <- if(n.strata>1){paste0(" and ",n.strata," strata")}else{""}
                  
              ## *** display
              if(print){
                  cat("       Generalized pairwise comparisons ",txt.endpoint,txt.strata,"\n\n", sep = "")

                  txt.statistic <- c(" - statistic       : ",
                                     ifelse(any(!is.na(object@restriction)),"restricted",""),
                                     switch(statistic,
                                            "winRatio" = ifelse(object@add.halfNeutral,"win odds","win ratio"),
                                            "netBenefit" = "net treatment benefit",
                                            "favorable" = ifelse(object@add.halfNeutral,"proportion in favor of treatment","proportion strictly in favor of treatment"),
                                            "unfavorable" = ifelse(object@add.halfNeutral,"proportion in favor of control","proportion strictly in favor of control")
                                            )
                                     )
                  cat(paste(txt.statistic, collapse = "")," (delta: endpoint specific, Delta: global) \n")

                  if(method.inference != "none"){
                      if(all(is.na(null))){
                          cat(" - null hypothesis : requires specification of the argument \'null\' \n", sep = "")
                          table.print$p.value <- NULL
                      }else{
                          cat(" - null hypothesis : Delta == ",stats::na.omit(null)[1]," \n", sep = "")
                      }
                  }

                  if(method.inference != "none"){
                      cat(" - confidence level: ",1-alpha," \n", sep = "")

                      if(attr(method.inference,"permutation")){
                          txt.method <- "permutation test"
                      }else if(attr(method.inference,"bootstrap")){
                          txt.method <- "bootstrap resampling"
                      }else if(attr(method.inference,"ustatistic")){
                          test.model.tte <- all(unlist(lapply(object@iidNuisance,dim))==0)
                          if(attr(object@scoring.rule,"test.match")){
                              txt.method <- paste0("variability of the estimate across strata")
                          }else{
                              txt.method <- paste0("H-projection of order ",attr(method.inference,"hprojection"))
                          }
                          if(transform != "id"){
                              txt.method <- paste0(txt.method," after ",transform," transformation \n")
                          }else{
                              txt.method <- paste0(txt.method," \n")
                          }
                          if(test.model.tte && (scoring.rule_like_peron || object@correction.uninf > 0)){
                              txt.method <- paste0(txt.method,"                     (ignoring the uncertainty of the nuisance parameters)")
                          }
                      }

                  if(attr(method.inference,"permutation") || attr(method.inference,"bootstrap") ){
                      
                      if(method.inference == "varexact permutation"){
                          txt.method <- paste0(txt.method, " with all possible samples \n")
                          table.print$n.resampling <- NULL
                      }else if(all(n.resampling[1]==n.resampling)){
                          txt.method <- paste0(txt.method, " with ",n.resampling[1]," samples \n")
                          table.print$n.resampling <- NULL
                      }else{
                          txt.method <- paste0(txt.method, " with [",min(n.resampling)," ; ",max(n.resampling),"] samples \n")
                      }

                      if(attr(method.inference,"permutation")){
                          txt.method.ci <- switch(method.ci.resampling,
                                                  "percentile" = "p-value computed using the permutation distribution",
                                                  "studentized" = "p-value computed using the studentized permutation distribution",
                                                  )
                      }else if(attr(method.inference,"bootstrap")){
                          txt.method.ci <- switch(method.ci.resampling,
                                                  "percentile" = "CI computed using the percentile method; p-value by test inversion",
                                                  "gaussian" = "CI/p-value computed assuming normality",
                                                  "studentized" = "CI computed using the studentized method; p-value by test inversion",
                                                  )
                      }
                          
                      txt.method <- paste0(txt.method,"                     ",txt.method.ci," \n")
                  }
                  cat(" - inference       : ",txt.method, sep = "")
              }
                  
              cat(" - treatment groups: ",object@level.treatment[2]," (treatment) vs. ",object@level.treatment[1]," (control) \n", sep = "")

                  if(n.strata>1){
                      txt.strataPC <- paste0(round(100*object@weightStrata, digit[1]),"%")
                      table.strataPC <- table(txt.strataPC)
                      if(any(table.strataPC>1)){
                          shorttxt.strataPC <- ifelse(table.strataPC>1,paste0(names(table.strataPC), " (x",table.strataPC,")"),names(table.strataPC))
                      }else{
                          shorttxt.strataPC <- txt.strataPC
                      }
                      cat(" - strata weights  : ",paste(shorttxt.strataPC, collapse = ", ")," \n", sep = "")
                  }else if(attr(scoring.rule,"test.match") & length(object@weightStrata)>1){
                      table.weightStrata <- table(paste0(round(100*object@weightStrata, digit[1]),"%"))
                      cat(" - pair weights    : ",paste(names(table.weightStrata), collapse = ", ")," (K=",paste(table.weightStrata, collapse=","),")\n", sep = "")
                  }else  if(any(object@type == "tte") && any(attr(scoring.rule,"test.censoring"))){

                          if(all(attr(scoring.rule,"method.score")[object@type=="tte"]=="CRPeron")){
                              txt.Peron <- "cif"
                          }else if(all(attr(scoring.rule,"method.score")[object@type=="tte"]=="SurvPeron")){
                              txt.Peron <- "survival"
                          }else{
                              txt.Peron <- "survival/cif"
                          }
                          if(scoring.rule_like_peron){
                              txt.scoring.rule <- paste0("probabilistic score based on the ",txt.Peron," curves")
                              if(scoring.rule == "efron"){
                                  txt.scoring.rule <- paste0(txt.scoring.rule, "\n \t\t     (set to 0 beyond available follow-up)")
                              }
                              if(scoring.rule == "latta"){
                                  txt.scoring.rule <- paste0(txt.scoring.rule, "\n \t\t     (common ",txt.Peron," across treatment groups)")
                              }
                          }else if(scoring.rule == "gehan"){
                              txt.scoring.rule <- "deterministic score or uninformative"
                          }
                          cat(" - censored pairs  : ",txt.scoring.rule,"\n", sep = "")
                      }
                  if(hierarchical && n.endpoint>1 && any(object@count.neutral>0)){
                      Uneutral.as.uninf <- unique(object@neutral.as.uninf)
                      if(identical(Uneutral.as.uninf,TRUE)){
                          txt.neutral <- "re-analyzed using lower priority endpoints"
                      }else if(identical(Uneutral.as.uninf,FALSE)){
                          txt.neutral <- "ignored at lower priority endpoints"
                      }else{
                          txt.neutral <- paste0("re-analyzed using lower priority endpoints for endpoint ",
                                                paste(which(object@neutral.as.uninf), collapse = ", "),
                                                " \n                     otherwise ignored at lower priority endpoints")
                      }
                      cat(" - neutral pairs   : ",txt.neutral,"\n", sep = "")
                  }
                  if(!( (object@correction.uninf == 0) && (all(object@count.uninf==0)) )){
                      txt.uninf <- switch(as.character(object@correction.uninf),
                                          "0" = if(n.endpoint==1){"no contribution"}else{"no contribution at the current endpoint, analyzed at later endpoints"},
                                          "1" = "score equals the averaged score of all informative pairs",
                                          "2" = "no contribution, their weight is passed to the informative pairs using IPCW"
                                          )
                      cat(" - uninformative pairs: ",txt.uninf,"\n", sep = "")
                  }

              cat(" - results\n")
              table.print2 <- table.print
              if("significance" %in% names(table.print)){
                  names(table.print2)[names(table.print2) == "significance"] <- ""
              }
              print(table.print2, row.names = FALSE)
              }
                  
              ## ** export
              return(invisible(table.print))
          }
          )




##----------------------------------------------------------------------
### summary.R ends here
