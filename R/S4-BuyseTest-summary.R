## * Documentation - summary
#' @docType methods
#' @name S4BuyseTest-summary
#' @title Summary Method for Class "S4BuyseTest"
#' @aliases summary,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param print [logical] Should the table be displayed?.
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' \code{"unfavorable"} displays the proportion in favor of the control.
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param type.display [numeric or character] the results/summary statistics to be displayed.
#' Either an integer indicating refering to a type of display in \code{BuyseTest.options()}
#' or the name of the column to be output (e.g. \code{c("strata","Delta","p.value")}).
#' @param strata [character vector] the name of the strata to be displayed. Can also be \code{"global"} to display the average over all strata.
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
#'   \item \code{total} lists the total number of pairs to be analyzed at the current priority.
#'   \item \code{total(\%)} lists the total percentage of pairs to be analyzed at the current priority.
#'   \item \code{favorable} lists the number of pairs classified in favor of the treatment at the current priority.
#'   \item \code{favorable(\%)} lists the number of pairs classified in favor of the treatment at the current priority.
#'   \item \code{unfavorable} lists the number of pairs classified in favor of the control at the current priority.
#'   \item \code{unfavorable(\%)} lists the percentage of pairs classified in favor of the control at the current priority.
#'   \item \code{neutral} lists the number of pairs classified as neutral at the current priority.
#'   \item \code{neutral(\%)} lists the percentage of pairs classified as neutral at the current priority.
#'   \item \code{uninf} lists the number of pairs that could not be classified at the current priority (due to missing values/censoring).
#'   \item \code{uninf(\%)} lists the percentage of pairs that could not be classified at the current priority (due to missing values/censoring).
#'   \item \code{delta} lists the value of the statistic (e.g. net benefit or win ratio) computed on the pairs analyzed at the current priority only.
#'   \item \code{Delta} lists the value of the statistic (e.g. net benefit or win ratio) computed on all the pairs analyzed up to the current priority.
#'   \item \code{Delta(\%)} lists the  net benefit or win ratio fraction (i.e. statistic up to the current priority divided by the final statistic).
#'   \item \code{information(\%)} lists the information fraction (i.e. number of favorable and unfavorable pairs up to the current priority divided by the final number of favorable and unfavorable pairs).
#'   \item \code{CI} Confidence interval for the value of \code{Delta} (performed independently at each priority, no adjustment for multiple comparison).
#'   \item \code{p.value} p-value for the test \code{Delta=0} (performed independently at each priority, no adjustment for multiple comparison).
#'   \item \code{resampling} number of samples used to compute the confidence intervals or p-values from permutations or bootstrap samples.
#' Only displayed if some bootstrap samples have been discarded, for example, they did not lead to sample any case or control.
#' }
#' Note: when using the Peron scoring rule or a correction for uninformative pairs, the columns \code{total}, \code{favorable}, \code{unfavorable}, \code{neutral}, and \code{uninf} are computing by summing the contribution of the pairs. This may lead to a decimal value.
#' 
#' \bold{statistic}: when considering a single endpoint and denoting
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
#'   \code{\link{S4BuyseTest-class}} for a presentation of the \code{S4BuyseTest} object. \cr
#'   \code{\link{S4BuyseTest-confint}} to output confidence interval and p-values in a matrix format.
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
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937 \
#' 
#' @keywords summary S4BuyseTest-method
#' @author Brice Ozenne

## * method - summary
#' @rdname S4BuyseTest-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "S4BuyseTest",
          definition = function(object, print = TRUE, percentage = TRUE, statistic = NULL,
                                conf.level = NULL,
                                strata = if(length(object@level.strata)==1){"global"}else{NULL},
                                type.display = 1,
                                digit = c(2,4,5), ...){

              ## ** normalize and check arguments
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }
              if(is.null(conf.level)){
                  conf.level <- option$conf.level
              }
              
              validLogical(print,
                           name1 = "print",
                           valid.length = 1,
                           method = "summary[S4BuyseTest]")
              
              validLogical(percentage,
                           name1 = "percentage",
                           valid.length = 1,
                           refuse.NA = FALSE, 
                           method = "summary[S4BuyseTest]")

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
                             method = "summary[S4BuyseTest]")

              validCharacter(strata,
                             name1 = "strata",
                             valid.length = NULL,
                             valid.values = c("global",object@level.strata),
                             refuse.NULL = FALSE,
                             method = "summary[S4BuyseTest]")

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
                               max = length(option$summary.display),
                               valid.length = 1)
                  type.display <- option$summary.display[[type.display]]
              }else{
                  validCharacter(type.display,
                                 name1 = "type.display",
                                 valid.values = c("endpoint","threshold","strata","weight","total","favorable","unfavorable","neutral","uninf","information(%)",
                                                  "delta","Delta","Delta(%)",
                                                  "p.value","CI","significance"),
                                 valid.length = NULL)
              }
              
              ## ** load info from object
              hierarchical <- object@hierarchical
              endpoint <- object@endpoint
              n.endpoint <- length(endpoint)
              n.strata <- length(object@level.strata)

              delta <- coef(object, statistic = statistic, cumulative = FALSE, stratified = TRUE)
              Delta <- coef(object, statistic = statistic, cumulative = TRUE)
              n.resampling <- object@n.resampling

              method.inference <- object@method.inference
              if(is.na(conf.level)){
                  method.inference[] <- "none" ## uses [] to not remove the attributees of method.inference
              }

              alpha <- 1-conf.level
              qInf <- round(100*alpha/2, digits = digit[2])
              qSup <- round(100*(1-alpha/2), digits = digit[2])
              name.ci <- paste0("CI [",qInf,"% ; ",qSup,"%]")

              ## ** update type.display
              ## update the name of the columns according to the request
              if("CI" %in% type.display){
                  type.display[type.display=="CI"] <- name.ci
              }
              vec.tfunu <- c("total","favorable","unfavorable","neutral","uninf")
              if(is.na(percentage)){
                  type.display <- setdiff(type.display,vec.tfunu)
              }

              ## remove columns not requested by the user
              rm.display <- NULL
              if(hierarchical){
                  rm.display <- c(rm.display,"weight")
              }
              if(is.na(percentage)){
                  rm.display <- c(rm.display,"total","favorable","unfavorable","neutral","uninf")
              }

              if(method.inference == "none"){
                  rm.display <- c(rm.display,name.ci,"p.value","significance","n.resampling")
              }else if(attr(method.inference,"ustatistic")){
                  rm.display <- c(rm.display,"n.resampling")
              }else if(attr(method.inference,"permutation")){
                  rm.display <- c(rm.display,name.ci)
              }
              if(identical(strata, "global")){
                  rm.display <- c(rm.display,"strata")
              }
              if("delta" %in% type.display && "Delta" %in% type.display && n.endpoint == 1){
                  rm.display <- union("delta", rm.display)
              }
              type.display <- setdiff(type.display,rm.display)

              ## ** compute confidence intervals and p-values
              outConfint  <- confint(object, conf.level = conf.level, statistic = statistic, ...)

              ## ** generate summary table
              ## *** prepare
              table <- data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=18),
                                  stringsAsFactors = FALSE)
              names(table) <- c("endpoint","threshold","weight","strata",
                                "total","favorable","unfavorable","neutral","uninf",
                                "delta","Delta","Delta(%)","information(%)",
                                "CIinf.Delta","CIsup.Delta","p.value","significance","n.resampling")
            
              index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1

              table[index.global,"favorable"] <- colSums(object@count.favorable)
              table[index.global,"unfavorable"] <- colSums(object@count.unfavorable)
              table[index.global,"neutral"] <- colSums(object@count.neutral)
              table[index.global,"uninf"] <- colSums(object@count.uninf)
              table[index.global,"total"] <- rowSums(table[index.global,c("favorable","unfavorable","neutral","uninf")])
            
              table[index.global,"endpoint"] <- object@endpoint
              table[index.global,"threshold"] <- object@threshold
              table[index.global,"weight"] <- object@weight
              table[index.global,"strata"] <- "global"

              if(statistic=="netBenefit"){ ##
                  table[index.global,"delta"] <- (colSums(object@count.favorable)-colSums(object@count.unfavorable))/sum(object@n.pairs)
              }else if(statistic == "winRatio"){
                  table[index.global,"delta"] <- colSums(object@count.favorable)/colSums(object@count.unfavorable)
              }else if(statistic == "favorable"){
                  table[index.global,"delta"] <- colSums(object@count.favorable)/sum(object@n.pairs)
              }else if(statistic == "unfavorable"){
                  table[index.global,"delta"] <- colSums(object@count.unfavorable)/sum(object@n.pairs)
              }
              table[index.global,"Delta"] <- Delta
              table[index.global,"Delta(%)"] <- 100*Delta/Delta[n.endpoint]
             
              for(iStrata in 1:n.strata){
                  index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iStrata
              
                  table[index.strata,"favorable"] <- object@count.favorable[iStrata,]
                  table[index.strata,"unfavorable"] <- object@count.unfavorable[iStrata,]
                  table[index.strata,"neutral"] <- object@count.neutral[iStrata,]
                  table[index.strata,"uninf"] <- object@count.uninf[iStrata,]
                  table[index.strata,"total"] <- rowSums(table[index.strata,c("favorable","unfavorable","neutral","uninf")])
              
                  table[index.strata,"strata"] <- object@level.strata[iStrata]
                  table[index.strata,"endpoint"] <- object@endpoint
                  table[index.strata,"threshold"] <- object@threshold
                  table[index.strata,"weight"] <- object@weight
                  table[index.strata,"delta"] <- delta[iStrata,]
              }

              ## *** information fraction and co
              table[index.global,"information(%)"] <- 100*cumsum(colSums(object@count.favorable+object@count.unfavorable)/sum(object@count.favorable+object@count.unfavorable))

              ## *** convert to percentage
              if(identical(percentage, TRUE)){
                  table[,"favorable"] <- 100*table[,"favorable"]/table[1,"total"]
                  table[,"unfavorable"] <- 100*table[,"unfavorable"]/table[1,"total"]
                  table[,"neutral"] <- 100*table[,"neutral"]/table[1,"total"]
                  table[,"uninf"] <- 100*table[,"uninf"]/table[1,"total"]
                  table[,"total"] <- 100*table[,"total"]/table[1,"total"]
              }
              
              ## *** compute CI and p-value
              if(!attr(method.inference,"permutation")){
                  table[index.global,"CIinf.Delta"] <- outConfint[,"lower.ci"]
                  table[index.global,"CIsup.Delta"] <- outConfint[,"upper.ci"]
              }
              table[index.global,"p.value"] <- outConfint[,"p.value"]
              table[index.global,"n.resampling"] <- attr(outConfint,"n.resampling")

              ## ** generate print table
              table.print <- table

              ## *** add column with stars
              if(method.inference != "none"){
                  colStars <- rep("",NROW(table.print))
                  colStars[index.global] <- sapply(table.print[index.global,"p.value"],function(x){
                      if(is.na(x)){""}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
                  })
                  table.print[,"significance"] <- colStars
              }

              ## *** restrict to strata
              if(!is.null(strata)){
                  table.print <- table.print[table.print$strata %in% strata,,drop = FALSE]                      
              }

              ## *** rounding
              ## counts
              if(!is.na(digit[1])){
                  param.signif <- c("total","favorable","unfavorable","neutral","uninf")
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digits = digit[1])
              }
              if(!is.na(digit[2])){
                  param.signif <- c("delta","Delta","CIinf.Delta","CIsup.Delta","Delta(%)","information(%)")
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digits = digit[2])
              }
              if(!is.na(digit[3])){
                  table.print[!is.na(table.print$p.value),"p.value"] <- format.pval(table.print[!is.na(table.print$p.value),"p.value"], digits = digit[3])                      
              }

              ## *** set Inf to NA in summary
              ## e.g. in the case of no unfavorable pairs the win ratio is Inf
              ##      this is not a valid estimate and it is set to NA
              if(any(is.infinite(table.print$delta))){
                  table.print[is.infinite(table.print$delta), "delta"] <- NA
              }
              if(any(is.nan(table.print$delta))){
                  table.print[is.nan(table.print$delta), "delta"] <- NA
              }
              if(any(is.infinite(table.print$Delta))){
                  table.print[is.infinite(table.print$Delta), "Delta"] <- NA
              }
              if(any(is.nan(table.print$Delta))){
                  table.print[is.nan(table.print$Delta), "Delta"] <- NA
              }
              
              ## *** convert NA to ""
              if(any(is.na(table.print$Delta))){
                  table.print[is.na(table.print$Delta), "Delta"] <- ""
              }
              if(any(is.na(table.print$CIinf.Delta))){
                  table.print[is.na(table.print$CIinf.Delta), "CIinf.Delta"] <- ""
              }
              if(any(is.na(table.print$CIsup.Delta))){
                  table.print[is.na(table.print$CIsup.Delta), "CIsup.Delta"] <- ""
              }
              if(any(is.na(table.print$p.value))){
                  table.print[is.na(table.print$p.value), "p.value"] <- ""
              }

              ## *** remove duplicated values in endpoint/threshold
              test.duplicated <- duplicated(interaction(table.print$endpoint,table.print$threshold,table.print$weight))
              table.print[which(test.duplicated),c("endpoint","threshold","weight")] <- ""

              ## *** merge CI inf and CI sup column
              if(method.inference != "none" && !attr(method.inference,"permutation")){
                  if("strata" %in% names(table.print)){
                      index.tempo <- which(table.print$strata == "global")                                       
                  }else{
                      index.tempo <- 1:NROW(table.print)
                  }
                  ## remove CI when the estimate is not defined
                  index.tempo <- intersect(index.tempo,
                                           intersect(which(!is.infinite(table.print$Delta)),
                                                     which(!is.na(table.print$Delta)))
                                           )

                  table.print$CIinf.Delta[index.tempo] <- paste0("[",table.print[index.tempo,"CIinf.Delta"],
                                                                 ";",table.print[index.tempo,"CIsup.Delta"],"]")

                  names(table.print)[names(table.print) == "CIinf.Delta"] <- name.ci
                  table.print$CIsup.Delta <- NULL
              }

              ## *** select relevant columns
              table.print <- table.print[,type.display]
              ##              type.display[type.display %in% names(table.print) == FALSE]
              if(identical(percentage,TRUE) & any(vec.tfunu %in% type.display)){
                  names(table.print)[names(table.print) %in% vec.tfunu] <- paste0(names(table.print)[names(table.print) %in% vec.tfunu],"(%)")
              }

              ## ** display
              if(print){
                  ## *** additional text
                  if(n.endpoint>1){
                      txt.endpoint <- paste0("with ",n.endpoint," ",ifelse(hierarchical, "prioritized ", ""),"endpoints", sep = "")
                  }else{
                      txt.endpoint <- paste0("with 1 endpoint")
                  }
                  txt.strata <- if(n.strata>1){paste0(" and ",n.strata," strata")}else{""}
                  
                  ## *** display
                  cat("       Generalized pairwise comparisons ",txt.endpoint,txt.strata,"\n\n", sep = "")
                  if(statistic == "winRatio"){
                      cat(" - statistic       : win ratio (delta: endpoint specific, Delta: global) \n",
                          " - null hypothesis : Delta == 1 \n", sep = "")
                  }else {
                      cat(" - statistic       : net benefit (delta: endpoint specific, Delta: global) \n",
                          " - null hypothesis : Delta == 0 \n", sep = "")
                  }
                  if(method.inference != "none"){
                      cat(" - confidence level: ",1-alpha," \n", sep = "")

                      if(attr(method.inference,"permutation")){
                          txt.method <- "permutation test"
                      }else if(attr(method.inference,"bootstrap")){
                          txt.method <- "bootstrap resampling"
                      }else if(attr(method.inference,"ustatistic")){
                          test.model.tte <- all(unlist(lapply(object@iidNuisance,dim))==0)
                          txt.method <- paste0("H-projection of order ",attr(method.inference,"hprojection"),"\n")
                          if(test.model.tte && (object@scoring.rule == "Peron" || object@correction.uninf > 0)){
                              txt.method <- paste0(txt.method,"                     (ignoring the uncertainty of the nuisance parameters) \n")
                          }
                      
                      }

                      if(attr(method.inference,"permutation") || attr(method.inference,"bootstrap") ){
                          ok.resampling <- all(n.resampling[1]==n.resampling)
                          if(ok.resampling){
                              txt.method <- paste0(txt.method, " with ",n.resampling[1]," samples \n")
                              table.print$n.resampling <- NULL
                          }else{
                              txt.method <- paste0(txt.method, " with [",min(n.resampling)," ; ",max(n.resampling),"] samples \n")
                          }

                          if(attr(method.inference,"permutation")){
                              txt.method.ci <- switch(attr(outConfint,"method.ci.resampling"),
                                                      "percentile" = "p-value computed using the permutation distribution",
                                                      "studentized" = "p-value computed using the studentized permutation distribution",
                                                      )
                          }else if(attr(method.inference,"bootstrap")){
                              txt.method.ci <- switch(attr(outConfint,"method.ci.resampling"),
                                                      "percentile" = "CI computed using the percentile method; p-value by test inversion",
                                                      "gaussian" = "CI/p-value computed assuming normality",
                                                      "studentized" = "CI computed using the studentized method; p-value by test inversion",
                                                      )
                          }
                          
                          txt.method <- paste0(txt.method,"                     ",txt.method.ci," \n")
                      }
                      cat(" - inference       : ",txt.method, sep = "")
                  }
                  
                  cat(" - treatment groups: ",object@level.treatment[1]," (control) vs. ",object@level.treatment[2]," (treatment) \n", sep = "")
                  if(any(object@type == "TimeToEvent") && any(attr(object@scoring.rule,"test.censoring"))){
                      
                      if(all(attr(object@scoring.rule,"method.score")[object@type=="TimeToEvent"]==5)){
                          txt.Peron <- "cif"
                      }else if(all(attr(object@scoring.rule,"method.score")[object@type=="TimeToEvent"]==4)){
                          txt.Peron <- "survival"
                      }else{
                          txt.Peron <- "survival/cif"
                      }

                      txt.scoring.rule <- switch(object@scoring.rule,
                                                 "Gehan" = "deterministic score or uninformative",
                                                 "Peron" = paste0("probabilistic score based on the ",txt.Peron," curves")
                                                 )
                      cat(" - censored pairs  : ",txt.scoring.rule,"\n", sep = "")
                  }
                  if(n.endpoint>1 && any(object@count.neutral>0)){
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
              return(invisible(list(table = table,
                                    table.print = table.print))
                     )
            
          }
)



