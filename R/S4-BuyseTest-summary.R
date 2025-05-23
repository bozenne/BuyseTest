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
#' @param print [logical] Should the results be displayed in the console?
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' \code{"unfavorable"} displays the proportion in favor of the control.
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
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937.
#' 
#' @keywords print
#' @author Brice Ozenne

## * method - summary
#' @rdname S4BuyseTest-summary
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
                  table.print$threshold[table.print$threshold<=1e-12] <- ""
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
                          if(test.model.tte && (scoring.rule == "Peron" || object@correction.uninf > 0)){
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
                  cat(" - strata weights  : ",paste(paste0(round(100*object@weightStrata, digit[1]),"%"), collapse = ", ")," \n", sep = "")
              }else if(attr(scoring.rule,"test.match") & length(object@weightStrata)>1){
                  table.weightStrata <- table(paste0(round(100*object@weightStrata, digit[1]),"%"))
                  cat(" - pair weights    : ",paste(names(table.weightStrata), collapse = ", ")," (K=",paste(table.weightStrata, collapse=","),")\n", sep = "")
              }else 
              if(any(object@type == "tte") && any(attr(scoring.rule,"test.censoring"))){

                  if(all(attr(scoring.rule,"method.score")[object@type=="tte"]=="CRPeron")){
                      txt.Peron <- "cif"
                  }else if(all(attr(scoring.rule,"method.score")[object@type=="tte"]=="SurvPeron")){
                      txt.Peron <- "survival"
                  }else{
                      txt.Peron <- "survival/cif"
                  }

                  txt.scoring.rule <- switch(scoring.rule,
                                             "Gehan" = "deterministic score or uninformative",
                                             "Peron" = paste0("probabilistic score based on the ",txt.Peron," curves")
                                             )
                  if(attr(scoring.rule,"efron")){
                      txt.scoring.rule <- paste0(txt.scoring.rule, "\n \t\t     (set to 0 beyond available follow-up)")
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



