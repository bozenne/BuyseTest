## * Documentation - summary
#' @docType methods
#' @name BuyseRes-summary
#' @title Summary Method for Class "BuyseRes"
#' @aliases summary,BuyseRes-method
#' @include BuyseRes-object.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param print [logical] Should the table be displayed?.
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' whereas \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016).
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param strata [character vector] the name of the strata to be displayed. Can also be \code{"global"} to display the average over all strata.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param ... arguments to be passed to \code{\link{BuyseRes-confint}}
#'
#' @details
#' \bold{Content of the output} \cr
#' The "results" table in the output show the result of the GPC at each endpoint, as well as its contribution to the global statistics.
#' More precisely, the column:
#' \itemize{
#'   \item \code{endpoint} lists the endpoints, by order of priority.
#'   \item \code{threshold} lists the threshold associated to each endpoint.
#'   \item \code{total} lists the total number of pairs (argument \code{percentage=FALSE}) or the percentage of pairs (argument \code{percentage=TRUE}) to be analyzed at the current priority.
#'   \item \code{favorable} lists the number of pairs (argument \code{percentage=FALSE}) or the percentage of pairs (argument \code{percentage=TRUE}) classified in favor of the treatment at the current priority.
#'   \item \code{unfavorable} lists the number of pairs (argument \code{percentage=FALSE}) or the percentage of pairs (argument \code{percentage=TRUE}) classified in favor of the control at the current priority.
#'   \item \code{neutral} lists the number of pairs (argument \code{percentage=FALSE}) or the percentage of pairs (argument \code{percentage=TRUE}) classified as neutral at the current priority.
#'   \item \code{uninf} lists the number of pairs (argument \code{percentage=FALSE}) or the percentage of pairs (argument \code{percentage=TRUE}) that could not be classified at the current priority (due to missing values/censoring).
#'   \item \code{delta} lists the value of the statistic (i.e. net benefit or win ratio) computed on the pairs analyzed at the current priority only.
#'   \item \code{Delta} lists the value of the statistic (i.e. net benefit or win ratio) computed on all the pairs analyzed up to the current priority.
#'   \item \code{CI} Confidence interval for the value of \code{Delta} (performed independently at each priority, no adjustment for multiple comparison).
#'   \item \code{p.value} p-value for the test \code{Delta=0} (performed independently at each priority, no adjustment for multiple comparison).
#'   \item \code{resampling} number of samples used to compute the confidence intervals or p-values from permutations or bootstrap samples.
#' Only displayed if some bootstrap samples have been discarded, for example, they did not lead to sample any case or control.
#' }
#' Note: when using the Peron scoring rule or a correction for uninformative pairs, the columns \code{total}, \code{favorable}, \code{unfavorable}, \code{neutral}, and \code{uninf} are computing by summing the contribution of the pairs. This may lead to a decimal value.
#' 
#' \bold{Statistical inference} \cr
#' When the interest is in obtaining p-values, we recommand the use of a permutation test.
#' However, when using a permutation test confidence intervals are not displayed in the summary.
#' This is because there is no (to the best of our knowledge) straightforward way to obtain good confidence intervals with permutations. 
#' An easy way consist in using the quantiles of the permutation distribution and then shift by the point estimate of the statistic.
#' This is what is output by \code{\link{BuyseRes-confint}}.
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
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#'   \code{\link{BuyseRes-confint}} to output confidence interval and p-values in a matrix format.
#' 
#' @examples
#' library(data.table)
#' 
#' dt <- simBuyseTest(1e2, n.strata = 3)
#' 
#'  \dontrun{
#'  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status) + Bin(toxicity), data=dt)
#'  }
#'  \dontshow{
#'  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status) + Bin(toxicity), data=dt, n.resampling = 10, trace = 0)
#'  }
#'  summary(BT)
#'  summary(BT, percentage = FALSE)
#'  summary(BT, statistic = "winRatio")
#' 
#' @keywords summary BuyseRes-method

## * method - summary
#' @rdname BuyseRes-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "BuyseRes",
          definition = function(object, print = TRUE, percentage = TRUE, statistic = NULL,
                                conf.level = NULL,
                                strata = if(length(object@level.strata)==1){"global"}else{NULL},
                                digit = c(2,4), ...){

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
                           method = "summary[BuyseRes]")
              
              validLogical(percentage,
                           name1 = "percentage",
                           valid.length = 1,
                           refuse.NA = FALSE, 
                           method = "summary[BuyseRes]")

              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio"),
                             valid.length = 1,
                             method = "summary[BuyseRes]")

              validCharacter(strata,
                             name1 = "strata",
                             valid.length = NULL,
                             valid.values = c("global",object@level.strata),
                             refuse.NULL = FALSE,
                             method = "summary[BuyseRes]")

              if(length(digit) == 1){digit <- rep(digit,2)}
              validInteger(digit,
                           name1 = "digit",
                           min = 0,
                           valid.length = 2,
                           method = "summary[BuyseRes]")

              ## ** load info from object
              hierarchical <- object@hierarchical
              endpoint <- object@endpoint
              n.endpoint <- length(endpoint)
              n.strata <- length(object@level.strata)

              delta <- slot(object, name = paste0("delta.",statistic))
              Delta <- slot(object, name = paste0("Delta.",statistic))
              n.resampling <- object@n.resampling

              method.inference <- object@method.inference
              if(is.na(conf.level)){
                  method.inference[] <- "none" ## uses [] to not remove the attributees of method.inference
              }

              alpha <- 1-conf.level
              
              ## ** compute confidence intervals and p-values
              outConfint  <- confint(object, conf.level = conf.level, statistic = statistic, ...)

              ## ** generate summary table
              ## *** prepare
              table <- data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=15))
              names(table) <- c("endpoint","threshold","weight","strata",
                                "n.total","n.favorable","n.unfavorable","n.neutral","n.uninf",
                                "delta","Delta","CIinf.Delta","CIsup.Delta","p.value","n.resampling")
            
              index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1
            
              table[index.global,"n.favorable"] <- colSums(object@count.favorable)
              table[index.global,"n.unfavorable"] <- colSums(object@count.unfavorable)
              table[index.global,"n.neutral"] <- colSums(object@count.neutral)
              table[index.global,"n.uninf"] <- colSums(object@count.uninf)
              table[index.global,"n.total"] <- rowSums(table[index.global,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
            
              table[index.global,"endpoint"] <- object@endpoint
              table[index.global,"threshold"] <- object@threshold
              table[index.global,"weight"] <- object@weight
              table[index.global,"strata"] <- "global"

              if(statistic=="netBenefit"){ ##
                  table[index.global,"delta"] <- (colSums(object@count.favorable)-colSums(object@count.unfavorable))/sum(object@n.pairs)
              }else{
                  table[index.global,"delta"] <- colSums(object@count.favorable)/colSums(object@count.unfavorable)
              }
              table[index.global,"Delta"] <- Delta
             
              for(iStrata in 1:n.strata){
                  index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iStrata
              
                  table[index.strata,"n.favorable"] <- object@count.favorable[iStrata,]
                  table[index.strata,"n.unfavorable"] <- object@count.unfavorable[iStrata,]
                  table[index.strata,"n.neutral"] <- object@count.neutral[iStrata,]
                  table[index.strata,"n.uninf"] <- object@count.uninf[iStrata,]
                  table[index.strata,"n.total"] <- rowSums(table[index.strata,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
              
                  table[index.strata,"strata"] <- object@level.strata[iStrata]
                  table[index.strata,"endpoint"] <- object@endpoint
                  table[index.strata,"threshold"] <- object@threshold
                  table[index.strata,"weight"] <- object@weight
                  table[index.strata,"delta"] <- delta[iStrata,]
              }

              ## *** convert to percentage
              if(identical(percentage, TRUE)){
                  table[,"n.favorable"] <- 100*table[,"n.favorable"]/table[1,"n.total"]
                  table[,"n.unfavorable"] <- 100*table[,"n.unfavorable"]/table[1,"n.total"]
                  table[,"n.neutral"] <- 100*table[,"n.neutral"]/table[1,"n.total"]
                  table[,"n.uninf"] <- 100*table[,"n.uninf"]/table[1,"n.total"]
                  table[,"n.total"] <- 100*table[,"n.total"]/table[1,"n.total"]
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
                  if(method.inference %in% c("u-statistic","u-statistic-bebu")){
                      table.print <- cbind(table.print[,setdiff(names(table.print), "n.resampling")],
                                           "significance" = colStars)
                  }else{
                      table.print <- cbind(table.print[,setdiff(names(table.print), "n.resampling")],
                                           "significance" = colStars,
                                           "n.resampling" = table.print[["n.resampling"]])
                  }
              }

              ## *** restrict to strata
              if(!is.null(strata)){
                  table.print <- table.print[table.print$strata %in% strata,,drop = FALSE]                      
              }

              ## *** remove useless columns
              if(is.na(percentage)){
                  table.print$n.favorable <- NULL
                  table.print$n.unfavorable <- NULL
                  table.print$n.neutral <- NULL
                  table.print$n.uninf <- NULL
                  table.print$n.total <- NULL
              }

              if(method.inference == "none"){
                  keep.cols <- setdiff(names(table.print),
                                       c("CIinf.Delta","CIsup.Delta","n.resampling","p.value"))
                  table.print <- table.print[,keep.cols, drop = FALSE]
              }else if(attr(method.inference,"ustatistic")){
                  keep.cols <- setdiff(names(table.print), "n.resampling")
                  table.print <- table.print[,keep.cols, drop = FALSE]
              }else if(attr(method.inference,"permutation")){
                  keep.cols <- setdiff(names(table.print), c("CIinf.Delta","CIsup.Delta"))
                  table.print <- table.print[,keep.cols, drop = FALSE]
              }

              if(identical(strata, "global")){
                  keep.cols <- which(names(table.print) %in% setdiff(names(table.print), "strata"))
                  table.print <- table.print[,keep.cols,drop = FALSE]
              }
            
              ## *** rounding
              ## counts
              if(!is.na(percentage) && !is.na(digit[1])){
                  param.signif <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digits = digit[1])
              }
              if(!is.na(digit[2])){
                  param.signif <- c("delta","Delta")
                  if(method.inference != "none"){
                      table.print[!is.na(table.print$p.value),"p.value"] <- format.pval(table.print[!is.na(table.print$p.value),"p.value"], digits = digit[2])                      
                  }
                  if(method.inference != "none" && !attr(method.inference,"permutation")){
                      param.signif <- c(param.signif, "CIinf.Delta","CIsup.Delta")
                  }
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digits = digit[2])
              }
              ## *** set names
              if(identical(percentage,TRUE)){
                  oldnames <- c("n.favorable","n.unfavorable","n.neutral","n.uninf","n.total")
                  newnames <- c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf","pc.total")
                  names(table)[match(oldnames,names(table))] <- newnames
                  names(table.print)[match(oldnames,names(table.print))] <- newnames
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
              if(!is.null(table.print$CIinf.Delta) && any(is.na(table.print$CIinf.Delta))){
                  table.print[is.na(table.print$CIinf.Delta), "CIinf.Delta"] <- ""
              }
              if(!is.null(table.print$CIsup.Delta) && any(is.na(table.print$CIsup.Delta))){
                  table.print[is.na(table.print$CIsup.Delta), "CIsup.Delta"] <- ""
              }
              if(!is.null(table.print$p.value) && any(is.na(table.print$p.value))){
                  table.print[is.na(table.print$p.value), "p.value"] <- ""
              }

              ## *** remove name significance
              if(method.inference != "none"){
                  colnames(table.print)[colnames(table.print)=="significance"] <- ""
              }
                  
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
                  qInf <- round(100*alpha/2, digits = digit[2])
                  qSup <- round(100*(1-alpha/2), digits = digit[2])

                  names(table.print)[names(table.print) == "CIinf.Delta"] <- paste0("CI [",qInf," ; ",qSup,"]")
                  table.print$CIsup.Delta <- NULL
              }

              ## *** simplify names
              if(identical(percentage,TRUE)){
                  names(table.print) <- gsub("^pc.","",names(table.print))
              }else if(identical(percentage,FALSE)){
                  names(table.print) <- gsub("^n.","",names(table.print))
              }

              ## *** remove duplicated values in endpoint/threshold
              test.duplicated <- duplicated(interaction(table.print$endpoint,table.print$threshold))
              table.print[which(test.duplicated),c("endpoint","threshold","weight")] <- ""

              ## *** rm weights
              if(hierarchical){
                  table.print$weight <- NULL
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
                      cat(" > statistic       : win ratio (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 1 \n", sep = "")
                  }else {
                      cat(" > statistic       : net benefit (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 0 \n", sep = "")
                  }
                  if(method.inference != "none"){
                      cat(" > confidence level: ",1-alpha," \n", sep = "")

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
                          if(attr(method.inference,"bootstrap")){
                              txt.method.ci <- switch(attr(outConfint,"method.ci.resampling"),
                                                      "gaussian" = "quantiles of a Gaussian distribution",
                                                      "student" = "quantiles of a Student's t-distribution",
                                                      "percentile" = "quantiles of the empirical distribution"
                                                      )
                          }else{
                              txt.method.ci <- "quantiles of the empirical distribution"
                          }
                          txt.method <- paste0(txt.method,"                     confidence intervals/p-values computed using the ",txt.method.ci," \n")
                      }
                      cat(" > inference       : ",txt.method, sep = "")
                  }
                  
                  cat(" > treatment groups: ",object@level.treatment[1]," (control) vs. ",object@level.treatment[2]," (treatment) \n", sep = "")
                  if(any(object@type == "TimeToEvent")){
                      txt.scoring.rule <- switch(object@scoring.rule,
                                               "Gehan" = "uninformative pairs",
                                               "Peron" = "use Kaplan Meier survival curves to compute the score"
                                               ) 

                      cat(" > censored pairs  : ",txt.scoring.rule,"\n", sep = "")
                  }
                  if(!( (object@correction.uninf == 0) && (all(object@count.uninf==0)) )){
                      txt.uninf <- switch(as.character(object@correction.uninf),
                                          "0" = "no contribution at the current endpoint, analyzed at later endpoints (if any)",
                                          "1" = "score equals the averaged score of all informative pairs",
                                          "2" = "no contribution, their weight is passed to the informative pairs using IPCW"
                                          )
                      cat(" > uninformative pairs: ",txt.uninf,"\n", sep = "")
                  }
                  
                  cat(" > results\n")
                  print(table.print, row.names = FALSE)
              }
              ## ** export
              return(invisible(list(table = table,
                                    table.print = table.print))
                     )
            
          }
)



