## * Documentation - summary
#' @docType methods
#' @name BuyseRes-summary
#' @title Summary Method for Class "BuyseRes"
#' @aliases summary summary,BuyseRes
#' @include BuyseRes-object.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param show [logical] Should the table be displayed?.
#' @param percentage [logical] Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netChance"} displays the net chance in favor of treatment, as described in Buyse (2010) and Peron et al. (2016)),
#' whereas \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016).
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param strata [character vector] the name of the strata to be displayed. Can also be \code{"global"} to display the average over all strata.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param ... arguments to be passed from the generic method to the class specific method [not relevant to the user]
#'
#' @details 
#' WARNING : when using a permutation test, the confidence interval is computed using quantiles of the distribution of the cumulative proportion in favor of the treatment
#' under the null hypothesis. 
#' It thus may not be valid if this hypothesis is rejected.
#' 
#' WARNING: For the win ratio, the proposed implementation enables the use of thresholds, endpoints that are not time to events,
#' and the correction proposed in Peron et al. (2016) to account for censoring. 
#' These development have not been examined by Wang et al. (2016), or in other papers (at out knowledge). They are only provided here by implementation convenience.
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object.
#' 
#' @examples
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

#' @rdname BuyseRes-summary
setGeneric(name = "summary", 
           def = function(object, ...){standardGeneric("summary")}
)

## * method - summary
#' @rdname BuyseRes-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "BuyseRes",
          definition = function(object, show = TRUE, percentage = TRUE,
                                statistic = BuyseTest.options()$statistic,
                                conf.level = 0.95, alternative = "two.sided",
                                strata = if(length(object@level.strata)==1){"global"}else{NULL},                                
                                digit = c(2,3)){

              ## ** normalize and check arguments
              validLogical(show,
                           name1 = "show",
                           valid.length = 1,
                           method = "summary[BuyseRes]")
              
              validLogical(percentage,
                           name1 = "percentage",
                           valid.length = 1,
                           refuse.NA = FALSE, 
                           method = "summary[BuyseRes]")

              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netchance" = "netChance",
                                  "winratio" = "winRatio",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netChance","winRatio"),
                             valid.length = 1,
                             method = "summary[BuyseRes]")

              validNumeric(conf.level,
                           name1 = "conf.level",
                           min = 0, max = 1,
                           refuse.NA = FALSE,
                           valid.length = 1,
                           method = "summary[BuyseRes]")
              if(is.na(conf.level)){
                  method.inference <- "none"
              }
              
              validCharacter(alternative,
                             name1 = "alternative",
                             valid.values = c("two.sided","less","greater"),
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
              endpoint <- object@endpoint
              n.endpoint <- length(endpoint)
              n.strata <- length(object@level.strata)

              delta <- slot(object, name = paste0("delta.",statistic))
              Delta <- slot(object, name = paste0("Delta.",statistic))
              n.resampling <- object@n.resampling
              method.inference <- object@method.inference
              alpha <- 1-conf.level

              ## ** generate summary table
              ## *** prepare
              table <- data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=14))
              names(table) <- c("endpoint","threshold","strata",
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
              table[index.global,"strata"] <- "global"

              if(statistic=="netChance"){ ##
                  table[index.global,"delta"] <- colSums(delta)
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
              if(method.inference %in% c("permutation","stratified permutation")){
                  outCI <- calcCIpermutation(Delta = Delta,
                                             Delta.permutation = slot(object, name = paste0("DeltaResampling.",statistic)),
                                             statistic = statistic,
                                             endpoint = endpoint,
                                             alternative = alternative,
                                             alpha =  alpha)
                  table[index.global,"CIinf.Delta"] <- outCI$Delta.CI[1,]
                  table[index.global,"CIsup.Delta"] <- outCI$Delta.CI[2,]
                  table[index.global,"p.value"] <- outCI$Delta.pvalue
                  table[index.global,"n.resampling"] <- outCI$n.resampling_real
                  
              }else if(method.inference %in% c("bootstrap","stratified bootstrap") ){
              }else if(method.inference == "asymptotic"){
              }
            
              ## ** generate print table
              table.print <- table
                  
              ## *** add column with stars
              if(method.inference != "none"){
                  colStars <- sapply(table.print[index.global,"p.value"],function(x){
                      if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
                  })

                  if(method.inference != "asymptotic"){
                      table.print <- cbind(table.print[,setdiff(names(table.print), "n.resampling")],
                                           "significance" = colStars)
                  }else{
                      table.print <- cbind(table.print[,setdiff(names(table.print), "n.resampling")],
                                           "significance" = colStars,
                                           table.print[,"n.resampling"])
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
                  keep.cols <- setdiff(names(table.print), c("CIinf.Delta","CIsup.Delta","n.resampling","p.value"))
                  table.print <- table.print[,keep.cols, drop = FALSE]
              }else if(method.inference == "asymptotic"){
                  keep.cols <- setdiff(names(table.print), "n.resampling")
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
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digit = digit[1])
              }
              
              if(!is.na(digit[2])){
                  param.signif <- c("delta","Delta")
                  if(method.inference != "none"){
                     param.signif <- c(param.signif, "CIinf.Delta","CIsup.Delta")
                  }
                  table.print[,param.signif] <- sapply(table.print[,param.signif], round, digit = digit[2])
              }

              ## *** set names
              if(identical(percentage,TRUE)){
                  oldnames <- c("n.favorable","n.unfavorable","n.neutral","n.uninf","n.total")
                  newnames <- c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf","pc.total")
                  names(table)[match(oldnames,names(table))] <- newnames
                  names(table.print)[match(oldnames,names(table.print))] <- newnames
              }

              ## *** convert NA to ""
              table.print[is.na(table.print)] <- ""

              ## *** remove name significance
              if(method.inference != "none"){
                  colnames(table.print)[colnames(table.print)=="significance"] <- ""
              }
                  
              ## *** merge CI inf and CI sup column
              if(method.inference != "none"){
                  if("strata" %in% names(table.print)){
                      index.tempo <- which(table.print$strata == "global")
                  }else{
                      index.tempo <- 1:NROW(table.print)
                  }
                  table.print$CIinf.Delta[index.tempo] <- paste0("[",table.print[index.tempo,"CIinf.Delta"],
                                                                 ";",table.print[index.tempo,"CIsup.Delta"],"]")
                  qInf <- round(100*alpha/2,digit[2])
                  qSup <- round(100*(1-alpha/2),digit[2])

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
              table.print[which(test.duplicated),c("endpoint","threshold")] <- ""
              
              ## ** display
              if(show){
                  ## *** additional text
                  txt.strata <- if(n.strata>1){paste0(" and ",n.strata," strata")}else{""}
                  txt.endpoint <- paste0("with ",n.endpoint," prioritized endpoint")
                  if(n.endpoint>1){txt.endpoint <- paste0(txt.endpoint,"s")}
                  
                  ## *** display                  
                  cat("        Generalized pairwise comparison ",txt.endpoint,txt.strata,"\n\n", sep = "")
                  if(statistic == "winRatio"){
                      cat(" > statistic       : win ratio (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 1 \n", sep = "")
                  }else {
                      cat(" > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 0 \n", sep = "")
                  }
                  if(method.inference %in% c("permutation","bootstrap", "stratified permutation", "stratified bootstrap")){
                      ok.permutation <- all(n.resampling[1]==n.resampling)
                      if(ok.permutation){
                          txt.permutation <- n.resampling[1]
                          table.print$n.resampling <- NULL
                      }else{
                          txt.permutation <- paste0("[",min(n.resampling)," ; ",max(n.resampling),"]")
                      }
                      txt.method <- switch(method.inference,
                                           "permutation" = "permutation test",
                                           "stratified permutation" = "stratified permutation test",
                                           "bootstrap" = "bootstrap resampling",
                                           "stratified bootstrap" = "stratified bootstrap resampling"
                                           )
                      cat(" > ",txt.method,": ",txt.permutation," samples, confidence level ",1-alpha," \n", sep = "")
                  }else if(method.inference == "asymptotic"){
                      cat(" >  test: ",txt.permutation," samples, confidence level ",1-alpha," \n", sep = "")
                  }
                  
                  cat(" > treatment groups: ",object@level.treatment[1]," (control) vs. ",object@level.treatment[2]," (treatment) \n", sep = "")
                  if(any(object@type == "TimeToEvent")){
                      txt.method.tte <- switch(object@method.tte,
                                               "Gehan" = "uninformative pairs \n",
                                               "Peto" = "imputation using Kaplan Meier \n",
                                               "Efron" = "imputation using Kaplan Meier stratified by treatment group \n",
                                               "Peron" = "imputation using Kaplan Meier stratified by treatment group \n"
                                               ) 

                      cat(" > censored pairs  : ",txt.method.tte,"\n", sep = "")
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



