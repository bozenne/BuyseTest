## * Documentation - summary
#' @docType methods
#' @name BuyseRes-summary
#' @title Summary Method for Class "BuyseRes"
#' @aliases summary summary,BuyseRes
#' @include BuyseRes-object.R
#' 
#' @description Summarize the results from the \code{\link{BuyseTest}} function.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param show Should the table be displayed? \emph{logical}.
#' @param percentage Should the percentage of pairs of each type be displayed ? Otherwise the number of pairs is displayed.
#' @param statistic the statistic summarizing the pairwise comparison: \code{"netChance"} displays the net chance in favor of treatment, as described in Buyse (2010) and Peron et al. (2016)), whereas \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016). 
#' @param strata the name of the strata to be displayed \emph{character vector}. Can also be \code{"global"} which displays the overall results.
#' @param digit the number of digit to use for printing the results. \emph{integer}. 
#' @param ... arguments to be passed from the generic method to the class specific method [not relevant to the user]
#' 
#' @details 
#' WARNING : the confidence interval is computed using quantiles of the distribution of the cumulative proportion in favor of the treatment under the null hypothesis. 
#' It thus may not be valid if this hypothesis is rejected. In particular, if the cumulative proportion in favor of the treatment is close to 1, the upper limit of the confidence interval may exceed 1. 
#' 
#' WARNING: For the win ratio, the proposed implementation enables the use of thresholds, endpoints that are not time to events, and the correction proposed in Peron et al. (2016) to account for censoring. 
#' These development have not been examined by Wang et al. (2016), or in other papers (at out knowledge). They are only provided here by implementation convenience.
#' 
#' @return 
#'   A matrix containing the endpoint (and the strata) in rows and the results of the pairwise comparison in columns:
#'     \itemize{
#'       \item\code{"n.favorable","n.unfavorable","n.neutral","n.uninformative"} : The favorable, unfavorable, neutral and uninformative pairs are reported in number of pairs classified in each category. \cr 
#'       \item\code{"delta"} indicates for each endpoint and each strata the summary statistic of the pairwise comparison. \cr
#'       \item\code{"Delta"} indicates for each endpoint the summary statistic of the pairwise comparison over all strata. \cr
#'       \item\code{"CIinf.Delta","CIsup.Delta"} the confidence interval of the summary statistic. \cr
#'       \item\code{"p.value"} the p.value relative to the null hypothesis. 
#'     }
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
#'  BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status) + Bin(toxicity), data=dt, n.permutation = 10, trace = 0)
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
                                strata = if(length(object@strata)==1){"global"}else{NULL},                                
                                digit = c(2,3)){

              ## make case insensitive and remove space
              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netchance" = "netChance",
                                  "winratio" = "winRatio",
                                  statistic)

### ** preparation
              validLogical(show,
                           name1 = "show",
                           valid.length = 1,
                           method = "summary[BuyseRes]")
              validLogical(percentage,
                           name1 = "percentage",
                           valid.length = 1,
                           refuse.NA = FALSE, 
                           method = "summary[BuyseRes]")
              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netChance","winRatio"),
                             valid.length = 1,
                             method = "summary[BuyseRes]")
            
### ** mise en forme
              n.endpoint <- length(object@endpoint)
              n.strata <- length(object@strata)
            
              delta <- object@delta[[statistic]]
              Delta <- object@Delta[[statistic]]
              Delta_quantile <- object@Delta_quantile[[statistic]]
              n_permutation <- object@n_permutation[[statistic]]
              p.value <- object@p.value[[statistic]]
              alpha <- 1-object@conf.level
                  
              table <- data.frame(matrix(NA,nrow=(n.strata+1)*n.endpoint,ncol=15))
              names(table) <- c("endpoint","threshold","strata","n.total","n.favorable","n.unfavorable","n.neutral","n.uninf","delta","Delta","CIinf.Delta","CIsup.Delta","n.permutation","p.value","")
              test.permutation <- !all(is.na(n_permutation))
            
              ### ** fill
              index.global <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1
            
              table[index.global,"n.favorable"] <- colSums(object@count_favorable)
              table[index.global,"n.unfavorable"] <- colSums(object@count_unfavorable)
              table[index.global,"n.neutral"] <- colSums(object@count_neutral)
              table[index.global,"n.uninf"] <- colSums(object@count_uninf)
              table[index.global,"n.total"] <- rowSums(table[index.global,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
            
              table[index.global,"endpoint"] <- object@endpoint
              table[index.global,"threshold"] <- object@threshold
              table[index.global,"strata"] <- "global"
              if(statistic=="netChance"){
                  table[index.global,"delta"] <- colSums(delta)
              }else{
                  table[index.global,"delta"] <- colSums(object@count_favorable)/colSums(object@count_unfavorable)
              }
              table[index.global,"Delta"] <- Delta
              table[index.global,"CIinf.Delta"] <- Delta + Delta_quantile[1,]
              table[index.global,"CIsup.Delta"] <- Delta + Delta_quantile[2,]
              table[index.global,"n.permutation"] <- n_permutation
              table[index.global,"p.value"] <- p.value
              table[index.global,ncol(table)] <- sapply(p.value,function(x){
                  if(is.na(x)){NA}else if(x<0.001){"***"}else if(x<0.01){"**"}else if(x<0.05){"*"}else if(x<0.1){"."}else{""}
              })
            
            for(iter_strata in 1:n.strata){
              index.strata <- seq(0,n.endpoint-1,by=1)*(n.strata+1)+1+iter_strata
              
              table[index.strata,"n.favorable"] <- object@count_favorable[iter_strata,]
              table[index.strata,"n.unfavorable"] <- object@count_unfavorable[iter_strata,]
              table[index.strata,"n.neutral"] <- object@count_neutral[iter_strata,]
              table[index.strata,"n.uninf"] <- object@count_uninf[iter_strata,]
              table[index.strata,"n.total"] <- rowSums(table[index.strata,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])
              
              table[index.strata,"strata"] <- object@strata[iter_strata]
              table[index.strata,"endpoint"] <- object@endpoint
              table[index.strata,"threshold"] <- object@threshold
              table[index.strata,"delta"] <- delta[iter_strata,]
            }
             
            ### ** posttreatment
            
              ## *** normalization of the counts
              if(is.na(percentage)){
                  table$n.favorable <- NULL
                  table$n.unfavorable <- NULL
                  table$n.neutral <- NULL
                  table$n.uninf <- NULL
                  table$n.total <- NULL
              }else if(percentage == TRUE){
                  table[,"n.favorable"] <- 100*table[,"n.favorable"]/table[1,"n.total"]
                  table[,"n.unfavorable"] <- 100*table[,"n.unfavorable"]/table[1,"n.total"]
                  table[,"n.neutral"] <- 100*table[,"n.neutral"]/table[1,"n.total"]
                  table[,"n.uninf"] <- 100*table[,"n.uninf"]/table[1,"n.total"]
                  table[,"n.total"] <- 100*table[,"n.total"]/table[1,"n.total"]
              }
            
              ## *** Permutation test
              if(test.permutation == FALSE){
                  keep.cols <- setdiff(names(table), c("CIinf.Delta","CIsup.Delta","n.permutation","p.value",""))
                  table <- table[,keep.cols, drop = FALSE]
              }
             
              ## *** strata
              if(!is.null(strata)){
                  table <- table[table$strata %in% strata,,drop = FALSE]
                  if(identical(strata, "global")){
                      keep.cols <- which(names(table) %in% setdiff(names(table), "strata"))
                      table <- table[,keep.cols,drop = FALSE]
                  }
              }
            
              ## *** rounding
              if(length(digit) == 1){digit <- rep(digit,2)}
            
              if(!is.na(percentage) && !is.na(digit[1]) && digit[1]>=0){
                  param.signif <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
                  table[,param.signif] <- sapply(table[,param.signif],round,digit=digit[1])
              }
              
              if(!is.na(digit[2]) && digit[2]>=0){
                  param.signif <- c("delta","Delta")
                  if(test.permutation){
                      param.signif <- c(param.signif, "CIinf.Delta","CIsup.Delta")
                  }
                  table[,param.signif] <- sapply(table[,param.signif],signif,digit=digit[2])
              }
            
              if(identical(percentage,TRUE)){
                  oldnames <- c("n.favorable","n.unfavorable","n.neutral","n.uninf","n.total")
                  newnames <- c("pc.favorable","pc.unfavorable","pc.neutral","pc.uninf","pc.total")
                  names(table)[match(oldnames,names(table))] <- newnames
              }
            
### ** display
              if(show){
                  table.print <- table
                  ## clean NA
                  emptyCol <- which(names(table.print) %in% c("Delta","CIinf.Delta","CIsup.Delta","n.permutation","p.value",""))
                  for(iterCol in emptyCol){
                      table.print[is.na(table.print[,iterCol]), iterCol] <- ""
                  }
                  table.print$threshold[is.na(table.print$threshold)] <- ""

                  ## merge CI inf and CI sup column
                  if(test.permutation){
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

                  ## simplify names
                  if(identical(percentage,TRUE)){
                      names(table.print) <- gsub("^pc.","",names(table.print))
                  }else if(identical(percentage,FALSE)){
                      names(table.print) <- gsub("^n.","",names(table.print))
                  }
                  
                  ## remove duplicated values in endpoint/threshold
                  test.duplicated <- duplicated(interaction(table.print$endpoint,table.print$threshold))
                  table.print[which(test.duplicated),c("endpoint","threshold")] <- ""

                  ## additional text
                  txt.strata <- if(n.strata>1){paste0(" and ",n.strata," strata")}else{""}
                  txt.endpoint <- paste0("with ",n.endpoint," prioritized endpoint")
                  if(n.endpoint>1){txt.endpoint <- paste0(txt.endpoint,"s")}
                  
                  ## display                  
                  cat("        Generalized pairwise comparison ",txt.endpoint,txt.strata,"\n\n", sep = "")
                  if(statistic == "winRatio"){
                      cat(" > statistic       : win ratio (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 1 \n", sep = "")
                  }else {
                      cat(" > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) \n",
                          " > null hypothesis : Delta == 0 \n", sep = "")
                  }
                  if(test.permutation){
                      ok.permutation <- all(n_permutation[1]==n_permutation)
                      if(ok.permutation){
                          txt.permutation <- n_permutation[1]
                          table.print$n.permutation <- NULL
                      }else{
                          txt.permutation <- paste0("[",min(n_permutation)," ; ",max(n_permutation),"]")
                      }
                      cat(" > permutation test: ",txt.permutation," samples, confidence level ",1-alpha," \n", sep = "")
                  }
                  cat(" > groups          : ",object@levels.treatment[1],"(control) vs. ",object@levels.treatment[2],"(treatment) \n", sep = "")
                  cat(" > results\n")
                  print(table.print, row.names = FALSE)         
              }
            
              ### ** export
              return(invisible(table))
            
          }
)



