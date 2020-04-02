## * Documentation - summary
#' @docType methods
#' @name BuyseSim-summary
#' @title Summary Method for Class "BuyseSim"
#' @aliases summary,BuyseSim-method
#' @include BuyseSim-object.R
#' 
#' @description Summarize the results from the \code{\link{powerBuyseTest}} function.
#' 
#' @param object output of \code{\link{powerBuyseTest}}
#' @param print [logical] Should the table be displayed?.
#' @param statistic [character] statistic relative to which the power should be computed:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"mannWhitney"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' Default value read from \code{BuyseTest.options()}.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param legend [logical] should explainations about the content of each column be displayed? 
#' @param col.rep [logical] should the number of successful simulations be displayed? 
#' @param method.inference [character vector] for which inference method the rejection rate should be displayed?
#'
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing a simulation study for generalized pairwise comparison. \cr
#'
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937 \
#' 
#' @keywords summary BuyseSim-method
#' @author Brice Ozenne

## * method - summary
#' @rdname BuyseSim-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "BuyseSim",
          definition = function(object, print = TRUE, statistic = NULL, legend = TRUE, method.inference = NULL,
                                col.rep = TRUE, digit = 4){

              dt.res <- slot(object, name = "results")
              alpha <- 1-slot(object, name = "conf.level")
              null <- slot(object, name = "null")
              
              ## ** normalize and check arguments
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- option$statistic
              }

              validLogical(print,
                           name1 = "print",
                           valid.length = 1,
                           method = "summary[BuyseSim]")

              statistic <- sapply(gsub("[[:blank:]]", "", tolower(statistic)),
                                  switch,
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                             valid.length = 1:2,
                             method = "summary[BuyseSim]")
              
              if(is.null(method.inference)){
                  method.inference <- unique(dt.res$method.inference)
              }else{
                  validCharacter(method.inference,
                                 name1 = "method.inference",
                                 valid.values = c("order=1 - transformation=FALSE","order=1 - transformation=TRUE",
                                                  "order=2 - transformation=FALSE","order=2 - transformation=TRUE"),
                                 valid.length = 1:4, refuse.NULL = FALSE,
                                 method = "summary[BuyseSim]")
              }
              
              ## ** extract information
              out <- setNames(vector(mode = "list", length = length(statistic)), statistic)
              outW <- setNames(vector(mode = "list", length = length(statistic)), statistic)
              for(iStatistic in statistic){
                  statistic.cols <- paste0(iStatistic, c("",".se",".p.value"))
                  
                  ## set infinite values to NA
                  dt.res[is.infinite(dt.res[[statistic.cols[[1]]]]), c(statistic.cols) := NA]
                  dt.res[is.infinite(dt.res[[statistic.cols[[2]]]]), c(statistic.cols[2:3]) := NA]

                  ##
                  out[[iStatistic]] <- dt.res[,list(rep.estimate = sum(!is.na(.SD[[1]])),
                                                    rep.se = sum(!is.na(.SD[[2]])),
                                                    mean.estimate = mean(.SD[[1]], na.rm = TRUE),
                                                    sd.estimate = stats::sd(.SD[[1]], na.rm = TRUE),
                                                    mean.se = mean(.SD[[2]], na.rm = TRUE),
                                                    rejection.rate = mean(.SD[[3]]<=alpha, na.rm = TRUE)),
                                              by = c("n.T","n.C","method.inference"), .SDcols = statistic.cols]
                  out[[iStatistic]] <- out[[iStatistic]][out[[iStatistic]]$method.inference %in% method.inference]
                  
                  if(object@method.inference == "u-statistic"){
                      out[[iStatistic]][, c("order") := grepl("order=2",.SD$method.inference)+1]
                      out[[iStatistic]][, c("transformation") := grepl("transformation=TRUE",.SD$method.inference)]

                      outW[[iStatistic]] <- dcast(out[[iStatistic]],
                                                  formula = n.T + n.C + rep.estimate + rep.se + mean.estimate + sd.estimate + order + mean.se ~ transformation,
                                                  value.var = "rejection.rate")
                      name.tempo <- c("n.T", "n.C", "rep.estimate", "rep.se", "mean.estimate", "sd.estimate", "order", "mean.se")
                      names(outW[[iStatistic]])[names(outW[[iStatistic]]) %in% name.tempo == FALSE] <- paste0("rejection (",setdiff(names(outW[[iStatistic]]),name.tempo),")")
                  }else{
                      outW[[iStatistic]] <- out[[iStatistic]][,.SD, .SDcols = c("n.T", "n.C", "rep.estimate", "rep.se", "mean.estimate", "sd.estimate")]
                  }
              }
              
              ## ** print              
              if(print){
                  cat("        Simulation study with Generalized pairwise comparison\n\n", sep = "")
                  rm.duplicate <- c("n.T", "n.C", "rep.estimate", "rep.se", "mean.estimate", "sd.estimate")
                  
                  if("netBenefit" %in% statistic){
                      cat(" > statistic   : net benefit (null hypothesis Delta=",null["netBenefit"],")\n", sep = "")
                      printNetBenefit <- as.data.frame(outW$netBenefit)
                      printNetBenefit <- round(printNetBenefit, digits = digit)
                      if(length(outW$netBenefit$order)>1){ ## remove duplicated values due to order = 1:2
                          printNetBenefit[printNetBenefit$order==2, rm.duplicate] <- ""
                      }
                      if(col.rep == FALSE){
                          printNetBenefit$rep.estimate <- NULL
                          printNetBenefit$rep.se <- NULL
                      }
                      print(printNetBenefit, row.names = FALSE)
                      cat("\n")
                  }
                  
                  if("winRatio" %in% statistic){
                      cat(" > statistic   : win ratio (null hypothesis Delta=",null["winRatio"],")\n", sep = "")
                      printWinRatio <- as.data.frame(outW$winRatio, stringsAsFactors = FALSE)
                      printWinRatio <- round(printWinRatio, digits = digit)
                      if(length(outW$winRatio$order)>1){ ## remove duplicated values due to order = 1:2
                          printWinRatio[printWinRatio$order==2, rm.duplicate] <- ""
                      }
                      if(col.rep == FALSE){
                          printWinRatio$rep.estimate <- NULL
                          printWinRatio$rep.se <- NULL
                      }
                      print(printWinRatio, row.names = FALSE)
                      cat("\n")
                  }
                  
                  if("favorable" %in% statistic){
                      cat(" > statistic   : proportion in favor of treatment (null hypothesis Delta=",null["favorable"],")\n", sep = "")
                      printFavorable <- as.data.frame(outW$favorable, stringsAsFactors = FALSE)
                      printFavorable <- round(printFavorable, digits = digit)
                      if(length(outW$favorable$order)>1){ ## remove duplicated values due to order = 1:2
                          printFavorable[printFavorable$order==2, rm.duplicate] <- ""
                      }
                      if(col.rep == FALSE){
                          printFavorable$rep.estimate <- NULL
                          printFavorable$rep.se <- NULL
                      }
                      print(printFavorable, row.names = FALSE)
                      cat("\n")
                  }

                  if("unfavorable" %in% statistic){
                      cat(" > statistic   : proportion in favor of control (null hypothesis Delta=",null["unfavorable"],")\n", sep = "")
                      printUnfavorable <- as.data.frame(outW$unfavorable, stringsAsFactors = FALSE)
                      printUnfavorable <- round(printUnfavorable, digits = digit)
                      if(length(outW$unfavorable$order)>1){ ## remove duplicated values due to order = 1:2
                          printUnfavorable[printUnfavorable$order==2, rm.duplicate] <- ""
                      }
                      if(col.rep == FALSE){
                          printUnfavorable$rep.estimate <- NULL
                          printUnfavorable$rep.se <- NULL
                      }
                      print(printUnfavorable, row.names = FALSE)
                      cat("\n")
                  }
                  
                  if(legend){
                      M <- rbind(c(" n.T",":","number of observations in the treatment group"),
                                 c(" n.C",":","number of observations in the control group"),
                                 c(" rep.estimate",":","number of sucessful simulations for the point estimation"),
                                 c(" rep.se",":","number of sucessful simulations for the estimation of the distribution of the estimate"),
                                 c(" mean.estimate",":","average estimate over simulations"),
                                 c(" sd.estimate",":","standard deviation of the estimate over simulations"))
                      if(object@method.inference == "u-statistic"){                          
                          M <- rbind(M,
                                     c(" order",":","order of the H-decomposition used to compute the asymptotic variance"),
                                     c(" mean.se",":","average estimated standard error of the estimate over simulations"),
                                     c(" rejection",":","frequency of the rejection of the null hypothesis over simulations"),
                                     c(" (TRUE/FALSE)",":","the parenthesis refers to whether transformation has been used to compute the p-values")
                                     )
                      }

                      nchar.1 <- sapply(M[,1],nchar)
                      M[,1] <- paste0(M[,1],
                                      sapply(max(nchar.1) - nchar.1, function(iX){paste0(rep(" ",time = iX),collapse = "")}))
                      txt.legend <- apply(M, 1, function(iRow){paste(iRow[1],iRow[2]," ",iRow[3],"\n",sep = "")})
                      cat(txt.legend,sep ="")
                      cat("\n")
                  }
              }
              
              ## ** export
              return(invisible(outW))
            
          }
)



