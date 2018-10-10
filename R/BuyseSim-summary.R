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
#' @param statistic [character] the statistic summarizing the pairwise comparison:
#' \code{"netBenefit"} displays the net benefit
#' whereas \code{"winRatio"} displays the win ratio.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param legend [logical] should explainations about the content of each column be displayed? 
#' @param ... arguments to be passed from the generic method to the class specific method [not relevant to the user]
#'
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing a simulation study for generalized pairwise comparison. \cr
#' 
#' @keywords summary BuyseSim-method

## * method - summary
#' @rdname BuyseSim-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "BuyseSim",
          definition = function(object, print = TRUE, statistic = NULL, legend = TRUE,
                                digit = 4){

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
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio"),
                             valid.length = 1:2,
                             method = "summary[BuyseSim]")

              ## ** extract information
              out <- setNames(vector(mode = "list", length = length(statistic)), statistic)
              for(iStatistic in statistic){
                  statistic.cols <- paste0(iStatistic, c("",".se",".p.value"))
                  dt.res <- slot(object, name = "results")
                  alpha <- 1-slot(object, name = "conf.level")

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
                                              by = c("n.T","n.C"), .SDcols = statistic.cols]
              }
              
              ## ** print              
              if(print){
                  cat("        Simulation study with Generalized pairwise comparison\n\n", sep = "")
                 
                  if("netBenefit" %in% statistic){
                      cat(" > statistic   : net benefit\n")
                      print(out$netBenefit, digits = digit)
                      cat("\n")
                  }
                  
                  if("winRatio" %in% statistic){
                      cat(" > statistic   : win ratio\n")
                      print(out$winRatio, digits = digit)
                      cat("\n")
                  }

                  if(legend){
                  M <- rbind(c(" n.T",":","number of observations in the treatment group"),
                             c(" n.C",":","number of observations in the control group"),
                             c(" rep.estimate",":","number of sucessful simulations for the point estimation"),
                             c(" rep.se",":","number of sucessful simulations for the estimation of the distribution of the estimate"),
                             c(" mean.estimate",":","average estimate over simulations"),
                             c(" sd.estimate",":","standard deviation of the estimate over simulations"),
                             c(" mean.se",":","average estimated standard error of the estimate over simulations"),
                             c(" rejection.rate",":","frequency of the rejection of the null hypothesis over simulations")
                             )

                  nchar.1 <- sapply(M[,1],nchar)
                  M[,1] <- paste0(M[,1],
                                  sapply(max(nchar.1) - nchar.1, function(iX){paste0(rep(" ",time = iX),collapse = "")}))
                  txt.legend <- apply(M, 1, function(iRow){paste(iRow[1],iRow[2]," ",iRow[3],"\n",sep = "")})
                  cat(txt.legend,sep ="")
                  cat("\n")
                  }
              }
              
              ## ** export
              return(invisible(out))
            
          }
)



