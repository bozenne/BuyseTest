## * Documentation - summary
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

## * method - summary
#' @rdname S4BuysePower-summary
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
