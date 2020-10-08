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
#' @param print [logical] Should the table be displayed?.
#' @param statistic [character] statistic relative to which the power should be computed:
#' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
#' \code{"mannWhitney"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
#' Default value read from \code{BuyseTest.options()}.
#' @param endpoint [character vector] the endpoints to be displayed: must be the name of the endpoint followed by an underscore and then by the threshold.
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#' @param digit [integer vector] the number of digit to use for printing the counts and the delta.  
#' @param legend [logical] should explainations about the content of each column be displayed? 
#' @param col.rep [logical] should the number of successful simulations be displayed? 
#'
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing a simulation study for generalized pairwise comparison. \cr
#'
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937 \
#' 
#' @keywords summary S4BuysePower-method
#' @author Brice Ozenne

## * method - summary
#' @rdname S4BuysePower-summary
#' @exportMethod summary
setMethod(f = "summary",
          signature = "S4BuysePower",
          definition = function(object, print = TRUE,
                                statistic = NULL, endpoint = NULL, order.Hprojection = NULL, transformation = NULL,
                                legend = TRUE, col.rep = FALSE, digit = 4){

              dt.res <- slot(object, name = "results")
              alpha <- 1-slot(object, name = "conf.level")
              null <- slot(object, name = "null")

              ## ** normalize and check arguments
              valid.endpoint <- paste0(object@endpoint,"_",object@threshold)
              valid.statistic <- unique(dt.res$statistic)
              valid.order <- unique(dt.res$order)
              valid.transformation <- unique(dt.res$transformation)
              option <- BuyseTest.options()
              if(is.null(statistic)){
                  statistic <- unique(dt.res$statistic)
              }
              if(is.null(endpoint)){
                  endpoint <- utils::tail(valid.endpoint, 1)
              }
              if(is.null(order.Hprojection)){
                  order.Hprojection <- max(dt.res$order.Hprojection)
              }
              if(is.null(transformation)){
                  transformation <- dt.res$transformation[which.max(dt.res$transformation)[1]]
              }
              validLogical(print,
                           name1 = "print",
                           valid.length = 1,
                           method = "summary[S4BuysePower]")

              statistic <- sapply(gsub("[[:blank:]]", "", tolower(statistic)),
                                  switch,
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = valid.statistic,
                             valid.length = 1:2,
                             method = "summary[S4BuysePower]")

              validCharacter(endpoint,
                             name1 = "endpoint",
                             valid.length = NULL,
                             valid.values = valid.endpoint,
                             refuse.duplicates = TRUE,
                             refuse.NULL = TRUE,
                             method = "summary[S4BuysePower]")

              validLogical(transformation,
                           name1 = "transformation",
                           valid.length = 1,
                           method = "summary[S4BuysePower]")

              validInteger(order.Hprojection,
                           name1 = "order.Hprojection",
                           valid.length = 1,
                           min = min(valid.order),
                           max = max(valid.order),
                           method = "summary[S4BuysePower]")

              ## ** subset
              index.subset <- which((dt.res$endpoint %in% endpoint) * (dt.res$order == order.Hprojection) * (dt.res$transformation == transformation) == 1)
              if(object@method.inference == "u-statistic"){                          
                  dtS.res <- dt.res[index.subset,list(rep.estimate = sum(!is.na(.SD$estimate)),
                                                      rep.se = sum(!is.na(.SD$se)),
                                                      mean.estimate = mean(.SD$estimate, na.rm = TRUE),
                                                      sd.estimate = stats::sd(.SD$estimate, na.rm = TRUE),
                                                      mean.se = mean(.SD$se, na.rm = TRUE),
                                                      rejection.rate = mean(.SD$p.value<=alpha, na.rm = TRUE)),
                                    by = c("n.T","n.C","endpoint","statistic"),]
                  col.value <- c("mean.estimate","sd.estimate","mean.se","rejection.rate","rep.estimate","rep.se")
              }else{
                  dtS.res <- dt.res[index.subset,list(rep.estimate = sum(!is.na(.SD$estimate)),
                                                      mean.estimate = mean(.SD$estimate, na.rm = TRUE)),
                                    by = c("n.T","n.C","endpoint","statistic"),]
                  col.value <- c("mean.estimate","rep.estimate")
              }
              index.endpoint <- match(dtS.res$endpoint, valid.endpoint)
              dtS.res$endpoint <- object@endpoint[index.endpoint]
              dtS.res$threshold <- object@threshold[index.endpoint]
              if(any(object@type[index.endpoint]==1)){
                  dtS.res$threshold[object@type[index.endpoint]==1] <- NA
              }
              data.table::setkeyv(dtS.res, c("endpoint","n.T"))
              data.table::setcolorder(dtS.res, neworder = c("statistic","endpoint","threshold","n.T","n.C",col.value))

              ## ** print              
              if(print){
                  cat("        Simulation study with Generalized pairwise comparison\n", sep = "")
                  cat("        with ",object@n.rep," samples\n\n", sep = "")
                  rm.duplicate <- c("n.T", "n.C", "rep.estimate", "rep.se", "mean.estimate", "sd.estimate")

                  
                  for(iStatistic in statistic){
                      name.statistic <- switch(iStatistic,
                                               "netBenefit" = "net benefit",
                                               "winRatio" = "win ratio",
                                               "favorable" = "proportion in favor of treatment",
                                               "unfavorable" = "proportion in favor of control"
                                               )
                      cat(" - statistic   : ",name.statistic," (null hypothesis Delta=",null[statistic],")\n", sep = "")

                      df.print <- as.data.frame(dtS.res[dtS.res$statistic == iStatistic])
                      df.print$statistic <- NULL
                      df.print[,col.value] <- round(df.print[,col.value], digits = digit)
                      if(col.rep == FALSE){
                          df.print$rep.estimate <- NULL
                          df.print$rep.se <- NULL
                      }
                      df.print[duplicated(df.print[,c("endpoint","threshold")]),c("endpoint","threshold")] <- as.character(NA)
                      df.print[] <- lapply(df.print, as.character)
                      df.print[is.na(df.print)] <- ""
                      print(df.print, row.names = FALSE, quote = FALSE)
                      cat("\n")
                  }
                  
                  if(legend){
                      M <- rbind(c(" n.T",":","number of observations in the treatment group"),
                                 c(" n.C",":","number of observations in the control group"),
                                 c(" mean.estimate",":","average estimate over simulations"),
                                 c(" sd.estimate",":","standard deviation of the estimate over simulations"))
                      if(object@method.inference == "u-statistic"){                          
                          M <- rbind(M,
                                     c(" mean.se",":","average estimated standard error of the estimate over simulations"),
                                     c(" rejection",":","frequency of the rejection of the null hypothesis over simulations")
                                     )
                          txt.note <- paste0("(standard error: H-projection of order ",order.Hprojection,"| p-value:")
                          if(transformation){
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
              return(invisible(dtS.res))
            
          }
)
