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
          definition = function(object, print = TRUE, statistic = NULL,
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

              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  statistic)

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio"),
                             valid.length = 1,
                             method = "summary[BuyseSim]")

              ## ** process
              statistic.cols <- paste0(statistic,c("",".se",".p.value"))
              dt.res <- slot(object, name = "results")
              alpha <- 1-slot(object, name = "conf.level")

              
              dtS.res <- dt.res[,list(repetitions = .N,
                                      mean.estimate = mean(.SD[[1]]),
                                      sd.estimate = stats::sd(.SD[[1]]),
                                      mean.se = mean(.SD[[2]]),
                                      rejection.rate = mean(.SD[[3]]<=alpha)),
                                by = c("n.T","n.C"), .SDcols = statistic.cols]

              ## ** print
              if(print){
                  cat("        Simulation study with Generalized pairwise comparison\n\n", sep = "")
                  if(statistic == "winRatio"){
                      cat(" > statistic       : win ratio\n")
                  }else {
                      cat(" > statistic       : net benefit\n")
                  }

                  print(dtS.res)
              }
              
              ## ** export
              return(invisible(dtS.res))
            
          }
)



