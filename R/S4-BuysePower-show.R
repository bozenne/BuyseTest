## * Documentation - show
#' @docType methods
#' @name S4BuysePower-show
#' @title Show Method for Class "S4BuysePower"
#' @aliases show,S4BuysePower-method
#' @include S4-BuysePower.R S4-BuysePower-summary.R
#' 
#' @description Display the main results stored in a \code{S4BuysePower} object.
#' 
#' @param object an \R object of class \code{S4BuysePower}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuysePower-summary}} for a more detailed presentation of the \code{S4BuysePower} object.
#'  
#' @keywords summary S4BuysePower-method
#' @author Brice Ozenne

## * Method - show
#' @rdname S4BuysePower-show
#' @exportMethod show
setMethod(f = "show",
          signature = "S4BuysePower",
          definition = function(object){
              outSummary <- as.data.frame(summary(object, print = FALSE,
                                                  statistic = unique(object@results$statistic)[1], 
                                                  col.rep = FALSE, legend = FALSE))
              outSummary$statistic <- NULL
              outSummary$rep.estimate <- NULL
              outSummary$rep.se <- NULL
              if(object@method.inference == "none"){     
                  col.value <- c("mean.estimate")
              }else{
                  col.value <- c("mean.estimate","sd.estimate","mean.se","rejection.rate")
              }
              outSummary[,col.value] <- round(outSummary[,col.value], digits = 3)
              outSummary[duplicated(outSummary[,c("endpoint","threshold")]),c("endpoint","threshold")] <- as.character(NA)
              outSummary[] <- lapply(outSummary, as.character)
              outSummary[is.na(outSummary)] <- ""
              print(outSummary, row.names = FALSE, quote = FALSE)

              return(invisible(NULL))
          }
          )
