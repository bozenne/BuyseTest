## * Documentation - show
#' @docType methods
#' @name BuyseSim-show
#' @title Show Method for Class "BuyseSim"
#' @aliases show,BuyseSim-method
#' @include BuyseSim-object.R BuyseSim-summary.R
#' 
#' @description Display the main results stored in a \code{\link{BuyseSim}} object.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseSim}}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseSim-summary}} for a more detailed presentation of the \code{BuyseSim} object.
#'  
#' @keywords summary BuyseSim-method

## * Method - show
#' @rdname BuyseSim-show
#' @exportMethod show
setMethod(f = "show",
          signature = "BuyseSim",
          definition = function(object){
              outSummary <- summary(object, print = FALSE)

              outSummary <- lapply(outSummary, function(iT){
                  iT[,c("rep.estimate","rep.se") := NULL]
              })

              print(outSummary)
           
              return(invisible(NULL))
          }
          )
