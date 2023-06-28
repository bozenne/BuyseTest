## * Documentation - show
#' @docType methods
#' @name S4BuysePower-show
#' @title Show Method for Class "S4BuysePower"
#' @aliases show,S4BuysePower-method
#' @include S4-BuysePower.R S4-BuysePower-print.R
#' 
#' @description Display the main results stored in a \code{S4BuysePower} object.
#' 
#' @param object an \R object of class \code{S4BuysePower}, i.e., output of \code{\link{powerBuyseTest}}
#' 
#' @seealso 
#'   \code{\link{powerBuyseTest}} for performing power calculation based on GPC. \cr
#'   \code{\link{S4BuysePower-summary}} for a more detailed presentation of the \code{S4BuysePower} object.
#'  
#' @return invisible \code{NULL}
#' @keywords print
#' 
#' @author Brice Ozenne

## * Method - show
#' @rdname S4BuyseTest-show
#' @exportMethod show
setMethod(f = "show",
          signature = "S4BuysePower",
          definition = function(object){

              print(object)

              return(invisible(NULL))
          }
          )
