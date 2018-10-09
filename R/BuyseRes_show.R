## * Documentation - show
#' @docType methods
#' @name BuyseRes-show
#' @title Show Method for Class "BuyseRes"
#' @aliases show,BuyseRes-method
#' @include BuyseRes-object.R BuyseRes-summary.R
#' 
#' @description Display the main results stored in a \code{\link{BuyseRes}} object.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseRes-summary}} for a more detailed presentation of the \code{BuyseRes} object.
#'  
#' @keywords summary BuyseRes-method

## * Method - show
#' @rdname BuyseRes-show
#' @exportMethod show
setMethod(f = "show",
          signature = "BuyseRes",
          definition = function(object){
            outSummary <- summary(object, conf.level = NA, print = FALSE, percentage = NA, strata = "global")

            table.print <- outSummary$table.print
            exclude.col <- c("CI [NA ; NA]","p.value","","n.resampling")
            table.print <- table.print[,setdiff(names(table.print), exclude.col)]
              
            print(table.print, row.names = FALSE)
           
            return(invisible(NULL))
          }
          )
