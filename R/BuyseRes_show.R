## * Documentation - show
#' @docType methods
#' @name BuyseRes-show
#' @title Show Method for Class "BuyseRes"
#' @aliases show show,BuyseRes-method
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
            
              table <- summary(object, show = FALSE, percentage = NA, strata="global")              
              table$threshold[is.na(table$threshold)] <- ""
              if("p.value" %in% names(table)){
                  table$CIinf.Delta <- paste0("[",table$CIinf.Delta,
                                              ";",table$CIsup.Delta,"]")
                  names(table)[names(table) == "CIinf.Delta"] <- "CI"
                  table$CIsup.Delta <- NULL
                  table$n.resampling <- NULL
              }
              
              print(table, row.names = FALSE)
           
              return(invisible(NULL))
          }
          )
