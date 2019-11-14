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
#' @author Brice Ozenne

## * Method - show
#' @rdname BuyseRes-show
#' @exportMethod show
setMethod(f = "show",
          signature = "BuyseRes",
          definition = function(object){

              ## compute summary statistics
              outSummary <- summary(object, print = FALSE, strata = "global")$table.print

              ## only keep certain columns
              type.display <- BuyseTest.options()$print.display
              vec.tfunu <- c("total","favorable","unfavorable","neutral","uninformative")
              if(any(vec.tfunu %in% type.display)){
                  type.display[type.display %in% vec.tfunu] <- paste0(type.display[type.display %in% vec.tfunu],"(%)")
              }
              if("CI" %in% type.display){
                  type.display <- c(setdiff(type.display,"CI"),grep("^CI",names(outSummary),value=TRUE))
              }
              type.display <- intersect(names(outSummary),type.display)
              
              ## display
              table.print <- outSummary[,type.display,drop=FALSE]
              if("significance" %in% names(table.print)){
                  names(table.print)[names(table.print) == "significance"] <- ""
              }
              print(table.print, row.names = FALSE)
           
              return(invisible(NULL))
          }
          )
