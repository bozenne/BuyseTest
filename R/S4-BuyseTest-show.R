## * Documentation - show
#' @docType methods
#' @name S4BuyseTest-show
#' @title Show Method for Class "S4BuyseTest"
#' @aliases show,S4BuyseTest-method
#' @include S4-BuyseTest.R S4-BuyseTest-summary.R
#' 
#' @description Display the main results stored in a \code{S4BuyseTest} object.
#' 
#' @param object an \R object of class \code{S4BuyseTest}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#'  
#' @keywords print
#' @author Brice Ozenne

## * Method - show
#' @rdname S4BuyseTest-show
#' @exportMethod show
setMethod(f = "show",
          signature = "S4BuyseTest",
          definition = function(object){

              ## compute summary statistics
              outSummary <- summary(object, print = FALSE, strata = "global")

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
