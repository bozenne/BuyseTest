## * Documentation BuyseSim
#' @name BuyseSim-class
#' @title Class "BuyseSim" (output of BuyseTest)
#' @aliases  BuyseSim BuyseSim-class
#' 
#' @description A \code{\link{powerBuyseTest}} output is reported in a \code{BuyseSim} object.
#' 
#' @seealso 
#' \code{\link{powerBuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{BuyseSim-summary}} for the summary of the BuyseTest function results
#' 
#' @keywords classes BuyseSim-class
#' 

## * Class BuyseSim
#' @rdname BuyseSim-class
#' @exportClass BuyseSim
setClass(
  
  Class = "BuyseSim",
  
  representation(
      alternative = "character",
      conf.level = "numeric",
      n.rep = "numeric",
      results = "data.table",
      transformation = "logical"
      )

)

## * Initialize BuyseSim objects
methods::setMethod(
             f = "initialize", 
             signature = "BuyseSim", 
             definition = function(.Object,
                                   alternative,
                                   conf.level,
                                   n.rep,
                                   results,
                                   transformation){

                 .Object@alternative <- alternative
                 .Object@conf.level <- conf.level
                 .Object@n.rep <- n.rep
                 .Object@results <- results
                 .Object@transformation <- transformation
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor BuyseSim objects
BuyseSim <- function(...) new("BuyseSim", ...) 
