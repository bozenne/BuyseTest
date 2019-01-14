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
      method.inference = "character",
      conf.level = "numeric",
      null = "numeric",
      n.rep = "numeric",
      results = "data.table"
      )

)

## * Initialize BuyseSim objects
methods::setMethod(
             f = "initialize", 
             signature = "BuyseSim", 
             definition = function(.Object,
                                   alternative,
                                   method.inference,
                                   conf.level,
                                   null,
                                   n.rep,
                                   results){

                 .Object@alternative <- alternative
                 .Object@method.inference <- method.inference
                 .Object@conf.level <- conf.level
                 .Object@null <- null
                 .Object@n.rep <- n.rep
                 .Object@results <- results
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor BuyseSim objects
BuyseSim <- function(...) new("BuyseSim", ...) 
