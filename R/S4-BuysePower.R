## * Documentation S4BuysePower
#' @name S4BuysePower-class
#' @title Class "S4BuysePower" (output of BuyseTest)
#' 
#' @description A \code{\link{powerBuyseTest}} output is reported in a \code{S4BuysePower} object.
#' 
#' @seealso 
#' \code{\link{powerBuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{S4BuysePower-summary}} for the summary of the BuyseTest function results
#' 
#' @keywords classes S4BuysePower-class
#' @author Brice Ozenne 

## * Class S4BuysePower
#' @rdname S4BuysePower-class
#' @exportClass S4BuysePower
setClass(
  
  Class = "S4BuysePower",
  
  representation(
      alternative = "character",
      method.inference = "character",
      conf.level = "numeric",
      endpoint = "character",
      null = "numeric",
      n.rep = "numeric",
      results = "data.table",
      threshold = "numeric",
      type = "character"
      )

)

## * Initialize S4BuysePower objects
methods::setMethod(
             f = "initialize", 
             signature = "S4BuysePower", 
             definition = function(.Object,
                                   alternative,
                                   method.inference,
                                   conf.level,
                                   endpoint,
                                   null,
                                   n.rep,
                                   results,
                                   threshold,
                                   type){

                 .Object@alternative <- alternative
                 .Object@method.inference <- method.inference
                 .Object@conf.level <- conf.level
                 .Object@endpoint <- endpoint
                 .Object@null <- null
                 .Object@n.rep <- n.rep
                 .Object@results <- results
                 .Object@threshold <- threshold
                 .Object@type <- type
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor S4BuysePower objects
S4BuysePower <- function(...) new("S4BuysePower", ...) 
