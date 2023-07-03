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
#' @keywords classes
#' @author Brice Ozenne 

## * Class S4BuysePower
#' @rdname S4BuysePower-class
#' @exportClass S4BuysePower
setClass(
  
  Class = "S4BuysePower",
  
  representation(
      args = "list",
      endpoint = "character",
      results = "data.table",
      sample.size = "matrix",
      seed = "numeric"
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
                                   power,
                                   n.rep,
                                   results,
                                   threshold,
                                   restriction,
                                   type,
                                   max.sample.size,
                                   sample.sizeC,
                                   sample.sizeT,
                                   seed){

                 ## ** store
                 .Object@args <- list(alternative = alternative,
                                      conf.level = conf.level,
                                      method.inference = method.inference,
                                      n.rep = n.rep,
                                      null = null,
                                      restriction = restriction,
                                      threshold = threshold,
                                      type = type
                                      )
                 
                 .Object@endpoint <- stats::setNames(endpoint,paste0(endpoint,ifelse(!is.na(restriction),paste0("_r",restriction),""),ifelse(threshold>1e-12,paste0("_t",threshold),"")))
                 .Object@results <- results
                 .Object@sample.size <- cbind("C" = sample.sizeC, "T" = sample.sizeT)
                 if(!is.null(power)){
                     .Object@args$power <- power
                     .Object@args$max.sample.size <- max.sample.size
                     attr(.Object@sample.size, "sample") <- cbind("C" = attr(sample.sizeC,"sample"), "T" = attr(sample.sizeT, "sample"))
                 }
                 if(!is.null(seed)){
                     .Object@seed <- seed
                 }

                 ## ** export
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor S4BuysePower objects
S4BuysePower <- function(...) new("S4BuysePower", ...) 
