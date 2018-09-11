## * Documentation BuyseRes
#' @name BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
#' @aliases  BuyseRes BuyseRes-class
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{BuyseRes-summary}} for the summary of the BuyseTest function results
#' 
#' @keywords classes BuyseRes-class
#' 

## * Class BuyseRes
#' @rdname BuyseRes-class
#' @exportClass BuyseRes
setClass(
  
  Class = "BuyseRes",
  
  representation(
      count.favorable = "matrix",      
      count.unfavorable = "matrix",
      count.neutral = "matrix",
      count.uninf = "matrix",
      n.pairs = "numeric",
      delta.netChance = "matrix",
      delta.winRatio = "matrix",
      Delta.netChance = "vector",
      Delta.winRatio = "vector",
      index.neutralT = "vector",
      index.neutralC = "vector",
      index.uninfT = "vector",
      index.uninfC = "vector",
      type = "vector",
      endpoint = "vector",
      level.treatment = "vector",
      level.strata = "vector",
      method.tte = "data.frame",
      method.inference = "character",
      strata = "vector",
      threshold = "numeric",
      n.resampling = "numeric",
      deltaResampling.netChance = "array",
      deltaResampling.winRatio = "array",
      DeltaResampling.netChance = "matrix",
      DeltaResampling.winRatio = "matrix",
      covariance = "matrix",
      tablePairScore = "list",
      tableSurvival = "list"
      )

)

## * Initialize BuyseRes objects
methods::setMethod(
             f = "initialize", 
             signature = "BuyseRes", 
             definition = function(.Object, 
                                   count.favorable,
                                   count.unfavorable,
                                   count.neutral,
                                   count.uninf,
                                   n.pairs,
                                   delta.netChance,
                                   delta.winRatio,
                                   Delta.netChance,
                                   Delta.winRatio,
                                   index.neutralT,
                                   index.neutralC,
                                   index.uninfT,
                                   index.uninfC,
                                   type,
                                   endpoint,
                                   level.strata,
                                   level.treatment,
                                   method.tte,
                                   method.inference,
                                   strata,
                                   threshold,
                                   n.resampling,
                                   deltaResampling.netChance,
                                   deltaResampling.winRatio,
                                   DeltaResampling.netChance,
                                   DeltaResampling.winRatio,
                                   covariance,
                                   tablePairScore,
                                   tableSurvival,
                                   args){

                 name.endpoint <- paste0(endpoint,"_",threshold)
                 
                 ## ** count
                 dimnames(count.favorable) <- list(level.strata, name.endpoint)
                 dimnames(count.unfavorable) <- list(level.strata, name.endpoint)
                 dimnames(count.neutral) <- list(level.strata, name.endpoint)
                 dimnames(count.uninf) <- list(level.strata, name.endpoint)

                 ## ** delta/Delta
                 dimnames(delta.netChance) <- list(level.strata, name.endpoint)
                 dimnames(delta.winRatio) <- list(level.strata, name.endpoint)
                 names(Delta.netChance) <- name.endpoint
                 names(Delta.winRatio) <- name.endpoint
                 
                 ## ** endpoint
                 D <- length(endpoint)
                 
                 ## ** strata
                 if(is.null(strata)){
                     strata <- as.character(NA)
                 }
                 n.strata <- length(level.strata)

                 ## ** store
                 .Object@count.favorable <- count.favorable      
                 .Object@count.unfavorable <- count.unfavorable
                 .Object@count.neutral <- count.neutral   
                 .Object@count.uninf <- count.uninf
                 .Object@n.pairs <- n.pairs
                 
                 .Object@delta.netChance <- delta.netChance
                 .Object@delta.winRatio <- delta.winRatio
                 .Object@Delta.netChance <- Delta.netChance
                 .Object@Delta.winRatio <- Delta.winRatio

                 .Object@index.neutralT <- index.neutralT
                 .Object@index.neutralC <- index.neutralC
                 .Object@index.uninfT <- index.uninfT
                 .Object@index.uninfC <- index.uninfC
                 
                 .Object@type <- type
                 .Object@endpoint <- endpoint
                 .Object@level.strata <- level.strata
                 .Object@level.treatment <- level.treatment
                 .Object@method.tte <- method.tte
                 .Object@method.inference <- method.inference
                 .Object@strata <- strata
                 .Object@threshold <- threshold
                 
                 .Object@n.resampling <- n.resampling

                 .Object@deltaResampling.netChance <- deltaResampling.netChance
                 .Object@deltaResampling.winRatio <- deltaResampling.winRatio
                 .Object@DeltaResampling.netChance <- DeltaResampling.netChance
                 .Object@DeltaResampling.winRatio <- DeltaResampling.winRatio

                 .Object@covariance <- covariance

                 .Object@tablePairScore <- tablePairScore
                 .Object@tableSurvival <- tableSurvival
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor BuyseRes objects
BuyseRes <- function(...) new("BuyseRes", ...) 
