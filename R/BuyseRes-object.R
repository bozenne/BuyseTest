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
#' @author Brice Ozenne

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
      delta.netBenefit = "matrix",
      delta.winRatio = "matrix",
      Delta.netBenefit = "vector",
      Delta.winRatio = "vector",
      type = "vector",
      endpoint = "vector",
      level.treatment = "vector",
      level.strata = "vector",
      scoring.rule = "character",
      hierarchical = "logical",
      neutral.as.uninf = "logical",
      correction.uninf = "numeric",
      method.inference = "character",
      strata = "vector",
      threshold = "numeric",
      n.resampling = "numeric",
      deltaResampling.netBenefit = "array",
      deltaResampling.winRatio = "array",
      DeltaResampling.netBenefit = "matrix",
      DeltaResampling.winRatio = "matrix",
      covariance = "matrix",
      covarianceResampling = "array",
      weight = "numeric",
      iidAverage = "list",
      iidNuisance = "list",
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
                                   delta.netBenefit,
                                   delta.winRatio,
                                   Delta.netBenefit,
                                   Delta.winRatio,
                                   type,
                                   endpoint,
                                   level.strata,
                                   level.treatment,
                                   scoring.rule,
                                   hierarchical,
                                   neutral.as.uninf,
                                   correction.uninf,
                                   method.inference,
                                   strata,
                                   threshold,
                                   n.resampling,
                                   deltaResampling.netBenefit,
                                   deltaResampling.winRatio,
                                   DeltaResampling.netBenefit,
                                   DeltaResampling.winRatio,
                                   covariance,
                                   covarianceResampling,
                                   weight,
                                   iidAverage_favorable,
                                   iidAverage_unfavorable,
                                   iidNuisance_favorable,
                                   iidNuisance_unfavorable,
                                   tablePairScore,
                                   tableSurvival,
                                   args){

                 name.endpoint <- paste0(endpoint,"_",threshold)
                 
                 ## ** count
                 dimnames(count.favorable) <- list(level.strata, name.endpoint)
                 dimnames(count.unfavorable) <- list(level.strata, name.endpoint)
                 dimnames(count.neutral) <- list(level.strata, name.endpoint)
                 dimnames(count.uninf) <- list(level.strata, name.endpoint)

                 names(weight) <- name.endpoint
                 
                 ## ** delta/Delta
                 dimnames(delta.netBenefit) <- list(level.strata, name.endpoint)
                 dimnames(delta.winRatio) <- list(level.strata, name.endpoint)
                 names(Delta.netBenefit) <- name.endpoint
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
                 
                 .Object@delta.netBenefit <- delta.netBenefit
                 .Object@delta.winRatio <- delta.winRatio
                 .Object@Delta.netBenefit <- Delta.netBenefit
                 .Object@Delta.winRatio <- Delta.winRatio

                 .Object@type <- type
                 .Object@endpoint <- endpoint
                 .Object@level.strata <- level.strata
                 .Object@level.treatment <- level.treatment
                 .Object@scoring.rule <- scoring.rule
                 .Object@hierarchical <- hierarchical
                 .Object@neutral.as.uninf <- neutral.as.uninf
                 .Object@correction.uninf <- correction.uninf
                 .Object@method.inference <- method.inference
                 .Object@strata <- strata
                 .Object@threshold <- threshold
                 
                 .Object@n.resampling <- n.resampling

                 .Object@deltaResampling.netBenefit <- deltaResampling.netBenefit
                 .Object@deltaResampling.winRatio <- deltaResampling.winRatio
                 .Object@DeltaResampling.netBenefit <- DeltaResampling.netBenefit
                 .Object@DeltaResampling.winRatio <- DeltaResampling.winRatio

                 .Object@covariance <- covariance
                 if(NCOL(covariance)>0){
                     dimnames(.Object@covariance) <- list(name.endpoint,
                                                          c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }
                 .Object@covarianceResampling <- covarianceResampling
                 if(dim(.Object@covarianceResampling)[1]>0){
                     dimnames(.Object@covarianceResampling) <- list(NULL,name.endpoint,
                                                            c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }
                 .Object@weight <- weight
                 
                 .Object@iidAverage <- list(favorable = iidAverage_favorable,
                                            unfavorable = iidAverage_unfavorable)
                 if(!is.null(.Object@iidAverage[[1]])){
                     colnames(.Object@iidAverage[[1]]) <- name.endpoint
                 }
                 if(!is.null(.Object@iidAverage[[2]])){
                     colnames(.Object@iidAverage[[2]]) <- name.endpoint
                 }
                 .Object@iidNuisance <- list(favorable = iidNuisance_favorable,
                                             unfavorable = iidNuisance_unfavorable)
                 .Object@tablePairScore <- tablePairScore
                 .Object@tableSurvival <- tableSurvival
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor BuyseRes objects
BuyseRes <- function(...) new("BuyseRes", ...) 
