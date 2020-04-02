## * Documentation BuyseRes
#' @name BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
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
      delta = "array",
      Delta = "matrix",
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
      deltaResampling = "array",
      DeltaResampling = "array",
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
                                   delta,
                                   Delta,
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
                                   deltaResampling,
                                   DeltaResampling,
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
                 dimnames(delta) <- list(level.strata,
                                         name.endpoint,
                                         c("favorable","unfavorable","netBenefit","winRatio"))
                 dimnames(Delta) <- list(name.endpoint,
                                         c("favorable","unfavorable","netBenefit","winRatio"))

                 if(length(covariance)>0){
                     dimnames(covariance) <- list(name.endpoint,
                                                  c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }
                 
                 ## ** delta/Delta resampling
                 if(length(deltaResampling)>0){
                     dimnames(deltaResampling) <- list(NULL,
                                                       name.endpoint,
                                                       c("favorable","unfavorable","netBenefit","winRatio"),
                                                       level.strata)
                 }
                 
                 if(length(DeltaResampling)>0){
                     dimnames(DeltaResampling) <- list(NULL,
                                                       name.endpoint,
                                                       c("favorable","unfavorable","netBenefit","winRatio"))
                 }
                 
                 if(length(covarianceResampling)>0){
                     dimnames(covarianceResampling) <- list(NULL,
                                                            name.endpoint,
                                                            c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }
                 
                 ## ** endpoint
                 D <- length(endpoint)
                 
                 ## ** strata
                 if(is.null(strata)){
                     strata <- as.character(NA)
                 }
                 n.strata <- length(level.strata)

                 ## ** scoring rule
                 if(!is.null(attr(scoring.rule,"method.score"))){
                     attr(scoring.rule,"method.score") <- setNames(attr(scoring.rule,"method.score"), name.endpoint)
                 }
                 
                 ## ** store
                 .Object@count.favorable <- count.favorable      
                 .Object@count.unfavorable <- count.unfavorable
                 .Object@count.neutral <- count.neutral   
                 .Object@count.uninf <- count.uninf
                 .Object@n.pairs <- n.pairs
                 
                 .Object@delta <- delta
                 .Object@Delta <- Delta

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

                 .Object@deltaResampling <- deltaResampling
                 .Object@DeltaResampling <- DeltaResampling

                 .Object@covariance <- covariance

                 .Object@covarianceResampling <- covarianceResampling

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
