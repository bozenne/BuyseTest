## * Documentation S4BuyseTest
#' @name S4BuyseTest-class
#' @title Class "S4BuyseTest" (output of BuyseTest)
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{S4BuyseTest} object.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{S4BuyseTest-summary}} for the summary of the BuyseTest function results
#' 
#' @keywords classes S4BuyseTest-class
#' @author Brice Ozenne

## * Class S4BuyseTest
#' @rdname S4BuyseTest-class
#' @exportClass S4BuyseTest
setClass(
  
  Class = "S4BuyseTest",
  
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

## * Initialize S4BuyseTest objects
methods::setMethod(
             f = "initialize", 
             signature = "S4BuyseTest", 
             definition = function(.Object, 
                                   count_favorable, ## from cpp object
                                   count_unfavorable, ## from cpp object
                                   count_neutral, ## from cpp object
                                   count_uninf, ## from cpp object
                                   delta, ## from cpp object
                                   Delta, ## from cpp object
                                   n_pairs, ## from cpp object
                                   iidAverage_favorable, ## from cpp object
                                   iidAverage_unfavorable, ## from cpp object
                                   iidNuisance_favorable, ## from cpp object
                                   iidNuisance_unfavorable, ## from cpp object
                                   covariance, ## from cpp object
                                   tableScore, ## from cpp object
                                   tableSurvival = NULL, ## added to the cpp object by .BuyseTest when requested by the user
                                   index.C,
                                   index.T,
                                   type,
                                   endpoint,
                                   level.strata,
                                   level.treatment,
                                   scoring.rule,
                                   hierarchical,
                                   neutral.as.uninf,
                                   correction.uninf,
                                   method.inference,
                                   method.score,
                                   strata,
                                   threshold,
                                   weight,
                                   n.resampling,
                                   deltaResampling = NULL, ## from inferenceResampling
                                   DeltaResampling = NULL, ## from inferenceResampling
                                   covarianceResampling = NULL, ## from inferenceResampling
                                   args){

                 name.endpoint <- paste0(endpoint,"_",threshold)
                 
                 ## ** count
                 dimnames(count_favorable) <- list(level.strata, name.endpoint)
                 dimnames(count_unfavorable) <- list(level.strata, name.endpoint)
                 dimnames(count_neutral) <- list(level.strata, name.endpoint)
                 dimnames(count_uninf) <- list(level.strata, name.endpoint)

                 ## ** delta/Delta
                 dimnames(delta) <- list(level.strata,
                                         name.endpoint,
                                         c("favorable","unfavorable","netBenefit","winRatio"))
                 dimnames(Delta) <- list(name.endpoint,
                                         c("favorable","unfavorable","netBenefit","winRatio"))

                 ## ** n_pairs
                 names(n_pairs) <- level.strata

                 ## ** iid and variance
                 if(!is.null(iidAverage_favorable) && NCOL(iidAverage_favorable)>0){
                     colnames(iidAverage_favorable) <- name.endpoint
                 }
                 
                 if(!is.null(iidAverage_unfavorable) && NCOL(iidAverage_unfavorable)>0){
                     colnames(iidAverage_unfavorable) <- name.endpoint
                 }

                 if(!is.null(iidNuisance_favorable) && NCOL(iidNuisance_favorable)>0){
                     colnames(iidNuisance_favorable) <- name.endpoint
                 }
                 
                 if(!is.null(iidNuisance_unfavorable) && NCOL(iidNuisance_unfavorable)>0){
                     colnames(iidNuisance_unfavorable) <- name.endpoint
                 }

                 if(!is.null(covariance) && length(covariance)>0){
                     dimnames(covariance) <- list(name.endpoint,
                                            c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }

                 ## ** tableScore
                 if(!is.null(tableScore) && length(tableScore)>0 && any(sapply(tableScore, data.table::is.data.table)==FALSE)){
                     tableScore <- pairScore2dt(tableScore,
                                                    level.treatment = level.treatment,
                                                    level.strata = level.strata,
                                                    n.strata = length(level.strata),
                                                    endpoint = endpoint,
                                                    threshold = threshold)
                 }
                 
                 ## ** tableSurvival

                 ## ** type
                 type <- stats::setNames(c("Binary","Continuous","TimeToEvent")[type], name.endpoint)

                 ## ** endpoint
                 names(endpoint) <- name.endpoint

                 ## ** level.strata
                 ## attr(outArgs$level.strata,"index") <- outArgs$index.strata
                 
                 ## ** level.treatment
                 attr(level.treatment,"indexC") <- index.C
                 attr(level.treatment,"indexT") <- index.T

                 ## ** scoring.rule
                 scoring.rule <- c("Gehan","Peron")[scoring.rule+1]
                 attr(scoring.rule,"method.score") <- stats::setNames(method.score, name.endpoint)

                 ## ** hierarchical
                 
                 ## ** neutral.as.uninf
                 
                 ## ** correction.uninf
                 
                 ## ** method.inference
                 
                 ## ** method.score

                 ## ** strata
                 if(is.null(strata)){
                     strata <- as.character(NA)
                 }
                 ## ** threshold
                 names(threshold) <- name.endpoint
                 
                 ## ** weight
                 names(weight) <- name.endpoint

                 ## ** n.resampling

                 ## ** Resampling
                 if(!is.null(deltaResampling) && length(deltaResampling)>0){
                     dimnames(deltaResampling) <- list(NULL,
                                                       name.endpoint,
                                                       c("favorable","unfavorable","netBenefit","winRatio"),
                                                       level.strata)
                 }
                 
                 if(!is.null(DeltaResampling) && length(DeltaResampling)>0){
                     dimnames(DeltaResampling) <- list(NULL,
                                                       name.endpoint,
                                                       c("favorable","unfavorable","netBenefit","winRatio"))
                 }
                 
                 if(!is.null(covarianceResampling) && length(covarianceResampling)>0){
                     dimnames(covarianceResampling) <- list(NULL,
                                                            name.endpoint,
                                                            c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }
                 
                 ## ** store
                 ## *** from c++ object
                 .Object@count.favorable <- count_favorable      
                 .Object@count.unfavorable <- count_unfavorable
                 .Object@count.neutral <- count_neutral   
                 .Object@count.uninf <- count_uninf
                 .Object@n.pairs <- n_pairs
                 .Object@delta <- delta
                 .Object@Delta <- Delta
                 .Object@iidAverage <- list(favorable = iidAverage_favorable,
                                            unfavorable = iidAverage_unfavorable)
                 .Object@iidNuisance <- list(favorable = iidNuisance_favorable,
                                             unfavorable = iidNuisance_unfavorable)

                 if(!is.null(covariance)){
                     .Object@covariance <- covariance
                 }
                 .Object@tablePairScore <- tableScore

                 ## *** required additional information
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
                 .Object@weight <- weight
                 .Object@n.resampling <- n.resampling
                 
                 ## *** optional information
                 ## resampling
                 if(!is.null(deltaResampling)){
                     .Object@deltaResampling <- deltaResampling
                 }
                 if(!is.null(DeltaResampling)){
                     .Object@DeltaResampling <- DeltaResampling
                 }
                 if(!is.null(DeltaResampling)){
                     .Object@covarianceResampling <- covarianceResampling
                 }

                 ## survival
                 if(!is.null(tableSurvival)){
                     .Object@tableSurvival <- tableSurvival
                 }

                 ## ** export
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor S4BuyseTest objects
S4BuyseTest <- function(...) new("S4BuyseTest", ...) 
