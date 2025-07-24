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
#' @keywords classes
#' @author Brice Ozenne

## * Class S4BuyseTest
#' @rdname S4BuyseTest-class
#' @exportClass S4BuyseTest
setClass(
  
  Class = "S4BuyseTest",
  
  representation(
      call = "list",
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
      add.halfNeutral = "logical",
      correction.uninf = "numeric",
      method.inference = "character",
      strata = "vector",
      threshold = "numeric",
      restriction = "numeric",
      nResampling = "array",
      deltaResampling = "array",
      DeltaResampling = "array",
      covariance = "matrix",
      covarianceResampling = "array",
      weightObs = "numeric",
      weightEndpoint = "numeric",
      weightStrata = "numeric",
      weightStrataResampling = "array",
      iidAverage = "list",
      iidNuisance = "list",
      seed = "numeric",
      tablePairScore = "list",
      tableSurvival = "list"
      )

)

## * Initialize S4BuyseTest objects
methods::setMethod(
             f = "initialize", 
             signature = "S4BuyseTest", 
             definition = function(.Object,
                                   call,
                                   count_favorable, ## from cpp object
                                   count_unfavorable, ## from cpp object
                                   count_neutral, ## from cpp object
                                   count_uninf, ## from cpp object
                                   delta, ## from cpp object
                                   Delta, ## from cpp object
                                   n_pairs, ## from cpp object
                                   iidAverage_favorable, ## from cpp object
                                   iidAverage_unfavorable, ## from cpp object
                                   iidAverage_neutral, ## from cpp object
                                   iidNuisance_favorable, ## from cpp object
                                   iidNuisance_unfavorable, ## from cpp object
                                   iidNuisance_neutral, ## from cpp object
                                   covariance, ## from cpp object
                                   tableScore, ## from cpp object
                                   tableSurvival = NULL, ## added to the cpp object by .BuyseTest when requested by the user
                                   index.C,
                                   index.T,
                                   index.strata,
                                   type,
                                   endpoint,
                                   level.strata,
                                   level.treatment,
                                   scoring.rule, 
                                   hierarchical,
                                   neutral.as.uninf,
                                   add.halfNeutral,
                                   correction.uninf,
                                   method.inference,
                                   method.score,
                                   seed,
                                   strata,
                                   threshold,
                                   restriction,
                                   weightObs,
                                   weightEndpoint,
                                   weightStrata,
                                   pool.strata,
                                   grid.strata,
                                   n.resampling,
                                   nResampling = NULL, ## from inferenceResampling
                                   deltaResampling = NULL, ## from inferenceResampling
                                   DeltaResampling = NULL, ## from inferenceResampling
                                   weightStrataResampling = NULL, ## from inferenceResampling
                                   covarianceResampling = NULL ## from inferenceResampling
                                   ){

                 name.endpoint <- paste0(endpoint,ifelse(!is.na(restriction),paste0("_r",restriction),""),ifelse(threshold>1e-12,paste0("_t",threshold),""))
                 level.strata2 <- rownames(grid.strata)

                 ## ** call
                 call <- call[-1]

                 ## ** count
                 dimnames(count_favorable) <- list(level.strata2, name.endpoint)
                 dimnames(count_unfavorable) <- list(level.strata2, name.endpoint)
                 dimnames(count_neutral) <- list(level.strata2, name.endpoint)
                 dimnames(count_uninf) <- list(level.strata2, name.endpoint)

                 ## ** delta/Delta
                 dimnames(delta) <- list(level.strata2,
                                         name.endpoint,
                                         c("favorable","unfavorable","neutral","uninf","netBenefit","winRatio"))
                 dimnames(Delta) <- list(name.endpoint,
                                         c("favorable","unfavorable","neutral","uninf","netBenefit","winRatio"))

                 ## ** n_pairs
                 names(n_pairs) <- level.strata2

                 ## ** iid and variance
                 if(!is.null(iidAverage_favorable) && NCOL(iidAverage_favorable)>0){
                     colnames(iidAverage_favorable) <- name.endpoint
                 }
                 
                 if(!is.null(iidAverage_unfavorable) && NCOL(iidAverage_unfavorable)>0){
                     colnames(iidAverage_unfavorable) <- name.endpoint
                 }

                 if(!is.null(iidAverage_neutral) && NCOL(iidAverage_neutral)>0){
                     colnames(iidAverage_neutral) <- name.endpoint
                 }

                 if(!is.null(iidNuisance_favorable) && NCOL(iidNuisance_favorable)>0){
                     colnames(iidNuisance_favorable) <- name.endpoint
                 }
                 
                 if(!is.null(iidNuisance_unfavorable) && NCOL(iidNuisance_unfavorable)>0){
                     colnames(iidNuisance_unfavorable) <- name.endpoint
                 }

                 if(!is.null(iidNuisance_neutral) && NCOL(iidNuisance_neutral)>0){
                     colnames(iidNuisance_neutral) <- name.endpoint
                 }

                 if(!is.null(covariance) && length(covariance)>0){
                     dimnames(covariance) <- list(name.endpoint,
                                            c("favorable","unfavorable","covariance","netBenefit","winRatio"))
                 }

                 ## ** tableScore
                 if(!is.null(tableScore) && length(tableScore)>0 && any(sapply(tableScore, data.table::is.data.table)==FALSE)){
                     tableScore <- pairScore2dt(tableScore,
                                                level.treatment = level.treatment,
                                                level.strata = level.strata2,
                                                n.strata = length(level.strata2),
                                                endpoint = endpoint,
                                                threshold = threshold,
                                                restriction = restriction)
                 }
                 
                 ## ** tableSurvival

                 ## ** type
                 type <- stats::setNames(type, name.endpoint)

                 ## ** endpoint
                 names(endpoint) <- name.endpoint

                 ## ** level.strata
                 attr(level.strata2,"index") <- index.strata
                 attr(level.strata2,"original") <- level.strata

                 ## ** level.treatment
                 attr(level.treatment,"indexC") <- index.C
                 attr(level.treatment,"indexT") <- index.T

                 ## ** scoring.rule
                 efron <- (scoring.rule==2)
                 scoring.rule <- c("Gehan","Peron")[(scoring.rule>0)+1]
                 attr(scoring.rule,"test.censoring") <- attr(method.score, "test.censoring")
                 attr(method.score,"test.censoring") <- NULL
                 attr(scoring.rule,"test.CR") <- attr(method.score, "test.CR")
                 attr(method.score,"test.CR") <- NULL
                 attr(scoring.rule,"test.match") <- !is.null(strata) && attr(strata,"match")
                 attr(scoring.rule,"method.score") <- stats::setNames(method.score, name.endpoint)
                 attr(scoring.rule,"efron") <- efron

                 ## ** hierarchical
                 
                 ## ** neutral.as.uninf

                 ## ** add.halfNeutral
                 
                 ## ** correction.uninf
                 
                 ## ** method.inference
                 
                 ## ** method.score

                 ## ** strata
                 if(is.null(strata)){
                     strata <- as.character(NA)
                 }
                 ## ** restriction
                 names(restriction) <- name.endpoint

                 ## ** threshold
                 names(threshold) <- name.endpoint
                 
                 ## ** weightEndpoint
                 names(weightEndpoint) <- name.endpoint

                 ## ** weightStrata
                 weightStrata <- as.double(weightStrata)
                 attr(weightStrata,"type") <- attr(pool.strata,"type")

                 ## ** prepare resampling object
                 if(!is.null(deltaResampling)){                     
                     dimnames(deltaResampling)[[3]] <- name.endpoint
                     dimnames(DeltaResampling)[[2]] <- name.endpoint
                     if(attr(method.inference,"studentized")){
                         dimnames(covarianceResampling)[[2]] <- name.endpoint
                     }
                 }

                 ## ** resampling

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
                                            unfavorable = iidAverage_unfavorable,
                                            neutral = iidAverage_neutral)
                 .Object@iidNuisance <- list(favorable = iidNuisance_favorable,
                                             unfavorable = iidNuisance_unfavorable,
                                             neutral = iidNuisance_neutral)

                 if(!is.null(covariance)){
                     .Object@covariance <- covariance
                 }
                 .Object@tablePairScore <- tableScore

                 ## *** required additional information
                 .Object@call <- call
                 .Object@type <- type
                 .Object@endpoint <- endpoint
                 .Object@level.strata <- level.strata2
                 .Object@level.treatment <- level.treatment
                 .Object@scoring.rule <- scoring.rule
                 .Object@hierarchical <- hierarchical
                 .Object@neutral.as.uninf <- neutral.as.uninf
                 .Object@add.halfNeutral <- add.halfNeutral
                 .Object@correction.uninf <- correction.uninf
                 .Object@method.inference <- method.inference
                 .Object@strata <- strata
                 .Object@threshold <- threshold
                 .Object@restriction <- restriction
                 .Object@weightObs <- weightObs
                 .Object@weightEndpoint <- weightEndpoint
                 .Object@weightStrata <- weightStrata
                 if(!missing(seed)){
                     .Object@seed <- seed
                 }

                 ## *** optional information
                 ## resampling
                 if(!is.null(deltaResampling)){
                     .Object@nResampling <- nResampling
                     .Object@deltaResampling <- deltaResampling
                     .Object@DeltaResampling <- DeltaResampling
                     .Object@weightStrataResampling <- weightStrataResampling
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
