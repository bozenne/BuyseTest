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
      tableComparison = "list")

  ## ### ** Check validity of the object
  ## , validity = function(object){

  ##     n.strata <- length(object@level.strata)
  ##     n.outcome <- length(object@endpoint)

  ##     ## *** count
  ##     validDimension(object@count.favorable,
  ##                    name1 = "@count.favorable",
  ##                    valid.dimension = c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@count.unfavorable,
  ##                    name1 = "@count.unfavorable",
  ##                    valid.dimension =  c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@count.neutral,
  ##                    name1 = "@count.neutral",
  ##                    valid.dimension = c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@count.uninf,
  ##                    name1 = "@count.uninf",
  ##                    valid.dimension = c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** delta
  ##     validDimension(object@delta.netChance,
  ##                    name1 = "@delta.netChance",
  ##                    valid.dimension = c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@delta.winRatio,
  ##                    name1 = "@delta.winRatio",
  ##                    valid.dimension = c(n.strata, n.outcome),
  ##                    type = c("NROW","NCOL"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** Delta
  ##     validDimension(object@Delta.netChance,
  ##                    name1 = "@Delta.netChance",
  ##                    valid.dimension = n.outcome,
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@Delta.winRatio,
  ##                    name1 = "@Delta.winRatio",
  ##                    valid.dimension = n.outcome,
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** index
  ##     validDimension(object@index.neutralT,
  ##                    name1 = "@index.neutralT",
  ##                    value2 = object@index.neutralC,
  ##                    name2 = "@index.neutralC",
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")
  ##     validDimension(object@index.uninfT,
  ##                    name1 = "@index.uninfT",
  ##                    value2 = object@index.uninfC,
  ##                    name2 = "@index.uninfC",
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** level.strata
  ##     validDimension(object@level.strata,
  ##                    name1 = "@level.strata",
  ##                    valid.dimension = n.strata,
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** level.treatment
  ##     validDimension(object@level.treatment,
  ##                    name1 = "@level.treatment",
  ##                    valid.dimension = 2,
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")

  ##     ## *** threshold
  ##     validDimension(object@threshold,
  ##                    name1 = "@threshold",
  ##                    valid.dimension = n.outcome,
  ##                    type = c("length"),
  ##                    method = "validity[BuyseTest]")
  
  ##     ## *** delta resampling
  ##     ## resampling can fail - should not be tested against object@n.resampling
      
  ##     ## *** Delta resampling
  ##     ## resampling can fail - should not be tested against object@n.resampling
     
  ##   return(TRUE)} 
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
                                   tableComparison,
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

                 ## ** no resampling
                 if(method.inference %in% c("none","asymptotic")){
                     deltaResampling.netChance <- array(dim=c(0,0,0))
                     deltaResampling.winRatio <- array(dim=c(0,0,0))
                     DeltaResampling.netChance <- matrix(NA, nrow = 0, ncol = 0)
                     DeltaResampling.winRatio <- matrix(NA, nrow = 0, ncol = 0)
                     n.resampling <- as.double(NA)
                 }

                 ## ** tableComparison
                 nComparison <- unlist(lapply(tableComparison,length))
                 if(any(nComparison>0)){
                     tableComparison <- tableComparison2dt(tableComparison,
                                                           level.treatment = level.treatment,
                                                           level.strata = level.strata,
                                                           n.strata = n.strata,
                                                           endpoint = endpoint,
                                                           threshold = threshold,
                                                           indexT = args$indexT,
                                                           indexC = args$indexC)
                 }
                     
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

                 .Object@tableComparison <- tableComparison
                 
                 ## validObject(.Object)
                 return(.Object)
                 
             })


## * Constructor BuyseRes objects
BuyseRes <- function(...) new("BuyseRes", ...) 
