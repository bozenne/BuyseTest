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
      endpoint = "vector",
      level.treatment = "vector",
      level.strata = "vector",
      method = "character",
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
                                   count.favorable, count.unfavorable, count.neutral, count.uninf, n.pairs,
                                   delta.netChance, delta.winRatio, Delta.netChance, Delta.winRatio,
                                   index.neutralT, index.neutralC, index.uninfT, index.uninfC, 
                                   endpoint, level.strata, level.treatment, method, method.inference, strata, threshold,
                                   n.resampling,
                                   deltaResampling.netChance, deltaResampling.winRatio, DeltaResampling.netChance, DeltaResampling.winRatio, 
                                   tableComparison,
                                   args){


                 ## ** count
                 dimnames(count.favorable) <- list(level.strata, endpoint)
                 dimnames(count.unfavorable) <- list(level.strata, endpoint)
                 dimnames(count.neutral) <- list(level.strata, endpoint)
                 dimnames(count.uninf) <- list(level.strata, endpoint)

                 ## ** delta/Delta
                 dimnames(delta.netChance) <- list(level.strata, endpoint)
                 dimnames(delta.winRatio) <- list(level.strata, endpoint)
                 names(Delta.netChance) <- endpoint
                 names(Delta.winRatio) <- endpoint
                 
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
                     ## Rcpp outputs vector: convert to matrix and rename
                     name.tempo <- c("strata",
                                     paste0("index.",level.treatment[2]), ## treated
                                     paste0("index.",level.treatment[1]), ## control
                                     paste0("indexWithinStrata.",level.treatment[2]), ## treated
                                     paste0("indexWithinStrata.",level.treatment[1]), ## control
                                     "favorable","unfavorable","neutral","uninformative")

                     tableComparison <- lapply(tableComparison, function(iC){
                         iM <- as.data.frame(matrix(iC, ncol = 9, byrow = FALSE,
                                                    dimnames = list(NULL,name.tempo)))
                         iM[,"strata"] <- factor(iM[,"strata"], levels = 0:(n.strata-1), labels = level.strata) ## indexes start at 1 in R and not at 0 as in C++
                         ## recall that indexes start at 1 in R and not at 0 as in C++
                         iM[,2] <- args$indexT[iM[,2]+1] ## restaure position in the original dataset, not the datasets relative to T and C
                         iM[,3] <- args$indexC[iM[,3]+1]
                         iM[,4:5] <- iM[,4:5] + 1 
                         return(iM)
                     })
                     names(tableComparison) <- paste0(endpoint,"_",threshold)
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
                 
                 .Object@endpoint <- endpoint
                 .Object@level.strata <- level.strata
                 .Object@level.treatment <- level.treatment
                 .Object@method <- method
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
