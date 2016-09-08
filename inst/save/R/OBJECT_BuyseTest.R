#' @name BuyseRes-class
#' @rdname BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
#' @aliases  BuyseRes BuyseRes-class
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
#' 
#' @slot levels.treatment the name of each group. \emph{character vector}.
#' @slot endpoint the name of the endpoint. \emph{character vector}.
#' @slot strata the name of the strata \emph{character vector}.
#' @slot threshold the threshold associated to each endpoint \emph{numeric vector}.
#' @slot n_pairs the total number of pairs. \emph{integer}.
#' @slot count_favorable the probability for a random pair to be favorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_unfavorable the probability for a random pair to be unfavorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_neutral the probability for a random pair to be neutral for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_uninf the probability for a random pair to be uninformative for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot index_neutralT the index in the dataset of the treatment observations from remaining neutral pairs. \emph{integer vector}.
#' @slot index_neutralC the index in the dataset of the control observations from remaining neutral pairs. \emph{integer vector}.
#' @slot index_uninfT the index in the dataset of the treatment observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot index_uninfC the index in the dataset of the control observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot delta the chance of a better outcome (net difference between the probability for a random pair to be favorable minus the probability to be unfavorable divided by the total number of pairs) for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot delta_boot the chance of a better outcome within each strata (first dimension), over the endpoint (second dimension) and for each bootstrap dataset (third dimension). \emph{array}.
#' @slot Delta_quantile the randomization test-based 2.5\% and the 97.5\% quantiles of the cumulative chance of a better outcome (in rows) for each endpoint (in columns). \emph{matrix}.
#' @slot p.value the p.value associated to the chance of a better outcome at each prioritized endpoint. \emph{numeric vector}.
#'   
#' @seealso 
#' \code{\link{BuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{BuyseRes-summary}} for the summary of the BuyseTest function results
#' 
#' @examples
#'   n.Treatment_testBin <- 500
#'   n.Control_testBin <- 500
#'   prob.Treatment_testBin <- c(0.5,0.75)
#'   prob.Control_testBin <- c(0.5,0.25)
#'   
#'   set.seed(10)
#'   data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
#'   data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
#'   data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
#'   data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)
#'   
#'   #### no strata, n.bootsrap=0
#'   BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                 treatment="treatment", type=c("bin","bin"))
#'   
#'   class(BuyseTest_object)
#'   
#' 
#' @keywords classes BuyseRes-class
#' 

#' @rdname BuyseRes-class
#' @exportClass BuyseRes
setClass(
  
  Class = "BuyseRes",
  
  representation(
    levels.treatment = "vector",
    endpoint = "vector",
    strata = "vector",
    threshold = "vector",
    n_pairs = "numeric",
    count_favorable = "matrix",      
    count_unfavorable = "matrix",
    count_neutral = "matrix",    
    count_uninf = "matrix",
    index_neutralT = "vector",
    index_neutralC = "vector",
    index_uninfT = "vector",
    index_uninfC = "vector",
    delta = "matrix", 
    delta_boot = "array", 
    Delta_quantile = "matrix",
    p.value = "vector"    
  ),
  
  validity = function(object){
    #cat("--- BuyseRes : checking --- ")
    
    n.strata <- length(object@strata)
    n.outcome <- length(object@endpoint)
    
    validDimension(object@delta, name1 = "@delta", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_favorable, name1 = "@count_favorable", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_unfavorable, name1 = "@count_unfavorable", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_neutral, name1 = "@count_neutral", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_uninf, name1 = "@count_uninf", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    
    validDimension(object@index_neutralT, name1 = "@index_neutralT", value2 = object@index_neutralC, name2 = "@index_neutralC", type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@index_uninfT, name1 = "@index_uninfT", value2 = object@index_uninfC, name2 = "@index_uninfC", type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@p.value, name1 = "@p.value", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@delta_boot, name1 = "@delta_boot", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@Delta_quantile, name1 = "@Delta_quantile", validDimension =  c(2, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    
    validDimension(object@threshold, name1 = "@threshold", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@levels.treatment, name1 = "@levels.treatment", validDimension = 2, type = c("length"), method = "validity[BuyseTest]")
    
    #cat(" : valid BuyseRes  \n")
    return(TRUE)} 
)

methods::setMethod(
  f = "initialize", 
  signature = "BuyseRes", 
  definition = function(.Object, 
                        levels.treatment,endpoint, strata, threshold, n_pairs, 
                        count_favorable, count_unfavorable, count_neutral, count_uninf, 
                        index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                        delta, delta_boot, Delta_quantile, p.value){
    
    n.strata <- length(strata)
    D <- length(endpoint)
    
    if(missing(delta_boot)){
      delta_boot <- array(NA,dim = c(n.strata,D,1))
    }
    if(missing(p.value)){
      p.value <- rep(NA,D)
    }
    if(missing(Delta_quantile)){
      Delta_quantile <- matrix(NA,nrow = 2, ncol = D, dimnames = list(c("2.5%","97.5%")))
    }
    
    .Object@delta <- delta
    .Object@count_favorable <- count_favorable      
    .Object@count_unfavorable <- count_unfavorable
    .Object@count_neutral <- count_neutral   
    .Object@count_uninf <- count_uninf
    .Object@index_neutralT <- index_neutralT
    .Object@index_neutralC <- index_neutralC
    .Object@index_uninfT <- index_uninfT
    .Object@index_uninfC <- index_uninfC
    .Object@n_pairs <- n_pairs
    .Object@delta_boot <- delta_boot
    .Object@p.value <- p.value    
    .Object@Delta_quantile <- Delta_quantile
    .Object@endpoint <- endpoint
    .Object@threshold <- threshold
    .Object@strata <- strata
    .Object@levels.treatment <- levels.treatment
    
    validObject(.Object)
    return(.Object)
    
})

BuyseRes <- function(...) new("BuyseRes", ...) 
