## * Documentation BuyseRes
#' @name BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
#' @aliases  BuyseRes BuyseRes-class
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
#' 
#' @slot count_favorable the probability for a random pair to be favorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_unfavorable the probability for a random pair to be unfavorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_neutral the probability for a random pair to be neutral for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot count_uninf the probability for a random pair to be uninformative for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot delta the chance of a better outcome (net difference between the probability for a random pair to be favorable minus the probability to be unfavorable divided by the total number of pairs) for each strata (in rows) and each endpoint (in columns). \emph{list of matrix}.
#' @slot Delta the chance of a better outcome over all strata for each endpoint. \emph{list of numeric vector}.
#' @slot delta_boot the chance of a better outcome within each strata (first dimension), over the endpoint (second dimension) and for each bootstrap dataset (third dimension). \emph{list of array}.
#' @slot Delta_quantile the randomization test-based 2.5\% and the 97.5\% quantiles of the cumulative chance of a better outcome (in rows) for each endpoint (in columns). \emph{list of matrix}.
#' @slot endpoint the name of the endpoint. \emph{character vector}.
#' @slot index_neutralT the index in the dataset of the treatment observations from remaining neutral pairs. \emph{integer vector}.
#' @slot index_neutralC the index in the dataset of the control observations from remaining neutral pairs. \emph{integer vector}.
#' @slot index_uninfT the index in the dataset of the treatment observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot index_uninfC the index in the dataset of the control observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot levels.treatment the name of each group. \emph{character vector}.
#' @slot n_pairs the total number of pairs. \emph{integer}.
#' @slot n_bootstrap the number of sucessful bootstrap samples of the resampling procedure. \emph{integer}.
#' @slot p.value the p.value associated to the chance of a better outcome at each prioritized endpoint. \emph{list of numeric vector}.
#' @slot strata the name of the strata \emph{character vector}.
#' @slot threshold the threshold associated to each endpoint \emph{numeric vector}.
#'   
#' @details slots 
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
#'   #### no strata,
#'   BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                 treatment="treatment", type=c("bin","bin"), n.bootstrap = 0)
#'   
#'   class(BuyseTest_object)
#'   
#' 
#' @keywords classes BuyseRes-class
#' 

## * Class BuyseRes
#' @rdname BuyseRes-class
#' @exportClass BuyseRes
setClass(
  
  Class = "BuyseRes",
  
  representation(
    count_favorable = "matrix",      
    count_unfavorable = "matrix",
    count_neutral = "matrix",    
    count_uninf = "matrix",
    delta = "list", 
    Delta = "list", 
    delta_boot = "list", 
    Delta_quantile = "list",
    endpoint = "vector",
    index_neutralT = "vector",
    index_neutralC = "vector",
    index_uninfT = "vector",
    index_uninfC = "vector",
    levels.treatment = "vector",
    n_bootstrap = "list",
    n_pairs = "numeric",
    p.value = "list",
    strata = "vector",
    threshold = "vector"
    ),

  ### ** Check validity of the object
  validity = function(object){
    
    n.strata <- length(object@strata)
    n.outcome <- length(object@endpoint)
    
    
    validDimension(object@count_favorable, name1 = "@count_favorable", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_unfavorable, name1 = "@count_unfavorable", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_neutral, name1 = "@count_neutral", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@count_uninf, name1 = "@count_uninf", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    
    validDimension(object@delta[[1]], name1 = "@delta[[1]]", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@delta[[2]], name1 = "@delta[[2]]", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    
    validDimension(object@Delta[[1]], name1 = "@Delta[[1]]", validDimension =  n.outcome, type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@Delta[[2]], name1 = "@Delta[[2]]", validDimension =  n.outcome, type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@delta_boot[[1]], name1 = "@delta_boot[[1]]", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]") # do not check the number of bootstrap
    validDimension(object@delta_boot[[2]], name1 = "@delta_boot[[2]]", validDimension =  c(n.strata, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]") # do not check the number of bootstrap
    
    validDimension(object@Delta_quantile[[1]], name1 = "@Delta_quantile[[1]]", validDimension =  c(2, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    validDimension(object@Delta_quantile[[2]], name1 = "@Delta_quantile[[2]]", validDimension =  c(2, n.outcome), type = c("NROW","NCOL"), method = "validity[BuyseTest]")
    
    validDimension(object@index_neutralT, name1 = "@index_neutralT", value2 = object@index_neutralC, name2 = "@index_neutralC", type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@index_uninfT, name1 = "@index_uninfT", value2 = object@index_uninfC, name2 = "@index_uninfC", type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@levels.treatment, name1 = "@levels.treatment", validDimension = 2, type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@n_bootstrap[[1]], name1 = "@n_bootstrap[[1]]", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@n_bootstrap[[2]], name1 = "@n_bootstrap[[2]]", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@p.value[[1]], name1 = "@p.value[[1]]", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    validDimension(object@p.value[[2]], name1 = "@p.value[[2]]", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
    
    validDimension(object@threshold, name1 = "@threshold", validDimension = n.outcome, type = c("length"), method = "validity[BuyseTest]")
  
    return(TRUE)} 
)

## * Initialiwe BuyseRes objects
methods::setMethod(
  f = "initialize", 
  signature = "BuyseRes", 
  definition = function(.Object, 
                        count_favorable, count_unfavorable, count_neutral, count_uninf, 
                        delta, Delta, delta_boot, Delta_quantile,
                        endpoint,
                        index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                        levels.treatment, n_bootstrap, n_pairs, p.value, strata, threshold){
    
    n.strata <- length(strata)
    D <- length(endpoint)
    
    if(missing(delta_boot)){
      delta_boot <- list(netChance = array(NA,dim = c(n.strata,D,1)),
                         winRatio = array(NA,dim = c(n.strata,D,1))
      )
    }
    if(missing(Delta_quantile)){
      Delta_quantile <- list(netChance = matrix(NA,nrow = 2, ncol = D, dimnames = list(c("lower%","upper%"))),
                             winRatio = matrix(NA,nrow = 2, ncol = D, dimnames = list(c("lower%","upper%")))
      )
    }
    if(missing(p.value)){
      p.value <- list(netChance = rep(NA,D),
                      winRatio = rep(NA,D)
      )
    }
    if(missing(n_bootstrap)){
      n_bootstrap <- list(netChance = rep(NA,D),
                          winRatio = rep(NA,D)
      )
    }
    
    .Object@count_favorable <- count_favorable      
    .Object@count_unfavorable <- count_unfavorable
    .Object@count_neutral <- count_neutral   
    .Object@count_uninf <- count_uninf
    .Object@delta <- delta
    .Object@Delta <- Delta
    .Object@delta_boot <- delta_boot
    .Object@Delta_quantile <- Delta_quantile
    .Object@endpoint <- endpoint
    .Object@index_neutralT <- index_neutralT
    .Object@index_neutralC <- index_neutralC
    .Object@index_uninfT <- index_uninfT
    .Object@index_uninfC <- index_uninfC
    .Object@levels.treatment <- levels.treatment
    .Object@n_bootstrap <- n_bootstrap
    .Object@n_pairs <- n_pairs
    .Object@p.value <- p.value    
    .Object@strata <- strata
    .Object@threshold <- threshold
    
    validObject(.Object)
    return(.Object)
    
})


## * Constructor BuyseRes objects
BuyseRes <- function(...) new("BuyseRes", ...) 
