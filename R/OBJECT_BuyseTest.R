#' @name BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
#' @aliases  BuyseRes BuyseRes-class
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
#' 
#' @slot @delta: the chance of a better outcome (net difference between the probability for a random pair to be favorable minus the probability to be unfavorable divided by the total number of pairs) for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot @count_favorable :  the probability for a random pair to be favorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot @count_unfavorable.group : the probability for a random pair to be unfavorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot @count_neutral :  the probability for a random pair to be neutral for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot @count_uninf : the probability for a random pair to be uninformative for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
#' @slot @index_neutralT :  the index in the dataset of the treatment observations from remaining neutral pairs. \emph{integer vector}.
#' @slot @index_neutralC : the index in the dataset of the control observations from remaining neutral pairs. \emph{integer vector}.
#' @slot @index_uninfT :  the index in the dataset of the treatment observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot @index_uninfC : the index in the dataset of the control observations from remaining uninformative pairs. \emph{integer vector}.
#' @slot @n_pairs :  the total number of pairs. \emph{integer}.
#' @slot @delta_boot :  the chance of a better outcome within each strata (first dimension), over the endpoint (second dimension) and for each bootstrap dataset (third dimension). \emph{array}.
#' @slot @p.value : the p.value associated to the chance of a better outcome at each prioritized endpoint. \emph{numeric vector}.
#' @slot @Delta_quantile :  the randomization test-based 2.5\% and the 97.5\% quantiles of the cumulative chance of a better outcome (in rows) for each endpoint (in columns). \emph{matrix}.
#' @slot @endpoint : the name of the endpoint. \emph{character vector}.
#' @slot @strata : the name of the strata \emph{character vector}.
#'   
#' @seealso 
#' \code{\link{BuyseTest}} for the function computing generalized pairwise comparisons. \cr
#' \code{\link{summary,BuyseRes-method}} for the summary of the BuyseTest function results
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
    delta = "matrix", 
    count_favorable = "matrix",      
    count_unfavorable = "matrix",
    count_neutral = "matrix",    
    count_uninf = "matrix",
    index_neutralT = "vector",
    index_neutralC = "vector",
    index_uninfT = "vector",
    index_uninfC = "vector",
    n_pairs = "numeric",
    delta_boot = "array", 
    p.value = "vector",    
    Delta_quantile = "matrix",
    endpoint = "vector",
    threshold = "vector",
    strata = "vector",
    levels.treatment = "vector"
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

#' Wrapper function BuyseRes
#'
#' @rdname BuyseRes-class
#' @export
BuyseRes <- function(...) new("BuyseRes", ...) 
