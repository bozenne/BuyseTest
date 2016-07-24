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
  
  Class="BuyseRes",
  
  representation(
    delta = "matrix", 
    count_favorable = "matrix",      
    count_unfavorable = "matrix",
    count_neutral = "matrix",    
    count_uninf = "matrix",
    index_neutralT="vector",
    index_neutralC = "vector",
    index_uninfT="vector",
    index_uninfC="vector",
    n_pairs = "numeric",
    delta_boot = "array", 
    p.value = "vector",    
    Delta_quantile = "matrix",
    endpoint = "vector",
    threshold = "vector",
    strata = "vector"
  ),
  
  validity = function(object){
    #cat("--- BuyseRes : checking --- ")
    
    n.strata <- nrow(object@delta)
    n.outcome <- ncol(object@delta)
    
    if(nrow(object@count_favorable)!=n.strata || ncol(object@count_favorable)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_favorable\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_favorable) : ",row(object@count_favorable)," ",ncol(object@count_favorable)," \n")
    }
    
    if(nrow(object@count_unfavorable)!=n.strata || ncol(object@count_unfavorable)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_unfavorable\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_favorable) : ",row(object@count_unfavorable)," ",ncol(object@count_unfavorable)," \n")
    }
    
    if(nrow(object@count_neutral)!=n.strata || ncol(object@count_neutral)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_neutral\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_neutral) : ",row(object@count_neutral)," ",ncol(object@count_neutral)," \n")
    }
    
    if(nrow(object@count_uninf)!=n.strata || ncol(object@count_uninf)!=n.outcome)
    {stop("validity[BuyseRes] : \'@count_uninf\' does not match  \'@delta\' dimensions \n",
          "dim(@delta) : ",n.strata," ",n.outcome," \n",
          "dim(@count_uninf) : ",row(object@count_uninf)," ",ncol(object@count_uninf)," \n")
    }
    
    if(length(object@index_neutralT)!=length(object@index_neutralC))
    {stop("validity[BuyseRes] : \'@index_neutralT\' does not match  \'@index_neutralC\' dimensions \n",
          "length(@index_neutralT) : ",length(object@index_neutralT)," \n",
          "length(@index_neutralC) : ",length(object@index_neutralC)," \n")
    }
    
    if(length(object@index_uninfT)!=length(object@index_uninfC))
    {stop("validity[BuyseRes] : \'@index_uninfT\' does not match  \'@index_uninfT\' dimensions \n",
          "length(@index_uninfT) : ",length(object@index_uninfT)," \n",
          "length(@index_uninfC) : ",length(object@index_uninfC)," \n")
    }
    
    #     if(object@n_pairs %% 1 != 0)
    #     {stop("validity[BuyseRes] : wrong specification of \'@n_pairs\' \n",
    #           "\'n_pairs\' must be an integer \n",
    #           "@n_pairs : ",object@n_pairs," \n")
    #     }
    
    if(length(object@p.value)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@n.outcome\' \n",
          "must have length n.outcome : ",n.outcome," \n",
          "length(@p.value) : ",length(object@p.value)," \n")
    }
    
    if(dim(object@delta_boot)[1]!=n.strata || dim(object@delta_boot)[2]!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@delta_boot\' \n",
          "must have dim[1]=",n.strata," and dim[2]=",n.outcome," \n",
          "dim(object@delta_boot) : ",paste(dim(object@delta_boot),collapse=" ")," \n")
    }
    
    if(nrow(object@Delta_quantile)!=2 || ncol(object@Delta_quantile)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@Delta_quantile\' \n",
          "must have dimensions : 2 ",n.outcome," \n",
          "dim(@Delta_quantile) : ",nrow(object@Delta_quantile)," ",ncol(object@Delta_quantile)," \n")
    }
    
    if(length(object@endpoint)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@endpoint\' \n",
          "must have length n.outcome : ",n.outcome," \n",
          "length(@endpoint) : ",length(object@endpoint)," \n")
    }
    
    if(length(object@threshold)!=n.outcome)
    {stop("validity[BuyseRes] : wrong specification of \'@threshold\' \n",
          "must have length n.outcome : ",n.outcome," \n",
          "length(@threshold) : ",length(object@threshold)," \n")
    }
    
    #cat(" : valid BuyseRes  \n")
    return(TRUE)} 
  
  
)

#' Wrapper function BuyseRes
#'
#' @rdname BuyseRes-class
#' @export
BuyseRes <- function(...) new("BuyseRes", ...) 
