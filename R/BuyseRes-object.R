## * Documentation BuyseRes
#' @name BuyseRes-class
#' @title Class "BuyseRes" (output of BuyseTest)
#' @aliases  BuyseRes BuyseRes-class
#' 
#' @description A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
#' 
#' @slot count_favorable [matrix] the probability for a random pair to be favorable for each strata (in rows) and each endpoint (in columns).
#' @slot count_unfavorable [matrix] the probability for a random pair to be unfavorable for each strata (in rows) and each endpoint (in columns).
#' @slot count_neutral [matrix] the probability for a random pair to be neutral for each strata (in rows) and each endpoint (in columns). 
#' @slot count_uninf [matrix] the probability for a random pair to be uninformative for each strata (in rows) and each endpoint (in columns). 
#' @slot delta [list of matrix] the chance of a better outcome
#' (net difference between the probability for a random pair to be favorable minus the probability to be unfavorable divided by the total number of pairs)
#' for each strata (in rows) and each endpoint (in columns).
#' @slot Delta [list of numeric vector] the chance of a better outcome over all strata for each endpoint.
#' @slot delta_permutation [list of array] the chance of a better outcome within each strata (first dimension),
#' over the endpoint (second dimension) and for each sample from the permutation test (third dimension).
#' @slot Delta_quantile [list of matrix] the randomization test-based 2.5\% and the 97.5\% quantiles of the cumulative chance of a better outcome (in rows) for each endpoint (in columns).
#' @slot endpoint [character vector] the name of the endpoint. 
#' @slot index_neutralT [integer vector] the index in the dataset of the treatment observations from remaining neutral pairs.
#' @slot index_neutralC [integer vector] the index in the dataset of the control observations from remaining neutral pairs.
#' @slot index_uninfT [integer vector] the index in the dataset of the treatment observations from remaining uninformative pairs.
#' @slot index_uninfC [integer vector] the index in the dataset of the control observations from remaining uninformative pairs.
#' @slot levels.treatment [character vector] the name of each group.
#' @slot n_pairs [integer] the total number of pairs.
#' @slot n_permutation [integer] the number of sucessful samples of the permutation test.
#' @slot p.value [list of numeric vector] the p.value associated to the chance of a better outcome at each prioritized endpoint.
#' @slot strata [character vector] the name of the strata.
#' @slot threshold [numeric vector] the threshold associated to each endpoint.
#' @slot conf.level [numeric] the confidence level.
#' @slot tableComparison [list] detail of the result of each pairwise comparison relative to each endpoint.
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
#'                                 treatment="treatment", type=c("bin","bin"), n.permutation = 0)
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
    delta_permutation = "list", 
    Delta_quantile = "list",
    endpoint = "vector",
    index_neutralT = "vector",
    index_neutralC = "vector",
    index_uninfT = "vector",
    index_uninfC = "vector",
    levels.treatment = "vector",
    n_permutation = "list",
    n_pairs = "numeric",
    p.value = "list",
    strata = "vector",
    threshold = "vector",
    conf.level = "numeric",
    tableComparison = "list"
    ),

  ### ** Check validity of the object
  validity = function(object){

      n.strata <- length(object@strata)
      n.outcome <- length(object@endpoint)
    
    
    validDimension(object@count_favorable,
                   name1 = "@count_favorable",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    validDimension(object@count_unfavorable,
                   name1 = "@count_unfavorable",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    validDimension(object@count_neutral,
                   name1 = "@count_neutral",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    validDimension(object@count_uninf,
                   name1 = "@count_uninf",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@delta[[1]],
                   name1 = "@delta[[1]]",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    validDimension(object@delta[[2]],
                   name1 = "@delta[[2]]",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@Delta[[1]],
                   name1 = "@Delta[[1]]",
                   valid.dimension =  n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    validDimension(object@Delta[[2]],
                   name1 = "@Delta[[2]]",
                   valid.dimension =  n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@delta_permutation[[1]],
                   name1 = "@delta_permutation[[1]]",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]") # do not check the number successful of sample for the permutation test
    validDimension(object@delta_permutation[[2]],
                   name1 = "@delta_permutation[[2]]",
                   valid.dimension =  c(n.strata, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]") # do not check the number successful of sample for the permutation test
    
    validDimension(object@Delta_quantile[[1]],
                   name1 = "@Delta_quantile[[1]]",
                   valid.dimension =  c(2, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    validDimension(object@Delta_quantile[[2]],
                   name1 = "@Delta_quantile[[2]]",
                   valid.dimension =  c(2, n.outcome),
                   type = c("NROW","NCOL"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@index_neutralT,
                   name1 = "@index_neutralT",
                   value2 = object@index_neutralC,
                   name2 = "@index_neutralC",
                   type = c("length"),
                   method = "validity[BuyseTest]")
    validDimension(object@index_uninfT,
                   name1 = "@index_uninfT",
                   value2 = object@index_uninfC,
                   name2 = "@index_uninfC",
                   type = c("length"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@levels.treatment,
                   name1 = "@levels.treatment",
                   valid.dimension = 2,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@n_permutation[[1]],
                   name1 = "@n_permutation[[1]]",
                   valid.dimension = n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    validDimension(object@n_permutation[[2]],
                   name1 = "@n_permutation[[2]]",
                   valid.dimension = n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    
    validDimension(object@p.value[[1]],
                   name1 = "@p.value[[1]]",
                   valid.dimension = n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    validDimension(object@p.value[[2]],
                   name1 = "@p.value[[2]]",
                   valid.dimension = n.outcome,
                   type = c("length"),
                   method = "validity[BuyseTest]")
    
      validDimension(object@threshold,
                     name1 = "@threshold",
                     valid.dimension = n.outcome,
                     type = c("length"),
                     method = "validity[BuyseTest]")
      validNumeric(object@conf.level,
                   name1 = "@conf.level",
                   valid.length = 1,
                   min = 0, max = 1, refuse.NA = FALSE,
                   method = "validity[BuyseTest]")
  
    return(TRUE)} 
)

## * Initialize BuyseRes objects
methods::setMethod(
  f = "initialize", 
  signature = "BuyseRes", 
  definition = function(.Object, 
                        count_favorable, count_unfavorable, count_neutral, count_uninf, 
                        delta, Delta, delta_permutation, Delta_quantile,
                        endpoint,
                        index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                        levels.treatment, n_permutation, n_pairs, p.value, strata, threshold, conf.level,
                        tableComparison, args){
    
    n.strata <- length(strata)
    D <- length(endpoint)
    
      if(missing(delta_permutation)){
          delta_permutation <- list(netChance = array(NA,dim = c(n.strata,D,1)),
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
      if(missing(n_permutation)){
          n_permutation <- list(netChance = rep(NA,D),
                                winRatio = rep(NA,D)
                                )
      }
      if(missing(tableComparison) || is.null(tableComparison)){
          tableComparison <- list()
      }else{
          ## Rcpp outputs vector: convert to matrix and rename
          name.tempo <- c("strata",
                          paste0("index.",levels.treatment[2]), ## treated
                          paste0("index.",levels.treatment[1]), ## control
                          paste0("indexWithinStrata.",levels.treatment[2]), ## treated
                          paste0("indexWithinStrata.",levels.treatment[1]), ## control
                          "favorable","unfavorable","neutral","uninformative")

          tableComparison <- lapply(tableComparison, function(iC){
              iM <- as.data.frame(matrix(iC, ncol = 9, byrow = FALSE,
                                         dimnames = list(NULL,name.tempo)))              
              iM[,"strata"] <- factor(iM[,"strata"], levels = 0:(n.strata-1), labels = strata) ## indexes start at 1 in R and not at 0 as in C++

              ## recall that indexes start at 1 in R and not at 0 as in C++
              iM[,2] <- args$indexT[iM[,2]+1] ## restaure position in the original dataset, not the datasets relative to T and C
              iM[,3] <- args$indexC[iM[,3]+1]
              iM[,4:5] <- iM[,4:5] + 1 
              return(iM)
          })
          names(tableComparison) <- paste0(endpoint,"_",threshold)
      }
      .Object@count_favorable <- count_favorable      
      .Object@count_unfavorable <- count_unfavorable
      .Object@count_neutral <- count_neutral   
      .Object@count_uninf <- count_uninf
      .Object@delta <- delta
      .Object@Delta <- Delta
      .Object@delta_permutation <- delta_permutation
      .Object@Delta_quantile <- Delta_quantile
      .Object@endpoint <- endpoint
      .Object@index_neutralT <- index_neutralT
      .Object@index_neutralC <- index_neutralC
      .Object@index_uninfT <- index_uninfT
      .Object@index_uninfC <- index_uninfC
      .Object@levels.treatment <- levels.treatment
      .Object@n_permutation <- n_permutation
      .Object@n_pairs <- n_pairs
      .Object@p.value <- p.value    
      .Object@strata <- strata
      .Object@threshold <- threshold
      .Object@conf.level <- conf.level
      .Object@tableComparison <- tableComparison
    
    validObject(.Object)
    return(.Object)
    
})


## * Constructor BuyseRes objects
BuyseRes <- function(...) new("BuyseRes", ...) 
