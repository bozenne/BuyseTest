### BuyseRes-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2019 (11:20) 
## Version: 
## Last-Updated: feb 20 2020 (13:56) 
##           By: Brice Ozenne
##     Update #: 102
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - iid
#' @docType methods
#' @name BuyseRes-iid
#' @title Extract the H-decomposition of the Estimator
#' @aliases iid iid,BuyseRes-method
#' @include BuyseRes-object.R
#' 
#' @description Extract the H-decomposition of the GPC estimator.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [character] for which endpoint(s) the H-decomposition should be output?
#' If \code{NULL} returns the sum of the H-decomposition over all endpoints.
#' @param type [character] type of iid to be output.
#' Can be only for the nuisance parameters (\code{"nuisance"}),
#' or for the u-statistic given the nuisance parameters (\code{"u-statistic"}),
#' or both.
#' @param normalize [logical] if \code{TRUE} the iid is centered and multiplied by the sample size.
#' Otherwise not.
#' @param cluster [numeric vector] return the H-decomposition aggregated by cluster.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-summary}} for a more detailed presentation of the \code{BuyseRes} object.
#' 
#' @keywords iid BuyseRes-method
#' @author Brice Ozenne

## * Method - iid
#' @rdname BuyseRes-iid
#' @exportMethod iid
setMethod(f = "iid",
          signature = "BuyseRes",
          definition = function(object, endpoint = NULL, normalize = TRUE, type = "all", cluster = NULL){

              n.obs <- NROW(object@iidAverage$favorable)
              valid.endpoint <- paste0(object@endpoint,"_",object@threshold)
              n.endpoint <- length(valid.endpoint)
              if(!is.null(cluster) && !is.numeric(cluster)){
                  cluster <- as.numeric(as.factor(cluster))
              }
              
              ## ** check arguments              
              if(is.numeric(endpoint)){
                  validInteger(endpoint,
                               name1 = "endpoint",
                               min = 1, max = length(valid.endpoint),
                               valid.length = NULL,
                               method = "iid[BuyseTest]")
                  endpoint <- valid.endpoint[endpoint]
              }   
              validCharacter(endpoint,
                             valid.length = 1:length(valid.endpoint),
                             valid.values = valid.endpoint,
                             refuse.NULL = FALSE)
              validCharacter(type,
                             valid.length = 1,
                             valid.values = c("all","nuisance","u-statistic"),
                             refuse.NULL = FALSE)
              if(object@method.inference != "u-statistic"){
                  stop("No H-decomposition in the object \n",
                       "Set the argument \'method.inference\' to \"u-statistic\" when calling BuyseTest \n")
              }
              validInteger(cluster, valid.length = n.obs, min = 1, max = n.obs, refuse.NA = TRUE, refuse.NULL = FALSE, refuse.duplicates = FALSE)

              ## ** extract H-decomposition
              if(type %in% c("all","u-statistic")){
                  object.iid <- object@iidAverage
              }else{
                  object.iid <- list(favorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                        dimnames = list(NULL, valid.endpoint)),
                                     unfavorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                          dimnames = list(NULL, valid.endpoint))
                                     )
              }
              if(type %in% c("all","nuisance") && (object@scoring.rule=="Peron")){
                  object.iid$favorable <- object.iid$favorable + object@iidNuisance$favorable
                  object.iid$unfavorable <- object.iid$unfavorable + object@iidNuisance$unfavorable
              }

              if(normalize==FALSE){
                  delta.favorable <- colSums(object@count.favorable)/sum(object@n.pairs)
                  delta.unfavorable <- colSums(object@count.unfavorable)/sum(object@n.pairs)
                  indexC <- attr(object@level.treatment,"indexC")
                  indexT <- attr(object@level.treatment,"indexT")

                  ## remove scaling 
                  object.iid$favorable[indexC,] <- length(indexC) * object.iid$favorable[indexC,]
                  object.iid$favorable[indexT,] <- length(indexT) * object.iid$favorable[indexT,]

                  object.iid$unfavorable[indexC,] <- length(indexC) * object.iid$unfavorable[indexC,]
                  object.iid$unfavorable[indexT,] <- length(indexT) * object.iid$unfavorable[indexT,]

                  ## remove centering
                  object.iid$unfavorable <- sweep(object.iid$unfavorable, MARGIN = 2, FUN = "+", STATS = cumsum(delta.unfavorable*object@weight))
                  object.iid$favorable <- sweep(object.iid$favorable, MARGIN = 2, FUN = "+", STATS = cumsum(delta.favorable*object@weight))
              }
              ## ** accumulate H-decomposition
              if(is.null(endpoint)){                  
                  ## iid decomposition over all endpoints
                  object.iid <- do.call(cbind,lapply(object.iid, function(iI){
                      iIID <- iI[, NCOL(iI)]
                      if(!is.null(cluster)){
                          iIID <- tapply(iIID,cluster,sum)
                      }
                      return(iIID)
                  }))
              }else{
                  ## iid decomposition for each endpoint
                  object.iid <- lapply(endpoint, function(iE){
                      iIID <- cbind(favorable = object.iid$favorable[,iE],
                                    unfavorable = object.iid$unfavorable[,iE])
                      if(!is.null(cluster)){
                          iIID <- cbind(favorable = tapply(iIID[,"favorable"],cluster,sum),
                                       unfavorable = tapply(iIID[,"unfavorable"],cluster,sum))
                      }
                      return(iIID)
                  })
                  names(object.iid) <- endpoint
              }

              ## ** output H-decomposition
              return(object.iid)
    
})

######################################################################
### BuyseRes-iid.R ends here
