### BuyseRes-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2019 (11:20) 
## Version: 
## Last-Updated: mar 29 2019 (11:46) 
##           By: Brice Ozenne
##     Update #: 36
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
#' @param order [integer 1,2] which order of the H-decomposition should be output?
#' Can be the first or second order. \code{NULL} output a list containing all terms.
#'  
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-summary}} for a more detailed presentation of the \code{BuyseRes} object.
#' 
#' @keywords iid BuyseRes-method

## * Method - iid
#' @rdname BuyseRes-iid
#' @exportMethod iid
setMethod(f = "iid",
          signature = "BuyseRes",
          definition = function(object,
                                endpoint = NULL){


              ## ** check arguments
              valid.endpoint <- object@endpoint
              validCharacter(endpoint, valid.length = 1:length(valid.endpoint), valid.values = valid.endpoint, refuse.NULL = FALSE)

              ## ** extract H-decomposition
              object.iid <- object@iid
              if(object@method.inference != "u-statistic"){
                  stop("No H-decomposition in the object \n",
                       "Set the argument \'method.inference\' to \"u-statistic\" when calling BuyseTest \n")
              }

              ## ** accumulate H-decomposition
              if(is.null(endpoint)){                  
                  ## iid decomposition over all endpoints
                  object.iid <- do.call(cbind,lapply(object.iid, function(iI){iI[, NCOL(iI)]}))
              }else{
                  ## iid decomposition for each endpoint
                  object.iid <- lapply(object.iid, function(iIID){iIID[,endpoint,drop=FALSE]})
              }

              ## ** output H-decomposition
              return(object.iid)
    
})

######################################################################
### BuyseRes-iid.R ends here
