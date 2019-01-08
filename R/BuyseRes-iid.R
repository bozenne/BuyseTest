### BuyseRes-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2019 (11:20) 
## Version: 
## Last-Updated: jan  8 2019 (10:16) 
##           By: Brice Ozenne
##     Update #: 24
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
                                endpoint = NULL,
                                order = 1){


              ## ** check arguments
              valid.endpoint <- object@endpoint
              validInteger(order, valid.length = 1, valid.values = 1:2, refuse.NULL = FALSE)
              validCharacter(endpoint, valid.length = 1, valid.values = valid.endpoint, refuse.NULL = FALSE)

              ## ** extract H-decomposition
              object.iid <- object@iid
              if(object@method.inference != "asymptotic"){
                  stop("No H-decomposition in the object \n",
                       "Set the argument \'method.inference\' to \"asymptotic\" when calling BuyseTest \n")
              }

              ## ** accumulate H-decomposition
              if(is.null(order) ||order == 1){
                  if(is.null(endpoint)){
                      ## iid decomposition over all endpoints
                      object.iid$iid1 <- apply(object.iid$iid1, MARGIN = c(1,3), sum)
                  }else{
                      ## iid decomposition for each endpoint
                      object.iid$iid1 <- object.iid$iid1[,endpoint,,drop=FALSE]
                  }
              }
              if(is.null(order) ||order == 2){
                  if(is.null(endpoint)){
                      ## iid decomposition over all endpoints
                      object.iid$iid2 <- apply(object.iid$iid2, MARGIN = c(1,3), sum)
                  }else{
                      ## iid decomposition for each endpoint
                      object.iid$iid2 <- object.iid$iid2[,endpoint,,drop=FALSE]
                  }
              }

              ## ** output H-decomposition
              if(is.null(order)){
                  return(object.iid)
              }else if(order == 1){
                  return(object.iid$iid1)
              }else if(order == 2){
                  return(object.iid$iid2)
              }              

    
})

######################################################################
### BuyseRes-iid.R ends here
