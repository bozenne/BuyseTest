#' @docType methods
#' @name BuyseRes-getCount
#' @title get Method for Class "BuyseRes"
#' @aliases getCount BuyseRes-get 
#' @include OBJECT_BuyseTest.R
#' 
#' @description Extract the number of pairs.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param type the type of pairs to be counted. Can be \code{"favorable"}, \code{"unfavorable"}, \code{neutral}, or \code{uninf}. Can also be \code{"all"} to select all of them.
#' 
#' @return 
#'   A \code{"vector"} containing the number of pairs
#' 
#' @examples
#' dt <- simulBT(1e2)
#' BT <- BuyseTest(data=dt,endpoint="Y_TTE1",treatment="Treatment",type="timeToEvent",censoring="event1", n.bootstrap = 0)
#' getCount(BT)
#' getCount(BT, type = "favorable")
#'
#' @keywords getCount BuyseRes-method

#' @rdname BuyseRes-getCount
setGeneric(name = "getCount", 
           def = function(object, type){standardGeneric("getCount")}
)

#' @rdname BuyseRes-getCount
#' @exportMethod getCount
setMethod(f = "getCount",
          signature = "BuyseRes",
          definition = function(object, type){
            
            if(missing(type)){
              type <- c("favorable","unfavorable","neutral","uninf")
            }
            
            validCharacter(type, validLength = NULL, validValues = c("favorable","unfavorable","neutral","uninf"), method = "getCount")
            
            out <- NULL
            if("favorable" %in% type){out <- c(out, favorable = object@count_favorable)}
            if("unfavorable" %in% type){out <- c(out, unfavorable = object@count_unfavorable)}
            if("neutral" %in% type){out <- c(out, neutral = object@count_neutral)}
            if("uninf" %in% type){out <- c(out, uninf = object@count_uninf)}
            
            return(out)
            }
)