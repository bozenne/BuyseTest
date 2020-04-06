## * Allocator alloc (for BuyseTest-options)
setGeneric(name = "alloc", 
           def = function(object, ...){standardGeneric("alloc")}
)

## * Selector select (for BuyseTest-options)
setGeneric(name = "select", 
           def = function(object, ...){standardGeneric("select")}
)

## * Selector getCount (for S4BuyseTest)
#' @rdname S4BuyseTest-getCount
#' @exportMethod getCount
setGeneric(name = "getCount",
           def = function(object, type){standardGeneric("getCount")}
)

## * Selector getPairScore (for S4BuyseTest)
#' @rdname S4BuyseTest-getPairScore
#' @exportMethod getPairScore
setGeneric(name = "getPairScore",
           def = function(object, endpoint = NULL, strata = NULL, sum = FALSE,
                          rm.withinStrata = TRUE, rm.strata = is.na(object@strata),
                          rm.indexPair = TRUE, rm.weight = FALSE, rm.corrected = (object@correction.uninf==0),
                          unlist = TRUE, trace = 1){
               standardGeneric("getPairScore")
           }
           )

## * Selector getSurvival (for S4BuyseTest)
#' @rdname S4BuyseTest-getSurvival
#' @exportMethod getSurvival
setGeneric(name = "getSurvival",
           def = function(object, type = NULL, endpoint = NULL, strata = NULL, unlist = TRUE, trace = TRUE){
               standardGeneric("getSurvival")
           }
)

## * Selector getIid (for S4BuyseTest)
#' @rdname S4BuyseTest-getIid
#' @exportMethod getIid
setGeneric(name = "getIid",
           def = function(object, endpoint = NULL, normalize = TRUE, type = "all", cluster = NULL){
               standardGeneric("getIid")
           }
)

