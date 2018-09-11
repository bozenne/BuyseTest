## * Allocator (for BuyseTest-options)
setGeneric(name = "alloc", 
           def = function(object, ...){standardGeneric("alloc")}
)

## * Selector (for BuyseTest-options)
setGeneric(name = "select", 
           def = function(object, ...){standardGeneric("select")}
)

## * Selector (for BuyseRes)
#' @rdname BuyseRes-getCount
#' @exportMethod getCount
setGeneric(name = "getCount",
           def = function(object, type){standardGeneric("getCount")}
)

## * Selector (for BuyseRes)
#' @rdname BuyseRes-getPairScore
#' @exportMethod getPairScore
setGeneric(name = "getPairScore",
           def = function(object, ...){standardGeneric("getPairScore")}
)

## * Selector (for BuyseRes)
#' @rdname BuyseRes-getSurvival
#' @exportMethod getSurvival
setGeneric(name = "getSurvival",
           def = function(object, ...){standardGeneric("getSurvival")}
)
