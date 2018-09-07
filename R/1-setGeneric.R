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
#' @rdname BuyseRes-getIndividualScore
#' @exportMethod getIndividualScore
setGeneric(name = "getIndividualScore",
           def = function(object, ...){standardGeneric("getIndividualScore")}
)

## * Selector (for BuyseRes)
#' @rdname BuyseRes-getSurvival
#' @exportMethod getSurvival
setGeneric(name = "getSurvival",
           def = function(object, ...){standardGeneric("getSurvival")}
)
