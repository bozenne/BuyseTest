## * Allocator (for BuyseTest-option)
setGeneric(name = "alloc", 
           def = function(object, ...){standardGeneric("alloc")}
)

## * Selector (for BuyseTest-option)
setGeneric(name = "select", 
           def = function(object, ...){standardGeneric("select")}
)

## * Selector (for BuyseRes)
#' @rdname BuyseRes-getCount
#' @exportMethod getCount
setGeneric(name = "getCount",
           def = function(object, type){standardGeneric("getCount")}
)
