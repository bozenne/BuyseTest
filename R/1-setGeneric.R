## * Allocator alloc (for BuyseTest-options)
setGeneric(name = "alloc", 
           def = function(object, ...){standardGeneric("alloc")}
)

## * Selector select (for BuyseTest-options)
setGeneric(name = "select", 
           def = function(object, ...){standardGeneric("select")}
)

## * Selector getCount (for S4BuyseTest)
setGeneric(name = "getCount",
           def = function(object, type){standardGeneric("getCount")}
)

## * Selector getPairScore (for S4BuyseTest)
setGeneric(name = "getPairScore",
           def = function(object, endpoint = NULL, strata = NULL, cumulative = FALSE,
                          rm.withinStrata = TRUE, rm.strata = is.na(object@strata),
                          rm.indexPair = TRUE, rm.weight = FALSE, rm.corrected = (object@correction.uninf==0),
                          unlist = TRUE, trace = 1){
               standardGeneric("getPairScore")
           }
           )

## * Selector getPseudovalue (for S4BuyseTest)
setGeneric(name = "getPseudovalue",
           def = function(object, statistic = NULL, endpoint = NULL){
               standardGeneric("getPseudovalue")
           }
)

## * Selector getSurvival (for S4BuyseTest)
setGeneric(name = "getSurvival",
           def = function(object, type = NULL, endpoint = NULL, strata = NULL, unlist = TRUE, trace = TRUE){
               standardGeneric("getSurvival")
           }
)

## * Selector getIid (for S4BuyseTest)
setGeneric(name = "getIid",
           def = function(object, endpoint = NULL, statistic = NULL, strata = FALSE, cumulative = TRUE, center = TRUE, scale = TRUE, type = "all", cluster = NULL, simplify = FALSE){
               standardGeneric("getIid")
           }
)

## * method sensitivity (for BuyseTest)
setGeneric(name = "sensitivity", 
           def = function(object, ...){standardGeneric("sensitivity")}
)
