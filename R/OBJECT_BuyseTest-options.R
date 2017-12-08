#' @name BuyseTest.options-class
#' @title Class "BuyseTest.options" (global setting for the BuyseTest package)
#' @aliases BuyseTest.options-class
#' @include METHOD_generic.R
#' 
#' @description Class defining the global settings for the BuyseTest package.
#' 
#' @inheritParams BuyseTest
#' @param conf.level the confidence level of the confidence interval
#' @param keep.bootstrap should the result of each boostrap sample be stored in BuyseRes objects?
#' 
#' @seealso 
#' \code{\link{BuyseTest.options}} to select or update global settings.#' 
#' 
#' @keywords classes options BuyseTest.options-class
setClass(
  Class = "BuyseTest.options",
  
  representation(
    conf.level = "numeric",
    cpus = "numeric",
    keep.bootstrap = "logical",
    method = "character",
    n.bootstrap = "numeric",
    neutralAsUninf = "logical",
    trace = "numeric",
    seed = "numeric",
    statistic = "character"
  ),
  
  validity = function(object){
    
    validNumeric(object@conf.level, name1 = "@conf.level", min = 0, max = 1, validLength = 1, method = "Class BuyseTest.options")
    validInteger(object@cpus, name1 = "@cpus", min = 1, validLength = 1, method = "Class BuyseTest.options")
    validLogical(object@keep.bootstrap, name1 = "@keep.bootstrap", validLength = 1, method = "Class BuyseTest.options")
    validCharacter(object@method, name1 = "@method", validValues = c("Peron","Efron","Peto","Gehan"), validLength = 1, method = "Class BuyseTest.options")
    validInteger(object@n.bootstrap, name1 = "@n.bootstrap", min = 0, validLength = 1, method = "Class BuyseTest.options")
    validLogical(object@neutralAsUninf, name1 = "@neutralAsUninf", validLength = 1, method = "Class BuyseTest.options")
    validInteger(object@trace, name1 = "@trace", min = 0, validLength = 1, method = "Class BuyseTest.options")
    validInteger(object@seed, name1 = "@seed", min = 1, validLength = 1, method = "Class BuyseTest.options")
    validCharacter(object@statistic, name1 = "@statistic", validValues = c("netChance","winRatio"), validLength = 1, method = "Class BuyseTest.options")
    return(TRUE)} 
)


#' @name BuyseTest.options-methods
#' @title Methods for the class "BuyseTest.options" 
#' @aliases BuyseTest.options-methods alloc,BuyseTest.options-method select,BuyseTest.options-method
#'
#' @description Methods to update or select global settings
#' 
#' @param object an object of class \code{BuyseTest.options}.
#' @param field a \code{list} named with the name of the fields to update and containing the values to assign to these fields
#' @param name.field a \code{character vector} containing the names of the field to be selected.

#' @rdname BuyseTest.options-methods
setMethod(f = "alloc",
          signature = "BuyseTest.options",
          definition = function(object, field){
            
            name.field <- names(field)
            n.field <- length(field)
            
            for (iter_field in 1:n.field) {
              slot(object, name.field[iter_field]) <- field[[iter_field]]
            }
            
            return(object)
          }
)

#' @rdname BuyseTest.options-methods
setMethod(f = "select",
          signature = "BuyseTest.options",
          definition = function(object, name.field){
            
            if (is.null(name.field)) {
              name.field <- slotNames(object)
            }
            n.field <- length(name.field)
            
            ls.slots <- setNames(vector(mode = "list", length = n.field), name.field)
            for (iter_field in 1:n.field) {
              ls.slots[[iter_field]] <- slot(object, name.field[iter_field])
            }
            
            return(ls.slots)
          }
)

