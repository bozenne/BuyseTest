## * Documentation BuyseTest.option
#' @title Class "BuyseTest.option" (global setting for the BuyseTest package)
#' @name BuyseTest.option-class
#' @aliases BuyseTest.option-class
#' @include 1-setGeneric.R
#' 
#' @description Class defining the global settings for the BuyseTest package.
#' 
#' @inheritParams BuyseTest
#' @param conf.level the confidence level of the confidence interval
#' @param keep.permutation should the result of each sample from the permutation test be stored in BuyseRes objects?
#' 
#' @seealso 
#' \code{\link{BuyseTest.option}} to select or update global settings.
#' 
#' @keywords classes options BuyseTest.option-class

## * Class BuyseTest.option
#' @rdname BuyseTest.option-class
setClass(
  Class = "BuyseTest.option",
  
  representation(
    conf.level = "numeric",
    cpus = "numeric",
    keep.permutation = "logical",
    method = "character",
    n.permutation = "numeric",
    neutralAsUninf = "logical",
    keepComparison = "logical",
    trace = "numeric",
    seed = "numeric",
    statistic = "character"
  ),

  ### ** Check validity of the object
  validity = function(object){
    validNumeric(object@conf.level,
                 name1 = "@conf.level",
                 min = 0,
                 max = 1,
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validInteger(object@cpus,
                 name1 = "@cpus",
                 min = 1,
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validLogical(object@keep.permutation,
                 name1 = "@keep.permutation",
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validCharacter(object@method,
                   name1 = "@method",
                   valid.values = c("Peron","Efron","Peto","Gehan"),
                   valid.length = 1,
                   method = "Class BuyseTest.option")
    validCharacter(object@method,
                   name1 = "@method",
                   valid.values = c("Peron","Efron","Peto","Gehan"),
                   valid.length = 1,
                   method = "Class BuyseTest.option")
    validInteger(object@n.permutation,
                 name1 = "@n.permutation",
                 min = 0,
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validLogical(object@neutralAsUninf,
                 name1 = "@neutralAsUninf",
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validLogical(object@keepComparison,
                 name1 = "@keepComparison",
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validInteger(object@trace,
                 name1 = "@trace",
                 min = 0,
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validInteger(object@seed,
                 name1 = "@seed",
                 min = 1,
                 valid.length = 1,
                 method = "Class BuyseTest.option")
    validCharacter(object@statistic,
                   name1 = "@statistic",
                   valid.values = c("netChance","winRatio"),
                   valid.length = 1,
                   method = "Class BuyseTest.option")
    return(TRUE)} 
)

#' @title Methods for the class "BuyseTest.option" 
#' @name BuyseTest.option-methods
#' @aliases BuyseTest.option-methods alloc,BuyseTest.option-method select,BuyseTest.option-method
#'
#' @description Methods to update or select global settings
#' 
#' @param object an object of class \code{BuyseTest.option}.
#' @param field a \code{list} named with the name of the fields to update and containing the values to assign to these fields
#' @param name.field a \code{character vector} containing the names of the field to be selected.

## * Alloc BuyseTest.option
#' @rdname BuyseTest.option-methods
setMethod(f = "alloc",
          signature = "BuyseTest.option",
          definition = function(object, field){
            
            name.field <- names(field)
            n.field <- length(field)
            
            for (iter_field in 1:n.field) {
              slot(object, name.field[iter_field]) <- field[[iter_field]]
            }
            
            return(object)
          }
)

## * Select BuyseTest.option
#' @rdname BuyseTest.option-methods
setMethod(f = "select",
          signature = "BuyseTest.option",
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

