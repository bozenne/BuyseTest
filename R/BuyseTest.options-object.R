## * Documentation BuyseTest.options
#' @title Class "BuyseTest.options" (global setting for the BuyseTest package)
#' @name BuyseTest.options-class
#' @aliases BuyseTest.options-class
#' @include 1-setGeneric.R
#' 
#' @description Class defining the global settings for the BuyseTest package.
#' 
#' @inheritParams BuyseTest
#' 
#' @seealso 
#' \code{\link{BuyseTest.options}} to select or update global settings.
#' 
#' @keywords classes options BuyseTest.options-class

## * Class BuyseTest.options
#' @rdname BuyseTest.options-class
setClass(
  Class = "BuyseTest.options",
  
  representation(
      alternative = "character",
      check = "logical",
      conf.level = "numeric",
      continuity.correction = "logical",
      correction.uninf = "logical",
      cpus = "numeric",
      keep.pairScore = "logical",
      keep.survival = "logical",
      method.inference = "character",
      method.tte = "character",
      n.resampling = "numeric",
      neutral.as.uninf = "logical",
      order.Hprojection = "numeric",
      statistic = "character",
      trace = "numeric",
      transformation = "logical"
  ),

### ** Check validity of the object
  validity = function(object){
      validCharacter(object@alternative,
                     name1 = "@alternative",
                     valid.values = c("two.sided","greater","less"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validLogical(object@check,
                   name1 = "@check",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validNumeric(object@conf.level,
                   name1 = "@conf.level",
                   min = 0,
                   max = 1,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@correction.uninf,
                   name1 = "@correction.uninf",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validInteger(object@cpus,
                   name1 = "@cpus",
                   min = 1,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@keep.pairScore,
                   name1 = "@keep.pairScore",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@keep.survival,
                   name1 = "@keep.survival",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@method.inference,
                     name1 = "@resampling",
                     valid.values = c("bootstrap", "stratified boostrata",
                                      "permutation", "stratified permutation",
                                      "asymptotic"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validCharacter(object@method.tte,
                     name1 = "@method.tte",
                     valid.values = c("Gehan","Peron"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validInteger(object@n.resampling,
                   name1 = "@n.resampling",
                   min = 0,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@neutral.as.uninf,
                   name1 = "@neutral.as.uninf",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validInteger(object@order.Hprojection,
                   name1 = "@order.Hprojection",
                   min = 1,
                   max = 2,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@statistic,
                     name1 = "@statistic",
                     valid.values = c("netBenefit","winRatio"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validInteger(object@trace,
                   name1 = "@trace",
                   min = 0, max = 2,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@transformation,
                   name1 = "@transformation",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      return(TRUE)} 
)

#' @title Methods for the class "BuyseTest.options" 
#' @name BuyseTest.options-methods
#' @aliases BuyseTest.options-methods alloc,BuyseTest.options-method select,BuyseTest.options-method
#'
#' @description Methods to update or select global settings
#' 
#' @param object an object of class \code{BuyseTest.options}.
#' @param field a \code{list} named with the name of the fields to update and containing the values to assign to these fields
#' @param name.field a \code{character vector} containing the names of the field to be selected.

## * Alloc BuyseTest.options
#' @rdname BuyseTest.options-methods
setMethod(f = "alloc",
          signature = "BuyseTest.options",
          definition = function(object, field){
            
            name.field <- names(field)
            n.field <- length(field)
            
            for (iField in 1:n.field) {
              slot(object, name.field[iField]) <- field[[iField]]
            }
            
            return(object)
          }
)

## * Select BuyseTest.options
#' @rdname BuyseTest.options-methods
setMethod(f = "select",
          signature = "BuyseTest.options",
          definition = function(object, name.field){
            
            if (is.null(name.field)) {
              name.field <- slotNames(object)
            }
            n.field <- length(name.field)
            
            ls.slots <- setNames(vector(mode = "list", length = n.field), name.field)
            for (iField in 1:n.field) {
              ls.slots[[iField]] <- slot(object, name.field[iField])
            }
            
            return(ls.slots)
          }
)

