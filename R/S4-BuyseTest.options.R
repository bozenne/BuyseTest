## * Documentation BuyseTest.options
#' @title Class "BuyseTest.options" (global setting for the BuyseTest package)
#' @name BuyseTest.options-class
#' @include 1-setGeneric.R
#' 
#' @description Class defining the global settings for the BuyseTest package.
#' 
#' @seealso 
#' \code{\link{BuyseTest.options}} to select or update global settings.
#' 
#' @keywords classes options BuyseTest.options-class
#' @author Brice Ozenne

## * Class BuyseTest.options
#' @rdname BuyseTest.options-class
setClass(
  Class = "BuyseTest.options",
  
  representation(
      alternative = "character",
      args.model.tte = "list",
      check = "logical",
      conf.level = "numeric",
      correction.uninf = "numeric",
      cpus = "numeric",
      debug = "numeric",
      engine = "character",
      fitter.model.tte = "character",
      hierarchical = "logical",
      keep.pairScore = "logical",
      keep.survival = "logical",      
      method.inference = "character",
      scoring.rule = "character",
      n.resampling = "numeric",
      strata.resampling = "character",
      neutral.as.uninf = "logical",
      add.halfNeutral = "logical",
      add.1.pperm = "logical",
      order.Hprojection = "numeric",
      pool.strata = "character",
      precompute = "logical",
      print.display = "character",
      statistic = "character",
      summary.display = "list",
      trace = "numeric",
      transformation = "logical",
      warning.correction = "numeric"
  ),

### ** Check validity of the object
  validity = function(object){
      validNames.summary <- c("endpoint","restriction","threshold","weight","strata","total","favorable","unfavorable","neutral","uninf","information(%)",
                              "delta","Delta","Delta(%)",
                              "p.value","CI","significance")
      
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
      validInteger(object@correction.uninf,
                   name1 = "@correction.uninf",
                   min = 0,
                   max = 2,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validInteger(object@cpus,
                   name1 = "@cpus",
                   min = 1,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validInteger(object@debug,
                   name1 = "@debug",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@engine,
                     name1 = "@engine",
                     valid.values = c("GPC_cpp","GPC2_cpp"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validCharacter(object@fitter.model.tte,
                     name1 = "@fitter.model.tte",
                     valid.values = c("prodlim", "survreg"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validLogical(object@hierarchical,
                   name1 = "@hierarchical",
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
                     valid.values = c("bootstrap", "stratified bootstrap", "studentized bootstrap", "studentized stratified bootstrap",
                                      "permutation", "stratified permutation",
                                      "none", "u-statistic", "u-statistic-bebu"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validCharacter(object@scoring.rule,
                     name1 = "@scoring.rule",
                     valid.values = c("Gehan","Peron"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validInteger(object@n.resampling,
                   name1 = "@n.resampling",
                   min = 0,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@strata.resampling,
                     name1 = "@n.resampling",
                     valid.values = c(as.character(NA),"treatment","strata"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validLogical(object@neutral.as.uninf,
                   name1 = "@neutral.as.uninf",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validLogical(object@add.halfNeutral,
                   name1 = "@add.halfNeutral",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validInteger(object@order.Hprojection,
                   name1 = "@order.Hprojection",
                   min = 1,
                   max = 2,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@pool.strata,
                     name1 = "@pool.strata",
                     valid.values = c("Buyse","CMH","equal","var-favorable","var-unfavorable","var-netBenefit","var-winRatio"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validLogical(object@precompute,
                   name1 = "@precompute",
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      validCharacter(object@print.display,
                     name1 = "@print.display",
                     valid.values = validNames.summary,
                     valid.length = NULL,
                     method = "Class BuyseTest.options")
      validCharacter(object@statistic,
                     name1 = "@statistic",
                     valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      lapply(object@summary.display,validCharacter,
             name1 = "@summary.display",
             valid.values = validNames.summary,
             valid.length = NULL,
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
      validNumeric(object@warning.correction,
                   name1 = "@warning.correction",
                   min = 0,
                   max = 1,
                   valid.length = 1,
                   method = "Class BuyseTest.options")
      return(TRUE)} 
)

#' @title Methods for the class "BuyseTest.options" 
#' @name BuyseTest.options-methods
#' @aliases alloc,BuyseTest.options-method
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
              validObject(object)
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
            
            ls.slots <- stats::setNames(vector(mode = "list", length = n.field), name.field)
            for (iField in 1:n.field) {
              ls.slots[[iField]] <- slot(object, name.field[iField])
            }
            
            return(ls.slots)
          }
)

