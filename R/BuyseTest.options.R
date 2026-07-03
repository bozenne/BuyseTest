## * Documentation - BuyseTest.options
#' @title Global options for BuyseTest package
#' @name BuyseTest.options
#' @include 0-onLoad.R
#'
#' @description Update or select global options for the BuyseTest package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'
#' @keywords utilities
#' 
#' @examples
#' library(data.table)
#' 
#' ## see all global parameters
#' BuyseTest.options()
#' 
#' ## see some of the global parameters
#' BuyseTest.options("n.resampling", "trace")
#' 
#' ## update some of the global parameters
#' BuyseTest.options(n.resampling = 10, trace = 1)
#' BuyseTest.options("n.resampling", "trace")
#' 
#' ## reinitialise all global parameters
#' BuyseTest.options(reinitialise = TRUE)

## * Function BuyseTest.options
#' @rdname BuyseTest.options
#' @export
BuyseTest.options <- function(..., reinitialise = FALSE){
  
    if (reinitialise == TRUE) {
        assign(".BuyseTest-options", 
               new("BuyseTest.options",
                   add.halfNeutral = FALSE, ## default value of argument add.halfNeutral in BuyseTest()
                   add.1.presample = TRUE, ## if TRUE p-value when using permutation is computed as (#more extreme+1)/(#sample+1) otherwise #more exterme/#sample.
                                           ## Another correction is used for the bootstrap to ensure strictly positive p-values
                   alternative = "two.sided", ## type of alternative hypothesis when doing hypothesis testing: less, greater, two.sided
                   args.model.tte = list(), ## additional argument passed to prodlim when fitting the survival model in BuyseTest() --> calcPeron
                   check = TRUE, ## should arguments be checked when running BuyseTest()
                   conf.level = 0.95, ## coverage of confidence intervals 
                   correction.uninf = 0, ## default value of argument correction.uninf in BuyseTest()               
                   cpus = 1, ## cpus used to performance inference via resampling in BuyseTest()
                   debug = -1, ## hidden argument in BuyseTest to display progress of the C++ code
                   engine = "GPC2_cpp", ## C++ function used to perform GPC calculation for BuyseTest()
                   hierarchical = TRUE, ## default value of argument hierarchical in BuyseTest()
                   keep.pairScore = FALSE, ## default value of argument keep.pairScore in BuyseTest()
                   keep.survival = FALSE, ## hidden argument to export survival values for the Peron Scoring rule in BuyseTest()
                   method.inference = "u-statistic", ## default value of argument method.inference in BuyseTest()
                   n.grid = 1000, ## default value of argument n.grid in BuyseTTEM.survreg() when discretising parametric survival curves into piecewise constant curves
                   n.resampling = 1000, ## default value of argument n.resampling in BuyseTest()
                   neutral.as.uninf = TRUE, ## default value of argument neutral.as.uninf in BuyseTest()
                   order.Hprojection = 1, ## hidden argument in BuyseTest() to control the type of H-projection when using method.inference="u-statistic". Can be 1 or 2
                   pool.strata = "CMH", ## default weighting scheme to pool estimates across strata in BuyseTest(). Can be "Buyse" (weighted average proportional to number of pairs per strata), "CMH" (Cochran-Mantel-Haenszel weights), "equal" (equal weights for all strata), "standardization" (average over all possible combinations of strata), "var-netbenefit" (optimal weight w.r.t. net benefit)                   
                   precompute = TRUE, ## hidden argument in BuyseTest() to pre-compute integrals over time before the C++ routine
                   print.display = c("endpoint","restriction","threshold","delta","Delta"), ## what to display when showing a S4BuyseTest object 
                   scoring.rule = "Peron", ## default value of argument scoring.rule in BuyseTest()               
                   statistic = "netBenefit",  ## what is the default statistic output by summary, confint ...
                   strata.resampling = as.character(NA), ## default value of argument strata.resampling in BuyseTest()
                   summary.display = list(c("endpoint","restriction","threshold","weight","strata","total","favorable","unfavorable","neutral","uninf","delta","Delta","CI","p.value","significance"),
                                          c("endpoint","restriction","threshold","weight","strata","favorable","unfavorable","delta","Delta","Delta(%)","information(%)")),
                   transformation = TRUE, ## should p-value/CI be computed after transformation (and appropriate backtransformation)
                   trace = 2, ## default value of argument trace in BuyseTest()
                   warning.correction = 0.25), ## display a warning when the correction lead to a very large change in estimate
               envir = BuyseTest.env)
    
    return(invisible(get(".BuyseTest-options", envir = BuyseTest.env)))
    
  }else{
    
    args <- list(...)
    object <- get(".BuyseTest-options", envir = BuyseTest.env)
    
    if (!is.null(names(args))) { # write
      validCharacter(names(args),
                     name1 = "...",
                     valid.length = NULL,
                     valid.values = slotNames(object),
                     refuse.duplicates = TRUE,
                     refuse.NULL = FALSE,
                     method = "BuyseTest.options")
      
      value <- alloc(object, field = args)
      
      assign(".BuyseTest-options", 
             value, 
             envir = BuyseTest.env)
      
      return(invisible(value))
      
    } else {# read
      
      validCharacter(args,
                     name1 = "...",
                     valid.length = NULL,
                     valid.values = slotNames(object),
                     refuse.duplicates = TRUE,
                     refuse.NULL = FALSE,
                     method = "BuyseTest.options")
      value <- select(object, name.field = unlist(args))
      return(value)
    }
    
  }
  
  
}

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
#' @keywords classes
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
      hierarchical = "logical",
      keep.pairScore = "logical",
      keep.survival = "logical",      
      method.inference = "character",
      scoring.rule = "character",
      n.grid = "numeric",
      n.resampling = "numeric",
      strata.resampling = "character",
      neutral.as.uninf = "logical",
      add.halfNeutral = "logical",
      add.1.presample = "logical",
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
      validCharacter(gsub("-"," ",tolower(object@method.inference), fixed = TRUE),
                     name1 = "@resampling",
                     valid.values = c("bootstrap", "stratified bootstrap", "studentized bootstrap", "studentized stratified bootstrap",
                                      "none",
                                      "permutation", "stratified permutation", "varexact permutation",
                                      "u statistic", "u statistic bebu"),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validCharacter(object@scoring.rule,
                     name1 = "@scoring.rule",
                     valid.values = c("Gehan","Peron","Efron","Latta",names(survival::survreg.distributions)),
                     valid.length = 1,
                     method = "Class BuyseTest.options")
      validInteger(object@n.grid,
                   name1 = "@n.grid",
                   min = 0,
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

