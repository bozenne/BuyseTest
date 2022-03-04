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
                   alternative = "two.sided", ## type of alternative hypothesis when doing hypothesis testing: less, greater, two.sided
                   args.model.tte = list(), ## additional argument passed to prodlim when fitting the survival model in BuyseTest() --> calcPeron
                   check = TRUE, ## should arguments be checked when running BuyseTest()
                   conf.level = 0.95, ## coverage of confidence intervals 
                   cpus = 1, ## cpus used to performance inference via resampling in BuyseTest()
                   debug = -1, ## hidden argument in BuyseTest to display progress of the C++ code
                   engine = "GPC2_cpp", ## C++ function used to perform GPC calculation for BuyseTest()
                   fitter.model.tte = "prodlim", ## survival model in BuyseTest() --> calcPeron
                   hierarchical = TRUE, ## default value of argument hierarchical in BuyseTest()
                   keep.pairScore = FALSE, ## default value of argument keep.pairScore in BuyseTest()
                   keep.survival = FALSE, ## hidden argument to export survival values for the Peron Scoring rule in BuyseTest()
                   method.inference = "u-statistic", ## default value of argument method.inference in BuyseTest()
                   scoring.rule = "Peron", ## default value of argument scoring.rule in BuyseTest()               
                   correction.uninf = 0, ## default value of argument correction.uninf in BuyseTest()               
                   n.resampling = 1000, ## default value of argument n.resampling in BuyseTest()
                   strata.resampling = as.character(NA), ## default value of argument strata.resampling in BuyseTest()
                   neutral.as.uninf = TRUE, ## default value of argument neutral.as.uninf in BuyseTest()
                   add.halfNeutral = FALSE, ## default value of argument add.halfNeutral in BuyseTest()
                   add.1.pperm = TRUE, ## if TRUE p-value are computed as (#more extreme+1)/(#perm + 1) otherwise #more exterme/#perm
                   order.Hprojection = 1, ## hidden argument in BuyseTest() to control the type of H-projection when using method.inference="u-statistic". Can be 1 or 2
                   precompute = TRUE, ## hidden argument in BuyseTest() to pre-compute integrals over time before the C++ routine
                   print.display = c("endpoint","restriction","threshold","delta","Delta"), ## what to display when showing a S4BuyseTest object 
                   statistic = "netBenefit",  ## what is the default statistic output by summary, confint ...
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
