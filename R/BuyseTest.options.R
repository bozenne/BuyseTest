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
#' @details It only affects the \code{\link{BuyseTest}} function
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
                   alternative = "two.sided",
                   args.model.tte = list(),
                   check = TRUE,
                   conf.level = 0.95,
                   cpus = 1,
                   debug = -1,
                   engine = "GPC2_cpp",
                   fitter.model.tte = "prodlim",
                   hierarchical = TRUE,
                   keep.pairScore = FALSE,
                   keep.survival = FALSE,
                   method.inference = "u-statistic",
                   scoring.rule = "Peron",               
                   correction.uninf = 0,               
                   n.resampling = 1000,
                   strata.resampling = as.character(NA),
                   neutral.as.uninf = TRUE,
                   order.Hprojection = 1,
                   precompute = TRUE,
                   print.display = c("endpoint","threshold","delta","Delta"),
                   statistic = "netBenefit",               
                   summary.display = list(c("endpoint","threshold","weight","strata","total","favorable","unfavorable","neutral","uninf","delta","Delta","CI","p.value","significance"),
                                          c("endpoint","threshold","weight","strata","favorable","unfavorable","delta","Delta","Delta(%)","information(%)")),
                   transformation = TRUE,
                   trace = 2), 
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
