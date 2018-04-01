## * Documentation - BuyseTest.option
#' @title Global options for BuyseTest package
#' @name BuyseTest.option
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
#' ## see all global parameters
#' BuyseTest.option()
#' 
#' ## see some of the global parameters
#' BuyseTest.option("n.permutation", "trace")
#' 
#' ## update some of the global parameters
#' BuyseTest.option(n.permutation = 10, trace = 1)
#' BuyseTest.option("n.permutation", "trace")
#' 
#' ## reinitialise all global parameters
#' BuyseTest.option(reinitialise = TRUE)

## * Function BuyseTest.option
#' @rdname BuyseTest.option
#' @export
BuyseTest.option <- function(..., reinitialise = FALSE){
  
  if (reinitialise == TRUE) {
    assign(".BuyseTest-option", 
           new("BuyseTest.option",
               conf.level = 0.95,
               cpus = 1,
               keep.permutation = TRUE,
               method = "Peron",
               n.permutation = 1000,
               neutralAsUninf = TRUE,
               keepComparison = FALSE,
               trace = 3,
               seed = 10,
               statistic = "netChance"), 
           envir = BuyseTest.env)
    
    return(invisible(get(".BuyseTest-option", envir = BuyseTest.env)))
    
  }else{
    
    args <- list(...)
    object <- get(".BuyseTest-option", envir = BuyseTest.env)
    
    if (!is.null(names(args))) { # write
      
      validCharacter(names(args),
                     name1 = "...",
                     valid.length = NULL,
                     valid.values = slotNames(object),
                     refuse.duplicates = TRUE,
                     refuse.NULL = FALSE,
                     method = "BuyseTest.option")
      
      value <- alloc(object, field = args)
      
      assign(".BuyseTest-option", 
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
                     method = "BuyseTest.option")
      value <- select(object, name.field = unlist(args))
      return(value)
    }
    
  }
  
  
}
