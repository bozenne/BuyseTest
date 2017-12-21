## * Documentation - BuyseTest.options
#' @title Global options for BuyseTest package
#' @name BuyseTest.options
#' @include onload.R
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
#' BuyseTest.options()
#' 
#' ## see some of the global parameters
#' BuyseTest.options("n.bootstrap", "trace")
#' 
#' ## update some of the global parameters
#' BuyseTest.options(n.bootstrap = 10, trace = 1)
#' BuyseTest.options("n.bootstrap", "trace")
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
               conf.level = 0.95,
               cpus = 1,
               keep.bootstrap = TRUE,
               method = "Peron",
               n.bootstrap = 1000,
               neutralAsUninf = TRUE,
               trace = 3,
               seed = 10,
               statistic = "netChance"), 
           envir = BuyseTest.env)
    
    return(invisible(get(".BuyseTest-options", envir = BuyseTest.env)))
    
  }else{
    
    args <- list(...)
    object <- get(".BuyseTest-options", envir = BuyseTest.env)
    
    if (!is.null(names(args))) { # write
      
      validCharacter(names(args), name1 = "...", validLength = NULL, validValues = slotNames(object), refuse.duplicates = TRUE, refuse.NULL = FALSE, method = "BuyseTest.options")
      
      value <- alloc(object, field = args)
      
      assign(".BuyseTest-options", 
             value, 
             envir = BuyseTest.env)
      
      return(invisible(value))
      
    } else {# read
      
      validCharacter(args, name1 = "...", validLength = NULL, validValues = slotNames(object), refuse.duplicates = TRUE, refuse.NULL = FALSE, method = "BuyseTest.options")
      value <- select(object, name.field = unlist(args))
      return(value)
    }
    
  }
  
  
}
