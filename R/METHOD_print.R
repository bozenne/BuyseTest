#' @docType methods
#' @name BuyseRes-show
#' @title Show Method for Class "BuyseRes"
#' @aliases show show,BuyseRes-method
#' @include OBJECT_BuyseTest.R METHOD_summary.R
#' 
#' @description Display the main results stored in a \code{\link{BuyseRes}} object.
#' 
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' 
#' @seealso 
#'   \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#'   \code{\link{BuyseRes-summary}} for a more detailed presentation of the \code{BuyseRes} object.
#' 
#' @examples
#'   n.Treatment_testBin <- 500
#'   n.Control_testBin <- 500
#'   prob.Treatment_testBin <- c(0.5,0.75)
#'   prob.Control_testBin <- c(0.5,0.25)
#'   
#'   set.seed(10)
#'   data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
#'   data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
#'   data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
#'                               rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
#'   data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)
#'   
#'   #### no strata
#'   \dontrun{
#'     BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                   treatment="treatment",type=c("bin","bin"),n.bootstrap=10000)
#'   }
#'   \dontshow{
#'     BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
#'                                   treatment="treatment",type=c("bin","bin"),
#'                                   n.bootstrap=10,trace=0)
#'   }
#'   
#'   BuyseTest_object
#' 
#' @keywords summary BuyseRes-method

#' @exportMethod show
setMethod(f = "show",
          signature = "BuyseRes",
          definition = function(object){
            
           table <- summary(object, strata="global", show = NULL)$pc
           table$threshold[is.na(table$threshold)] <- ""
           print(table, row.names = FALSE)
           
           return(invisible(NULL))
          }
)
