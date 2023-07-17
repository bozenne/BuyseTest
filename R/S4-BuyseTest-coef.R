### S4BuyseTest-coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2019 (10:45) 
## Version: 
## Last-Updated: jul 14 2023 (23:58) 
##           By: Brice Ozenne
##     Update #: 299
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - coef
#' @docType methods
#' @name S4BuyseTest-coef
#' @title Extract Summary Statistics from GPC
#' @aliases coef,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'  
#' @description Extract summary statistics (net benefit, win ratio, ...) from GPC.
#' 
#' @param object a \code{S4BuyseTest} object, output of \code{\link{BuyseTest}}.
#' @param statistic [character] the type of summary statistic. See the detail section.
#' @param endpoint [character] for which endpoint(s) the summary statistic should be output?
#' If \code{NULL} returns the summary statistic for all endpoints.
#' @param strata [character vector] the name of the strata to be displayed.
#' Can also be \code{"global"} or \code{FALSE} to display the statistic pooled over all strata,
#' or \code{TRUE} to display each strata-specific statistic.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param resampling [logical] should the summary statistic obtained by resampling be output?
#' @param simplify [logical] should the result be coerced to the lowest possible dimension?
#' @param ... ignored.
#'
#' @details
#' One of the following statistic can be specified:
#' \itemize{
#' \item \code{"netBenefit"}: returns the net benefit.
#' \item \code{"winRatio"}: returns the win ratio.
#' \item \code{"favorable"}: returns the proportion in favor of the treatment (also called Mann-Whitney parameter).
#' \item \code{"unfavorable"}: returns the proportion in favor of the control.
#' \item \code{"unfavorable"}: returns the proportion of neutral pairs.
#' \item \code{"unfavorable"}: returns the proportion of uninformative pairs.
#' \item \code{"count.favorable"}: returns the number of pairs in favor of the treatment.
#' \item \code{"count.unfavorable"}: returns the number of pairs in favor of the control.
#' \item \code{"count.neutral"}: returns the number of neutral pairs.
#' \item \code{"count.uninf"}: returns the number of uninformative pairs.
#' }
#' 
#' @return When \code{resampling=FALSE} and \code{strata=FALSE}, return a numeric vector.
#' When \code{resampling=FALSE} and \code{strata=TRUE}, return a matrix.
#' When \code{resampling=TRUE} and \code{strata=FALSE}, return a matrix.
#' When \code{resampling=TRUE} and \code{strata=TRUE}, return an 3-dimensional array.
#' 
#' @keywords method
#' @author Brice Ozenne

## * method - coef
#' @rdname S4BuyseTest-coef
#' @exportMethod coef
setMethod(f = "coef",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                strata = NULL,
                                cumulative = NULL,
                                resampling = FALSE,
                                simplify = TRUE,
                                ...){
              
              ## ** normalize arguments
              option <- BuyseTest.options()
              mycall <- match.call()
              dots <- list(...)
              if(length(dots)>0){
                  stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
              }
    
              ## statistic
              if(is.null(statistic)){
                  statistic <- option$statistic
              }

              statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                  "netbenefit" = "netBenefit",
                                  "winratio" = "winRatio",
                                  "favorable" = "favorable",
                                  "unfavorable" = "unfavorable",
                                  "uninformative" = "uninf",
                                  statistic)
              
              type.count <- c("count.favorable","count.unfavorable","count.neutral","count.uninf")

              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio","favorable","unfavorable","neutral","uninf",type.count),
                             valid.length = 1,
                             method = "coef[S4BuyseTest]")

              if(is.null(cumulative)){
                  if(statistic %in% c("count.neutral","neutral","count.uninf","uninf")){
                      cumulative <- FALSE
                  }else{
                      cumulative <- TRUE
                  }
              }else if(statistic %in% c("count.neutral","neutral","count.uninf","uninf") && cumulative){
                  stop("Argument \'cumulative\' must be FALSE when argument \'statistic\' is set to \"neutral\" or \"uninf\". \n")
              }

              ## endpoint
              valid.endpoint <- names(object@endpoint)
              n.endpoint <- length(valid.endpoint)
              if(is.null(endpoint)){
                  endpoint <- valid.endpoint
              }else if(!is.null(endpoint)){
                  if(is.numeric(endpoint)){
                      validInteger(endpoint,
                                   name1 = "endpoint",
                                   min = 1, max = length(valid.endpoint),
                                   valid.length = NULL,
                                   method = "iid[BuyseTest]")
                      endpoint <- valid.endpoint[endpoint]
                  }else{
                      validCharacter(endpoint,
                                     valid.length = 1:length(valid.endpoint),
                                     valid.values = valid.endpoint,
                                     refuse.NULL = FALSE)
                  }
              }
              weightEndpoint <- slot(object, "weightEndpoint")
                  
              ## strata
              level.strata <- object@level.strata
              n.strata <- length(level.strata)
              weightStrata <- object@weightStrata
              if(is.null(strata)){
                  if(length(object@level.strata)==1){
                      strata <- "global"                      
                  }else{
                      strata <- c("global", object@level.strata)
                  }
              }else if(identical(strata,FALSE)){
                  strata <- "global"
              }else if(identical(strata,TRUE)){
                  strata <- object@level.strata
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",object@level.strata),
                                 refuse.NULL = FALSE,
                                 method = "coef[S4BuyseTest]")
              }
              
              ## resampling
              if(resampling){

                  if(!attr(slot(object, "method.inference"),"permutation") && !attr(slot(object, "method.inference"),"bootstrap")){
                      stop("No resampling procedure was performed so cannot output the corresponding coefficients. \n")
                  }

                  if(statistic %in% type.count){
                      stop("The number of ",gsub("count.","",statistic)," pairs when performing resampling is not saved. \n")
                  }

                  n.resampling <- slot(object, "n.resampling")
                  weightStrataResampling <- slot(object, "weightStrataResampling")

              }
              
                            
              ## ** win ratio as a special case
              if(statistic == "winRatio"){
                  type.weightStrata <- attr(weightStrata,"type")

                  if(identical(strata,"global") && cumulative){
                      ## avoid any operation, rely on C++ code
                      if(resampling){
                          if(simplify){
                              out <- slot(object, "DeltaResampling")[,endpoint,statistic]
                          }else{
                              out <- matrix(slot(object, "DeltaResampling")[,endpoint,statistic],
                                            ncol = length(endpoint), dimnames = list(NULL, endpoint))
                          }                      
                      }else{
                          if(simplify){
                              out <- slot(object, "Delta")[endpoint,statistic]
                          }else{
                              out <- matrix(slot(object, "Delta")[endpoint,statistic],
                                            nrow = 1, ncol = length(endpoint), dimnames = list("global", endpoint))
                          }
                      }
                      
                  }else if(type.weightStrata != "var-winratio" || all(strata %in% level.strata)){
                      ## from C++ code
                      ## Delta.col(5) = arma::trans(arma::sum(wcumdelta_favorable,0)/arma::sum(wcumdelta_unfavorable,0));

                      out.fav <- coef(object, statistic = "favorable",
                                      endpoint = endpoint, strata = strata, cumulative = cumulative,
                                      resampling = resampling, simplify = simplify)
                      out.unfav <- coef(object, statistic = "unfavorable",
                                        endpoint = endpoint, strata = strata, cumulative = cumulative,
                                        resampling = resampling, simplify = simplify)

                      out <- out.fav/out.unfav

                  }else if(type.weightStrata == "var-winratio" && "global" %in% strata){
                      ## from C++ code
                      ## arma::mat winRatioStrata = cumdelta_favorable/cumdelta_unfavorable;
                      ## Delta.col(5) = arma::trans(arma::sum(winRatioStrata.each_col() % weightPool,0));
                      
                      out.fav <- coef(object, statistic = "favorable",
                                      endpoint = valid.endpoint, strata = level.strata, cumulative = cumulative,
                                      resampling = resampling, simplify = FALSE)
                      out.unfav <- coef(object, statistic = "unfavorable",
                                        endpoint = valid.endpoint, strata = level.strata, cumulative = cumulative,
                                        resampling = resampling, simplify = FALSE)

                      out.ratio <- out.fav/out.unfav

                      if(resampling){
                          out <- array(NA, dim = c(n.resampling, n.strata+1, n.endpoint),
                                       dimnames = list(NULL, c("global",level.strata), valid.endpoint))
                          out[,level.strata,valid.endpoint] <- out.ratio
                              
                          for(iE in 1:n.endpoint){ ## iE <- 1
                              if(n.strata == 1){
                                  deltaResampling[,"global",iE] <- out.ratio[,1,iE]
                              }else{
                                  deltaResampling[,"global",iE] <- rowSums(out.ratio[,,iE]*weightStrataResampling)
                              }
                          }
                          out <- out[,strata,endpoint,drop=simplify]
                          
                      }else{
                          out <- rbind(global = colSums(.colMultiply_cpp(out.ratio, weightStrata)),
                                       object.delta)[strata,endpoint,drop=simplify]
                      }
                  }

                  return(out)
                  
              }

              ## ** normalize element in object (add global or stratified result)
              if(statistic %in% type.count){

                  object.statistic <- slot(object, statistic)                  

                  delta <- rbind(global = colSums(object.statistic),
                                 object.statistic)
                  Delta <- matrix(.rowCumSum_cpp(delta),
                                  nrow = n.strata+1, ncol = n.endpoint,
                                  dimnames = list(c("global",level.strata), valid.endpoint))

              }else if(if(statistic != "winRatio" && resampling){
                  
                  object.deltaResampling <- array(slot(object, "deltaResampling")[,,,statistic],
                                                  dim = c(n.resampling, n.strata, n.endpoint),
                                                  dimnames = list(NULL, level.strata, valid.endpoint))
                  object.DeltaResampling <- matrix(slot(object, "DeltaResampling")[,,statistic],
                                                   ncol = n.endpoint, dimnames = list(NULL, valid.endpoint))
                  
                  deltaResampling <- array(NA, dim = c(n.resampling, n.strata+1, n.endpoint),
                                           dimnames = list(NULL, c("global",level.strata), valid.endpoint))
                  deltaResampling[,level.strata,valid.endpoint] <- object.deltaResampling
                  for(iE in 1:n.endpoint){ ## iE <- 1
                      if(n.strata == 1){
                          deltaResampling[,"global",iE] <- object.deltaResampling[,1,iE]
                      }else{
                          deltaResampling[,"global",iE] <- rowSums(object.deltaResampling[,,iE]*weightStrataResampling)
                      }
                  }

                  DeltaResampling <- array(NA, dim = c(n.resampling, n.strata+1, n.endpoint),
                                           dimnames = list(NULL, c("global",level.strata), valid.endpoint))
                  DeltaResampling[,"global",valid.endpoint] <- object.DeltaResampling 
                  for(iS in 1:n.strata){ ## iS <- 1
                      if(n.endpoint == 1){
                          DeltaResampling[,iS+1,1] <- object.deltaResampling[,iS,1]*weightEndpoint
                      }else{
                          DeltaResampling[,iS+1,] <- .rowCumSum_cpp(.rowMultiply_cpp(object.deltaResampling[,iS,], weightEndpoint))
                      }
                  }

              }else if(statistic != "winRatio"){

                  object.delta <- matrix(slot(object, "delta")[,,statistic], 
                                         nrow = n.strata, ncol = n.endpoint, dimnames = list(level.strata, valid.endpoint))
                  object.Delta <- matrix(slot(object, "Delta")[,statistic],
                                         nrow = 1, ncol = n.endpoint, dimnames = list("global", valid.endpoint))                  

                  delta <- rbind(global = colSums(.colMultiply_cpp(object.delta, weightStrata)),
                                 object.delta)
                  Delta <- rbind(object.Delta,
                                 .rowCumSum_cpp(.rowMultiply_cpp(object.delta, weightEndpoint)))

              }
              
              ## ** extract information
              if(resampling == FALSE){

                  if(cumulative==TRUE){
                      out <- Delta[strata,endpoint,drop=simplify]
                  }else if(cumulative == FALSE){
                      out <- delta[strata,endpoint,drop=simplify]
                  }
                  
              }else if(resampling){

                  if(cumulative==TRUE){
                      out <- DeltaResampling[,strata,endpoint,drop=simplify]
                  }else if(cumulative == FALSE){
                      out <- deltaResampling[,strata,endpoint,drop=simplify]
                  }
              }              

              ## ** export
              return(out)

          })

######################################################################
### S4BuyseTest-coef.R ends here
