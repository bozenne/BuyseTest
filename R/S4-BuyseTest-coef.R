### S4BuyseTest-coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2019 (10:45) 
## Version: 
## Last-Updated: Mar 13 2023 (09:32) 
##           By: Brice Ozenne
##     Update #: 245
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
#' @title Coef Method for Class "S4BuyseTest"
#' @aliases coef,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'  
#' @description Extract summary statistics from the result of a \code{\link{BuyseTest}} function.
#' 
#' @param object output of \code{\link{BuyseTest}}
#' @param statistic [character] the type of summary statistic. See the detail section.
#' @param endpoint [character] for which endpoint(s) the summary statistic should be output?
#' If \code{NULL} returns the summary statistic for all endpoints.
#' @param stratified [logical] should the summary statistic be strata-specific?
#' Otherwise a summary statistic over all strata is returned.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param resampling [logical] should the summary statistic obtained by resampling be output?
#' 
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
#' @return When \code{resampling=FALSE} and \code{stratified=FALSE}, return a numeric vector.
#' When \code{resampling=FALSE} and \code{stratified=TRUE}, return a matrix.
#' When \code{resampling=TRUE} and \code{stratified=FALSE}, return a matrix.
#' When \code{resampling=TRUE} and \code{stratified=TRUE}, return an 3-dimensional array.
#' 
#' @keywords coef S4BuyseTest-method
#' @author Brice Ozenne

## * method - coef
#' @rdname S4BuyseTest-coef
#' @exportMethod coef
setMethod(f = "coef",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                stratified = FALSE,
                                cumulative = NULL,
                                resampling = FALSE,
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
              if(!is.null(endpoint)){
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

              ## resampling
              if(resampling){
                  if(!attr(slot(object, "method.inference"),"permutation") && !attr(slot(object, "method.inference"),"bootstrap")){
                      stop("Not resampling procedure was performed so cannot output the corresponding coefficients. \n")
                  }
                  n.resampling <- slot(object, "n.resampling")
                  DeltaResampling <- slot(object, "DeltaResampling")
                  deltaResampling <- slot(object, "deltaResampling")
                  weightStrataResampling <- slot(object, "weightStrataResampling")
              }else if(statistic %in% type.count){
                  Acount <- array(c(slot(object, "count.favorable"),
                                    slot(object, "count.unfavorable"),
                                    slot(object, "count.neutral"),
                                    slot(object, "count.uninf")),
                                  dim = c(n.strata,n.endpoint,4),
                                  dimnames = list(level.strata, valid.endpoint, c("count.favorable","count.unfavorable","count.neutral","count.uninf"))
                                  )
                  Mcount <- matrix(Acount[,,statistic],
                                   nrow = n.strata, ncol = n.endpoint,
                                   dimnames = list(level.strata, valid.endpoint))
              }else{
                  Delta <- slot(object, "Delta")
                  delta <- slot(object, "delta")
              }
              if(statistic %in% type.count && resampling){
                  stop("The number of ",statistic," pairs when performing resampling is not saved. \n")
              }
              
              

              ## ** specific cases: redirect to the least calculation when having to choose between two equivalent possibilities
              if(n.strata == 1 && (statistic %in% type.count == FALSE) && ("stratified" %in% names(mycall) == FALSE)){
                  rm.strata <- TRUE
                  if(cumulative){stratified <- FALSE}else{stratified <- TRUE}
              }else if(n.endpoint == 1 && ("cumulative" %in% names(mycall) == FALSE)){
                  if(stratified){cumulative <- FALSE}else{cumulative <- TRUE}
                  rm.strata <- FALSE
              }else{
                  rm.strata <- FALSE
              }
                              
              ## ** extract information
              if(cumulative==TRUE && stratified==FALSE){

                  if(statistic %in% type.count){
                      out <- cumsum(colSums(Mcount))
                  }else if(resampling==FALSE){
                      out <- stats::setNames(Delta[,statistic],
                                             rownames(Delta))
                  }else{
                      out <- matrix(DeltaResampling[,,statistic],
                                    nrow = dim(DeltaResampling)[1], ncol = dim(DeltaResampling)[2],
                                    dimnames = dimnames(DeltaResampling)[1:2])
                  }

              }else if(cumulative==FALSE && stratified==TRUE){

                  if(statistic %in% type.count){
                      out <- Mcount
                  }else if(resampling==FALSE){
                      out <- matrix(delta[,,statistic],
                                    nrow = dim(delta)[1], ncol = dim(delta)[2],
                                    dimnames = dimnames(delta)[1:2])
                      if(rm.strata){
                          out <- stats::setNames(out[1,], dimnames(delta)[[2]])
                      }
                  }else{
                      out <- array(deltaResampling[,,,statistic],
                                   dim = dim(deltaResampling)[1:3],
                                   dimnames = dimnames(deltaResampling)[1:3])
                      if(rm.strata){
                          out <- matrix(out[,1,],
                                        nrow = dim(deltaResampling)[1], ncol = dim(deltaResampling)[3],
                                        dimnames = dimnames(deltaResampling)[c(1,3)])
                      }
                  }

              }else if(statistic %in% type.count){

                  if(cumulative==FALSE && stratified==FALSE){
                      out <- colSums(Mcount)
                  }else if(cumulative==TRUE && stratified==TRUE){
                      out <- matrix(.rowCumSum_cpp(Mcount),
                                    nrow = n.strata, ncol = n.endpoint,
                                    dimnames = list(level.strata, valid.endpoint))
                  }

              }else if(statistic != "winRatio"){
                  if(resampling==FALSE){
                      M.stat <- matrix(delta[,,statistic],
                                       nrow = dim(delta)[1], ncol = dim(delta)[2],
                                       dimnames = dimnames(delta)[1:2])

                      if(cumulative==FALSE && stratified==FALSE){

                          out <- stats::setNames(colSums(.colMultiply_cpp(M.stat, object@weightStrata)),
                                                 colnames(delta))
                   
                      }else if(cumulative==TRUE && stratified==TRUE){
                          
                          out <- matrix(.rowCumSum_cpp(.rowMultiply_cpp(M.stat, weightEndpoint)),
                                        nrow = dim(delta)[1], ncol = dim(delta)[2],
                                        dimnames = dimnames(delta)[1:2])
                          
                      }

                  }else{
                      M.stat <- array(deltaResampling[,,,statistic],
                                      dim = dim(deltaResampling)[1:3],
                                      dimnames = dimnames(deltaResampling)[1:3])

                      if(cumulative==FALSE && stratified==FALSE){

                          out <- matrix(apply(M.stat, MARGIN = 3, FUN = function(iM){rowSums(iM*weightStrataResampling)}),
                                        nrow = dim(deltaResampling)[1], ncol = dim(deltaResampling)[3],
                                        dimnames = dimnames(deltaResampling)[c(1,3)])                          

                      }else if(cumulative==TRUE && stratified==TRUE){

                          ls.out <- apply(M.stat, MARGIN = 2, FUN = function(iM){.rowCumSum_cpp(.rowMultiply_cpp(iM, weightEndpoint))}, simplify = FALSE)
                          out <- aperm(array(unlist(ls.out),
                                             dim = dim(deltaResampling)[c(1,3,2)],
                                             dimnames = dimnames(deltaResampling)[c(1,3,2)]),
                                       perm = c(1,3,2))                         

                      }
                  }                      
              }else{

                  if(attr(object@weightStrata,"type")=="var-winratio" && cumulative == FALSE && stratified == FALSE){
                      A.fav <- coef(object, cumulative = cumulative, stratified = TRUE, resampling = resampling, statistic = "favorable")
                      A.unfav <- coef(object, cumulative = cumulative, stratified = TRUE, resampling = resampling, statistic = "unfavorable")
                      A.win <- A.fav/A.unfav
                      if(resampling){
                          out <- matrix(apply(A.win, MARGIN = 3, FUN = function(iM){rowSums(iM*weightStrataResampling)}),
                                        nrow = dim(deltaResampling)[1], ncol = dim(deltaResampling)[3],
                                        dimnames = dimnames(deltaResampling)[c(1,3)])                          
                      }else{
                          out <- stats::setNames(colSums(.colMultiply_cpp(A.win, object@weightStrata)),
                                                 colnames(A.win))
                      }
                  }else{
                      M.fav <- coef(object, cumulative = cumulative, stratified = stratified, resampling = resampling, statistic = "favorable")
                      M.unfav <- coef(object, cumulative = cumulative, stratified = stratified, resampling = resampling, statistic = "unfavorable")
                      out <- M.fav/M.unfav
                  }
              }

              ## ** export
              if(!is.null(endpoint)){
                  if((!stratified || rm.strata) && !resampling){
                      out <- out[endpoint]
                  }else if((stratified && !resampling) || (!stratified && resampling)){
                      out <- out[,endpoint,drop=FALSE]                      
                  }else{
                      out <- out[,,endpoint,drop=FALSE]
                  }
              }
              return(out)

          })

######################################################################
### S4BuyseTest-coef.R ends here
