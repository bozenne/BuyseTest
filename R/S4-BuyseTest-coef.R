### S4BuyseTest-coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2019 (10:45) 
## Version: 
## Last-Updated: jul 23 2025 (17:16) 
##           By: Brice Ozenne
##     Update #: 426
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
#' @param statistic [character] the statistic summarizing the pairwise comparison: \itemize{
#' \item \code{"netBenefit"}: displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
#' \item \code{"winRatio"}: win ratio or win odds.
#' \item \code{"favorable"}: proportion strictly in favor of the treatment or Mann-Whitney parameter.
#' \item \code{"unfavorable"}: proportion in favor of the control.
#' \item \code{"neutral"}: proportion of neutral pairs.
#' \item \code{"uninf"}: proportion of uninformative pairs.
#' \item \code{"count.favorable"}: number of pairs in favor of the treatment.
#' \item \code{"count.unfavorable"}: number of pairs in favor of the control.
#' \item \code{"count.neutral"}: number of neutral pairs.
#' \item \code{"count.uninf"}: number of uninformative pairs.
#' }
#' Default value read from \code{BuyseTest.options()}.
#' @param endpoint [character] for which endpoint(s) the summary statistic should be output?
#' If \code{NULL} returns the summary statistic for all endpoints.
#' @param strata [character vector] the strata relative to which the statistic should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the statistic pooled over all strata,
#' or \code{TRUE} to output each strata-specific statistic.
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param resampling [logical] should the summary statistic obtained by resampling be output?
#' @param simplify [logical] should the result be coerced to the lowest possible dimension?
#' @param ... ignored.
#'
#' @details
#' \bold{statistic}: with a single endpoint denoted \eqn{Y} and \eqn{X} in the treatment and control group and a threshold of clinical relevance \eqn{\tau}: \itemize{
#' \item \code{"netBenefit"}: \eqn{P[Y \ge X + \tau] - P[X \ge Y + \tau]}. See Buyse (2010).
#' \item \code{"winRatio"}: the win ratio \eqn{\frac{P[Y \ge X + \tau]}{P[X \ge Y + \tau]}} or the win odds \eqn{\frac{P[Y \ge X + \tau]+0.5P[|Y - X|<\tau]}{P[X \ge Y + \tau]+0.5P[|Y - X|<\tau]}}. see Wang (2016) and Dong (2019).
#' \item \code{"favorable"}: \eqn{P[Y \ge X + \tau]} or the Mann-Whitney parameter \eqn{P[Y \ge X + \tau]+0.5P[|Y - X|<\tau]}. See Fay (2018).
#' \item \code{"unfavorable"}: \eqn{P[Y \le X + \tau]} or \eqn{P[Y \le X + \tau]+0.5P[|Y - X|<\tau]}.
#' }
#' The value of the argument \code{add.halfNeutral} used when running \code{\link{BuyseTest}} decides whether \eqn{0.5P[|Y - X|<\tau]} is considered, e.g. whether the win ratio or win odds is output.
#' 
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the Mann-Whitney parameter: Fay, Michael P. et al (2018). \bold{Causal estimands and confidence intervals asscoaited with Wilcoxon-Mann-Whitney tests in randomized experiments}. \emph{Statistics in Medicine} 37:2923-2937 \cr
#' On the win odds: Dong, G., Hoaglin, D. C., Qiu, J., Matsouaka, R. A., Chang, Y. W., Wang, J., & Vandemeulebroecke, M. (2019). \bold{The Win Ratio: On Interpretation and Handling of Ties}. \emph{Statistics in Biopharmaceutical Research}, 12(1), 99–106. https://doi.org/10.1080/19466315.2019.1575279 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245
#' 
#' @return When \code{resampling=FALSE} and \code{simplify=FALSE}, a matrix (strata, endpoint).
#' When \code{resampling=FALSE} and \code{simplify=FALSE}, an array (sample, strata, endpoint).
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
                                strata = FALSE,
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
    
              ## *** statistic
              add.halfNeutral <- slot(object,"add.halfNeutral")
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

              ## *** endpoint
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
                  
              ## *** strata
                  level.strata <- object@level.strata
              n.strata <- length(level.strata)
              weightStrata <- object@weightStrata
              type.weightStrata <- attr(weightStrata,"type")

              if(is.null(strata)){
                  if(length(level.strata)==1){
                      strata <- "global"                      
                  }else{
                      strata <- c("global", level.strata)
                  }
              }else if(identical(strata,FALSE)){
                  strata <- "global"
              }else if(identical(strata,TRUE)){
                  strata <- level.strata
              }else if(is.numeric(strata)){
                  validInteger(strata,
                               name1 = "strata",
                               valid.length = NULL,
                               min = 1,
                               max = length(level.strata),
                               refuse.NULL = TRUE,
                               refuse.duplicates = TRUE,
                               method = "autoplot[S4BuyseTest]")
                  strata <- level.strata[strata]
              }else{
                  validCharacter(strata,
                                 name1 = "strata",
                                 valid.length = NULL,
                                 valid.values = c("global",level.strata),
                                 refuse.NULL = FALSE,
                                 method = "coef[S4BuyseTest]")
              }
              

              ## *** resampling
              if(resampling){

                  if(!attr(slot(object, "method.inference"),"permutation") && !attr(slot(object, "method.inference"),"bootstrap")){
                      stop("No resampling procedure was performed so cannot output the corresponding coefficients. \n")
                  }

                  if(statistic %in% type.count){
                      stop("The number of ",gsub("count.","",statistic)," pairs when performing resampling is not saved. \n")
                  }
                  
                  n.resampling <- dim(slot(object, "nResampling"))[1]
                  weightStrataResampling <- slot(object, "weightStrataResampling")

              }
        
              ## ** normalize element in object (add global or stratified result)
              if(statistic %in% type.count){

                  object.statistic <- slot(object, statistic)                  

                  delta <- rbind(global = colSums(object.statistic),
                                 object.statistic)
                  Delta <- matrix(.rowCumSum_cpp(delta),
                                  nrow = n.strata+1, ncol = n.endpoint,
                                  dimnames = list(c("global",level.strata), valid.endpoint))

              }else if(resampling == FALSE){

                  ## *** outcome specific
                  delta.strata <- matrix(slot(object, "delta")[,,statistic], 
                                         nrow = n.strata, ncol = n.endpoint, dimnames = list(level.strata, valid.endpoint))
                  
                  if(statistic == "winRatio" && type.weightStrata != "var-winratio"){
                      out.fav <- coef(object, statistic = "favorable",
                                      endpoint = valid.endpoint, strata = "global", cumulative = FALSE,
                                      resampling = FALSE, simplify = FALSE)
                      out.unfav <- coef(object, statistic = "unfavorable",
                                        endpoint = valid.endpoint, strata = "global", cumulative = FALSE,
                                        resampling = FALSE, simplify = FALSE)
                      delta.global <- out.fav/out.unfav
                  }else{
                      delta.global <- colSums(.colMultiply_cpp(delta.strata, weightStrata))
                  }
                  delta <- rbind(global = delta.global, delta.strata)

                  ## *** cumulated over outcomes
                  Delta.global <- matrix(slot(object, "Delta")[,statistic],
                                         nrow = 1, ncol = n.endpoint, dimnames = list("global", valid.endpoint))                  

                  if(statistic == "winRatio"){
                      out.cumFav <- coef(object, statistic = "favorable",
                                         endpoint = valid.endpoint, strata = level.strata, cumulative = TRUE,
                                         resampling = FALSE, simplify = FALSE)
                      out.cumUnfav <- coef(object, statistic = "unfavorable",
                                           endpoint = valid.endpoint, strata = level.strata, cumulative = TRUE,
                                           resampling = FALSE, simplify = FALSE)
                      Delta.strata <- out.cumFav/out.cumUnfav
                  }else{
                      Delta.strata <- .rowCumSum_cpp(.rowMultiply_cpp(delta.strata, weightEndpoint))
                      if(add.halfNeutral && n.endpoint>1 && statistic %in% c("favorable","unfavorable")){
                          ## do not cumulate neutral over endpoints, only keep the proportion w.r.t. the last endpoint
                          neutral.strata <- matrix(slot(object, "delta")[,,"neutral"], 
                                                   nrow = n.strata, ncol = n.endpoint, dimnames = list(level.strata, valid.endpoint))
                          Delta.strata <- Delta.strata - 0.5 * .rowCumSum_cpp(cbind(0,.rowMultiply_cpp(neutral.strata, weightEndpoint)[,1:(n.endpoint-1),drop=FALSE]))
                      }
                  }
                  rownames(Delta.strata) <- rownames(delta.strata)
                  Delta <- rbind(global = Delta.global, Delta.strata)

              }else if(resampling){
                  
                  object.deltaResampling <- array(slot(object, "deltaResampling")[,,,statistic],
                                                  dim = c(n.resampling, n.strata, n.endpoint),
                                                  dimnames = list(NULL, level.strata, valid.endpoint))
                  object.DeltaResampling <- matrix(slot(object, "DeltaResampling")[,,statistic],
                                                   ncol = n.endpoint, dimnames = list(NULL, valid.endpoint))

                  deltaResampling <- array(NA, dim = c(n.resampling, n.strata+1, n.endpoint),
                                           dimnames = list(NULL, c("global",level.strata), valid.endpoint))
                  deltaResampling[,level.strata,valid.endpoint] <- object.deltaResampling

                  DeltaResampling <- array(NA, dim = c(n.resampling, n.strata+1, n.endpoint),
                                           dimnames = list(NULL, c("global",level.strata), valid.endpoint))
                  DeltaResampling[,"global",valid.endpoint] <- object.DeltaResampling 

                  if(statistic == "winRatio"){
                      if(statistic == "winRatio" && type.weightStrata != "var-winratio"){
                          favorableResampling <- coef(object, statistic = "favorable",
                                                      endpoint = valid.endpoint, strata = "global", cumulative = FALSE,
                                                      resampling = TRUE, simplify = FALSE)
                          unfavorableResampling <- coef(object, statistic = "unfavorable",
                                                        endpoint = valid.endpoint, strata = "global", cumulative = FALSE,
                                                        resampling = TRUE, simplify = FALSE)
                      }

                      cumFavorableResampling <- coef(object, statistic = "favorable",
                                                     endpoint = valid.endpoint, strata = level.strata, cumulative = TRUE,
                                                     resampling = TRUE, simplify = FALSE)
                      cumUnfavorableResampling <- coef(object, statistic = "unfavorable",
                                                       endpoint = valid.endpoint, strata = level.strata, cumulative = TRUE,
                                                       resampling = TRUE, simplify = FALSE)
                  }else if(add.halfNeutral && n.endpoint>1 && statistic %in% c("favorable","unfavorable")){
                      object.neutralResampling <- array(slot(object, "deltaResampling")[,,,"neutral"],
                                                        dim = c(n.resampling, n.strata, n.endpoint),
                                                        dimnames = list(NULL, level.strata, valid.endpoint))                              
                                  
                  }

                  for(iE in 1:n.endpoint){ ## iE <- 1
                      if(n.strata == 1){
                          deltaResampling[,"global",iE] <- object.deltaResampling[,1,iE]
                      }else{
                          if(statistic == "winRatio" && type.weightStrata != "var-winratio"){
                              deltaResampling[,"global",iE] <- favorableResampling[,"global",iE]/unfavorableResampling[,"global",iE]
                          }else{
                              deltaResampling[,"global",iE] <- rowSums(object.deltaResampling[,,iE]*weightStrataResampling)
                          }                  
                      }
                  }

                  for(iS in 1:n.strata){ ## iS <- 1
                      if(n.endpoint == 1){
                          DeltaResampling[,iS+1,1] <- object.deltaResampling[,iS,]
                      }else{
                          if(statistic == "winRatio"){
                              DeltaResampling[,iS+1,] <- cumFavorableResampling[,iS,]/cumUnfavorableResampling[,iS,]
                          }else{
                              DeltaResampling[,iS+1,] <- .rowCumSum_cpp(.rowMultiply_cpp(object.deltaResampling[,iS,], weightEndpoint))

                              if(add.halfNeutral && n.endpoint>1 && statistic %in% c("favorable","unfavorable")){
                                  ## do not cumulate neutral over endpoints, only keep the proportion w.r.t. the last endpoint
                                  DeltaResampling[,iS+1,] <- DeltaResampling[,iS+1,] - 0.5 * .rowCumSum_cpp(cbind(0,.rowMultiply_cpp(object.neutralResampling[,iS,], weightEndpoint)[,1:(n.endpoint-1),drop=FALSE]))
                                  
                              }
                              
                          }
                      }
                  }

              }

              ## ** extract information
              if(resampling == FALSE){

                  if(cumulative==TRUE){
                      out <- Delta[strata,endpoint,drop=simplify]
                      if(length(strata)==1 && is.null(names(out))){
                          names(out) <- endpoint
                      }
                  }else if(cumulative == FALSE){
                      out <- delta[strata,endpoint,drop=simplify]
                      if(length(strata)==1 && is.null(names(out))){
                          names(out) <- endpoint
                      }
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
