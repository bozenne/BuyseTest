## * getCount (documentation)
#' @docType methods
#' @name getCount
#' @title Extract the Number of Favorable, Unfavorable, Neutral, Uninformative pairs
#' @aliases getCount,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'
#' @description Extract the number of favorable, unfavorable, neutral, uninformative pairs.
#'
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param type the type of pairs to be counted. Can be \code{"favorable"}, \code{"unfavorable"}, \code{neutral}, or \code{uninf}. Can also be \code{"all"} to select all of them.
#'
#' @return
#'   A \code{"vector"} containing the number of pairs
#'
#' @keywords get S4BuyseTest-method
#' @author Brice Ozenne

## * getCount (code)
#' @rdname getCount
#' @exportMethod getCount
setMethod(f = "getCount",
          signature = "S4BuyseTest",
          definition = function(object, type){

            if (missing(type)) {
              type <- c("favorable","unfavorable","neutral","uninf")
            }

            validCharacter(type,
                           valid.length = NULL,
                           valid.values = c("favorable","unfavorable","neutral","uninf"),
                           method = "getCount")

            out <- NULL
            if ("favorable" %in% type) {out <- c(out, favorable = object@count.favorable)}
            if ("unfavorable" %in% type) {out <- c(out, unfavorable = object@count.unfavorable)}
            if ("neutral" %in% type) {out <- c(out, neutral = object@count.neutral)}
            if ("uninf" %in% type) {out <- c(out, uninf = object@count.uninf)}

            return(out)
          }
)


## * getIid (documentation)
#' @docType methods
#' @name getIid
#' @title Extract the H-decomposition of the Estimator
#' @aliases getIid,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Extract the H-decomposition of the GPC estimator.
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [character] for which endpoint(s) the H-decomposition should be output?
#' If \code{NULL} returns the sum of the H-decomposition over all endpoints.
#' @param type [character] type of H-decomposition to be output.
#' Can be only for the nuisance parameters (\code{"nuisance"}),
#' or for the u-statistic given the nuisance parameters (\code{"u-statistic"}),
#' or both.
#' @param center [logical] if \code{TRUE} the H-decomposition is centered around 0 (estimated statistic is substracted).
#' @param scale [logical] if \code{TRUE} the H-decomposition is rescaled (by the sample size in the corresponding arm) such that its sums of squares approximate the variance of the estimator.
#' @param cluster [numeric vector] return the H-decomposition aggregated by cluster.
#' @param statistic [character] the statistic relative to which the H-decomposition should be output: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @param cumulative [logical] should the H-decomposition be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param strata [character vector] the strata relative to which the H-decomposition of the statistic should be output.
#' Can also be \code{"global"} or \code{FALSE} to output the H-decompositon of the pooled statistic.
#' or \code{TRUE} to output the H-decompositon of each strata-specific statistic.
#' @param simplify [logical] should the result be coerced to the lowest possible dimension?
#' 
#' @details WARNING: argument \code{scale} and \code{center} should be used with care as when set to \code{FALSE} they may not lead to a meaningful decomposition.
#'  
#' @return A list of matrices, each element of the list correspond to a statistic (global or strata-specific) and each matrix has as many columns as endpoints and rows as observations.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#' 
#' @keywords S4BuyseTest-method
#' @author Brice Ozenne

## * getIid (code)
#' @rdname getIid
#' @exportMethod getIid
setMethod(f = "getIid",
          signature = "S4BuyseTest",
          definition = function(object,
                                endpoint = NULL,
                                statistic = NULL,
                                strata = FALSE,
                                cumulative = TRUE,
                                center = TRUE,
                                scale = TRUE,
                                type = "all",
                                cluster = NULL,
                                simplify = TRUE){

              ## ** normalize user input
              option <- BuyseTest.options()
              n.obs <- NROW(object@iidAverage$favorable)

              ## iid has been stored in object
              if(object@method.inference != "u statistic"){
                  stop("No H-decomposition in the object \n",
                       "Set the argument \'method.inference\' to \"u-statistic\" when calling BuyseTest \n")
              }

              ## index group
              indexC <- attr(object@level.treatment,"indexC")
              indexT <- attr(object@level.treatment,"indexT")

              ## strata
              level.strata <- object@level.strata
              weightStrata <- object@weightStrata
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
              if(attr(weightStrata,"type")=="standardization"){
                  if(any(strata %in% level.strata)){
                      stop("Cannot output the strata-specific H-decomposition when using standardization. \n")
                  }
                  n.strata <- 1
                  indexStrata <- list(union(indexC,indexT))
              }else{
                  n.strata <- length(level.strata)
                  indexStrata <- attr(level.strata,"index")                  
              }
              attr(level.strata,"index") <- NULL

              ## cluster
              if(!is.null(cluster)){
                  if(length(cluster) != n.obs){
                      stop("Incorrect length for argument \'cluster\'. Should have length ",n.obs ,".\n")
                  }
                  if(!is.factor(cluster)){
                      cluster <- as.factor(cluster)
                  }
                  Ucluster <- levels(cluster)
                  cluster <- as.numeric(cluster)                  
              }
              
              ## endpoint
              valid.endpoint <- names(object@endpoint)
              n.endpoint <- length(valid.endpoint)
              weightEndpoint <- object@weightEndpoint
              if(is.numeric(endpoint)){
                  validInteger(endpoint,
                               name1 = "endpoint",
                               min = 1, max = length(valid.endpoint),
                               valid.length = NULL,
                               method = "iid[BuyseTest]")
                  endpoint <- valid.endpoint[endpoint]
              }else if(is.null(endpoint)){
                  endpoint <- valid.endpoint
              }else{
              }
              validCharacter(endpoint,
                             valid.length = 1:length(valid.endpoint),
                             valid.values = valid.endpoint)
              
              ## pairs
              n.pairs <- object@n.pairs
              ntot.pair <- sum(n.pairs)
              
              ## type
              validCharacter(gsub("-"," ",tolower(type), fixed = TRUE),
                             valid.length = 1,
                             valid.values = c("all","nuisance","u statistic"),
                             refuse.NULL = FALSE)


              ##  statistic
              add.halfNeutral <- slot(object,"add.halfNeutral")
              if(is.null(statistic)){
                  statistic <- option$statistic
              }else{
                  statistic <- sapply(tolower(statistic), function(iChar){
                      switch(gsub("[[:blank:]]", "", iChar),
                             "netbenefit" = "netBenefit",
                             "winratio" = "winRatio",
                             "favorable" = "favorable",
                             "unfavorable" = "unfavorable",
                             statistic)})
                  
                  validCharacter(statistic,
                                 name1 = "statistic",
                                 valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                                 valid.length = 1,
                                 refuse.duplicates = TRUE,
                                 method = "getIid[S4BuyseTest]")
              }
              ## if(length(statistic)>1 && any(strata %in% level.strata)){
              ##     stop("Argument \'statistic\' must be of length one when asking for strata-specific H-decomposition. \n")
              ## }

              ## ** extract H-decomposition
              if(type %in% c("all","u statistic")){
                  object.iid <- object@iidAverage[c("favorable","unfavorable","neutral")]
              }else{
                  object.iid <- list(favorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                        dimnames = list(NULL, valid.endpoint)),
                                     unfavorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                          dimnames = list(NULL, valid.endpoint)),
                                     neutral = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                      dimnames = list(NULL, valid.endpoint))
                                     )
              }

              if(type %in% c("all","nuisance") && (object@scoring.rule=="Peron")){
                  if(length(object@iidNuisance$favorable)>0){
                      object.iid$favorable <- object.iid$favorable + object@iidNuisance$favorable
                  } ## otherwise model.tte has been passed as argument and there is no uncertainty regarding nuisance
                  if(length(object@iidNuisance$unfavorable)>0){
                      object.iid$unfavorable <- object.iid$unfavorable + object@iidNuisance$unfavorable
                  } ## otherwise model.tte has been passed as argument and there is no uncertainty regarding nuisance
                  if(length(object@iidNuisance$neutral)>0){
                      object.iid$neutral <- object.iid$neutral + object@iidNuisance$neutral
                  } ## otherwise model.tte has been passed as argument and there is no uncertainty regarding nuisance
              }

              ## ** remove normalization
              if((center==FALSE) || (scale==FALSE)){

                  if(center==FALSE){
                      delta.favorable <- coef(object, endpoint = valid.endpoint, cumulative = FALSE, statistic = "favorable", strata = level.strata, resampling = FALSE, simplify = FALSE)
                      delta.unfavorable <- coef(object, endpoint = valid.endpoint, cumulative = FALSE, statistic = "unfavorable", strata = level.strata, resampling = FALSE, simplify = FALSE)
                      delta.neutral <- coef(object, endpoint = valid.endpoint, cumulative = FALSE, statistic = "neutral", strata = level.strata, resampling = FALSE, simplify = FALSE)
                  }

                  for(iter_strata in 1:n.strata){  ## iter_strata <- 1
                      iStrataC <- intersect(indexC,indexStrata[[iter_strata]])
                      iStrataT <- intersect(indexT,indexStrata[[iter_strata]])
                      ## remove scaling (by the sample size: \sum_i IF_i^2 -> 1/n^2 \sum_i IF_i^2)
                      object.iid$favorable[iStrataC,] <- length(iStrataC) * object.iid$favorable[iStrataC,,drop=FALSE]
                      object.iid$favorable[iStrataT,] <- length(iStrataT) * object.iid$favorable[iStrataT,,drop=FALSE]

                      object.iid$unfavorable[iStrataC,] <- length(iStrataC) * object.iid$unfavorable[iStrataC,,drop=FALSE]
                      object.iid$unfavorable[iStrataT,] <- length(iStrataT) * object.iid$unfavorable[iStrataT,,drop=FALSE]

                      object.iid$neutral[iStrataC,] <- length(iStrataC) * object.iid$neutral[iStrataC,,drop=FALSE]
                      object.iid$neutral[iStrataT,] <- length(iStrataT) * object.iid$neutral[iStrataT,,drop=FALSE]
                  
                      ## remove centeringn
                      if(center==FALSE){
                          object.iid$favorable[c(iStrataC,iStrataT),] <- .rowCenter_cpp(object.iid$favorable[c(iStrataC,iStrataT),,drop=FALSE], - delta.favorable[iter_strata,,drop=FALSE])
                          object.iid$unfavorable[c(iStrataC,iStrataT),] <- .rowCenter_cpp(object.iid$unfavorable[c(iStrataC,iStrataT),,drop=FALSE], -delta.unfavorable[iter_strata,,drop=FALSE])
                          object.iid$neutral[c(iStrataC,iStrataT),] <- .rowCenter_cpp(object.iid$neutral[c(iStrataC,iStrataT),,drop=FALSE], -delta.neutral[iter_strata,,drop=FALSE])
                          if(scale){ ## restaure scaling
                              object.iid$favorable[iStrataC,] <- object.iid$favorable[iStrataC,,drop=FALSE]/length(iStrataC)
                              object.iid$favorable[iStrataT,] <- object.iid$favorable[iStrataT,,drop=FALSE]/length(iStrataT)

                              object.iid$unfavorable[iStrataC,] <- object.iid$unfavorable[iStrataC,,drop=FALSE]/length(iStrataC)
                              object.iid$unfavorable[iStrataT,] <- object.iid$unfavorable[iStrataT,,drop=FALSE]/length(iStrataT)

                              object.iid$neutral[iStrataC,] <- object.iid$neutral[iStrataC,,drop=FALSE]/length(iStrataC)
                              object.iid$neutral[iStrataT,] <- object.iid$neutral[iStrataT,,drop=FALSE]/length(iStrataT)
                          }
                      }
                  }

              }

              ## ** cumulate over endpoints
              if(cumulative && n.endpoint > 1){
                  keep.names <- list(favorable = colnames(object.iid$favorable),
                                     unfavorable = colnames(object.iid$unfavorable))

                  object.iid$favorable <- .rowCumSum_cpp(.rowMultiply_cpp(object.iid$favorable, weightEndpoint))
                  colnames(object.iid$favorable) <- keep.names$favorable
                  object.iid$unfavorable <- .rowCumSum_cpp(.rowMultiply_cpp(object.iid$unfavorable, weightEndpoint))
                  colnames(object.iid$unfavorable) <- keep.names$unfavorable
                  if(add.halfNeutral && statistic %in% c("favorable","unfavorable","winRatio")){
                      ## do not cumulate neutral over endpoints, only keep the proportion w.r.t. the last endpoint
                      object.iid$favorable <- object.iid$favorable - 0.5 * .rowCumSum_cpp(cbind(0,.rowMultiply_cpp(object.iid$neutral, weightEndpoint)[,1:(n.endpoint-1),drop=FALSE]))
                      object.iid$unfavorable <- object.iid$unfavorable - 0.5 * .rowCumSum_cpp(cbind(0,.rowMultiply_cpp(object.iid$neutral, weightEndpoint)[,1:(n.endpoint-1),drop=FALSE]))
                  }

              }
              
              ## ** add strata weights
              if(attr(weightStrata,"type")=="standardization"){
                  ls.iid.favorable <- list(global = object.iid$favorable)
                  ls.iid.unfavorable <- list(global = object.iid$unfavorable)
              }else{
                  MweightStrata <- matrix(NA, nrow = n.obs, ncol = n.endpoint, dimnames = list(NULL,valid.endpoint))
                  MweightStrata[unlist(indexStrata),] <- do.call(rbind,lapply(1:n.strata, function(iS){
                      matrix(weightStrata[iS], nrow = length(indexStrata[[iS]]), ncol = n.endpoint)
                  }))

                  ls.iid.favorable <- c(list(global = object.iid$favorable*MweightStrata),
                                        stats::setNames(lapply(1:n.strata, function(iStrata){ ## iStrata <- 1
                                            iM <- object.iid$favorable
                                            iM[-indexStrata[[iStrata]],] <- 0
                                            return(iM)
                                        }), level.strata))

                  ls.iid.unfavorable <- c(list(global = object.iid$unfavorable*MweightStrata),
                                          stats::setNames(lapply(1:n.strata, function(iStrata){ ## iStrata <- 1
                                              iM <- object.iid$unfavorable
                                              iM[-indexStrata[[iStrata]],] <- 0
                                              return(iM)
                                          }), level.strata))
              }

              ## ** aggregate at a cluster level
              if(!is.null(cluster)){
                  ls.iid.favorable <- lapply(ls.iid.favorable, function(iIID){ ## iIID <- ls.iid.favorable[[1]]
                      do.call(rbind,by(iIID,cluster,colSums, simplify = FALSE))
                  })
                  ls.iid.unfavorable <- lapply(ls.iid.unfavorable, function(iIID){
                      do.call(rbind,by(iIID,cluster,colSums, simplify = FALSE))
                  })
              }

              ## ** iid decomposition of the chosen statistic for each endpoint
              if("winRatio" %in% statistic){
                  Delta.favorable <- coef(object, endpoint = endpoint, cumulative = cumulative, statistic = "favorable", strata = strata, resampling = FALSE, simplify = FALSE)
                  Delta.unfavorable <- coef(object, endpoint = endpoint, cumulative = cumulative, statistic = "unfavorable", strata = strata, resampling = FALSE, simplify = FALSE)
              }
              
              out <- stats::setNames(lapply(strata, function(iS){ ## iS <- strata[1]

                  iIID.fav <- ls.iid.favorable[[iS]][,endpoint,drop=FALSE]
                  iIID.unfav <- ls.iid.unfavorable[[iS]][,endpoint,drop=FALSE]

                  if(statistic == "favorable"){
                      iOut <- iIID.fav
                  }else if(statistic == "unfavorable"){
                      iOut <- iIID.unfav
                  }else if(statistic == "netBenefit"){
                      iOut <- iIID.fav - iIID.unfav
                  }else if(statistic == "winRatio"){
                      iOut <- .rowScale_cpp(iIID.fav,Delta.unfavorable[iS,]) - .rowMultiply_cpp(iIID.unfav, Delta.favorable[iS,]/Delta.unfavorable[iS,]^2)
                      colnames(iOut) <- colnames(iIID.fav)
                  }                                        
                  return(iOut)
              }), strata)

              ## ** output H-decomposition
              if(simplify && length(strata)==1){
                  out <- out[[1]]
              }
              if(simplify && length(endpoint)==1){
                  if(length(strata)==1){
                      out <- out[,1]
                  }else{
                      out <- do.call(cbind,out)
                      colnames(out) <- strata
                  }
              }
              return(out)
              
          })

## * getPairScore (documentation)
#' @docType methods
#' @name getPairScore
#' @title Extract the Score of Each Pair
#' @aliases getPairScore,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'
#' @description Extract the score of each pair.
#'
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [integer/character vector] the endpoint for which the scores should be output.
#' @param strata [character vector] the strata relative to which the score should be output.
#' @param rm.withinStrata [logical] should the columns indicating the position of each member of the pair
#' within each treatment group be removed?
#' @param rm.strata [logical] should the column containing the level of the strata variable be removed from the output?
#' @param rm.indexPair [logical] should the column containing the number associated to each pair be removed from the output?
#' @param rm.weight [logical] should the column weight be removed from the output?
#' @param rm.corrected [logical] should the columns corresponding to the scores after weighting be removed from the output?
#' @param cumulative [logical] should the scores be cumulated over endpoints?
#' @param unlist [logical] should the structure of the output be simplified when possible?
#' @param trace [logical] should a message be printed to explain what happened
#' when the function returned \code{NULL}?
#'
#' @details The maximal output (i.e. with all columns) contains for each endpoint, a data.table with:
#' \itemize{
#' \item \code{"strata"}: the name of the strata to which the pair belongs.
#' \item \code{"index.T"}: the index of the treatment observation in the pair relative to the original dataset.
#' \item \code{"index.C"}: the index of the control observation in the pair relative to the original dataset.
#' \item \code{"indexWithinStrata.T"}: the index of the treatment observation in the pair relative to the treatment group and the strata.
#' \item \code{"indexWithinStrata.C"}: the index of the control observation in the pair relative to the control group and the strata.
#' \item \code{"favorable"}: the probability that the endpoint is better in the treatment arm vs. in the control arm.
#' \item \code{"unfavorable"}: the probability that the endpoint is worse in the treatment arm vs. in the control arm.
#' \item \code{"neutral"}: the probability that the endpoint is no different in the treatment arm vs. in the control arm.
#' \item \code{"uninformative"}: the weight of the pair that cannot be attributed to favorable/unfavorable/neutral.
#' \item \code{"weight"}: the residual weight of the pair to be analyzed at the current outcome. Each pair starts with a weight of 1.
#' \item \code{"favorable.corrected"}: same as \code{"favorable"}  after weighting.
#' \item \code{"unfavorable.corrected"}: same as \code{"favorable"} after weighting.
#' \item \code{"neutral.corrected"}: same as \code{"favorable"} after weighting.
#' \item \code{"uninformative.corrected"}: same as \code{"favorable"} after weighting.
#' }
#' Note that the \code{.T} and \code{.C} may change since they correspond of the label of the treatment and control arms.
#' The first weighting consists in multiplying the probability by the residual weight of the pair
#' (i.e. the weight of the pair that was not informative at the previous endpoint). This is always performed.
#' For time to event endpoint an additional weighting may be performed to avoid a possible bias in presence of censoring.
#' @keywords get S4BuyseTest-method
#' @author Brice Ozenne
#' @examples
#' library(data.table)
#' library(prodlim)
#' 
#' ## run BuyseTest
#' library(survival) ## import veteran
#'
#' BT.keep <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status") + cont(karno),
#'                      data = veteran, keep.pairScore = TRUE, 
#'                      trace = 0, method.inference = "none")
#'
#' ## Extract scores
#' pScore <- getPairScore(BT.keep, endpoint = 1)
#'
#' ## look at one pair
#' indexPair <- intersect(which(pScore$index.1 == 22),
#'                        which(pScore$index.2 == 71))
#' pScore[indexPair]
#'
#' ## retrive pair in the original dataset
#' pVeteran <- veteran[pScore[indexPair,c(index.1,index.2)],]
#' pVeteran
#' 
#' ## the observation from the control group is censored at 97
#' ## the observation from the treatment group has an event at 112
#' ## since the threshold is 20, and (112-20)<97
#' ## we know that the pair is not in favor of the treatment
#'
#' ## the formula for probability in favor of the control is
#' ## Sc(97)/Sc(112+20)
#' ## where Sc(t) is the survival at time t in the control arm.
#' 
#' ## we first estimate the survival in each arm
#' e.KM <- prodlim(Hist(time,status)~trt, data = veteran)
#'
#' ## and compute the survival
#' iSurv <- predict(e.KM, times =  c(97,112+20),
#'                  newdata = data.frame(trt = 1, stringsAsFactors = FALSE))[[1]]
#'
#' ## the probability in favor of the control is then
#' pUF <- iSurv[2]/iSurv[1]
#' pUF
#' ## and the complement to one of that is the probability of being neutral
#' pN <- 1 - pUF
#' pN
#' 
#' if(require(testthat)){
#'    testthat::expect_equal(pUF, pScore[indexPair, unfavorable])
#'    testthat::expect_equal(pN, pScore[indexPair, neutral])
#' }

## * getPairScore (code)
#' @rdname getPairScore
#' @exportMethod getPairScore
setMethod(f = "getPairScore",
          signature = "S4BuyseTest",
          definition = function(object, endpoint, strata, cumulative,
                                rm.withinStrata, rm.strata, rm.indexPair, rm.weight, rm.corrected,
                                unlist, trace){

              ## ** extract table
              if(length(object@tablePairScore)==0){
                  if(trace){
                      cat("pairScore was not exported from the object \n",
                          "Consider setting the argument \'keep.pairScore\' to \"TRUE\" when calling the \"BuyseTest\" function \n", sep = "")
                  }
                  return(invisible(NULL))
              }else{
                  out <- data.table::copy(object@tablePairScore)

                  endpoint.names <- object@endpoint
                  strata.names <- object@level.strata
                  
                  if(!is.null(endpoint)){
                      if(is.numeric(endpoint)){
                          validInteger(endpoint, min = 1, max = length(endpoint.names), valid.length = NULL,
                                       refuse.duplicates = TRUE)
                      }else if(is.character(endpoint)){
                          validCharacter(endpoint, valid.length = NULL, valid.values = endpoint.names,
                                         refuse.duplicates = TRUE)
                          endpoint <- match(endpoint, endpoint.names)
                      }else{
                          stop("Argument \'endpoint\' must be a numeric of character vector \n")
                      }

                      out <- out[endpoint] 
                  }

                  ## ** restrict to strata strata
                  if(!is.null(strata)){
                      if(is.numeric(strata)){
                          validInteger(strata, min = 1, max = length(strata.names), valid.length = NULL,
                                       refuse.duplicates = TRUE)
                          strata <- strata.names[strata]
                      }else if(is.character(strata)){
                          validCharacter(strata, valid.length = NULL, valid.values = strata.names,
                                         refuse.duplicates = TRUE)
                      }else{
                          stop("Argument \'endpoint\' must be a numeric of character vector \n")
                      }

                      if(length(strata.names)>1){
                          for(iEndpoint in 1:length(out)){ ## iEndpoint <- 1
                              index.strata <- which(out[[iEndpoint]]$strata %in% strata)
                              out[[iEndpoint]] <- out[[iEndpoint]][index.strata]
                          }
                      }

                  }

                  old.names <- c("index.C", "index.T", "indexWithinStrata.C", "indexWithinStrata.T")
                  new.names <- c(paste0("index.",object@level.treatment), paste0("indexWithinStrata.",object@level.treatment))

                  ## ** cumulate endpoints
                  if(cumulative){

                      ## *** add weights
                      weightEndpoint <- object@weightEndpoint[endpoint]
                      if(any(weightEndpoint!=1)){
                          test <- lapply(1:length(out), function(iE){
                              out[[iE]][, c("favorable","unfavorable","neutral","uninf") := list(.SD$favorable*weightEndpoint[iE],
                                                                                                 .SD$unfavorable*weightEndpoint[iE],
                                                                                                 .SD$neutral*weightEndpoint[iE],
                                                                                                 .SD$uninf*weightEndpoint[iE])]
                              out[[iE]][, c("favorableC","unfavorableC","neutralC","uninfC") := list(.SD$favorableC*weightEndpoint[iE],
                                                                                                     .SD$unfavorableC*weightEndpoint[iE],
                                                                                                     .SD$neutralC*weightEndpoint[iE],
                                                                                                     .SD$uninfC*weightEndpoint[iE])]
                              return(NULL)
                          })
                      }

                      ## *** cumulate
                      if(length(out)>1){
                          out <- lapply(out, function(iOut){data.table::setkeyv(iOut,"index.pair")})

                          for(iEndpoint in 2:length(out)){ ## iEndpoint <- 2
                              iOut.save <- out[[iEndpoint]]
                              out[[iEndpoint]] <- data.table::copy(out[[iEndpoint-1]])
                              out[[iEndpoint]][iOut.save$index.pair, c("favorable","unfavorable","neutral","uninf") := list(.SD$favorable + iOut.save$favorable,
                                                                                                                            .SD$unfavorable + iOut.save$unfavorable,
                                                                                                                            iOut.save$neutral,
                                                                                                                            iOut.save$uninf)]

                              out[[iEndpoint]][iOut.save$index.pair, c("favorableC","unfavorableC","neutralC","uninfC") := list(.SD$favorableC + iOut.save$favorableC,
                                                                                                                                .SD$unfavorableC + iOut.save$unfavorableC,
                                                                                                                                iOut.save$neutralC,
                                                                                                                                iOut.save$uninfC)]
                          }
                          out <- lapply(out, function(iOut){data.table::setkeyv(iOut,c("index.T","index.C"))})

                      }
                  }
 
                  for(iEndpoint in 1:length(out)){ ## iEndpoint <- 2
                      if(rm.withinStrata){
                          out[[iEndpoint]][,c("indexWithinStrata.T","indexWithinStrata.C") := NULL]
                          data.table::setnames(out[[iEndpoint]], old = old.names[1:2], new = new.names[1:2])
                      }else{
                          data.table::setnames(out[[iEndpoint]], old = old.names, new = new.names)
                      }
                      if(rm.indexPair){
                          out[[iEndpoint]][,c("index.pair") := NULL]
                      }
                      if(rm.strata){
                          out[[iEndpoint]][,c("strata") := NULL]
                      }
                      if(rm.weight){
                          out[[iEndpoint]][,c("weight") := NULL]
                      }
                      if(rm.corrected){
                          out[[iEndpoint]][,c("favorableC","unfavorableC","neutralC","uninfC") := NULL]
                      }
                  }
                  
                  if(length(out) == 1 && unlist == TRUE){
                      out <- out[[1]] 
                  }

                  return(out[])
              }
          })

## * getPseudovalue (documentation)
#' @docType methods
#' @name getPseudovalue
#' @title Extract the pseudovalues of the Estimator
#' @aliases getPseudovalue,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 

#' @description Extract the pseudovalues of the estimator.
#' The average of the pseudovalues is the estimate and their standard deviation the standard error of the estimate times a factor n
#' (i.e. a t-test on their mean will give asymptotically valid confidence intervals and p-values).
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [character] for which endpoint(s) the pseudovalues should be output?
#' If \code{NULL} returns the sum of the H-decomposition over all endpoints.
#' @param statistic [character] the statistic relative to which the pseudovalues should be computed: \code{"netBenefit"}, \code{"winRatio"}, \code{"favorable"}, \code{"unfavorable"}.
#' See the documentation of the \code{coef} method for further details.
#' Default value read from \code{BuyseTest.options()}. 
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#' 
#' @keywords method
#' @author Brice Ozenne
#' @examples
#' set.seed(10)
#' n <- 250
#' d <- simBuyseTest(n)
#'
#' e.BT <- BuyseTest(treatment ~ tte(eventtime,status,2) + bin(toxicity),
#'                  data = d, trace = 0)
#'
#' #### net Benefit
#' pseudo <- getPseudovalue(e.BT)
#' summary(lm(pseudo~1))$coef
#' ## asymptotically equivalent to
#' confint(e.BT, transformation = TRUE)
#' ## (small differences: small sample corrections)
#' 
#' summary(lm(getPseudovalue(e.BT, endpoint = 1)~1))$coef
#' 
#' #### win Ratio
#' pseudo <- getPseudovalue(e.BT, statistic = "winRatio")
#' summary(lm(pseudo~1))$coef ## wrong p-value (should compare to 1 instead of 0)
#' ## asymptotically equivalent to
#' confint(e.BT, statistic = "winRatio", transformation = TRUE)
#' 
#' #### favorable
#' pseudo <- getPseudovalue(e.BT, statistic = "favorable")
#' summary(lm(pseudo~1))$coef ## wrong p-value (should compare to 1/2 instead of 0)
#' ## asymptotically equivalent to
#' confint(e.BT, statistic = "favorable", transformation = TRUE)
#' 
#' #### unfavorable
#' pseudo <- getPseudovalue(e.BT, statistic = "unfavorable")
#' summary(lm(pseudo~1))$coef ## wrong p-value (should compare to 1/2 instead of 0)
#' ## asymptotically equivalent to
#' confint(e.BT, statistic = "unfavorable", transformation = TRUE)

## * getPseudovalue (code)
#' @rdname getPseudovalue
#' @exportMethod getPseudovalue
setMethod(f = "getPseudovalue",
          signature = "S4BuyseTest",
          definition = function(object, statistic = NULL, endpoint = NULL){

              option <- BuyseTest.options()

              ## ** normalize arguments

              ## endpoint
              valid.endpoint <- names(object@endpoint)

              if(is.null(endpoint)){
                  endpoint <- utils::tail(valid.endpoint,1)
              }else if(is.numeric(endpoint)){
                  validInteger(endpoint,
                               name1 = "endpoint",
                               min = 1, max = length(valid.endpoint),
                               valid.length = 1,
                               method = "iid[BuyseTest]")
                  endpoint <- valid.endpoint[endpoint]
              }else{
                  validCharacter(endpoint,
                                 valid.length = 1,
                                 valid.values = valid.endpoint,
                                 refuse.NULL = FALSE)
              }

              ##  statistics
              if(is.null(statistic)){
                  statistic <- option$statistic
              }else{
                  statistic <- switch(gsub("[[:blank:]]", "", tolower(statistic)),
                                      "netbenefit" = "netBenefit",
                                      "winratio" = "winRatio",
                                      "favorable" = "favorable",
                                      "unfavorable" = "unfavorable",
                                      statistic)
              }
              
              validCharacter(statistic,
                             name1 = "statistic",
                             valid.values = c("netBenefit","winRatio","favorable","unfavorable"),
                             valid.length = 1,
                             method = "getPseudovalue[S4BuyseTest]")
              
              ## ** compute pseudovalue
              object.delta <- coef(object, statistic = statistic)[endpoint]
              count.favorable <- coef(object, statistic = "favorable")[endpoint]
              count.unfavorable <- coef(object, statistic = "unfavorable")[endpoint]
              object.iid <- data.frame(favorable = unname(getIid(object, endpoint = endpoint, statistic = "favorable", strata = "global", simplify = FALSE)[["global"]]),
                                       unfavorable = unname(getIid(object, endpoint = endpoint, statistic = "unfavorable", strata = "global", simplify = FALSE)[["global"]]))
              n.obs <- NROW(object.iid)

              out <- switch(statistic,
                            "favorable" = n.obs * object.iid$favorable + object.delta,
                            "unfavorable" = n.obs * object.iid$unfavorable + object.delta,
                            "netBenefit" = n.obs * (object.iid$favorable - object.iid$unfavorable) + object.delta,
                            "winRatio" = n.obs * (object.iid$favorable / count.unfavorable - object.iid$unfavorable * (count.favorable/count.unfavorable^2)) + object.delta,
                            )

              ## ** export
              return(out)
          })

## * getSurvival (documentation)
#' @docType methods
#' @name getSurvival
#' @title Extract the Survival and Survival Jumps
#' @aliases getSurvival,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'
#' @description Extract the survival and survival jumps.
#'
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param type [character vector] the type of survival to be output. See details.
#' @param endpoint [integer/character vector] the endpoint for which the survival should be output.
#' @param strata [integer/character vector] the strata relative to which the survival should be output.
#' @param unlist [logical] should the structure of the output be simplified when possible.
#' @param trace [logical] should a message be printed to explain what happened
#' when the function returned \code{NULL}.
#' 
#' @details The argument \code{type} can take any of the following values:
#' \itemize{
#' \item \code{"survTimeC"}: survival at the event times for the observations of the control arm.
#' \item \code{"survTimeT"}: survival at the event times for the observations of the treatment arm.
#' \item \code{"survJumpC"}: survival at the jump times for the survival model in the control arm.
#' \item \code{"survJumpT"}: survival at the time times for the survival model in the treatment arm.
#' \item \code{"lastSurv"}: survival at the last event time.
#' }
#'
#' @keywords get S4BuyseTest-method
#' @author Brice Ozenne

## * getSurvival (code)
#' @rdname getSurvival
#' @exportMethod getSurvival
setMethod(f = "getSurvival",
          signature = "S4BuyseTest",
          definition = function(object, type, endpoint, strata, unlist, trace){

              if(length(object@tableSurvival)==0){
                  
                  if(trace>0){
                      if(all(tolower(object@type)!="timetoevent")){
                          add.txt <- "No endpoint of type time to event \n"
                      }else if(tolower(object@scoring.rule)!="peron"){
                          add.txt <- "Consider setting the argument \'scoring.rule\' to \"Peron\" when calling BuyseTest \n"
                      }else{
                          add.txt <- "Consider setting the argument \'keep.survival\' to TRUE in BuyseTest.options \n"
                      }
                      cat("Survival was not exported from the object \n",
                          add.txt, sep = "")    
                  }
                  return(invisible(NULL))
              }else{

                  if(is.null(type)){
                      type <- c("survTimeC","survTimeT","survJumpC","survJumpT","lastSurv")
                  }else{
                      validCharacter(type, valid.length = NULL, refuse.duplicates = TRUE,
                                     valid.values = c("survTimeC","survTimeT","survJumpC","survJumpT","lastSurv"))
                  }
                  if(!is.null(type)){
                      out <- object@tableSurvival[type]
                  }else{
                      out <- data.table::copy(object@tableSurvival)
                  }
                  
                  endpoint.names <- object@endpoint
                  strata.names <- object@level.strata
                  
                  if(!is.null(endpoint)){
                      if(is.numeric(endpoint)){
                          validInteger(endpoint, min = 1, max = length(endpoint.names), valid.length = NULL,
                                       refuse.duplicates = TRUE)
                      }else if(is.character(endpoint)){
                          validCharacter(endpoint, valid.length = NULL, valid.values = endpoint.names,
                                         refuse.duplicates = TRUE)
                          endpoint <- match(endpoint, endpoint.names)
                      }else{
                          stop("Argument \'endpoint\' must be a numeric of character vector \n")
                      }

                      if("survTimeC" %in% type){ out$survTimeC <- out$survTimeC[endpoint] }
                      if("survTimeT" %in% type){ out$survTimeT <- out$survTimeT[endpoint] } 
                      if("survJumpC" %in% type){ out$survJumpC <- out$survJumpC[endpoint] }
                      if("survJumpT" %in% type){ out$survJumpT <- out$survJumpT[endpoint] }
                      if("lastSurv" %in% type){ out$lastSurv <- out$lastSurv[endpoint] }
                  }
                  
                  if(!is.null(strata)){
                      if(is.numeric(strata)){
                          validInteger(strata, min = 1, max = length(strata.names), valid.length = NULL,
                                       refuse.duplicates = TRUE)
                      }else if(is.character(strata)){
                          validCharacter(strata, valid.length = NULL, valid.values = strata.names,
                                         refuse.duplicates = TRUE)
                          endpoint <- match(strata, strata.names)
                      }else{
                          stop("Argument \'endpoint\' must be a numeric of character vector \n")
                      }

                      for(iEndpoint in 1:length(out[[1]])){
                          if(length(strata)==1 && unlist == TRUE){
                              if("survTimeC" %in% type){ out$survTimeC[[iEndpoint]] <- out$survTimeC[[iEndpoint]][[1]] }
                              if("survTimeT" %in% type){ out$survTimeT[[iEndpoint]] <- out$survTimeT[[iEndpoint]][[1]] }
                              if("survJumpC" %in% type){ out$survJumpC[[iEndpoint]] <- out$survJumpC[[iEndpoint]][[1]] }
                              if("survJumpT" %in% type){ out$survJumpT[[iEndpoint]] <- out$survJumpT[[iEndpoint]][[1]] }
                              if("lastSurv" %in% type){ out$lastSurv[[iEndpoint]] <- out$lastSurv[[iEndpoint]][1,] }
                          }else{
                              if("survTimeC" %in% type){ out$survTimeC[[iEndpoint]] <- out$survTimeC[[iEndpoint]][strata] }
                              if("survTimeT" %in% type){ out$survTimeT[[iEndpoint]] <- out$survTimeT[[iEndpoint]][strata] }
                              if("survJumpC" %in% type){ out$survJumpC[[iEndpoint]] <- out$survJumpC[[iEndpoint]][strata] }
                              if("survJumpT" %in% type){ out$survJumpT[[iEndpoint]] <- out$survJumpT[[iEndpoint]][strata] }
                              if("lastSurv" %in% type){ out$lastSurv[[iEndpoint]] <- out$lastSurv[[iEndpoint]][strata,,drop=FALSE] }
                          }
                      }

                  }

                  if(length(endpoint) == 1 && unlist == TRUE){
                      if("survTimeC" %in% type){ out$survTimeC <- out$survTimeC[[1]] }
                      if("survTimeT" %in% type){ out$survTimeT <- out$survTimeT[[1]] }
                      if("survJumpC" %in% type){ out$survJumpC <- out$survJumpC[[1]] }
                      if("survJumpT" %in% type){ out$survJumpT <- out$survJumpT[[1]] }
                      if("lastSurv" %in% type){ out$lastSurv <- out$lastSurv[[1]] }
                  }

                  if(length(type) == 1 && unlist == TRUE){
                      out <- out[[type]]
                  }
                  return(out)
              }

              
          })
