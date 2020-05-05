## * getCount (documentation)
#' @name S4BuyseTest-getCount
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
#' @rdname S4BuyseTest-getCount
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
#' @name S4BuyseTest-getIid
#' @title Extract the H-decomposition of the Estimator
#' @aliases getIid,S4BuyseTest-method
#' @include S4-BuyseTest.R
#' 
#' @description Extract the H-decomposition of the GPC estimator.
#' 
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [character] for which endpoint(s) the H-decomposition should be output?
#' If \code{NULL} returns the sum of the H-decomposition over all endpoints.
#' @param type [character] type of iid to be output.
#' Can be only for the nuisance parameters (\code{"nuisance"}),
#' or for the u-statistic given the nuisance parameters (\code{"u-statistic"}),
#' or both.
#' @param normalize [logical] if \code{TRUE} the iid is centered and multiplied by the sample size.
#' Otherwise not.
#' @param cluster [numeric vector] return the H-decomposition aggregated by cluster.
#' 
#' @seealso 
#' \code{\link{BuyseTest}} for performing a generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-summary}} for a more detailed presentation of the \code{S4BuyseTest} object.
#' 
#' @keywords S4BuyseTest-method
#' @author Brice Ozenne

## * getIid (code)
#' @rdname S4BuyseTest-getIid
setMethod(f = "getIid",
          signature = "S4BuyseTest",
          definition = function(object, endpoint = NULL, normalize = TRUE, type = "all", cluster = NULL){

              n.obs <- NROW(object@iidAverage$favorable)
              valid.endpoint <- paste0(object@endpoint,"_",object@threshold)
              n.endpoint <- length(valid.endpoint)
              if(!is.null(cluster) && !is.numeric(cluster)){
                  cluster <- as.numeric(as.factor(cluster))
              }
              
              ## ** check arguments              
              if(is.numeric(endpoint)){
                  validInteger(endpoint,
                               name1 = "endpoint",
                               min = 1, max = length(valid.endpoint),
                               valid.length = NULL,
                               method = "iid[BuyseTest]")
                  endpoint <- valid.endpoint[endpoint]
              }   
              validCharacter(endpoint,
                             valid.length = 1:length(valid.endpoint),
                             valid.values = valid.endpoint,
                             refuse.NULL = FALSE)
              validCharacter(type,
                             valid.length = 1,
                             valid.values = c("all","nuisance","u-statistic"),
                             refuse.NULL = FALSE)
              if(object@method.inference != "u-statistic"){
                  stop("No H-decomposition in the object \n",
                       "Set the argument \'method.inference\' to \"u-statistic\" when calling BuyseTest \n")
              }
              validInteger(cluster, valid.length = n.obs, min = 1, max = n.obs, refuse.NA = TRUE, refuse.NULL = FALSE, refuse.duplicates = FALSE)

              ## ** extract H-decomposition
              if(type %in% c("all","u-statistic")){
                  object.iid <- object@iidAverage
              }else{
                  object.iid <- list(favorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                        dimnames = list(NULL, valid.endpoint)),
                                     unfavorable = matrix(0, nrow = n.obs, ncol = n.endpoint,
                                                          dimnames = list(NULL, valid.endpoint))
                                     )
              }
              if(type %in% c("all","nuisance") && (object@scoring.rule=="Peron")){
                  object.iid$favorable <- object.iid$favorable + object@iidNuisance$favorable
                  object.iid$unfavorable <- object.iid$unfavorable + object@iidNuisance$unfavorable
              }

              if(normalize==FALSE){
                  delta.favorable <- colSums(object@count.favorable)/sum(object@n.pairs)
                  delta.unfavorable <- colSums(object@count.unfavorable)/sum(object@n.pairs)
                  indexC <- attr(object@level.treatment,"indexC")
                  indexT <- attr(object@level.treatment,"indexT")

                  ## remove scaling 
                  object.iid$favorable[indexC,] <- length(indexC) * object.iid$favorable[indexC,]
                  object.iid$favorable[indexT,] <- length(indexT) * object.iid$favorable[indexT,]

                  object.iid$unfavorable[indexC,] <- length(indexC) * object.iid$unfavorable[indexC,]
                  object.iid$unfavorable[indexT,] <- length(indexT) * object.iid$unfavorable[indexT,]

                  ## remove centering
                  object.iid$unfavorable <- sweep(object.iid$unfavorable, MARGIN = 2, FUN = "+", STATS = cumsum(delta.unfavorable*object@weight))
                  object.iid$favorable <- sweep(object.iid$favorable, MARGIN = 2, FUN = "+", STATS = cumsum(delta.favorable*object@weight))
              }
              ## ** accumulate H-decomposition
              if(is.null(endpoint)){                  
                  ## iid decomposition over all endpoints
                  object.iid <- do.call(cbind,lapply(object.iid, function(iI){
                      iIID <- iI[, NCOL(iI)]
                      if(!is.null(cluster)){
                          iIID <- tapply(iIID,cluster,sum)
                      }
                      return(iIID)
                  }))
              }else{
                  ## iid decomposition for each endpoint
                  object.iid <- lapply(endpoint, function(iE){
                      iIID <- cbind(favorable = object.iid$favorable[,iE],
                                    unfavorable = object.iid$unfavorable[,iE])
                      if(!is.null(cluster)){
                          iIID <- cbind(favorable = tapply(iIID[,"favorable"],cluster,sum),
                                       unfavorable = tapply(iIID[,"unfavorable"],cluster,sum))
                      }
                      return(iIID)
                  })
                  names(object.iid) <- endpoint
              }

              ## ** output H-decomposition
              return(object.iid)
    
})

## * getPairScore (documentation)
#' @name S4BuyseTest-getPairScore
#' @title Extract the Score of Each Pair
#' @aliases getPairScore,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'
#' @description Extract the score of each pair.
#'
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [integer/character vector] the endpoint for which the scores should be output.
#' @param strata [integer/character vector] the strata for which the scores should be output.
#' @param rm.withinStrata [logical] should the columns indicating the position of each member of the pair
#' within each treatment group be removed?
#' @param rm.strata [logical] should the column containing the level of the strata variable be removed from the output?
#' @param rm.indexPair [logical] should the column containing the number associated to each pair be removed from the output?
#' @param rm.weight [logical] should the column weight be removed from the output?
#' @param rm.corrected [logical] should the columns corresponding to the scores after weighting be removed from the output?
#' @param sum [logical] should the scores be cumulated over endpoints?
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

## * getPairScore (examples)
#' @rdname S4BuyseTest-getPairScore
#' @examples
#' library(data.table)
#' library(prodlim)
#' 
#' ## run BuyseTest
#' data(veteran,package="survival")
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
setMethod(f = "getPairScore",
          signature = "S4BuyseTest",
          definition = function(object, endpoint, strata, sum,
                                rm.withinStrata, rm.strata, rm.indexPair, rm.weight, rm.corrected,
                                unlist, trace){

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
                  
                  if(!is.null(strata)){
                      if(is.numeric(strata)){
                          validInteger(strata, min = 1, max = length(strata.names), valid.length = NULL,
                                       refuse.duplicates = TRUE)
                      }else if(is.character(strata)){
                          validCharacter(strata, valid.length = NULL, valid.values = strata.names,
                                         refuse.duplicates = TRUE)
                          strata <- match(strata, strata.names)
                      }else{
                          stop("Argument \'endpoint\' must be a numeric of character vector \n")
                      }
                      
                      for(iEndpoint in 1:length(out)){ ## iEndpoint <- 1
                          index.strata <- which(out[[iEndpoint]]$strata %in% strata)
                          out[[iEndpoint]][, c("strata") := factor(.SD$strata, levels = 1:length(strata.names), labels = strata.names)]
                          out[[iEndpoint]] <- out[[iEndpoint]][index.strata]
                      }

                  }

                  old.names <- c("index.C", "index.T", "indexWithinStrata.C", "indexWithinStrata.T")
                  new.names <- c(paste0("index.",object@level.treatment), paste0("indexWithinStrata.",object@level.treatment))

                  if(sum && length(out)>1){
                      out.save <- out
                      out <- out.save[1]                     
                      for(iEndpoint in 2:length(out.save)){
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("favorable") := .SD$favorable + out.save[[iEndpoint]]$favorable]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("unfavorable") := .SD$unfavorable + out.save[[iEndpoint]]$unfavorable]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("neutral") := .SD$neutral + out.save[[iEndpoint]]$neutral]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("uninf") := .SD$uninf + out.save[[iEndpoint]]$uninf]

                          out[[1]][out.save[[iEndpoint]]$index.pair, c("weight") := .SD$weight + out.save[[iEndpoint]]$weight]

                          out[[1]][out.save[[iEndpoint]]$index.pair, c("favorableC") := .SD$favorableC + out.save[[iEndpoint]]$favorableC]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("unfavorableC") := .SD$unfavorableC + out.save[[iEndpoint]]$unfavorableC]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("neutralC") := .SD$neutralC + out.save[[iEndpoint]]$neutralC]
                          out[[1]][out.save[[iEndpoint]]$index.pair, c("uninfC") := .SD$uninfC + out.save[[iEndpoint]]$uninfC]
                      }
                  }
                  
                  for(iEndpoint in 1:length(out)){ ## iEndpoint <- 1
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

## * getSurvival (documentation)
#' @name S4BuyseTest-getSurvival
#' @title Extract the Survival and Survival Jumps
#' @aliases getSurvival,S4BuyseTest-method
#' @include S4-BuyseTest.R
#'
#' @description Extract the survival and survival jumps.
#'
#' @param object an \R object of class \code{\linkS4class{S4BuyseTest}}, i.e., output of \code{\link{BuyseTest}}
#' @param type [character vector] the type of survival to be output. See details.
#' @param endpoint [integer/character vector] the endpoint for which the survival should be output.
#' @param strata [integer/character vector] the strata for which the survival should be output.
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
#' @rdname S4BuyseTest-getSurvival
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
