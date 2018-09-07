## * Documentation - getCount
#' @name BuyseRes-getCount
#' @title Extract the Number of Favorable, Unfavorable, Neutral, Uninformative pairs
#' @include BuyseRes-object.R
#' @aliases BuyseRes-getCount getCount getCount,BuyseRes-method
#'
#' @description Extract the number of favorable, unfavorable, neutral, uninformative pairs.
#'
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param type the type of pairs to be counted. Can be \code{"favorable"}, \code{"unfavorable"}, \code{neutral}, or \code{uninf}. Can also be \code{"all"} to select all of them.
#'
#' @return
#'   A \code{"vector"} containing the number of pairs
#'
#' @keywords get BuyseRes-method

## * Method - getCount
#' @rdname BuyseRes-getCount
setMethod(f = "getCount",
          signature = "BuyseRes",
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


## * Documentation - getIndividualScore
#' @name BuyseRes-getIndividualScore
#' @title Extract the Score of Each Pair
#' @include BuyseRes-object.R
#' @aliases BuyseRes-getIndividualScore getIndividualScore getIndividualScore,BuyseRes-method
#'
#' @description Extract the score of each pair.
#'
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
#' @param endpoint [integer/character vector] the endpoint for which the scores should be output.
#' @param strata [integer/character vector] the strata for which the scores should be output.
#' @param rm.withinStrata [logical] should the columns indicating the position of each member of the pair
#' within each treatment group be removed?
#' @param rm.weight [logical] should the column weight be remove from the output?
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
#' \item \code{"weight"}: the residual weight of the pair to be analysed at the current outcome. Each pair starts with a weight of 1.
#' }
#' Note that the \code{.T} and \code{.C} may change since they correspond of the label of the treatment and control arms.
#' @keywords get BuyseRes-method

## * Method - getIndividualScore
#' @rdname BuyseRes-getIndividualScore
setMethod(f = "getIndividualScore",
          signature = "BuyseRes",
          definition = function(object, endpoint = NULL, strata = NULL,
                                rm.withinStrata = TRUE, rm.weight = FALSE,
                                unlist = TRUE, trace = 1){

              if(length(object@tableIndividualScore)==0){
                  if(trace){
                      cat("Survival was not exported from the object \n",
                          "Consider setting the argument \'keep.survival\' to \"TRUE\" in BuyseTest.options \n", sep = "")
                  }
                  return(invisible(NULL))
              }else{
                  out <- data.table::copy(object@tableIndividualScore)

                  endpoint.names <- names(object@Delta.netChance)
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
                          if(length(strata)==1 && unlist == TRUE){
                              out[[iEndpoint]][,c("strata") := NULL]
                          }
                      }

                  }

                  old.names <- c("index.T", "index.C", "indexWithinStrata.T","indexWithinStrata.C")
                  new.names <- c(paste0("index.",object@level.treatment), paste0("indexWithinStrata.",object@level.treatment))

                  for(iEndpoint in 1:length(out)){ ## iEndpoint <- 1
                      if(rm.withinStrata){
                          out[[iEndpoint]][,c("indexWithinStrata.T","indexWithinStrata.C") := NULL]
                          setnames(out[[iEndpoint]], old = old.names[1:2], new = new.names[1:2])
                      }else{
                          setnames(out[[iEndpoint]], old = old.names, new = new.names)
                      }
                      if(rm.weight){
                          out[[iEndpoint]][,c("weight") := NULL]
                      }
                  }

                  
                  if(length(endpoint) == 1 && unlist == TRUE){
                      out <- out[[1]] 
                  }

                  return(out[])
              }
          })

## * Documentation - getSurvival
#' @name BuyseRes-getSurvival
#' @title Extract the Survival and Survival Jumps
#' @include BuyseRes-object.R
#' @aliases BuyseRes-getSurvival getSurvival getSurvival,BuyseRes-method
#'
#' @description Extract the survival and survival jumps.
#'
#' @param object an \R object of class \code{\linkS4class{BuyseRes}}, i.e., output of \code{\link{BuyseTest}}
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
#' }
#'
#' @keywords get BuyseRes-method

## * Method - getSurvival
#' @rdname BuyseRes-getSurvival
setMethod(f = "getSurvival",
          signature = "BuyseRes",
          definition = function(object, type = NULL, endpoint = NULL, strata = NULL, unlist = TRUE, trace = TRUE){
              
              if(length(object@tableSurvival)==0){
                  
                  if(trace>0){
                      if(all(tolower(object@type)!="timetoevent")){
                          add.txt <- "No endpoint of type time to event \n"
                      }else if(tolower(object@method.tte$method)!="peron"){
                          add.txt <- "Consider setting the argument \'method.tte\' to \"Peron\" when calling BuyseTest \n"
                      }else{
                          add.txt <- "Consider setting the argument \'keep.survival\' to \"TRUE\" in BuyseTest.options \n"
                      }
                      cat("Survival was not exported from the object \n",
                          add.txt, sep = "")    
                  }
                  return(invisible(NULL))
              }else{

                  if(is.null(type)){
                      type <- c("survTimeC","survTimeT","survJumpC","survJumpT")
                  }else{
                      validCharacter(type, valid.length = NULL, refuse.duplicates = TRUE,
                                     valid.values = c("survTimeC","survTimeT","survJumpC","survJumpT"))
                  }
                  if(!is.null(type)){
                      out <- object@tableSurvival[type]
                  }else{
                      out <- data.table::copy(object@tableSurvival)
                  }
                  
                  endpoint.names <- names(object@Delta.netChance)
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
                          }else{
                              if("survTimeC" %in% type){ out$survTimeC[[iEndpoint]] <- out$survTimeC[[iEndpoint]][strata] }
                              if("survTimeT" %in% type){ out$survTimeT[[iEndpoint]] <- out$survTimeT[[iEndpoint]][strata] }
                              if("survJumpC" %in% type){ out$survJumpC[[iEndpoint]] <- out$survJumpC[[iEndpoint]][strata] }
                              if("survJumpT" %in% type){ out$survJumpT[[iEndpoint]] <- out$survJumpT[[iEndpoint]][strata] }
                          }
                      }

                  }

                  if(length(endpoint) == 1 && unlist == TRUE){
                      if("survTimeC" %in% type){ out$survTimeC <- out$survTimeC[[1]] }
                      if("survTimeT" %in% type){ out$survTimeT <- out$survTimeT[[1]] }
                      if("survJumpC" %in% type){ out$survJumpC <- out$survJumpC[[1]] }
                      if("survJumpT" %in% type){ out$survJumpT <- out$survJumpT[[1]] }
                  }

                  if(length(type) == 1 && unlist == TRUE){
                      out <- out[[type]]
                  }
                  return(out)
              }

              
          })
