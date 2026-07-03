### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  1 2019 (23:06) 
## Version: 
## Last-Updated: jul  3 2026 (14:42) 
##           By: Brice Ozenne
##     Update #: 138
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * iid.auc
#' @title Extract the idd Decomposition for the AUC
#'
#' @description Extract the iid decompotion relative to AUC estimate.
#' 
#' @param x object of class \code{BuyseTestAUC} (output of the \code{auc} function).
#' @param ... not used. For compatibility with the generic function.
#'
#' @return A column vector.
#' 
#' @method iid BuyseTestAuc
#' @export
iid.BuyseTestAuc <- function(x,...){
    object <- x
    return(attr(object,"iid")[,"global",drop=FALSE])
}

## * iid.BuyseTestBrier
##' @title Extract the idd Decomposition for the Brier Score
##'
##' @description Extract the iid decompotion relative to Brier score estimate.
##' 
##' @param x object of class \code{BuyseTestBrier} (output of the \code{brier} function).
##' @param ... not used. For compatibility with the generic function.
##'
##' @return A column vector.
##' 
##' @method iid BuyseTestBrier
##' @keywords methods
##' 
##' @export
iid.BuyseTestBrier <- function(x,...){
    object <- x
    return(attr(object,"iid"))
}


## * iid.BuyseTTEM
#' @export
iid.BuyseTTEM <- function(x, treatment, strata, cause = 1, ...){

    ## ** normalize input
    object <- x

    ## *** treatment
    if(length(treatment)!=1){
        stop("Argument \'treatment\' should have length 1. \n")
    }
    if(is.numeric(treatment)){
        if(treatment %in% 1:2 == FALSE){
            stop("Argument \'treatment\' should be 1 or 2 when numeric. \n")
        }
    }else if(is.character(treatment) || is.factor(treatment)){
        if(treatment %in% object$peron$level.treatment == FALSE){
            stop("Argument \'treatment\' should be \"",object$peron$level.treatment[1],"\" or \"",object$peron$level.treatment[2],"\" when character/factor. \n")
        }
        treatment <- which(treatment == object$peron$level.treatment)
    }

    ## *** strata
    if(missing(strata)){
        if(object$peron$n.strata == 1){
            strata <- 1
        }else{
            stop("Argument \'strata\' argument is missing. \n",
                 "Needed when the time to event model contains covariates. \n")
        }
    }else{
        if(inherits(strata,"data.frame")){

            if(NROW(strata)!=1){
                stop("Argument \'strata\' should have a single row when being a data.frame. \n")
            }
            if(any(attr(object$peron$strata.var,"original") %in% names(strata) == FALSE)){
                stop("Missing column(s) \"",paste(setdiff(attr(object$peron$strata.var,"original"), names(strata)), collapse = "\", \""),"\" in argument \'strata\'. \n")
            }            
            if(length(attr(object$peron$strata.var,"original"))==1){
                strata <- strata[[attr(object$peron$strata.var,"original")]]
            }else{
                strata <- interaction(as.data.frame(strata)[attr(object$peron$strata.var,"original")], sep = ".")
            }

            if(is.numeric(strata)){
                if(strata %in% 1:object$peron$n.strata == FALSE){
                    stop("Unknown strata defined by argument \'strata\'. \n",
                         "Should be one of ",paste(1:object$peron$n.strata, collapse = ", ")," \n")
                }
            }else if(strata %in% object$peron$level.strata == FALSE){
                stop("Unknown strata defined by argument \'strata\'. \n",
                     "Should be one of ",paste(object$peron$level.strata, collapse = ", ")," \n")
            }else{
                strata <- which(object$peron$level.strata==strata)
            }
        }else{

            if(length(strata)!=1){
                stop("Argument \'strata\' should either be a data.frame specifying the strata variables \n",
                     "or have length 1, e.g. be an integer indexing the strata. \n")
            }
            
            if(is.numeric(strata)){
                if(strata %in% 1:object$peron$n.strata == FALSE){
                    stop("When numeric argument \'strata\' should be one of ",paste(1:object$peron$n.strata, collapse = ", ")," \n")
                }
            }else if(is.character(strata) || is.factor(strata)){
                if(strata %in% object$peron$level.strata == FALSE){
                    stop("When character/factor argument \'strata\' should be one of ",paste(object$peron$level.strata, collapse = ", ")," \n")
                }
                strata <- which(object$peron$level.strata==strata)
            }else{
                stop("Unknown type for argument \'strata\': should be numeric/character/factor of length 1 or a data.frame. \n")
            }
        }
    }

    ## *** cause
    if(length(cause)!=1){
        stop("Argument \'cause\' should have length 1. \n")
    }

    if(is.numeric(cause)){
        if(cause %in% 1:object$peron$n.CR == FALSE){
            stop("Argument \'cause\' should be ",paste(1:object$peron$n.CR, collapse = ", "),". \n")
        }
    }else{
        stop("Unknown type for argument \'cause\': should be numeric. \n")
    }

    ## ** extract iid
    out <- object$peron$iid.cif[[strata]][[treatment]][[cause]]
    if(is.null(out)){
        stop("iid decomposition not available - consider setting the argument \'iid\' to TRUE when calling BuyseTTEM. \n")
    }
    type.CR <- ifelse(object$peron$n.CR==1,"survival","competing.risks")

    ## ** export
    if(type.CR=="survival"){
        return(-out)
    }else{
        return(out)
    }
}

## * iid.prodlim - documentation
#' @title Extract i.i.d. decomposition from a prodlim model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid.prodlim
#' 
#' @param x A prodlim object.
#' @param add0 [logical] add the 0 to vector of relevant times.
#' @param ... not used. For compatibility with the generic method.
#' 
#' @details
#' This function is a simplified version of the iidCox function of the riskRegression package.
#' Formula for the influence function can be found in (Ozenne et al., 2017).
#' 
#' @return A list containing:
#' \itemize{
#'  \item IFbeta: Influence function for the regression coefficient.
#'  \item IFhazard: Time differential of the influence function of the hazard.
#'  \item IFcumhazard: Influence function of the cumulative hazard.
#'  \item time: Times at which the influence function has been evaluated.
#'  \item etime.max: Last observation time (i.e. jump or censoring) in each strata.
#'  \item label.strata: Strata to which each observation belong.
#'  \item X: Design matrix.
#'  \item table: Hazard at each time for each strata.
#' }
#'
#' @keywords methods
#' 
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#'
#' @author Brice Ozenne

## * iid.prodlim - examples
#' @rdname iid.prodlim
#' @examples
#' library(data.table)
#' library(prodlim)
#' 
#' set.seed(10)
#' dt <- simBuyseTest(10)
#' setkeyv(dt, "treatment")
#' 
#' e.KM <- prodlim(Hist(eventtime,status)~treatment, data = dt)
#' lava::iid(e.KM)

## * iid.prodlim - code
#' @export
iid.prodlim <- function(x, add0 = FALSE, ...){
    object <- x
    
    if(object$type!="surv"){
        stop("Influence function only available for survival models \n")
    }
    
    ## ** extract elements from object
    X <- object$X
    is.strata <- !is.null(X)
    strataVar <- names(X)
    n.strataVar <- NCOL(X)
    
    if(is.strata){
        n.strata <- NROW(X)
    }else{
        n.strata <- 1
    }
    level.strata <- as.character(interaction(X))
    vec.strata <- factor(interaction(object$model.matrix[object$originalDataOrder,,drop=FALSE]), levels = level.strata)

    vec.strataNum <- as.numeric(vec.strata)
    vec.eventtime <- object$model.response[object$originalDataOrder,1]
    vec.status <- object$model.response[object$originalDataOrder,2]
    
    ## ** Extract baseline hazard + number at risk
    ## baseline hazard
    df.strata <- do.call(rbind,lapply(1:n.strata, function(iS){
        M <- matrix(X[iS,], ncol = n.strataVar, nrow = object$size.strata[iS], byrow = TRUE,
                    dimnames = list(NULL, strataVar))
        cbind("strata.index" = iS, data.frame(M, stringsAsFactors = FALSE))
    }))
    tableHazard <- data.table::data.table(df.strata,
                                          hazard = object$hazard,
                                          survival = object$surv,
                                          time = object$time,
                                          event = object$n.event,
                                          atrisk = object$n.risk)
    tableHazard.red <- tableHazard[tableHazard$event>0]
    if(add0){
        tableHazard0 <- tableHazard[,.SD[1], by = "strata.index"]
        tableHazard0[,c("hazard","survival","time","event","atrisk") := list(rep(0,n.strata),
                                                                             rep(1,n.strata),
                                                                             rep(0,n.strata),
                                                                             rep(0,n.strata),
                                                                             as.double(table(vec.strata))
                                                                             )]
        tableHazard.red <- rbind(tableHazard0,tableHazard.red)
        data.table::setkeyv(tableHazard.red, c("strata.index","time"))
    }    
    n.times <- NROW(tableHazard.red)
    n.obs <- NROW(object$model.matrix)
    
    ## ** Computation of the influence function
    ## -\Ind[strata] \int(\lambda0/S0) - jump/S0)
    IFhazard <- vector(mode = "list", length = n.strata)
    IFcumhazard <- vector(mode = "list", length = n.strata)
    IFsurvival <- vector(mode = "list", length = n.strata)
    IFcif <- vector(mode = "list", length = n.strata)
    ls.Utime1 <- vector(mode = "list", length = n.strata)
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        iTableHazard <- tableHazard.red[tableHazard.red$strata.index == iStrata]
        ls.Utime1[[iStrata]] <- iTableHazard$time        
        iN.time <- length(ls.Utime1[[iStrata]])
        iHazard <- iTableHazard$hazard
        iSurvival <- iTableHazard$survival

        ## prepare
        IFhazard[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)
        IFcumhazard[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)
        IFsurvival[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)
        IFcif[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)

        ## only keep observation in the strata and with eventtime at or after the first jump
        iSubsetObs <- intersect(which(vec.strataNum==iStrata),
                                which(vec.eventtime>=min(ls.Utime1[[iStrata]])))        
        iVec.eventtime <- vec.eventtime[iSubsetObs]
        iVec.status <- vec.status[iSubsetObs]

        iIndexJump <- prodlim::sindex(jump.times = ls.Utime1[[iStrata]], eval.times = iVec.eventtime)
        iDelta_iS0 <- iVec.status / iTableHazard$atrisk[iIndexJump]
        
        ## hazard
        iHazard_iS0 <- iHazard/iTableHazard$atrisk

        iIndEvent <- do.call(cbind, lapply(ls.Utime1[[iStrata]], function(iT){
            (abs(iT - iVec.eventtime ) < 1e-12) * iDelta_iS0
        }))
        iRatio <- do.call(cbind, lapply(1:iN.time, function(iT){
            (iT <= iIndexJump) * iHazard_iS0[iT]
        }))
        IFhazard[[iStrata]][iSubsetObs,] <- - iRatio + iIndEvent
        
        ## cumulative hazard
        IFcumhazard[[iStrata]][iSubsetObs,] <- .rowCumSum_cpp(IFhazard[[iStrata]][iSubsetObs,,drop=FALSE])

        ## survival
        ## note use exp(-surv) instead of product limit for consistency with riskRegression
        ## IFsurvival[[iStrata]][iSubsetObs,] <- .rowMultiply_cpp(-IFcumhazard[[iStrata]][iSubsetObs,,drop=FALSE], scale = exp(-cumsum(iTableHazard$hazard)))
        IFsurvival[[iStrata]][iSubsetObs,] <- .rowMultiply_cpp(-IFcumhazard[[iStrata]][iSubsetObs,,drop=FALSE], scale = iTableHazard$survival)
        IFcif[[iStrata]][iSubsetObs,] <- 1 - IFsurvival[[iStrata]][iSubsetObs,]
    }
    
    ## ** Modification used by BuyseTest to enable the user to easily specify model.tte
    if(!is.null(object$XX) && !identical(object$X,object$XX)){
        
        oldlevel.strata <- level.strata
        level.strata <- as.character(interaction(object$XX))

        X <- object$XX

        index.strata <- match(tableHazard.red[,interaction(.SD),.SDcols = names(object$X)],
                              oldlevel.strata)
        tableHazard.red[, c(names(object$X)) := NULL]        
        tableHazard.red <- cbind(tableHazard.red[,.SD, .SDcols = "strata.index"],
                                 object$XX[index.strata,,drop=FALSE],
                                 tableHazard.red[,.SD, .SDcols = c("hazard","survival","time","event","atrisk")])
        
    }
    
    ## ** Export
    return(list(IFhazard = IFhazard,
                IFcumhazard = IFcumhazard,
                IFsurvival = IFsurvival,
                time = ls.Utime1, 
                etime.max = tableHazard[,max(.SD$time),by = "strata.index"][[2]],
                label.strata = level.strata,
                X = X,
                table = tableHazard.red
                ))
}

## * iid.S3sensitivity
#' @export
iid.S3sensitivity <- function(x, ...){
    return(attr(x,"iid"))
}

## * getIid S4BuyseTest
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

              if(type %in% c("all","nuisance") && (tolower(object@scoring.rule)!="gehan")){
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


##----------------------------------------------------------------------
### iid.R ends here
