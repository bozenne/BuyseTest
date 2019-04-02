### prodlim-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  1 2019 (23:06) 
## Version: 
## Last-Updated: apr  2 2019 (16:08) 
##           By: Brice Ozenne
##     Update #: 66
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iidProdlim - documentation
#' @title Extract i.i.d. decomposition from a prodlim model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iidProdlim
#' 
#' @param object A prodlim object
#' 
#' @details
#' This function is a simplified version of the iidCox function of the riskRegression package.
#' Formula for the influence function can be found in (Ozenne et al., 2017).
#' 
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#'
#'

## * iidProdlim - examples
#' @rdname iidProdlim
#' @examples
#' library(data.table)
#' library(prodlim)
#' 
#' set.seed(10)
#' dt <- simBuyseTest(10)
#' setkeyv(dt, "Treatment")
#' 
#' e.KM <- prodlim(Hist(eventtime,status)~Treatment, data = dt)
#' iidProdlim(e.KM)
#'
#'
#' 
#' 
#'

## * iidProdlim - code
#' @rdname iidProdlim
#' @export
iidProdlim <- function(object){

    if(!inherits(object,"prodlim")){
        stop("Argument \'object\' must inherit from prodlim \n")
    }
    if(object$type!="surv"){
        stop("Influence function only available for survival models \n")
    }
    
    ## ** extract elements from object
    is.strata <- !is.null(object$X)
    strataVar <- names(object$X)
    n.strataVar <- NCOL(object$X)
    
    if(is.strata){
        n.strata <- NROW(object$X)
    }else{
        n.strata <- 1
    }
    level.strata <- as.character(interaction(object$X))
    vec.strata <- factor(interaction(object$model.matrix[object$originalDataOrder,,drop=FALSE]), levels = level.strata)
    vec.strataNum <- as.numeric(vec.strata)

    vec.eventtime <- object$model.response[object$originalDataOrder,1]
    vec.status <- object$model.response[object$originalDataOrder,2]
    
    ## ** Extract baseline hazard + number at risk
    ## baseline hazard
    df.strata <- do.call(rbind,lapply(1:n.strata, function(iS){
        M <- matrix(object$X[iS,], ncol = n.strataVar, nrow = object$size.strata[iS], byrow = TRUE,
                    dimnames = list(NULL, strataVar))
        cbind("strata.index" = iS, data.frame(M, stringsAsFactors = FALSE))
    }))
    tableHazard <- cbind(df.strata, hazard = object$hazard, survival = object$surv, time = object$time,
                         event = object$n.event,
                         atrisk = object$n.risk)
    
    tableHazard.red <- tableHazard[tableHazard$event>0,]
    n.times <- NROW(tableHazard.red)
    n.obs <- NROW(object$model.matrix)
    
    ## ** Computation of the influence function
    ## -\Ind[strata] \int(\lambda0/S0) - jump/S0)
    IFhazard <- vector(mode = "list", length = n.strata)
    IFcumhazard <- vector(mode = "list", length = n.strata)
    IFsurvival <- vector(mode = "list", length = n.strata)
    ls.Utime1 <- vector(mode = "list", length = n.strata)
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        iTableHazard <- tableHazard.red[tableHazard.red$strata.index == iStrata,]
        ls.Utime1[[iStrata]] <- iTableHazard$time
        iN.time <- length(ls.Utime1[[iStrata]])
        iHazard <- iTableHazard$hazard
        iSurvival <- iTableHazard$survival

        ## prepare
        IFhazard[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)
        IFcumhazard[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)
        IFsurvival[[iStrata]] <- matrix(0, nrow = n.obs, ncol = iN.time)

        ## only keep observation in the strata and with eventtime at or after the first jump
        iSubsetObs <- intersect(which(vec.strataNum==iStrata),
                                which(vec.eventtime>=min(ls.Utime1[[iStrata]])))
        iVec.eventtime <- vec.eventtime[iSubsetObs]
        iVec.status <- vec.status[iSubsetObs]
        
        iIndexJump <- prodlim::sindex(ls.Utime1[[iStrata]], iVec.eventtime)
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
        IFcumhazard[[iStrata]][iSubsetObs,] <- t(apply(IFhazard[[iStrata]][iSubsetObs,],1,cumsum))

        ## survival
        ## note use exp(-surv) instead of product limit for consistency with riskRegression
        IFsurvival[[iStrata]][iSubsetObs,] <- sweep(-IFcumhazard[[iStrata]][iSubsetObs,], FUN = "*", STATS = exp(-cumsum(iTableHazard$hazard)), MARGIN = 2)
        ## IFsurvival[[iStrata]][iSubsetObs,] <- sweep(-IFcumhazard[[iStrata]][iSubsetObs,], FUN = "*", STATS = iTableHazard$survival, MARGIN = 2)
    }

    
    ## ** Export
    return(list(IFhazard = IFhazard,
                IFcumhazard = IFcumhazard,
                IFsurvival = IFsurvival,
                time = ls.Utime1, 
                etime.max = as.double(tapply(tableHazard$time,tableHazard$strata.index, max)),
                label.strata = level.strata,
                X = object$X
                ))
}


##----------------------------------------------------------------------
### prodlim-iid.R ends here
