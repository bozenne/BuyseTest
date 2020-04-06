### prodlim-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  1 2019 (23:06) 
## Version: 
## Last-Updated: apr  6 2020 (10:58) 
##           By: Brice Ozenne
##     Update #: 104
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iid.prodlim - documentation
#' @title Extract i.i.d. decomposition from a prodlim model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid.prodlim
#' 
#' @param object A prodlim object.
#' @param add0 [logical] add the 0 to vector of relevant times.
#' @param ... not used. For compatibility with the generic method.
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
#' @rdname iid.prodlim
#' @export
iid.prodlim <- function(object, add0 = FALSE, ...){

    if(!inherits(object,"prodlim")){
        stop("Argument \'object\' must inherit from prodlim \n")
    }
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
    tableHazard <- data.table::data.table(df.strata, hazard = object$hazard, survival = object$surv, time = object$time,
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
        IFcumhazard[[iStrata]][iSubsetObs,] <- t(apply(IFhazard[[iStrata]][iSubsetObs,,drop=FALSE],1,cumsum))

        ## survival
        ## note use exp(-surv) instead of product limit for consistency with riskRegression
        IFsurvival[[iStrata]][iSubsetObs,] <- sweep(-IFcumhazard[[iStrata]][iSubsetObs,,drop=FALSE], FUN = "*", STATS = exp(-cumsum(iTableHazard$hazard)), MARGIN = 2)
        ## IFsurvival[[iStrata]][iSubsetObs,] <- sweep(-IFcumhazard[[iStrata]][iSubsetObs,], FUN = "*", STATS = iTableHazard$survival, MARGIN = 2)
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


##----------------------------------------------------------------------
### prodlim-iid.R ends here
