### prodlim-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  1 2019 (23:06) 
## Version: 
## Last-Updated: apr  1 2019 (23:14) 
##           By: Brice Ozenne
##     Update #: 12
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
#' @param object object A prodlim object
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
#' 
#' e.KM <- prodlim(Hist(eventtime,status)~Treatment, data = dt)
#' ## iidProdlim(e.KM)
#' 
#'

## * iidProdlim - code
#' @rdname iidProdlim
#' @export
iidProdlim <- function(object){

    ## ** extract elements from object
    is.strata <- !is.null(object$X)

    if(is.strata){
        n.strata <- NROW(object$X)
    }else{
        n.strata <- 1
    }
    
    ## ** Extract new observations
    vec.time <-  names(object)
    object$time
  
    ## ** Compute quantities of interest
    ## baseline hazard
    object$hazard
    lambda0 <- predictCox(object,
                          type = "hazard",
                          centered = FALSE,
                          keep.strata = TRUE)
    etime1.min <- rep(NA, nStrata)
    
    ## S0, E, jump times
    object.index_strata <- list() 
    object.order_strata <- list()
  
    object.eXb_strata <- list()
    object.LPdata_strata <- list()
    object.status_strata <- list()
    object.time_strata <- list()
  
    new.index_strata <- list()
    new.order_strata <- list()
  
    new.eXb_strata <- list()
    new.LPdata_strata <- list()
    new.status_strata <- list()
    new.time_strata <- list()
  
    Ecpp <- list()
    new.indexJump <- list()
    new.order <- NULL
  
    for(iStrata in 1:nStrata){

        ## reorder object data
        object.index_strata[[iStrata]] <- which(object.strata == object.levelStrata[iStrata])
        object.order_strata[[iStrata]] <- order(object.modelFrame[object.index_strata[[iStrata]], .SD$stop])

        indexTempo <- object.index_strata[[iStrata]][object.order_strata[[iStrata]]]
        object.eXb_strata[[iStrata]] <- object.eXb[indexTempo]
        if(nVar>0){
            object.LPdata_strata[[iStrata]] <- object.LPdata[indexTempo,,drop = FALSE]
        }else{
            object.LPdata_strata[[iStrata]] <- matrix(nrow = 0, ncol = 0)
        }
        object.status_strata[[iStrata]] <- object.modelFrame[indexTempo, .SD$status]
        object.time_strata[[iStrata]] <- object.modelFrame[indexTempo, .SD$stop]
    
        ## reorder new data
        if(!is.null(newdata)){
            new.index_strata[[iStrata]] <- which(new.strata == object.levelStrata[iStrata])
            new.order_strata[[iStrata]] <- order(new.time[new.index_strata[[iStrata]]])
      
            new.eXb_strata[[iStrata]] <- new.eXb[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
            if(nVar>0){
                new.LPdata_strata[[iStrata]] <- new.LPdata[new.index_strata[[iStrata]][new.order_strata[[iStrata]]],,drop = FALSE]
            }else{
                new.LPdata_strata[[iStrata]] <- matrix(nrow = 0, ncol = 0)
            }
            new.status_strata[[iStrata]] <- new.status[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
            new.time_strata[[iStrata]] <- new.time[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
        }else{
            new.index_strata[[iStrata]] <- object.index_strata[[iStrata]]
            new.order_strata[[iStrata]] <- object.order_strata[[iStrata]]
      
            new.eXb_strata[[iStrata]] <- object.eXb_strata[[iStrata]]
            new.LPdata_strata[[iStrata]] <- object.LPdata_strata[[iStrata]]
            new.status_strata[[iStrata]] <- object.status_strata[[iStrata]]
            new.time_strata[[iStrata]] <- object.time_strata[[iStrata]]
        }
    
        ## E
        Ecpp[[iStrata]] <-  calcE_cpp(status = object.status_strata[[iStrata]], 
                                      eventtime = object.time_strata[[iStrata]],
                                      eXb = object.eXb_strata[[iStrata]],
                                      X = object.LPdata_strata[[iStrata]],
                                      p = nVar, add0 = TRUE)
    
        new.indexJump[[iStrata]] <- prodlim::sindex(Ecpp[[iStrata]]$Utime1, new.time) - 1
        # if event/censoring is before the first event in the training dataset 
        # then sindex return 0 thus indexJump is -1
        # the following 3 lines convert -1 to 0
        if(any(new.indexJump[[iStrata]]<0)){
            new.indexJump[[iStrata]][new.indexJump[[iStrata]]<0] <- 0
        }
    
        ## store order
        if(length(new.order>0)){
            new.order <- c(new.order, new.index_strata[[iStrata]][new.order_strata[[iStrata]]])
        }else{
            new.order <- new.index_strata[[iStrata]][new.order_strata[[iStrata]]]
        }
    
    }
    
    ## ** Computation of the influence function (coefficients)
    IFbeta <- NULL
    IFcumhazard <- NULL
    IFhazard <- NULL
    calcIFhazard <- list(delta_iS0 = NULL,
                         Elambda0 = NULL,
                         cumElambda0 = NULL,
                         lambda0_iS0= NULL,
                         cumLambda0_iS0= NULL,
                         time1 = NULL)
    ls.Utime1 <- NULL
  
    #### beta
    for(iStrata in 1:nStrata){
    
        new.indexJump_strata <- new.indexJump[[iStrata]][new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
    
        ## IF
        if(nVar > 0){
            if(store.iid != "approx"){
                IFbeta_tempo <- IFbeta_cpp(newT = new.time_strata[[iStrata]],
                                           neweXb = new.eXb_strata[[iStrata]],
                                           newX = new.LPdata_strata[[iStrata]],
                                           newStatus = new.status_strata[[iStrata]], 
                                           newIndexJump = new.indexJump_strata, 
                                           S01 = Ecpp[[iStrata]]$S0,
                                           E1 = Ecpp[[iStrata]]$E,
                                           time1 = Ecpp[[iStrata]]$Utime1,
                                           iInfo = iInfo,
                                           p = nVar)
            }else{
                IFbeta_tempo <- IFbetaApprox_cpp(newX = new.LPdata_strata[[iStrata]],
                                                 newStatus = new.status_strata[[iStrata]],
                                                 newIndexJump = new.indexJump_strata,  
                                                 E1 = Ecpp[[iStrata]]$E,
                                                 iInfo = iInfo,
                                                 p = nVar)
            }
        }else{
            IFbeta_tempo <- matrix(NA, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
        }
    
        ## output
        IFbeta <- rbind(IFbeta, IFbeta_tempo)
    
    }
    
  
    ## set original order
    IFbeta <- IFbeta[order(new.order),,drop=FALSE]
    ## name coefficients
    colnames(IFbeta) <- infoVar$lpvars
    # }}}
    
    # {{{ Computation of the influence function (baseline hazard)
    if(baseline.iid){
        for(iStrata in 1:nStrata){
            ## hazard
            if(nStrata==1){ # select only the time,lambda corresponding to the events and not censored observations
                timeStrata <- lambda0$time[lambda0$time %in% Ecpp[[1]]$Utime1]
                lambda0Strata <- lambda0$hazard[lambda0$time %in% Ecpp[[1]]$Utime1]
            }else{ # same within the strata
                index.strata <- which(lambda0$strata == object.levelStrata[iStrata])
                index.keep <- index.strata[lambda0$time[index.strata] %in% Ecpp[[iStrata]]$Utime1]
      
                timeStrata <- lambda0$time[index.keep]
                lambda0Strata <- lambda0$hazard[index.keep]
            }
    
            etime1.min[iStrata] <- timeStrata[1]
    
            ## tau.hazard
            if(is.null(tau.hazard)){
                tau.hazard_strata <- object.time_strata[[iStrata]][object.status_strata[[iStrata]] == 1]
            }else if(is.list(tau.hazard)){
                tau.hazard_strata <- tau.hazard[[nStrata]]
            }else{
                tau.hazard_strata <- tau.hazard
            }
    
            ## E
            nUtime1_strata <- length(Ecpp[[iStrata]]$Utime1)
            if(nVar > 0){
                Etempo <- Ecpp[[iStrata]]$E[-NROW(Ecpp[[iStrata]]$E),,drop = FALSE]
            }else{
                Etempo <- matrix(0, ncol = 1, nrow = nUtime1_strata-1)
            }
    
            ## IF
            if(any(new.status_strata[[iStrata]]>0)){
                IFlambda_res <- IFlambda0_cpp(tau = tau.hazard_strata,
                                              IFbeta = IFbeta,
                                              newT = new.time, neweXb = new.eXb, newStatus = new.status, newIndexJump = new.indexJump[[iStrata]], newStrata = as.numeric(new.strata),
                                              S01 = Ecpp[[iStrata]]$S0,
                                              E1 = Etempo,
                                              time1 = timeStrata, lastTime1 = Ecpp[[iStrata]]$Utime1[nUtime1_strata], # here lastTime1 will not correspond to timeStrata[length(timeStrata)] when there are censored observations
                                              lambda0 = lambda0Strata,
                                              p = nVar, strata = iStrata,
                                              exact = (store.iid!="approx"), minimalExport = (store.iid=="minimal")
                                              )      
            }else{
                if(length(tau.hazard_strata)==0){tau.hazard_strata <- max(object.time_strata[[iStrata]])}
                IFlambda_res <- list(hazard = matrix(0, ncol = length(tau.hazard_strata), nrow = NROW(IFbeta)),
                                     cumhazard = matrix(0, ncol = length(tau.hazard_strata), nrow = NROW(IFbeta))
                                     )
                if(length(tau.hazard_strata)==0){tau.hazard_strata <- NA}
            }
    
            # output 
            ls.Utime1 <- c(ls.Utime1, list(tau.hazard_strata))
            if(store.iid=="minimal"){
                calcIFhazard$delta_iS0 <- c(calcIFhazard$delta_iS0, list(IFlambda_res$delta_iS0))
                calcIFhazard$Elambda0 <- c(calcIFhazard$Elambda0, list(IFlambda_res$Elambda0))
                calcIFhazard$cumElambda0 <- c(calcIFhazard$cumElambda0, list(IFlambda_res$cumElambda0))
                calcIFhazard$lambda0_iS0 <- c(calcIFhazard$lambda0_iS0, list(IFlambda_res$lambda0_iS0))
                calcIFhazard$cumLambda0_iS0 <- c(calcIFhazard$cumLambda0_iS0, list(IFlambda_res$cumLambda0_iS0))
                calcIFhazard$time1 <- c(calcIFhazard$time1, list(timeStrata)) # event time by strata
            }else{
                if(keep.times){
                    colnames(IFlambda_res$hazard) <- tau.hazard_strata
                    colnames(IFlambda_res$cumhazard) <- tau.hazard_strata
                }
                IFhazard <- c(IFhazard, list(IFlambda_res$hazard))
                IFcumhazard <- c(IFcumhazard, list(IFlambda_res$cumhazard))
            }
        }
    }
    # }}}
  
    # {{{ export
    return(list(IFbeta = IFbeta,  
                IFhazard = IFhazard,
                IFcumhazard = IFcumhazard,
                calcIFhazard = calcIFhazard,
                time = ls.Utime1,  # time at which the IF is assessed
                etime1.min = etime1.min,
                etime.max = lambda0$lastEventTime,
                indexObs = new.order,
                store.iid = store.iid
                ))
    # }}}
}


##----------------------------------------------------------------------
### prodlim-iid.R ends here
