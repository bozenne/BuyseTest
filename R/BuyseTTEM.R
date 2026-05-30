### BuyseTTEM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2020 (12:15) 
## Version: 
## Last-Updated: May 31 2026 (00:17) 
##           By: Brice Ozenne
##     Update #: 1018
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * BuyseTTEM (documentation)
#' @title Time to Event Model
#' @name BuyseTTEM
#' 
#' @description Pre-compute quantities of a time to event model useful for predictions.
#' Only does something for prodlim objects.
#' 
#' @param object time to event model. 
#' @param treatment [list or character] List containing the value of the treatment variable for each individual or name of the treatment variable.
#' The former is necessary when the survival model do not dependent on treatment.
#' @param iid [logical] Should the iid decomposition of the predictions be output. 
#' @param iid.surv [character] Estimator of the survival used when computing the influence function.
#' Can be the product limit estimator (\code{"prodlim"}) or an exponential approximation (\code{"exp"}, same as in \code{riskRegression::predictCoxPL}).
#' @param n.grid [integer, >0] Number of timepoints used to discretize the time scale. Not relevant for prodlim objects.
#' @param ... For compatibility with the generic method.
#'
#' @return An S3 object of class \code{BuyseTTEM}.
#' @keywords models
#' 
#' @examples
#' library(prodlim)
#' library(survival)
#' library(data.table)
#'
#' tau <- seq(0,3,length.out=10)
#' 
#' #######################
#' #### survival case ####
#' #######################
#' 
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#' df.Ta <- data.frame(treatment = "T", strata = "a")
#' df.Ca <- data.frame(treatment = "C", strata = "a")
#' 
#' #### Non-parametric model (stratified Kapan-Meier) ####
#'
#' ## estimate
#' eKM.prodlim <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.data)
#' ## plot(eKM.prodlim)
#'
#' ## convert to TTEM object
#' eKM.TTEM <- BuyseTTEM(eKM.prodlim, treatment = "treatment", iid = TRUE)
#'
#' ## extract jump times
#' jumpKM.Ta <- predict(eKM.TTEM, type = "jump", treatment = "T", strata = "a")
#' jumpKM.Ta
#' ## extract last observation time
#' last.timeKM.Ta <- predict(eKM.TTEM, type = "last.time", treatment = "T", strata = "a")
#' last.timeKM.Ta
#' ## extract survival at last observation time
#' last.survKM.Ta <- predict(eKM.TTEM, type = "last", treatment = "T", strata = "a")
#' last.survKM.Ta
#' ## extract jump size
#' changeKM.Ta <- predict(eKM.TTEM, type = "change", treatment = "T", strata = "a")
#' last.survKM.Ta-sum(changeKM.Ta$survival) ## sum to 1
#' ## prediction
#' predict(eKM.TTEM, time = tau, treatment = "T", strata = "a")
#' ## same as prodlim: predict(eKM.prodlim, times = tau, newdata = df.Ta)
#' predict(eKM.TTEM, time = tau, treatment = "C", strata = "a")
#' ## fix bug in prodlim: predict(eKM.prodlim, times = tau, newdata = df.Ca)
#' ## survival should be 0 if the last observation is an event
#' 
#' #### Parametric model ####
#' dist <- "lognormal"
#' 
#' ## estimate
#' eWei.survreg <- survreg(Surv(eventtime,status)~treatment+strata, data = df.data,
#'                         dist = dist)
#' ## gridCDF <- seq(.01,.99,by=.01)
#' ## plot(x = predict(eWei.survreg, newdata = df.Ta, type = "quantile", p = gridCDF,
#' ##                  y = rev(gridCDF), col = "red", type = "l")
#' ## points(x = predict(eWei.survreg, newdata = df.Ca, type = "quantile", p = gridCDF,
#' ##                  y = rev(gridCDF), col = "blue", type = "l")
#'
#' ## conver to TTEM object
#' eWei.TTEM <- BuyseTTEM(eWei.survreg, treatment = "treatment", iid = TRUE, n.grid = 1e3)
#'
#' ## extract jump times
#' jumpWei.Ta <- predict(eWei.TTEM, type = "jump", treatment = "T", strata = "a")
#' jumpWei.Ta
#' ## extract last observation time
#' last.timeWei.Ta <- predict(eWei.TTEM, type = "last.time", treatment = "T", strata = "a")
#' last.timeWei.Ta
#' ## extract survival at last observation time
#' last.survWei.Ta <- predict(eWei.TTEM, type = "last", treatment = "T", strata = "a")
#' last.survWei.Ta
#' ## extract jump size
#' changeWei.Ta <- predict(eWei.TTEM, type = "change", treatment = "T", strata = "a")
#' last.survWei.Ta-sum(changeWei.Ta$survival) ## sum to 1
#' ## prediction
#' predWei.Ta <- predict(eWei.TTEM, time = tau, treatment = "T", strata = "a")
#' predWei.Ta
#' ## survreg: predict(eWei.survreg, newdata = df.Ta, type = "quantile", p = 1-predWei.Ta$survival)-tau
#' ## approximate retrieve tau, up to the grid resolution
#'
#' #### competing risk case ####
#' df.dataCR <- copy(df.data)
#' df.dataCR$status <- rbinom(NROW(df.dataCR), prob = 0.5, size = 2)
#' 
#' eAJ.prodlim <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.dataCR)
#' ## plot(eAJ.prodlim)
#' 
#' eAJ.TTEM <- BuyseTTEM(eAJ.prodlim, treatment = "treatment", iid = TRUE)
#' 
#' ## extract jump times
#' jumpAJ.Ta <- predict(eAJ.TTEM, type = "jump", treatment = "T", strata = "a")
#' jumpAJ.Ta
#' ## extract last observation time
#' last.timeAJ.Ta <- predict(eAJ.TTEM, type = "last.time", treatment = "T", strata = "a")
#' last.timeAJ.Ta
#' ## extract survival at last observation time
#' last.cifAJ.Ta <- predict(eAJ.TTEM, type = "last", treatment = "T", strata = "a")
#' last.cifAJ.Ta
#' ## extract jump size
#' changeAJ.Ta <- predict(eAJ.TTEM, type = "change", treatment = "T", strata = "a")
#' last.cifAJ.Ta-sum(changeAJ.Ta$cif) ## sum to 0
#' ## prediction
#' predict(eAJ.TTEM, time = tau, treatment = "T", strata = "a")
#' ## same as prodlim: predict(eAJ.prodlim, times = tau, newdata = df.Ta, cause = 1)
#' 
#' predict(eAJ.TTEM, time = tau, treatment = "C", strata = "a")
#' ## same as prodlim: predict(eAJ.prodlim, times = tau, newdata = df.Ca, cause = 1)
#' @export
`BuyseTTEM` <-
    function(object,...) UseMethod("BuyseTTEM")


## * BuyseTTEM.default
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.formula <- function(object, treatment, iid, iid.surv = "exp", ...){
    e.prodlim <- prodlim(object, ...)
    return(BuyseTTEM(e.prodlim, treatment = treatment, iid = iid, iid.surv = iid.surv))
}

## * BuyseTTEM.prodlim
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.prodlim <- function(object, treatment, iid, iid.surv = "exp", ...){

    tol12 <- 1e-12
    tol11 <- 1e-11
    dots <- list(...)

    ## ** check arguments
    ## *** object
    if(any(object$time<=0)){
        stop("Only handles strictly positive event times \n")
    }

    ## *** iid surv
    iid.surv <- match.arg(iid.surv, c("prodlim","exp"))

    ## *** ...
    ## levels are converted to numeric by BuyseTest
    ## this makes sure to keep track of the original labels
    level.treatment <- dots$level.treatment
    level.strata <- dots$level.strata

    ## ** prepare output
    test.CR <- switch(object$type,
                      "surv" = FALSE,
                      "risk" = TRUE)
    if(test.CR){
        n.CR <- length(object$cuminc)
    }else{
        n.CR <- 1
        object$cuminc <- list(1-object$surv)
        object$cause.hazard <- list(object$hazard)
    }

    ## Note: object$xlevels is not in the right order
    object$peron <- .initPeron(object = object,
                               X = object$X,
                               treatment = treatment,
                               level.treatment = dots$level.treatment,
                               level.strata = dots$level.strata,
                               xlevels = NULL,
                               n.CR = n.CR)

    treatment <- object$peron$treatment
    level.treatment <- object$peron$level.treatment

    level.strata <- object$peron$level.strata
    n.strata <- object$peron$n.strata
    strata.var <- object$peron$strata.var
    
    ## ** compute start/stop indexes per strata
    ## position in the results of the pair (treatment,strata)
    index.start <- matrix(NA, nrow = n.strata, ncol = 2, dimnames = list(NULL,level.treatment))
    index.stop <- matrix(NA, nrow = n.strata, ncol = 2, dimnames = list(NULL,level.treatment))

    if(is.null(object$X)){ ## no covariate at all in the survival model

        index.start[] <- object$first.strata
        index.stop[] <- object$first.strata + object$size.strata - 1

    }else if(treatment %in% names(object$X)){ ## survival model includes treatment

        iIndexX.C <- which(object$peron$X[[treatment]]==0)
        iIndexX.T <- which(object$peron$X[[treatment]]==1)

        iMindexX <- do.call(rbind,lapply(1:n.strata, function(iStrata){ ## iStrata <- 1
            if(n.strata>1){
                iIndexS <- which(object$peron$X[["..strata.."]] == iStrata)
                return(c(C = intersect(iIndexX.C,iIndexS),
                         T = intersect(iIndexX.T,iIndexS)))
            }else{
                return(c(C = iIndexX.C,
                         T = iIndexX.T))
            }        
        }))

        index.start[] <- object$first.strata[iMindexX]
        index.stop[] <- object$first.strata[iMindexX] + object$size.strata[iMindexX] - 1

    }else{ ## stratified survival model but not on treatment

        index.start[] <- cbind(object$first.strata, object$first.strata)
        index.stop[] <- cbind(object$first.strata + object$size.strata - 1,
                              object$first.strata + object$size.strata - 1)

    }

    ## ** find last CIF value
    object$peron$last.time[,1] <- object$time[index.stop[,1]]
    object$peron$last.time[,2] <- object$time[index.stop[,2]]
    for(iCause in 1:n.CR){ ## iCause <- 1
        object$peron$last.estimate[,1,iCause] <- object$cuminc[[iCause]][index.stop[,1]]
        object$peron$last.estimate[,2,iCause] <- object$cuminc[[iCause]][index.stop[,2]]
    }
    
    ## ** table CIF 
    index.allJump <- sort(unlist(lapply(object$cause.hazard, function(iVec){which(iVec>0)})))
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        object$peron$jumpSurvHaz[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)

        for(iTreat in 1:2){ ## iTreat <- 1
            iMissingCIF <- sum(object$peron$last.estimate[iStrata,iTreat,])<(1-tol12)

            ## jump time (for all causes) in this strata
            iIndex.jumpObs <- intersect(index.allJump, index.start[iStrata,iTreat]:index.stop[iStrata,iTreat])
            if(length(iIndex.jumpObs)==0){
                object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- data.frame(index.obs = NA,
                                                                            index.jump = 1,
                                                                            time.jump = 0,
                                                                            survival = 1,
                                                                            setNames(as.list(rep(0,n.CR)),paste0("hazard",1:n.CR))
                                                                            )
                                                                            
            }else{
                object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- data.frame(index.obs = iIndex.jumpObs,
                                                                            index.jump = 1:length(iIndex.jumpObs),
                                                                            time.jump = object$time[iIndex.jumpObs],
                                                                            survival = object$surv[iIndex.jumpObs],
                                                                            do.call(cbind,setNames(lapply(object$cause.hazard, function(iVec){iVec[iIndex.jumpObs]}), paste0("hazard",1:n.CR))))
            }
            object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]][order(object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$time.jump),]
            iStrata.nJump <- length(iIndex.jumpObs)
            for(iEvent in 1:n.CR){ ## iEvent <- 1
                
                ## CIF at each jump (if any) and add time 0
                if(iStrata.nJump>0){
                    iStrata.extcuminc <- c(0,object$cuminc[[iEvent]][iIndex.jumpObs])
                    iStrata.exttime.jump <- c(-tol12,object$time[iIndex.jumpObs])
                }else{
                    iStrata.extcuminc <- 0
                    iStrata.exttime.jump <- -tol12
                }

                ## CIF just after the last observations
                if(iMissingCIF){
                    iStrata.extcuminc <- c(iStrata.extcuminc,NA)
                    iStrata.exttime.jump <- c(iStrata.exttime.jump, as.double(object$peron$last.time[iStrata,iTreat]) + tol11)
                }

                ##  index of the jumps
                if(iStrata.nJump>0){
                    iStrata.index.afterJump <- c(1,prodlim::sindex(jump.times = iStrata.exttime.jump, eval.times = object$time[iIndex.jumpObs] + tol12))
                    iStrata.obs <- c(NA,iIndex.jumpObs)
                }else{
                    iStrata.index.afterJump <- 1
                    iStrata.obs <- NA
                }
                if(iMissingCIF){
                    iStrata.index.afterJump <- c(iStrata.index.afterJump, NA)
                    iStrata.obs <- c(iStrata.obs, NA)
                }

                ## *** store
                object$peron$cif[[iStrata]][[iTreat]][[iEvent]] <- data.frame(time = iStrata.exttime.jump,
                                                                              cif = iStrata.extcuminc,
                                                                              index = iStrata.index.afterJump - 1,
                                                                              obs = iStrata.obs) ## move to C++ indexing
            }
        }
    }

    ## cif <- object$peron$cif[[1]][[1]][[1]]
    ## (cif[cif[,"index.cif.after"]+1,"cif"] - cif[cif[,"index.cif.before"]+1,"cif"]) - cif[,"dcif"]
    ## tail(object$peron$cif[[1]]$C[[1]])
    ## tail(object$peron$cif[[1]]$T[[1]])

    ## ** table iid
    if(iid){
        n.obs <- NROW(object$model.response)

        vec.eventtime <- object$model.response[object$originalDataOrder,"time"]
        if(test.CR){
            vec.status <- object$model.response[object$originalDataOrder,"event"] ## 1:n.CR event, n.CR+1 censoring
        }else{
            vec.status <- object$model.response[object$originalDataOrder,"status"] ## 0 censoring, 1 event
        }

        if(treatment %in% names(object$model.matrix)){
            model.matrix <- object$model.matrix[object$originalDataOrder,,drop=FALSE]
            model.matrix[[treatment]] <- factor(model.matrix[[treatment]], labels = level.treatment)
        }else{
            model.matrix <- NULL
        }
        if(is.null(attr(strata.var,"original"))){
            model.matrix <- cbind(model.matrix, "..strata.." = rep(1,n.obs))
        }else if(!identical(attr(strata.var,"original"),"..strata..")){
            model.matrix <- cbind(model.matrix, "..strata.." = as.numeric(factor(interaction(model.matrix[,strata.var,drop=FALSE]), levels = level.strata)))
        }

        object$peron$iid.cif <- setNames(vector(mode = "list", length=n.strata), level.strata)
        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            object$peron$iid.cif[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)

            for(iTreat in 1:2){ ## iTreat <- 1

                iHazard <- vector(mode = "list", length=n.CR)
                iHazard.iid <- vector(mode = "list", length=n.CR)
                object$peron$iid.cif[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)

                ## position of the subjects from the treatment group and strata in the dataset
                if(treatment %in% names(object$model.matrix)){
                    iIndexObs.strata <- which(model.matrix[[treatment]]==level.treatment[iTreat])
                }else{
                    iIndexObs.strata <- 1:n.obs
                }
                if(n.strata>1){
                    iIndexObs.strata <- intersect(iIndexObs.strata, which(model.matrix[["..strata.."]]==iStrata))
                }
                ## event times (censored or jump) of the subjects from the treatment group and strata in the dataset
                iTimeObs.strata <- vec.eventtime[iIndexObs.strata]
                ## type of event of the subjects from the treatment group and strata in the dataset
                iStatus.strata <- vec.status[iIndexObs.strata]
                ## jump time in the treatment group and strata
                iJump.time <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$time.jump
                iN.jump <- length(iJump.time)
                
                if(iN.jump==0){ ## no events: iid = 0
                    for(iCause in 1:n.CR){ ## iEvent <- 1
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iCause]] <- matrix(0, nrow = n.obs, ncol = 1)
                    }
                    next
                }
                ## matrix of at risk subjects over time
                iAtRisk <- do.call(cbind,lapply(iJump.time, function(iTime){
                    iTimeObs.strata >= iTime
                }))

                ## *** influence function for each cause-specific hazard
                for(iCause in 1:n.CR){ ## iCause <- 1                    
                    iHazard[[iCause]] <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]][[paste0("hazard",iCause)]]
                    iStatus <- do.call(cbind,lapply(iJump.time, function(iTime){
                        (abs(iTimeObs.strata-iTime)<tol11)*(iStatus.strata==iCause)
                    }))                                        
                    iHazard.iid[[iCause]] <- .rowScale_cpp(iStatus - .rowMultiply_cpp(iAtRisk, scale=iHazard[[iCause]]), scale = colSums(iAtRisk))
                }

                ## *** influence function for the overall survival
                if(iid.surv == "exp"){
                    iSurvival <- exp(-cumsum(rowSums(object$peron$jumpSurvHaz[[iStrata]][[iTreat]][,paste0("hazard",1:n.CR),drop=FALSE])))
                }else if(iid.surv == "prodlim"){
                    iSurvival <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$survival    
                }
                iSurvival.iid <- -.rowMultiply_cpp(.rowCumSum_cpp(Reduce("+",iHazard.iid)), iSurvival)

                ## *** influence function for the cumulative incidence
                for(iCause in 1:n.CR){ ## iCause <- 1
                    iParam <- na.omit(unique(object$peron$cif[[iStrata]][[iTreat]][[iCause]]$index))

                    object$peron$iid.cif[[iStrata]][[iTreat]][[iCause]] <- matrix(0, nrow = n.obs, ncol = length(iParam))
                    ## length(iParam)=1 corresponds to only time 0 which has IF 0, i.e. the initialization is enough
                    if(test.CR & length(iParam)>1){  
                        ## survival at t-
                        iSurvivalM1 <- c(1,iSurvival[1:(iN.jump-1)])
                        iSurvivalM1.iid <- cbind(0,iSurvival.iid[,1:(iN.jump-1),drop=FALSE])
                        ## add iid at time 0 and (if censoring) NA after the last event
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iCause]][iIndexObs.strata,] <- cbind(0,.rowCumSum_cpp(.rowMultiply_cpp(iSurvivalM1.iid,iHazard[[iCause]]) + .rowMultiply_cpp(iHazard.iid[[iCause]],iSurvivalM1)))
                    }else if(length(iParam)>1){
                        ## add iid at time 0 and (if censoring) NA after the last event
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iCause]][iIndexObs.strata,] <- cbind(0,-iSurvival.iid)
                    }
                }
            }}
    }

    ## ** export
    class(object) <- append("BuyseTTEM",class(object))
    return(object)
}

## * BuyseTTEM.survreg
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.survreg <- function(object, treatment, n.grid = NULL, iid, ...){

    tol12 <- 1e-12
    tol11 <- 1e-11    
    dots <- list(...)    
    if(is.null(n.grid)){ ## n.grid will sometimes be passed to BuyseTTEM.survreg as NULL when no value is provided by the user, this avoids to over-write the default
        n.grid <- BuyseTest.options()$n.grid
    }

    ## ** check arguments

    ## *** object
    if(any(object$y[,"time"]<=0)){
        stop("Only handles strictly positive event times \n")
    }

    ## *** ...
    level.treatment <- dots$level.treatment
    level.strata <- dots$level.strata

    ## ** prepare output
    mf <- stats::model.frame(object)
    object$peron <- .initPeron(object = object,
                               X = mf[,-1,drop=FALSE], ## first column for the Surv object,
                               treatment = treatment,
                               level.treatment = level.treatment,
                               level.strata = level.strata,
                               xlevels = object$xlevels,
                               n.CR = 1)

    n.strata <- object$peron$n.strata
    treatment <- object$peron$treatment
    level.treatment <- object$peron$level.treatment
    level.strata <- object$peron$level.strata

    if(treatment %in% names(mf)){
        if(is.null(object$xlevels[[treatment]])){
            mf[[treatment]] <- factor(mf[[treatment]], levels = sort(unique(mf[[treatment]])), labels = level.treatment)
        }else{
            mf[[treatment]] <- factor(mf[[treatment]], levels = object$xlevels[[treatment]], labels = level.treatment)
        }
    }
            
    ## ** handling competing risks
    object$peron$n.CR <- 1 ## survival case (one type of event)

    ## ** prepare for iid
    if(iid){

        ## *** extract information
        n.obs <- stats::nobs(object)
        X <- stats::model.matrix(stats::formula(object), mf)
        beta <- stats::coef(object)
        sigma <- object$scale
        object.iid <- lava::iid(object)
        object.iid.beta <- object.iid[,setdiff(colnames(object.iid), "logsigma"),drop=FALSE]
        if("logsigma" %in% colnames(object.iid)){
            object.iid.sigma <- object.iid[,"logsigma",drop=FALSE]
        }else{
            object.iid.sigma <- matrix(0, nrow = n.obs, ncol = 1, dimnames = list(1:n.obs,"logsigma"))
        }
        
        ## *** extract link and derivative
        object.dist <- survival::survreg.distributions[[object$dist]]
        if(is.null(object.dist$dist)){
            object.dist$dist <- object$dist
        }

        object.dist$quantileM1 <- switch(object.dist$dist,
                                      "t" = function(q,df){stats::pt(q,df)},
                                      "gaussian" = function(q,parms){stats::pnorm(q)},
                                      "logistic" = function(q,parms){1/(1+exp(-q))},
                                      "extreme" = function(q,parms){1-exp(-exp(q))},
                                      NULL)
        object.dist$dquantileM1 <- switch(object.dist$dist,
                                       "t" = function(q,df){stats::dt(q,df)},
                                       "gaussian" = function(q,parms){stats::dnorm(q)},
                                       "logistic" = function(q,parms){exp(-q)/(1+exp(-q))^2},
                                       "extreme" = function(q,parms){exp(q)*exp(-exp(q))},
                                       NULL)

        if("quantile" %in% names(object.dist) == FALSE){
            object.dist$quantile <- survival::survreg.distributions[[object.dist$dist]]$quantile
        }

        if("quantileM1" %in% names(object.dist) == FALSE){
            object.dist$quantileM1 <- switch(object.dist$dist,
                                             "t" = function(q,df){stats::pt(q,df)},
                                             "gaussian" = function(q,parms){stats::pnorm(q)},
                                             "logistic" = function(q,parms){1/(1+exp(-q))},
                                             "extreme" = function(q,parms){1-exp(-exp(q))},
                                             NULL)
            object.dist$dquantileM1 <- switch(object.dist$dist,
                                           "t" = function(q,df){stats::dt(q,df)},
                                           "gaussian" = function(q,parms){stats::dnorm(q)},
                                           "logistic" = function(q,parms){exp(-q)/(1+exp(-q))^2},
                                           "extreme" = function(q,parms){exp(q)*exp(-exp(q))},
                                           NULL)
        }
        
        if(is.null(object.dist$itrans)){
            object.dist$trans <- function(x){x}
            object.dist$itrans <- function(x){x}
            object.dist$dtrans <- function(x){1}
        }
    }

    ## ** table CIF and iid
    grid.quantile <- seq(from=0,to=1-1/n.grid,length.out=n.grid)
    object$peron$last.estimate[] <- utils::tail(grid.quantile,1) ## final modeled survival value close to 0 i.e. CIF close to 1
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        for(iTreat in level.treatment){ ## iTreat <- 1

            ## *** find a subject from the treatment and strata group
            if(treatment %in% names(mf)){
                ## restrict to cases? which(mf[,1][,2]==1)
                iIndex.obs <- which(mf[[treatment]]==iTreat)
            }else{
                iIndex.obs <- 1:n.obs
            }
            if(n.strata>1){
                iIndex.obs <- intersect(iIndex.obs,which(object$peron$X[,"..strata.."]==iStrata))
            }
            iNewdata <- mf[iIndex.obs[1],,drop=FALSE]
            ## *** get the timepoint corresponding to the grid of survival values
            iJump <- predict(object, newdata = iNewdata, p = grid.quantile, type = "quantile")

            object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- data.frame(index.obs = NA,
                                                                        index.jump = 1:length(grid.quantile),
                                                                        time.jump = iJump,
                                                                        survival = 1-grid.quantile,
                                                                        hazard1 = NA)
            
            object$peron$cif[[iStrata]][[iTreat]][[1]] <- data.frame(time = c(-tol12,iJump),
                                                                     cif = c(0,grid.quantile),
                                                                     index = 0:length(grid.quantile),
                                                                     obs = NA)

            object$peron$last.time[iStrata,iTreat] <- utils::tail(iJump,1)

            if(iid){
                iLP <- drop(X[iIndex.obs[1],,drop=FALSE] %*% beta)
                object.iid.iLP <- object.iid.beta %*% t(X[iIndex.obs[1],,drop=FALSE])
                
                ## *** compute time with se [not used]
                ## fit.trans <- iLP + object.dist$quantile(grid.quantile) * sigma
                ## range(object.dist$itrans(fit.trans) - iJump)
                ## fit.iid <- .rowScale_cpp(.colCenter_cpp(object.iid.sigma %*% object.dist$quantile(grid.quantile), center = - object.iid.iLP), object.dist$dtrans(object.dist$itrans(fit.trans)))
                ## fit.se <- sqrt(colSums(fit.iid^2))
                ## quantile(predict(object, type = "quantile", p = grid.quantile, newdata = mf[iIndex.obs[1],,drop=FALSE],se = TRUE)$se - fit.se, na.rm = TRUE)

                ## *** compute influence function
                iPred <- (object.dist$trans(iJump[iJump>0]) - iLP)/sigma ## range(grid.quantile - object.dist$quantileM1(iPred))

                object$peron$iid.cif[[iStrata]][[iTreat]][[1]] <- cbind(matrix(0, nrow = n.obs, ncol = 1+sum(iJump==0)),
                                                                        .rowMultiply_cpp(.colCenter_cpp(- object.iid.sigma %*% (object.dist$trans(iJump[iJump>0]) - iLP) / sigma^2,
                                                                                                        center = object.iid.iLP / sigma), scale = object.dist$dquantileM1(iPred)))
                ## range(object$peron$iid.cif[[iStrata]][[iTreat]][[1]][,10] - (- object.iid.sigma * (object.dist$trans(iJump[10]) - iLP) /sigma^2 - object.iid.iLP / sigma) * object.dist$dquantileM1(iPred[10]))
            }
        }
    }

    ## ** export
    class(object) <- append("BuyseTTEM",class(object))
    return(object)
}


## * BuyseTTEM.BuyseTTEM
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.BuyseTTEM <- function(object, ...){
    return(object)
}

## * .initPeron
##' @param object survival model
##' @param X [data.frame] containing the stratification variables
##' @param treatment [character] variable identifying the treatment variable in X
##' @param n.CR [integer] number of competing risks
##' @param levels.treatment [character vector] possible values of the treatment variable (optional)
##' @param levels.strata [character vector] possible values of the strata variable (optional)
##' @param xlevels [list] order of the levels of each factor variable (optional)
##'
##' @noRd
.initPeron <- function(object, data, X, n.CR, treatment = NULL, level.treatment = NULL, level.strata = NULL, xlevels = NULL){

    ## ** recover treatment if possible
    if(is.null(treatment)){
        if(NCOL(X)==1){
            treatment <- names(X)[1]
        }else if(NCOL(object$X)==0){
            stop("Argument \'treatment\' is missing (cannot guess it when there is no covariate in the survival model). \n")
        }else if(NCOL(object$X)>1){
            stop("Argument \'treatment\' is missing (cannot guess it when there is more than one covariate in the survival model). \n")
        }
    }else if(length(treatment)!=1){
        stop("Argument \'treatment\' should have length 1. \n")
    }

    ## ** recover level.treatment if possible
    if(is.null(level.treatment)){
        if(!is.null(xlevels)){
            level.treatment <- xlevels[[treatment]]
        }else if(!is.null(X) && treatment %in% names(X)){
            if(is.factor(X[[treatment]])){
                level.treatment <- levels(X[[treatment]])
            }else{
                level.treatment <- unique(X[[treatment]])
            }          
        }else{
            data <- eval(object$call$data)
            if(treatment %in% names(data) == FALSE){
                stop("Argument \'treatment\' does not match any of the column in the original dataset. \n",
                     "Argument value: ",treatment,".\n",
                     "Possible columns: \"",paste(names(data), collapse = "\", \""),"\"\n")
            }
            if(is.factor(data[[treatment]])){
                level.treatment <- levels(data[[treatment]])
            }else{
                level.treatment <- sort(unique(data[[treatment]]))
            } 
        }
    }

    ## ** recover level.treatment if possible
    if(is.null(X)){
        X <- as.data.frame(stats::setNames(list(level.treatment),treatment))
    }else if(treatment %in% names(X) == FALSE){
        X <- expand.grid(c(as.list(X), stats::setNames(list(level.treatment),treatment)))
    }

    ## ** normalize treatment and strata variable to numeric
    ## treatment variable (convert to numeric)
    if(!is.numeric(X[[treatment]]) || any(X[[treatment]] %in% 0:1 == FALSE)){
        X[[treatment]] <- as.numeric(factor(X[[treatment]], levels = level.treatment))-1
    }

    ##  strata variable (convert to factor with the right order)
    strata.var <- setdiff(names(X),treatment)
    if(length(strata.var)==0){ ## no existing strata variable
        X <- cbind(X, "..strata.." = 1) ## set all observation to strata 1
        strata.var <- "..strata.."
        attr(strata.var,"original") <- NULL
        if(is.null(level.strata)){level.strata <- "REF"}
    }else if(identical(strata.var,"..strata..")){ ## unique strata variable already in the right format
        attr(strata.var,"original") <- "..strata.."
        if(is.null(level.strata)){level.strata <- unique(X[["..strata.."]])}
    }else{

        ## check name of the strata variables
        ## if(any(strata.var == "..strata..")){
        ##     stop("Incorrect strata variable \n",
        ##          "Cannot use \"..strata..\" as it will be used internally \n.")
        ## }
        attr(strata.var,"original") <- strata.var
        
        ## update the design matrix with the right ordering of the factors
        for(iVar in setdiff(names(xlevels),treatment)){
            if(!is.null(xlevels)){
                X[[iVar]] <- factor(X[[iVar]], levels = xlevels[[iVar]])
            }else{
                X[[iVar]] <- factor(X[[iVar]])
            }
        }

        ## create the unique strata variable
        UX.strata <- interaction(X[,strata.var,drop=FALSE])
        if(is.null(level.strata)){level.strata <- levels(UX.strata)}
        
        X <- cbind(X, "..strata.." = as.numeric(factor(UX.strata, levels = level.strata)))
    }

    n.strata <- length(level.strata)

    ## ** collect elements
    out <- list(cif = setNames(lapply(1:n.strata, function(iS){setNames(lapply(1:2, function(iT){vector(mode = "list", length = n.CR)}), level.treatment)}), level.strata), ## cif at each obsevation time
                iid.cif = setNames(lapply(1:n.strata, function(iS){setNames(lapply(1:2, function(iT){vector(mode = "list", length = n.CR)}), level.treatment)}), level.strata), ## iid of the cif over time
                n.CR = n.CR, ## number of competing risks
                X = X, ## design matrix
                treatment = treatment, ## name of the treatment variable
                level.treatment = level.treatment,## levels of the treatment variable
                strata.var = strata.var, ## name of the original strata values (outside the treatment)
                level.strata = level.strata, ## vector contain all possible strata values (outside the treatment)
                n.strata = n.strata, ## number of strata (outside the treatment)
                last.estimate = array(NA, dim = c(n.strata, 2, n.CR), dimnames = list(level.strata,level.treatment,NULL)), ## estimate of each cumulative incidence at the last observed event (array strata x treatment x cause)
                last.time = matrix(NA, nrow = n.strata, ncol = 2, dimnames = list(level.strata,level.treatment)), ## time of the last observed event (matrix strata x treatment)
                jumpSurvHaz = setNames(lapply(1:n.strata, function(iS){setNames(vector(mode = "list", length = 2), level.treatment)}), level.strata) ## index of the jump times/suvival/cause-specific hazard
                ## for non-censored events (all causes) in the strata (list strata x treatment)
                )

    ## ** export
    return(out)

}

## * predict.BuyseTTEM
#' @title Prediction with Time to Event Model
#' @name predict.BuyseTTEM
#' 
#' @description Evaluate the cumulative incidence function (cif) / survival in one of the treatment groups. 
#' 
#' @param object time to event model. 
#' @param time [numeric vector] time at which to evaluate the cif/survival.
#' @param treatment [character/integer] Treatment or index of the treatment group.
#' @param strata  [character/integer] Strata or index of the strata.
#' @param cause [integer] The cause relative to which the cif will be evaluated.
#' @param se [logical] Should the standard errors associated with the cif/survival be output? 
#' @param iid [logical] Should the influence function associated with the cif/survival be output? 
#' @param type [character] by default (\code{NULL}), output the cif (competing risks) or survival (no competing risks) at the timepoints specified by argument times.
#' Otherwise, when argument times is missing, output the discretisation times for the survival (\code{"jump"}),
#' or the survival at the last observation time (\code{"last"}),
#' or the last observation time (\code{"last.time"}),
#' or the change in cif or survival between discretisation times (\code{"change"}).
#' @param ... not used, for compatibility with the generic method.
#'
#' @return a list containing the survival (element \code{survival}) or the cumulative incidence function (element \code{cif}),
#' and possible standard errors (element \code{.se}) and influence function (element \code{.iid}).
#' 
#' @keywords methods
#' 
#' @export
predict.BuyseTTEM <- function(object, time, treatment, strata, cause = 1, type = NULL,
                              se = NULL, iid = FALSE, ...){

    ## ** normalize input
    ## *** treatment
    if(!is.numeric(treatment)){
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

    ## *** type
    if(!is.null(type)){
        type <- match.arg(type, c("jump","last","change","last.time"))
        if(!missing(time)){
            stop("Argument \'type\' cannot be \"",type,"\" unless argument \'time\' is missing. \n")
        }
        if(iid){
            stop("Argument \'iid\' cannot be TRUE unless argument \'type\' is NULL. \n")
        }
        if(!is.null(se) && se){
            stop("Argument \'se\' cannot be TRUE unless argument \'type\' is NULL. \n")
        }
        se <- FALSE
    }else if(is.null(se)){
        se <- any(lengths(object$peron$iid.cif[[1]][[1]])>0)
    }
    
    ## ** output last cif estimate
    type.CR <- ifelse(object$peron$n.CR==1,"survival","competing.risks")
    if(!is.null(type) && type == "last.time"){        
        return(object$peron$last.time[strata,treatment])
    }else if(!is.null(type) && type == "last"){
        if(type.CR=="survival"){
            return(1-object$peron$last.estimate[strata,treatment,cause])
        }else{
            return(object$peron$last.estimate[strata,treatment,cause])
        }
    }else if(!is.null(type) && type == "jump"){
        return(object$peron$jumpSurvHaz[[strata]][[treatment]]$time.jump)
    }

    ## ** evaluate the cif at the requested times
    if(is.null(type)){
        table.peron <- object$peron$cif[[strata]][[treatment]][[cause]]
        index.table <- pmax(1,prodlim::sindex(jump.time = table.peron$time, eval.time = time)) ## pmin since 1 is taking care of negative times
        if(iid || se){
            tableIID.peron <- lava::iid(object, strata = strata, treatment = treatment, cause = cause)
        }
        
    }else if(type == "change"){
        if(type.CR=="survival"){ ## cif  = 1-survival
            table.peron <- object$peron$jumpSurvHaz[[strata]][[treatment]]
            table.peron$indexBefore.jump <- c(0,table.peron$index.jump[-NROW(table.peron)])
            table.peron$cif <- (1-table.peron$survival[table.peron$index.jump]) - (1-c(1,table.peron$survival)[table.peron$indexBefore.jump+1])
        }else{ ## competing risk case: can only be prodlim object
            table.peron <- object$peron$cif[[strata]][[treatment]][[cause]][!is.na(object$peron$cif[[strata]][[treatment]][[cause]]$obs),c("obs","index","time","cif"),drop=FALSE]
            names(table.peron)[1:3] <- c("index.obs", "index.jump", "time.jump")
            table.peron$indexBefore.jump <- c(0,table.peron$index.jump[-NROW(table.peron)])
            table.peron$cif <- table.peron$cif[table.peron$index.jump] - c(0,table.peron$cif)[table.peron$indexBefore.jump+1]
        }
        index.table <- 1:NROW(table.peron)
    }

    ## ** export
    if(is.null(type)){
        out <- list(time = table.peron[index.table,"time"],
                    index = table.peron[index.table,"index"])
    }else if(type == "change"){
        out <- list(time = table.peron[index.table,"time.jump"],
                    indexBefore = table.peron[index.table,"indexBefore.jump"],
                    indexAfter = table.peron[index.table,"index.jump"])
    }
    if(type.CR=="survival"){
        if(is.null(type)){
            out$survival <- 1-table.peron[index.table,"cif"]
        }else if(type == "change"){
            out$survival <- -table.peron[index.table,"cif"]
        }
        if(se){ ## only when type is NULL
            out$survival.se <- sqrt(colSums(tableIID.peron[,out$index+1,drop=FALSE]^2))
        }
        if(iid){ ## only when type is NULL
            out$survival.iid <- tableIID.peron[,out$index+1,drop=FALSE]
        }
        
    }else{
        out$cif <- table.peron[index.table,"cif"]
        if(se){## only when type is NULL
            out$cif.se <- sqrt(colSums(tableIID.peron[,out$index+1,drop=FALSE]^2))
        }
        if(iid){ ## only when type is NULL
            out$cif.iid <- tableIID.peron[,out$index+1,drop=FALSE]
        }        
    }
    return(out)
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

######################################################################
### BuyseTTEM.R ends here
