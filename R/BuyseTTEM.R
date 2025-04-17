### BuyseTTEM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2020 (12:15) 
## Version: 
## Last-Updated: apr  8 2025 (11:07) 
##           By: Brice Ozenne
##     Update #: 836
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
#' library(data.table)
#'
#' tau <- seq(0,3,length.out=10)
#'
#' #### survival case ####
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' e.prodlim <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.data)
#' ## plot(e.prodlim)
#' 
#' e.prodlim2 <- BuyseTTEM(e.prodlim, treatment = "treatment", iid = TRUE)
#' 
#' predict(e.prodlim2, time = tau, treatment = "T", strata = "a")
#' predict(e.prodlim, times = tau, newdata = data.frame(treatment = "T", strata = "a"))
#' 
#' predict(e.prodlim2, time = tau, treatment = "C", strata = "a")
#' predict(e.prodlim, times = tau, newdata = data.frame(treatment = "C", strata = "a"))
#' 
#' #### competing risk case ####
#' df.dataCR <- copy(df.data)
#' df.dataCR$status <- rbinom(NROW(df.dataCR), prob = 0.5, size = 2)
#' 
#' e.prodlimCR <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.dataCR)
#' ## plot(e.prodlimCR)
#' 
#' e.prodlimCR2 <- BuyseTTEM(e.prodlimCR, treatment = "treatment", iid = TRUE)
#' 
#' predict(e.prodlimCR2, time = tau, treatment = "T", strata = "a")
#' predict(e.prodlimCR, times = tau, newdata = data.frame(treatment = "T", strata = "a"), cause = 1)
#' 
#' predict(e.prodlimCR2, time = tau, treatment = "C", strata = "a")
#' predict(e.prodlimCR, times = tau, newdata = data.frame(treatment = "C", strata = "a"), cause = 1)
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
    object.stratified <- NCOL(object$X)>1
    
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
            if(object.stratified){
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
            iIndex.jump <- intersect(index.allJump, index.start[iStrata,iTreat]:index.stop[iStrata,iTreat])
    
            object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- data.frame(index.jump = iIndex.jump,
                                                                        time.jump = object$time[iIndex.jump],
                                                                        survival = object$surv[iIndex.jump],
                                                                        do.call(cbind,setNames(lapply(object$cause.hazard, function(iVec){iVec[iIndex.jump]}), paste0("hazard",1:n.CR))))
            object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]][order(object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$time.jump),]
            iIndex.jump <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$index.jump
            iStrata.nJump <- length(iIndex.jump)
            for(iEvent in 1:n.CR){ ## iEvent <- 1
                
                ## CIF at each jump (if any) and add time 0
                if(iStrata.nJump>0){
                    iStrata.extcuminc <- c(0,object$cuminc[[iEvent]][iIndex.jump])
                    iStrata.exttime.jump <- c(-tol12,object$time[iIndex.jump])
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
                    iStrata.index.afterJump <- c(1,prodlim::sindex(jump.times = iStrata.exttime.jump, eval.times = object$time[iIndex.jump] + tol12))
                }else{
                    iStrata.index.afterJump <- 1
                }
                if(iMissingCIF){
                    iStrata.index.afterJump <- c(iStrata.index.afterJump, NA)
                }

                ## *** store
                object$peron$cif[[iStrata]][[iTreat]][[iEvent]] <- data.frame(time = iStrata.exttime.jump,
                                                                              cif = iStrata.extcuminc,
                                                                              index = iStrata.index.afterJump - 1) ## move to C++ indexing
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

        object$peron$iid.hazard <- setNames(vector(mode = "list", length=n.strata), level.strata)
        object$peron$iid.survival <- setNames(vector(mode = "list", length=n.strata), level.strata)
        object$peron$iid.cif <- setNames(vector(mode = "list", length=n.strata), level.strata)
        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            object$peron$iid.hazard[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)
            object$peron$iid.survival[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)
            object$peron$iid.cif[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)

            for(iTreat in 1:2){ ## iTreat <- 1

                object$peron$iid.hazard[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)
                object$peron$iid.cif[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)

                if(treatment %in% names(object$model.matrix)){
                    iIndStrata <- which(model.matrix[[treatment]]==level.treatment[iTreat])
                }else{
                    iIndStrata <- 1:n.obs
                }
                if(object.stratified){
                    iIndStrata <- intersect(iIndStrata, which(model.matrix[["..strata.."]]==iStrata))
                }
                iIndex.jump <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$index.jump
                iN.jump <- length(iIndex.jump)
                
                if(length(iIndex.jump)==0){ ## no event: iid = 0
                    for(iEvent in 1:n.CR){ ## iEvent <- 1
                        object$peron$iid.hazard[[iStrata]][[iTreat]][[iEvent]] <- matrix(0, nrow = n.obs, ncol = 1)
                        object$peron$iid.survival[[iStrata]][[iTreat]][[iEvent]] <- matrix(0, nrow = n.obs, ncol = 1)
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iEvent]] <- matrix(0, nrow = n.obs, ncol = 1)
                    }
                    next
                }

                ## *** influence function for each cause-specific hazard
                for(iEvent in 1:n.CR){ ## iEvent <- 1
                    iJump.time <- object$time[iIndex.jump]
                    iHazard <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]][[paste0("hazard",iEvent)]]
                    iStatus <- do.call(cbind,lapply(iJump.time, function(iTime){
                        (abs(vec.eventtime[iIndStrata]-iTime)<tol11)*(vec.status[iIndStrata]==iEvent)
                    }))
                    
                    iAtRisk <- do.call(cbind,lapply(iJump.time, function(iTime){
                        vec.eventtime[iIndStrata] >= iTime
                    }))

                    object$peron$iid.hazard[[iStrata]][[iTreat]][[iEvent]] <- .rowScale_cpp(iStatus - .rowMultiply_cpp(iAtRisk, scale=iHazard), scale = colSums(iAtRisk))
                }

                ## *** influence function for the overall survival
                if(iid.surv == "exp"){
                    iSurv <- exp(-cumsum(rowSums(object$peron$jumpSurvHaz[[iStrata]][[iTreat]][,paste0("hazard",1:n.CR),drop=FALSE])))
                }else if(iid.surv == "prodlim"){
                    iSurv <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$survival    
                }
                object$peron$iid.survival[[iStrata]][[iTreat]] <- -.rowMultiply_cpp(.rowCumSum_cpp(Reduce("+",object$peron$iid.hazard[[iStrata]][[iTreat]])), iSurv)

                ## *** influence function for the cumulative incidence
                for(iCR in 1:n.CR){ ## iCR <- 1
                    iParam <- na.omit(unique(object$peron$cif[[iStrata]][[iTreat]][[iCR]]$index))

                    object$peron$iid.cif[[iStrata]][[iTreat]][[iCR]] <- matrix(0, nrow = n.obs, ncol = length(iParam))
                    if(test.CR){
                        iHazard <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]][[paste0("hazard",iCR)]] ## at t
                        iHazard.iid <- object$peron$iid.hazard[[iStrata]][[iTreat]][[iCR]] ## at t

                        iSurvival <- c(1,object$peron$jumpSurvHaz[[iStrata]][[iTreat]][1:(iN.jump-1),"survival"]) ## at t-
                        iSurvival.iid <- cbind(0,object$peron$iid.survival[[iStrata]][[iTreat]][,1:(iN.jump-1),drop=FALSE]) ## at t-

                        ## add iid at time 0 and (if censoring) NA after the last event
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iCR]][iIndStrata,] <- cbind(0,.rowCumSum_cpp(.rowMultiply_cpp(iSurvival.iid,iHazard) + .rowMultiply_cpp(iHazard.iid,iSurvival)))
                    }else{
                        ## add iid at time 0 and (if censoring) NA after the last event
                        object$peron$iid.cif[[iStrata]][[iTreat]][[iCR]][iIndStrata,] <- cbind(0,-object$peron$iid.survival[[iStrata]][[iTreat]])
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
BuyseTTEM.survreg <- function(object, treatment, n.grid = 1e3, iid, ...){

    tol12 <- 1e-12
    tol11 <- 1e-11
    dots <- list(...)

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
    object.stratified <- NCOL(object$X)>1

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

            if(treatment %in% names(mf)){
                ## restrict to cases? which(mf[,1][,2]==1)
                iIndex.obs <- which(mf[[treatment]]==iTreat)
            }else{
                iIndex.obs <- 1:n.obs
            }
            if(object.stratified){
                iIndex.obs <- intersect(iIndex.obs,which(object$peron$X[,"..strata.."]==iStrata))
            }

            ## jump time  in this strata
            iNewdata <- mf[iIndex.obs[1],,drop=FALSE]
            iJump <- predict(object, newdata = iNewdata, p = grid.quantile, type = "quantile")
            
            object$peron$jumpSurvHaz[[iStrata]][[iTreat]] <- data.frame(index.jump = NA,
                                                                        time.jump = iJump,
                                                                        survival = 1-grid.quantile)

            iTime <- sort(mf[iIndex.obs,1][,1])
            iIndex.param <- prodlim::sindex(jump.times = iJump, eval.times = iTime + tol12)
            iSurv <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$survival[iIndex.param]
            object$peron$cif[[iStrata]][[iTreat]][[1]] <- data.frame(time = c(-tol12,iTime),
                                                                     cif = c(0,1-iSurv),
                                                                     index = c(0,iIndex.param-1))

            object$peron$last.time[iStrata,iTreat] <- utils::tail(iJump,1)

            if(iid){
                iP <- NROW(object$peron$jumpSurvHaz[[iStrata]][[iTreat]])
                iLP <- drop(X[iIndex.obs[1],,drop=FALSE] %*% beta)
                object.iid.iLP <- object.iid.beta %*% t(X[iIndex.obs[1],,drop=FALSE])
                
                ## *** compute time with se [not used]
                ## fit.trans <- iLP + object.dist$quantile(grid.quantile) * sigma
                ## range(object.dist$itrans(fit.trans) - iJump)
                ## fit.iid <- .rowScale_cpp(.colCenter_cpp(object.iid.sigma %*% object.dist$quantile(grid.quantile), center = - object.iid.iLP), object.dist$dtrans(object.dist$itrans(fit.trans)))
                ## fit.se <- sqrt(colSums(fit.iid^2))
                ## quantile(predict(object, type = "quantile", p = grid.quantile, newdata = mf[iIndex.obs[1],,drop=FALSE],se = TRUE)$se - fit.se, na.rm = TRUE)

                ## *** compute survival with se
                iPred <- (object.dist$trans(iJump) - iLP)/sigma ## range(grid.quantile - object.dist$quantileM1(iPred))
                
                object$peron$iid.cif[[iStrata]][[iTreat]][[1]] <- .rowMultiply_cpp(.colCenter_cpp(- object.iid.sigma %*% (object.dist$trans(iJump) - iLP) / sigma^2,
                                                                                                  center = object.iid.iLP / sigma), scale = object.dist$dquantileM1(iPred))
                object$peron$iid.cif[[iStrata]][[iTreat]][[1]][,iJump==0] <- 0

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
        if(any(strata.var == "..strata..")){
            stop("Incorrect strata variable \n",
                 "Cannot use \"..strata..\" as it will be used internally \n.")
        }
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
#' @param iid [logical] Should the influence function associated with the cif/survival be output? 
#' @param ... not used, for compatibility with the generic method.
#'
#' @return a list containing the survival (element \code{survival}) or the cumulative incidence function (element \code{cif}),
#' and possible standard errors (element \code{.se}) and influence function (element \code{.iid}).
#' 
#' @keywords methods
#' 
#' @export
predict.BuyseTTEM <- function(object, time, treatment, strata, cause = 1, iid = FALSE, ...){

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
    
    type <- ifelse(object$peron$n.CR==1,"survival","competing.risks")

    ## ** output last cif estimate    
    if(identical(time,"last")){        
        if(type=="survival"){
            return(1-object$peron$last.estimate[strata,treatment,cause])
        }else{
            return(object$peron$last.estimate[strata,treatment,cause])
        }
    }else if(identical(time,"jump")){
        return(object$peron$jumpSurvHaz[[strata]][[treatment]]$time.jump)
    }

    ## ** extract information
    out <- list()

    table.cif <- object$peron$cif[[strata]][[treatment]][[cause]]
    index.table <- pmax(1,prodlim::sindex(jump.time = table.cif$time, eval.time = time)) ## pmin since 1 is taking care of negative times
    out$index <- table.cif[index.table,"index"]
    
    if(type=="survival"){
        out$survival <- 1-table.cif[index.table,"cif"]
        if(iid){
            out$survival.iid <- lava::iid(object, strata = strata, treatment = treatment, cause = cause)[,out$index+1,drop=FALSE]
            out$survival.se <- sqrt(colSums(out$survival.iid^2))
        }
    }else{
        out$cif <- table.cif[index.table,"cif"]
        if(iid){
            out$cif.iid <- lava::iid(object, strata = strata, treatment = treatment, cause = cause)[,out$index+1,drop=FALSE]
            out$cif.se <- sqrt(colSums(out$cif.iid^2))
        }
    }
    
    ## ** export
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
    type <- ifelse(object$peron$n.CR==1,"survival","competing.risks")

    ## ** export
    if(type=="survival"){
        return(-out)
    }else{
        return(out)
    }
}

######################################################################
### BuyseTTEM.R ends here
