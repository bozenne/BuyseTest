### BuyseTTEM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2020 (12:15) 
## Version: 
## Last-Updated: jan  5 2021 (12:16) 
##           By: Brice Ozenne
##     Update #: 310
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
#' @param treatment [character] Name of the treatment variable.
#' @param level.treatment [character/numeric vector] Unique values of the treatment variable.
#' @param iid [logical] Should the iid decomposition of the predictions be output. 
#' @param iid.surv [character] Estimator of the survival used when computing the influence function.
#' Can be the product limit estimator (\code{"prodlim"}) or an exponential approximation (\code{"exp"}, same as in \code{riskRegression::predictCoxPL}).
#' @param ... additional arguments passed to lower lever methods.
#'
#' @examples
#' library(prodlim)
#' library(data.table)
#'
#' #### survival case ####
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' e.prodlim <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.data)
#' ## plot(e.prodlim)
#' 
#' e.prodlim2 <- BuyseTTEM(e.prodlim,
#'           treatment = "treatment",
#'           level.treatment = unique(df.data$treatment),
#'           iid = TRUE)
#' predict(e.prodlim2, time = 1:10, treatment = "T", strata = "a")
#' predict(e.prodlim, times = 1:10, newdata = data.frame(treatment = "T", strata = "a"))
#' 
#' predict(e.prodlim2, time = 1:10, treatment = "C", strata = "a")
#' predict(e.prodlim, times = 1:10, newdata = data.frame(treatment = "C", strata = "a"))
#' 
#' #### competing risk case ####
#' df.dataCR <- copy(df.data)
#' df.dataCR$status <- rbinom(NROW(df.dataCR), prob = 0.5, size = 2)
#' 
#' e.prodlimCR <- prodlim(Hist(eventtime,status)~treatment+strata, data = df.dataCR)
#' ## plot(e.prodlimCR)
#' 
#' e.prodlimCR <- BuyseTTEM(e.prodlimCR,
#'           treatment = "treatment",
#'           level.treatment = unique(df.dataCR$treatment),
#'           iid = TRUE)
#' @export
`BuyseTTEM` <-
    function(object,...) UseMethod("BuyseTTEM")


## * BuyseTTEM.default
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.formula <- function(object, treatment, level.treatment = NULL, iid, iid.surv = "exp", ...){
    e.prodlim <- prodlim(object, ...)
    return(BuyseTTEM(e.prodlim, treatment = treatment, level.treatment = level.treatment, iid = iid, iid.surv = iid.surv))
}

## * BuyseTTEM.prodlim
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.prodlim <- function(object, treatment, level.treatment = NULL, iid, iid.surv = "exp", ...){

    X <- object$X
    p <- NCOL(X)
    tol12 <- 1e-12
    tol11 <- 1e-11
    
    ## ** check arguments
    iid.surv <- match.arg(iid.surv, c("prodlim","exp"))
    if(missing(treatment) && NCOL(object$X)==1){
        treatment <- names(object$X)[1]
    }
    if(treatment %in%  names(X) == FALSE){
        stop("Wrong specification of the argument \'treatment\' \n",
             "Could not be found in the prodlim object, e.g. in object$X. \n")
    }
    if(is.null(level.treatment)){
        level.treatment <- unique(X[[treatment]])
    }else if(length(level.treatment) != length(unique(X[[treatment]]))){
        stop("Wrong specification of the argument \'level.treatment\' \n",
             "Does not match the number of possible values in object$X. \n")
    }

    name.strata <- setdiff(names(X),treatment)
    if(length(name.strata)==0){
        X <- cbind(X, "..strata.." = 1)
        X.strata <- NULL
        n.strata <- 1
        name.strata <- "..strata.."
        Uname.strata <- "REF"
    }else{
        X.strata <- unique(X[,name.strata,drop=FALSE])
        n.strata <- NROW(X.strata)
        name.strata <- names(X.strata)
        Uname.strata <- interaction(X.strata)
    }
    
    if(NCOL(X.strata)>1 && any(name.strata == "..strata..")){
        stop("Incorrect strata variable \n",
             "Cannot use \"..strata..\" as it will be used internally \n.")
    }
    if(any(object$time<=0)){
        stop("Only handles strictly positive event times \n")
    }
    
    ## ** normalize prodlim for competing risk / survival
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
    
    ## ** prepare output
    object$peron <- list(cif = setNames(vector(mode = "list", length = n.strata), Uname.strata),
                         X = X, ## strata
                         level.treatment = level.treatment,
                         Ustrata = levels(Uname.strata),
                         n.CR = n.CR,
                         index.start = NULL, ## index (among all estimates) of the first estimate relative to a given strata (list)
                         index.stop = NULL, ## index (among all estimates) of the last estimate relative to a given strata (list)
                         last.estimate = NULL, ## estimate of the cumulative incidence at the last observed event (matrix strata x treatment)
                         last.time = NULL, ## time of the last observed event (matrix strata x treatment)
                         jumpSurvHaz = setNames(vector(mode = "list", length = n.strata), Uname.strata) ## index of the jump times/suvival/cause-specific hazard
                         ## for non-censored events (all causes) in the strata (list strata x treatment)
                         )
    
    ## ** update design matrix
    ## normalize treatment variable
    if(!is.numeric(X[[treatment]]) || any(X[[treatment]] %in% 0:1 == FALSE)){
        object$peron$X[[treatment]] <- as.numeric(factor(object$X[[treatment]], levels = level.treatment))-1
    }

    ##  normalize strata variable
    if(!identical(name.strata,"..strata..")){
        value.strata <- interaction(object$X[,name.strata,drop=FALSE])
        object$peron$X <- cbind(object$peron$X,
                                "..strata.." = as.numeric(as.factor(value.strata)))
    }

    ## ** compute start/stop indexes per strata
    iIndexX.C <- which(object$peron$X[[treatment]]==0)
    iIndexX.T <- which(object$peron$X[[treatment]]==1)

    ## position in the results of the pair (treatment,strata)
    iMindexX <- do.call(rbind,lapply(1:n.strata, function(iStrata){
        iIndexS <- which(object$peron$X[["..strata.."]]==iStrata)
        return(c(C = intersect(iIndexX.C,iIndexS),
                 T = intersect(iIndexX.T,iIndexS)))
    }))

    object$peron$index.start <- matrix(NA, nrow = n.strata, ncol = 2, dimnames = list(NULL,level.treatment))
    object$peron$index.start[] <- object$first.strata[iMindexX]

    object$peron$index.stop <- matrix(NA, nrow = n.strata, ncol = 2, dimnames = list(NULL,level.treatment))
    object$peron$index.stop[] <- object$first.strata[iMindexX] + object$size.strata[iMindexX] - 1

    ## ** find last CIF value
    object$peron$last.time <- cbind(object$time[object$peron$index.stop[,level.treatment[1]]],
                                    object$time[object$peron$index.stop[,level.treatment[2]]])
    colnames(object$peron$last.time) <- level.treatment

    object$peron$last.estimate <- array(unlist(lapply(object$cuminc, function(iVec){
        cbind(iVec[object$peron$index.stop[,level.treatment[1]]],
              iVec[object$peron$index.stop[,level.treatment[2]]])
    })), dim = c(n.strata,2,n.CR),
    dimnames = list(levels(Uname.strata),level.treatment,NULL))

    ## ** table CIF and dCIF
    index.allJump <- sort(unlist(lapply(object$cause.hazard, function(iVec){which(iVec>0)})))

    for(iStrata in 1:n.strata){ ## iStrata <- 1
        object$peron$value[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)
        object$peron$jumpSurvHaz[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)

        for(iTreat in level.treatment){ ## iTreat <- 1

            object$peron$value[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)

            iMissingCIF <- sum(object$peron$last.estimate[iStrata,iTreat,])<(1-tol12)
            
            ## jump time (for all causes) in this strata
            iIndex.jump <- intersect(index.allJump, object$peron$index.start[iStrata,iTreat]:object$peron$index.stop[iStrata,iTreat])
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

                ## dCIF at each jump
                if(iStrata.nJump>0){
                    ## find index of the CIF parameter before and after the jump (NA when after last observation point)
                    ## also add time 0
                    iStrata.index.beforeJump <- c(1,prodlim::sindex(jump.times = iStrata.exttime.jump, eval.times = object$time[iIndex.jump] - tol12))
                    iStrata.index.afterJump <- c(1,prodlim::sindex(jump.times = iStrata.exttime.jump, eval.times = object$time[iIndex.jump] + tol12))
                    iStrata.dextcuminc <- iStrata.extcuminc[iStrata.index.afterJump] - iStrata.extcuminc[iStrata.index.beforeJump]
                }else{
                    iStrata.index.beforeJump <- 1
                    iStrata.index.afterJump <- 1
                    iStrata.dextcuminc <- 0
                }

                ## complete dCIF after the last event
                if(iMissingCIF){
                    iStrata.index.beforeJump <- c(iStrata.index.beforeJump, iStrata.nJump+1)
                    iStrata.index.afterJump <- c(iStrata.index.afterJump, NA)
                    iStrata.dextcuminc <- c(iStrata.dextcuminc,NA)
                }

                ## *** store
                object$peron$cif[[iStrata]][[iTreat]][[iEvent]] <- data.frame(time = iStrata.exttime.jump,
                                                                              cif = iStrata.extcuminc,
                                                                              dcif = iStrata.dextcuminc,
                                                                              index.cif.before = iStrata.index.beforeJump - 1, ## move to C++ indexing
                                                                              index.cif.after = iStrata.index.afterJump - 1) ## move to C++ indexing
            }
        }
    }

    ## cif <- object$peron$cif[[1]][[1]][[1]]
    ## (cif[cif[,"index.cif.after"]+1,"cif"] - cif[cif[,"index.cif.before"]+1,"cif"]) - cif[,"dcif"]


    ## ** table iid
    if(iid){
        vec.eventtime <- object$model.response[object$originalDataOrder,"time"]
        if(test.CR){
            vec.status <- object$model.response[object$originalDataOrder,"event"] ## 1:n.CR event, n.CR+1 censoring
        }else{
            vec.status <- object$model.response[object$originalDataOrder,"status"] ## 0 censoring, 1 event
        }
        name.allStrata <- interaction(object$peron$X[,c(treatment,name.strata)])

        if(is.null(X.strata)){
            name.allStrataOriginal <- levels(object$X$treatment)
        }else{
            name.allStrataOriginal <- levels(interaction(object$X[,c(treatment,name.strata)]))
        }
        vec.allStrata <- as.numeric(factor(interaction(object$model.matrix[object$originalDataOrder,,drop=FALSE]), levels = name.allStrataOriginal))
        n.obs <- length(vec.allStrata)
        
        object$peron$iid.hazard <- setNames(vector(mode = "list", length=n.strata), Uname.strata)
        object$peron$iid.survival <- setNames(vector(mode = "list", length=n.strata), Uname.strata)
        object$peron$iid.cif <- setNames(vector(mode = "list", length=n.strata), Uname.strata)
        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            object$peron$iid.hazard[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)
            object$peron$iid.survival[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)
            object$peron$iid.cif[[iStrata]] <- setNames(vector(mode = "list", length=2), level.treatment)

            for(iTreat in level.treatment){ ## iTreat <- 1

                object$peron$iid.hazard[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)
                object$peron$iid.cif[[iStrata]][[iTreat]] <- vector(mode = "list", length=n.CR)

                iIndex.allStrata <- intersect(which(level.treatment[object$peron$X[,treatment]+1]==iTreat),which(object$peron$X[,"..strata.."]==iStrata))
                iIndStrata <- which(vec.allStrata==iIndex.allStrata)
                iIndex.jump <- object$peron$jumpSurvHaz[[iStrata]][[iTreat]]$index.jump

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
                        vec.eventtime[iIndStrata]>=iTime
                    }))
                    iN.jump <- length(iIndex.jump)
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
                for(iCR in 1:n.CR){
                    iParam <- na.omit(unique(c(object$peron$cif[[iStrata]][[iTreat]][[iCR]]$index.cif.before,object$peron$cif[[iStrata]][[iTreat]][[iCR]]$index.cif.after)))

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

## * BuyseTTEM.BuyseTTEM
#' @rdname BuyseTTEM
#' @export
BuyseTTEM.BuyseTTEM <- function(object, ...){
    return(object)
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
#' @export
predict.BuyseTTEM <- function(object, time, treatment, strata, cause = 1, iid = FALSE, ...){

    ## ** normalize input
    if(!is.numeric(treatment)){
        treatment <- which(treatment == object$peron$level.treatment)
    }
    if(missing(strata)){
        strata <- 1
    }
    
    ## ** output last cif estimate
    if(identical(time,"last")){
        if(object$type=="surv"){
            return(1-object$peron$last.estimate[strata,treatment,cause])
        }else{
            return(object$peron$last.estimate[strata,treatment,cause])
        }
    }else if(identical(time,"jump")){
        return(object$peron$jumpSurvHaz[[strata]][[treatment]]$time.jump)
    }

    ## ** extract information
    out <- list(index = NULL, cif = NULL)
    if(inherits(object,"prodlim")){
        if(object$type=="surv"){
            table.cif <- object$peron$cif[[strata]][[treatment]][[1]]
        }else{
            table.cif <- object$peron$cif[[strata]][[treatment]][[cause]]
        }
        index.table <- pmax(1,prodlim::sindex(jump.time = table.cif$time, eval.time = time)) ## pmin since 1 is taking care of negative times
        out$index <- table.cif[index.table,"index.cif.after"]

        if(object$type=="surv"){
            out$survival <- 1-table.cif[index.table,"cif"]
        }else{
            out$cif <- table.cif[index.table,"cif"]
        }
    
        if(iid){
            if(object$type=="surv"){
                out$survival.iid <- iid(object, strata = strata, treatment = treatment, cause = cause)[,out$index+1,drop=FALSE]
                out$survival.se <- sqrt(colSums(out$survival.iid^2))
            }else{
                out$cif.iid <- iid(object, strata = strata, treatment = treatment, cause = cause)[,out$index+1,drop=FALSE]
                out$cif.se <- sqrt(colSums(out$cif.iid^2))
            }
        }
    }

    ## ** export
    return(out)

}

## * iid.BuyseTTEM
#' @export
iid.BuyseTTEM <- function(x, treatment, strata, cause = 1){
    if(is.null(x$peron$iid.cif[[strata]][[treatment]][[cause]])){
        stop("iid decomposition not available - consider setting the argument \'iid\' to TRUE when calling BuyseTTEM. \n")
    }
    if(x$type=="surv"){
        return(-x$peron$iid.cif[[strata]][[treatment]][[cause]])
    }else{
        return(x$peron$iid.cif[[strata]][[treatment]][[cause]])
    }
}

######################################################################
### BuyseTTEM.R ends here
