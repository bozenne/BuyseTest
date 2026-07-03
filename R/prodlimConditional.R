### linyinglim.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 30 2026 (16:04) 
## Version: 
## Last-Updated: jul  3 2026 (15:01) 
##           By: Brice Ozenne
##     Update #: 176
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - prodlimConditional
##' @examples
##' seqT <- seq(0,0.75,length.out = 5)
##' 
##' ########################################
##' ##  Two time to events (no censoring) ##
##' ########################################
##' set.seed(1)
##' argsTTE <- list(scale.T = c(1,1), scale.C = c(0.5,2), scale.censoring.T = FALSE)
##' d2TTE <- simBuyseTest(100, argsTTE = argsTTE)
##'
##' eLYlim <- prodlimConditional(eventtime = c("eventtime1","eventtime2"),
##'                              status = c("status1","status2"), strata = "treatment",
##'                              data = d2TTE, restriction = rep(0.75,2), times = seqT, method = "IPCW")
##' 
##' 
##' ##########################################
##' ##  Two time to events (with censoring) ##
##' ##########################################

## * Function mprodlim
prodlimConditional <- function(eventtime, status, strata, data, restriction, times,
                               method = "minevent", iid = TRUE){


    ## ** normalize user input
    n.endpoint <- length(eventtime)
    data <- as.data.frame(data)
    
    ## *** eventtime
    if(n.endpoint==1){
        stop("Argument \'eventtime\' should have length at least two. \n")
    }
    if(any(duplicated(eventtime))){
        stop("Argument \'eventtime\' should not contain duplicated elements. \n")
    }
    if(any(is.na(eventtime))){
        stop("Argument \'eventtime\' should not contain NAs. \n")
    }
    if(any(eventtime %in% names(data) == FALSE)){
        stop("Mismatch between argument \'eventtime\' and argument \'data\'. \n",
             "Cannot find element(s) \"",paste(eventtime[eventtime %in% names(data) == FALSE], collapse = "\", \""),"\" from \'eventtime\' in \'data\'. \n")
    }
    if(any(data[eventtime]<0)){
        stop("Argument \'times\' should correspond to columns with non-negative values. \n")
    }

    ## *** status
    if(length(status)!=n.endpoint){
        stop("Mismatch between argument \'eventtime\' and argument \'status\'. \n",
             "They should have the same length. \n")
    }
    if(any(duplicated(status))){
        stop("Argument \'status\' should not contain duplicated elements. \n")
    }
    if(any(is.na(status))){
        stop("Argument \'status\' should not contain NAs. \n")
    }
    if(any(status %in% names(data) == FALSE)){
        stop("Mismatch between argument \'status\' and argument \'data\'. \n",
             "Cannot find element(s) \"",paste(status[status %in% names(data) == FALSE], collapse = "\", \""),"\" from \'status\' in \'data\'. \n")
    }
    if(any(unlist(data[status]) %in% 0:1 == FALSE)){
        stop("Argument \'status\' should only take value 0 or 1. \n")
    }

    ## *** strata
    if(any(duplicated(strata))){
        stop("Argument \'strata\' should not contain duplicated elements. \n")
    }
    if(any(is.na(strata))){
        stop("Argument \'strata\' should not contain NAs. \n")
    }
    if(any(strata %in% names(data) == FALSE)){
        stop("Mismatch between argument \'strata\' and argument \'data\'. \n",
             "Cannot find element(s) \"",paste(strata[strata %in% names(data) == FALSE], collapse = "\", \""),"\" from \'strata\' in \'data\'. \n")
    }

    ## *** restriction
    if(length(restriction)==n.endpoint){
        restriction <- restriction[1:(n.endpoint-1)]
    }else if(length(restriction)!=(n.endpoint-1)){
        stop("Mismatch between argument \'eventtime\' and argument \'restriction\'. \n",
             "The length of argument \'restriction\' should the same or exactly one less than the length of argument \'eventtime\'. \n")
    }
    restriction <- c(restriction, times[length(times)])

    if(any(is.na(restriction))){
        stop("Argument \'restriction\' should not contain NAs. \n")
    }
    if(!is.numeric(restriction)){
        stop("Argument \'restriction\' should be a numeric vector. \n")
    }
    if(any(restriction<=0)){
        stop("Argument \'restriction\' should contain strictly positive values. \n")
    }

    ## *** times
    Utimes <- sort(unique(times))
    n.times <- length(times)
    if(!is.numeric(times)){
        stop("Argument \'times\' should be a numeric vector. \n")
    }
    if(any(times<0)){
        stop("Argument \'times\' should contain non-negative values. \n")
    }

    ## *** method
    method <- match.arg(method, c("IPCW","minevent"))

    ## ** create a normalized dataset
    data2 <- cbind(data[eventtime],data[status],data[strata])
    names(data2) <- c(paste0("eventtime",1:n.endpoint),
                      paste0("status",1:n.endpoint),
                      strata)
    for(iStrata in strata){
        if(!is.factor(data2[[iStrata]])){
            data2[[iStrata]] <- as.factor(data2[[iStrata]])
        }
    }
    
    data2[paste0("eventtimeN",1:n.endpoint)] <- NA
    data2[paste0("eventtimeMax",1:n.endpoint)] <- NA
    data2[paste0("statusMax",1:n.endpoint)] <- NA

    ## ** marginal
    ## a) Prob[X_1 > \tau_1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2 | strata] = Prob[X_1/\tau_1 > 1, X_2/\tau_2 > 1 | strata] = Prob[min(X_1/\tau_1,X_2/\tau_2) > 1 | strata]

    ## normalize: X_1/\tau_1, X_2/tau_2, ...
    data2[paste0("eventtimeN",1:(n.endpoint-1))] <- sweep(data2[paste0("eventtime",1:(n.endpoint-1))], MARGIN = 1, FUN = "/", STATS = restriction)

    ## find the earlier time and its status
    if(n.endpoint==2){
        data2$eventtime.minM1 <- data2$eventtimeN1
        data2$status.minM1 <- data2$status1
    }else{
        rowindex.minM1 <- apply(data2[paste0("eventtimeN",1:(n.endpoint-1))], MARGIN = 1, which.min)
        posindex.minM1 <- 1:NROW(data) + (rowindex.minM1-1) * NROW(data)
        data2$eventtime.minM1 <- as.matrix(data2[paste0("eventtimeN",1:(n.endpoint-1))])[posindex.minM1]
        data2$status.minM1 <- as.matrix(data2[paste0("status",1:(n.endpoint-1))])[posindex.minM1]
    }
    
    ## fit Kaplan Meier
    ff.marginal <- stats::reformulate(termlabels = strata, response = "Hist(eventtime.minM1, status.minM1)")
    eKM.marginal <- prodlim::prodlim(ff.marginal, data = data2)

    ## extract survival probabilities
    eTTEM.marginal <- BuyseTTEM(eKM.marginal, treatment = strata[1], iid = iid)
    eKMsurv.marginal <- summary(eKM.marginal, times = 1, iid = iid)
    strata.names <- as.character(interaction(as.data.frame(eKMsurv.marginal)[,strata,drop=FALSE],sep="."))
    n.strata <- length(strata.names)

    ## ** joint
    eKMsurv.joint <- vector(length = n.times, mode = "list")
    
    ## *** special case: times = 0
    ## a) Prob[X_1 > \tau_1, X_2 > 0 | strata]                 = Prob[X_1/tau_1 > 1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > 0 | strata]   = Prob[min(X_1/tau_1, X_2/tau_2) > 1 | strata]
    if(any(Utimes==0)){
        eKMsurv.joint[which(Utimes==0)] <- replicate(sum(Utimes==0),list(eKMsurv.marginal))
    }else{
        eKMsurv.joint <- NULL
    }

    ## *** 'normal' case: times > 0
    ## a) Prob[X_1 > \tau_1, X_2 > t | strata]                 = Prob[min(X_1/tau_1, X_2 / t) > 1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > t | strata]   = Prob[min(X_1/tau_1, X_2/tau_2, X_3 / t) > 1 | strata]

    if(any(Utimes>0)){

        ## prepare
        if(method == "minevent"){
            ff.joint <- stats::reformulate(termlabels = strata, response = "Hist(eventtime.min, status.min)")
        }

        ## loop over times
        for(iT in which(Utimes>0)){ ## iT <- 6

            if(method == "minevent"){
                data2$eventtime.min <- pmin(data2$eventtime.minM1, data2[[paste0("eventtime",n.endpoint)]]/Utimes[iT])
                data2$status.min <- ifelse(data2$eventtime.minM1 <= data2[[paste0("eventtime",n.endpoint)]]/Utimes[iT],
                                           data2$status.minM1, data2[[paste0("status",n.endpoint)]])
                iKM.joint <- prodlim(ff.joint, data = data2)
                iTTEM.joint <- BuyseTTEM(iKM.joint, treatment = strata[1], iid = iid)
                eKMsurv.joint[[iT]] <- summary(iTTEM.joint, times = 1, iid = iid)

                eKMsurv.conditional <- cbind(eKMsurv.joint[,c(strata,"times")],
                                             surv = eKMsurv.joint$surv/eKMsurv.joint$marginal)
                
                if(iid){
                    term1.iid <- sweep(attr(eKMsurv.joint[[iT]],"iid"), margin = 1, FUN = "/", stats = eKMsurv.joint$marginal)
                    term2.iid <- - sweep(attr(eKMsurv.marginal[[iT]],"iid"), margin = 1, FUN = "/", stats = eKMsurv.joint$conditional/eKMsurv.joint$marginal^2)
                    attr(eKMsurv.conditional,"iid") <- term1.iid + term2.iid
                    eKMsurv.conditional$se.surv <- colSums(attr(eKMsurv.conditional,"iid")^2)
                    eKMsurv.conditional$lower <- eKMsurv.conditional$surv + stats::qnorm(0.025) * eKMsurv.conditional$se.surv
                    eKMsurv.conditional$upper <- eKMsurv.conditional$surv + stats::qnorm(0.975) * eKMsurv.conditional$se.surv
                }
            }else if(method == "IPCW"){
                ## ## number of events in the observed process
                ## Ncount <- do.call(rbind,lapply(times, function(iTime){
                ##     with(data2,tapply((eventtime1>tau)*(eventtime2>iTime),group,mean))
                ## }))

                ## ## denominator
                ## data2$eventtime12 <- pmax(data2$eventtime1, data2$eventtime2)
                ## data2$status12 <- 1-data2$status1*(data2$status2>0)
                ## eKM12 <- prodlim(Hist(eventtime12, status12) ~ group, data = data2)
                ## eIPCW <- do.call(rbind,lapply(times, function(iTime){
                ##     summary(eKM12, times = max(iTime,tau))$surv
                ## }))

                ## ## assemble
                ## eSurv <- cbind(Ncount/eIPCW, Ncount/(eIPCW*eCon))
            }
        }
    }

    ## ** assemble
    out <- list(conditional = eKMsurv.joint,
                joint = eKMsurv.joint,
                marginal = eKMsurv.marginal)
    out$conditional$n.risk <- NA
    out$conditional$n.event <- NA
    out$conditional$n.lost <- NA
    out$conditional$surv <- eKMsurv.joint$surv/eKMsurv.marginal$surv
    out$conditional$se.surv <- NA
    out$conditional$lower <- NA
    out$conditional$upper <- NA

    ## ** export    
    return(out)
}


##----------------------------------------------------------------------
### linyinglim.R ends here
