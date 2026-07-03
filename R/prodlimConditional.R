### prodlimConditional.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 30 2026 (16:04) 
## Version: 
## Last-Updated: jul  3 2026 (17:23) 
##           By: Brice Ozenne
##     Update #: 275
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * prodlimConditional (documentation)

## * prodlimConditional (examples)
##' @examples
##' seqT <- seq(0,0.75,length.out = 5)
##' 
##' ###################################################
##' ## Two time to events - Two censoring mechanisms ##
##' ###################################################
##' set.seed(1)
##' argsTTE <- list(scale.T = c(1,1), scale.C = c(0.5,2),
##'                 frailty.C = 0, frailty.T = 2)
##' d2TTE <- simBuyseTest(1000, argsTTE = argsTTE, latent = TRUE)
##' d2TTE$status1Uncensored <- 1
##' d2TTE$status2Uncensored <- 1
##'
##' ## Independent TTE in one group, positively correlated in the other
##' by(d2TTE,d2TTE$treatment, function(iData){
##' with(iData,cor(eventtime1Uncensored,eventtime2Uncensored))
##' })
##'
##' #### no censoring ####
##'
##' ## empirical average estimator
##' d2TTEsurv <- d2TTE[d2TTE$eventtime1Uncensored>0.75,]
##' do.call(cbind,lapply(seqT, function(iT){
##'  tapply(d2TTEsurv$eventtime2Uncensored>iT,d2TTEsurv$treatment,mean)
##' }))
##' ##   [,1] [,2]      [,3]      [,4]      [,5]
##' ## C    1 0.8883721 0.8093023 0.7255814 0.6883721
##' ## T    1 0.8993506 0.7824675 0.6948052 0.6363636
##'
##' ## minevent estimator
##' eProdlim2 <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                              status = paste0("status",1:2,"Uncensored"),
##'                              strata = "treatment", restriction = rep(0.75,2),
##'                              data = d2TTE, times = seqT)
##' eProdlim2
##' ##    treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.8883721 0.04622803 0.7977668 0.9789774
##' ## 4          T 0.1875 0.8993506 0.05350114 0.7944903 1.0042110
##' ## 5          C 0.3750 0.8093023 0.04413190 0.7228054 0.8957993
##' ## 6          T 0.3750 0.7824675 0.05089930 0.6827067 0.8822283
##' ## 7          C 0.5625 0.7255814 0.04198691 0.6432886 0.8078742
##' ## 8          T 0.5625 0.6948052 0.04849389 0.5997589 0.7898515
##' ## 9          C 0.7500 0.6883721 0.04086851 0.6082713 0.7684729
##' ## 10         T 0.7500 0.6363636 0.04685719 0.5445252 0.7282020
##' 
##' ## Peron (ignores conditionning)
##' e.peron <- prodlim(Hist(eventtime2Uncensored, status2Uncensored) ~ treatment,
##'                    data = d2TTE)
##' eS.peron <- summary(e.peron, times = seqT)
##' eS.peron[order(eS.peron$time),]
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000   1000       0      0 1.000 0.00000 1.000 1.000
##' ## 2          T 0.000   1000       0      0 1.000 0.00000 1.000 1.000
##' ## 3          C 0.188    912       0      0 0.912 0.00896 0.894 0.930
##' ## 4          T 0.188    569       0      0 0.569 0.01566 0.538 0.600
##' ## 5          C 0.375    836       0      0 0.836 0.01171 0.813 0.859
##' ## 6          T 0.375    432       0      0 0.432 0.01566 0.401 0.463
##' ## 7          C 0.562    755       0      0 0.755 0.01360 0.728 0.782
##' ## 8          T 0.562    359       0      0 0.359 0.01517 0.329 0.389
##' ## 9          C 0.750    688       0      0 0.688 0.01465 0.659 0.717
##' ## 10         T 0.750    311       0      0 0.311 0.01464 0.282 0.340
##' ## ok for C but too low for T
##' 
##' #### with censoring ####
##' 
##' ## empirical average estimator (inefficient)
##' d2TTEsurvC <- d2TTE[d2TTE$eventtime1>0.75,]
##' do.call(cbind,lapply(seqT, function(iT){
##'  tapply(d2TTEsurvC$eventtime2Uncensored>iT,d2TTEsurvC$treatment,mean)
##' }))
##' ##    [,1]      [,2]      [,3]      [,4]      [,5]
##' ## C    1 0.8878505 0.7943925 0.6915888 0.6728972
##' ## T    1 0.8972603 0.7739726 0.7054795 0.6712329
##' 
##' ## minevent estimator
##' eProdlim2C <- prodlimConditional(eventtime = paste0("eventtime",1:2),
##'                              status = paste0("status",1:2),
##'                              strata = "treatment", restriction = rep(0.75,2),
##'                              data = d2TTE, times = seqT)
##' eProdlim2C
##' ##    treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.9656701 0.06838695 0.8316342 1.0997061
##' ## 4          T 0.1875 0.8867831 0.06303572 0.7632354 1.0103308
##' ## 5          C 0.3750 0.8592152 0.06862112 0.7247203 0.9937101
##' ## 6          T 0.3750 0.7859406 0.06518644 0.6581776 0.9137037
##' ## 7          C 0.5625 0.7669355 0.07075982 0.6282488 0.9056223
##' ## 8          T 0.5625 0.7267865 0.06558515 0.5982419 0.8553310
##' ## 9          C 0.7500 0.7454226 0.07246591 0.6033920 0.8874532
##' ## 10         T 0.7500 0.6589528 0.06519537 0.5311722 0.7867333
##'
##' ##################################################
##' ## Two time to events - One censoring mechanism ##
##' ##################################################
##' 
##' set.seed(1)
##' argsTTE.1C <- list(scale.T = c(1,1), scale.C = c(0.5,2),
##'                    scale.censoring.T = 1,
##'                    frailty.C = 0, frailty.T = 2)
##' d2TTE.1C <- simBuyseTest(1000, argsTTE = argsTTE.1C, latent = TRUE)
##' d2TTE.1C$status1Uncensored <- 1
##' d2TTE.1C$status2Uncensored <- 1
##'
##' #### no censoring ####
##' d2TTEsurv.1C <- d2TTE.1C[d2TTE.1C$eventtime1Uncensored>0.75,]
##' do.call(cbind,lapply(seqT, function(iT){
##'  tapply(d2TTEsurv.1C$eventtime2Uncensored>iT,d2TTEsurv.1C$treatment,mean)
##' }))
##' ##    [,1]      [,2]      [,3]      [,4]      [,5]
##' ## C    1 0.8983051 0.8093220 0.7415254 0.6864407
##' ## T    1 0.8978873 0.7957746 0.7183099 0.6478873
##' 
##' ## minevent estimator
##' eProdlim2.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                                    status = paste0("status",1:2,"Uncensored"),
##'                                    strata = "treatment", restriction = rep(0.75,2),
##'                                    data = d2TTE.1C, times = seqT)
##' eProdlim2.1C
##' ##   treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.8983051 0.04688525 0.8064117 0.9901985
##' ## 4          T 0.1875 0.8978873 0.04965632 0.8005627 0.9952119
##' ## 5          C 0.3750 0.8093220 0.04485533 0.7214072 0.8972369
##' ## 6          T 0.3750 0.7957746 0.04730722 0.7030542 0.8884951
##' ## 7          C 0.5625 0.7415254 0.04309208 0.6570665 0.8259844
##' ## 8          T 0.5625 0.7183099 0.04530907 0.6295057 0.8071140
##' ## 9          C 0.7500 0.6864407 0.04148595 0.6051297 0.7677516
##' ## 10         T 0.7500 0.6478873 0.04334595 0.5629308 0.7328438
##' 
##' ## IPCW estimator
##' eIPCW2.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                                 status = paste0("status",1:2,"Uncensored"),
##'                                 strata = "treatment", restriction = rep(0.75,2),
##'                                 data = d2TTE.1C, times = seqT, method = "IPCW")
##' eIPCW2.1C
##' ##    treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.8983051 0.05527272 0.7899725 1.0066376
##' ## 4          T 0.1875 0.8978873 0.04870827 0.8024209 0.9933538
##' ## 5          C 0.3750 0.8093220 0.05250603 0.7064121 0.9122320
##' ## 6          T 0.3750 0.7957746 0.04587985 0.7058518 0.8856975
##' ## 7          C 0.5625 0.7415254 0.05025044 0.6430364 0.8400145
##' ## 8          T 0.5625 0.7183099 0.04355132 0.6329508 0.8036689
##' ## 9          C 0.7500 0.6864407 0.04831220 0.5917505 0.7811309
##' ## 10         T 0.7500 0.6478873 0.04128108 0.5669779 0.7287967

## * prodlimConditional (code)
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
    if(method == "IPCW"){
        index.censoring <- which(apply(data[status], MARGIN = 1, prod)==0)
        if(length(index.censoring)>0){
            maxEventTime <- apply(data[index.censoring,eventtime]*data[index.censoring,status], MARGIN = 1, FUN = max)
            maxCensoringTime <- apply(data[index.censoring,eventtime]*(1-data[index.censoring,status]), MARGIN = 1, FUN = max)
            if(any(maxCensoringTime < maxEventTime)){
                stop("Can only apply method = \"IPCW\" when the censoring mechanism is the same for all eventtimes. \n",
                     "(events observed after censoring of another type of event for some subjects) \n")
            }
            test.notIdentical <- apply(data[index.censoring,eventtime]*(1-data[index.censoring,status]), MARGIN = 1, FUN = function(iRow){any(diff(iRow[iRow>0])!=0)})
            if(any(test.notIdentical)){
                stop("Can only apply method = \"IPCW\" when the censoring mechanism is the same for all eventtimes. \n",
                     "(censoring times differs across type of event for some subjects) \n")
            }
            
        }
    }

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
    eKMsurv.marginal <- summary(eTTEM.marginal, times = 1, iid = iid)
    strata.names <- as.character(interaction(as.data.frame(eKMsurv.marginal)[,strata,drop=FALSE],sep="."))
    n.strata <- length(strata.names)

    ## ** joint
    eKMsurv.joint <- vector(length = n.times, mode = "list")
    eKMsurv.conditional <- vector(length = n.times, mode = "list")
    
    ## *** special case: times = 0
    ## a) Prob[X_1 > \tau_1, X_2 > 0 | strata]                 = Prob[X_1/tau_1 > 1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > 0 | strata]   = Prob[min(X_1/tau_1, X_2/tau_2) > 1 | strata]
    if(any(Utimes==0)){
        eKMsurv.joint[[which(Utimes==0)]] <- eKMsurv.marginal
        eKMsurv.conditional[[which(Utimes==0)]] <- cbind(eKMsurv.marginal[strata], time = 0, surv = 1)
        if(iid){
            eKMsurv.conditional[[which(Utimes==0)]]$se.surv <- 0
            eKMsurv.conditional[[which(Utimes==0)]]$lower <- 1
            eKMsurv.conditional[[which(Utimes==0)]]$upper <- 1
            attr(eKMsurv.conditional[[which(Utimes==0)]],"iid") <- matrix(0, nrow = NROW(data), ncol = NROW(eKMsurv.conditional[[which(Utimes==0)]]))
        }
    }

    ## *** 'normal' case: times > 0
    ## a) Prob[X_1 > \tau_1, X_2 > t | strata]                 = Prob[min(X_1/tau_1, X_2 / t) > 1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > t | strata]   = Prob[min(X_1/tau_1, X_2/tau_2, X_3 / t) > 1 | strata]

    if(any(Utimes>0)){

        ## prepare
        if(method == "minevent"){
            ff.joint <- stats::reformulate(termlabels = strata, response = "Hist(eventtime.min, status.min)")
        }else if(method == "IPCW"){
            data2$eventtimeIPCW <- apply(data2[paste0("eventtime",1:n.endpoint)],1,max)
            data2$statusIPCW <- 1-apply(data2[paste0("status",1:n.endpoint)],1,prod)
            data2$index <- 1:NROW(data2)

            ff.IPCW <- stats::reformulate(termlabels = strata, response = "Hist(eventtimeIPCW, statusIPCW)")
            eKM.IPCW <- prodlim(ff.IPCW, data = data2)
            eTTEM.IPCW <- BuyseTTEM(eKM.IPCW, treatment = strata[1], iid = iid)
            eKMsurv.IPCW <- summary(eTTEM.IPCW, times = pmax(times,max(restriction)), iid = iid)

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
                eKMsurv.joint[[iT]]$time <- Utimes[iT]
                eKMsurv.conditional[[iT]] <- cbind(eKMsurv.joint[[iT]][strata],
                                                   time = Utimes[iT],
                                                   surv = eKMsurv.joint[[iT]]$surv/eKMsurv.marginal$surv)
                if(iid){
                    term1.iid <- sweep(attr(eKMsurv.joint[[iT]],"iid"), MARGIN = 1, FUN = "/", STATS = eKMsurv.marginal$surv)
                    term2.iid <- - sweep(attr(eKMsurv.marginal,"iid"), MARGIN = 1, FUN = "/", STATS = eKMsurv.joint[[iT]]$surv/eKMsurv.marginal$surv^2)
                }
            }else if(method == "IPCW"){
                
                ## number of events in the observed process
                Ncount <- tapply((data2$eventtime.minM1>1) * (data2[[paste0("eventtime",n.endpoint)]]>Utimes[iT]), interaction(data2[strata], drop=TRUE), mean)
                
                ## assemble
                eKMsurv.conditional[[iT]] <- data.frame(eKMsurv.marginal[strata],
                                                        time = Utimes[iT],
                                                        surv = Ncount/(eKMsurv.IPCW$surv*eKMsurv.marginal$surv))
                if(iid){
                    denom.iid <- sweep(attr(eKMsurv.IPCW,"iid"), MARGIN = 1, FUN = "*", STATS = eKMsurv.marginal$surv) + sweep(attr(eKMsurv.marginal,"iid"), MARGIN = 1, FUN = "*", STATS = eKMsurv.IPCW$surv)
                    
                    Xstrata <- model.matrix(reformulate(termlabel = strata, intercept = FALSE), data2)
                    Xevent <- (data2$eventtime.minM1>1) * (data2[[paste0("eventtime",n.endpoint)]]>Utimes[iT])
                    Ncount.iid <- sweep(sweep(sweep(Xstrata, FUN = "*", MARGIN = 1, STATS = Xevent), FUN = "/", MARGIN = 2, STATS = colSums(Xstrata)), FUN = "-", MARGIN = 2, STATS = Ncount/NROW(data2))

                    term1.iid <- sweep(Ncount.iid, MARGIN = 2, FUN = "/", STATS = eKMsurv.IPCW$surv*eKMsurv.marginal$surv)
                    term2.iid <- -sweep(denom.iid, MARGIN = 2, FUN = "/", STATS = Ncount/(eKMsurv.IPCW$surv*eKMsurv.marginal$surv)^2)
                }
            }

            attr(eKMsurv.conditional[[iT]],"iid") <- term1.iid + term2.iid
            eKMsurv.conditional[[iT]]$se.surv <- sqrt(colSums(attr(eKMsurv.conditional[[iT]],"iid")^2))
            eKMsurv.conditional[[iT]]$lower <- eKMsurv.conditional[[iT]]$surv + stats::qnorm(0.025) * eKMsurv.conditional[[iT]]$se.surv
            eKMsurv.conditional[[iT]]$upper <- eKMsurv.conditional[[iT]]$surv + stats::qnorm(0.975) * eKMsurv.conditional[[iT]]$se.surv

        }
    }

    ## ** assemble
    out <- do.call(rbind,eKMsurv.conditional)
    rownames(out) <- NULL
    if(iid){
        attr(out,"iid") <- do.call(cbind,lapply(eKMsurv.conditional, function(iO){attr(iO,"iid")}))
        colnames(attr(out,"iid")) <- paste0(interaction(out[strata]),"(",out$time,")")
    }
    if(method == "IPCW"){
        attr(out,"IPCW") <- eKMsurv.IPCW
    }else if(method == "minevent"){
        attr(out,"joint") <- as.data.frame(as.list(do.call(rbind,eKMsurv.joint))) ## as.data.frame(as.list ... remove attributes like iid
    }
    attr(out,"marginal") <- as.data.frame(as.list(eKMsurv.marginal))
    
    ## ** export    
    return(out)
}


##----------------------------------------------------------------------
### prodlimConditional.R ends here
