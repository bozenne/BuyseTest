### prodlimConditional.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 30 2026 (16:04) 
## Version: 
## Last-Updated: jul  9 2026 (16:54) 
##           By: Brice Ozenne
##     Update #: 735
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * prodlimConditional (documentation)
##' @title Conditional survival estimator
##' @description Estimator of the survival of one time to event conditional on other time to events: \eqn{\mathbb{P}\left[T_m>t|T_1>\tau_1,...,T_{m-1}>\tau_{m-1}\right]}.
##' The last time to event (say \eqn{T_m}) is the one for which the survival is evaluated while the other ones (\eqn{T_1,...,T_{m-1}}) are use to condition.
##'
##' @param eventtime [character vector] name of the columns in argument \code{data} containing the time to events.
##' @param status [character vector] name of the columns in argument \code{data} containing the type of event (0: censoring, 1: event of interest).
##' @param strata [character vector] name of the columns in argument \code{data} indicating the groups for which a separate estimator should be used.
##' @param data [data.frame] dataset.
##' @param restriction [numeric vector] timepoints (\eqn{\tau_1,...,\tau_{m-1}}) used to condition all but the last time to event.
##' @param times [numeric vector] timepoints at which to evaluate the survival. No length requirement.
##' @param method [character] type of estimator for the joint survival distribution \eqn{\mathbb{P}\left[T_1>\tau_1,...,T_{m-1}>\tau_{m-1},T_m>t\right]}: \itemize{
##' \item \code{"prodlim"} Kaplan Meier estimators for the minimum normalized time to events \eqn{\frac{T_1}{\tau_1},...,\frac{T_{m-1}}{\tau_{m-1}},\frac{T_m}{t}}.
##' \item \code{"IPCW"} Inverse Probability of Censoring Weights assuming same censoring process for all timepoints - see Parner et al. 2025.
##' }
##' @param iid [logical] should the influence function of the estimator and standard errors for the survival be exported?
##' @param iid.surv [character] Estimator of the survival used when computing the influence function.
##' Can be the product limit estimator (\code{"prodlim"}) or an exponential approximation (\code{"exp"}, same as in \code{riskRegression::predictCox}).
##' @param type [character] should the survival that is output be the \code{"joint"} (\eqn{\mathbb{P}\left[T_1>\tau_1,...,T_{m-1}>\tau_{m-1},T_m\right]}) or the \code{"conditional"} survival (\eqn{\mathbb{P}\left[T_m>t|T_1>\tau_1,...,T_{m-1}>\tau_{m-1}\right]})?
##' 
##' @details The conditional distribution is computed as a ratio between the full joint survival \eqn{\mathbb{P}\left[T_1>\tau_1,...,T_{m-1}>\tau_{m-1},T_m>t\right]}
##' and the partial joint survival \eqn{\mathbb{P}\left[T_1>\tau_1,...,T_{m-1}>\tau_{m-1}\right]}.
##' The later is computed as \eqn{\mathbb{P}\left[\min\left(\frac{T_1}{\tau_1},...,\frac{T_{m-1}}{\tau_{m-1}}\right)>1\right]},
##' i.e., using a Kaplan-Meier estimator for the minimum normalized time to event.
##'
##' @references Parner, E. T., & Overgaard, M. (2025). Estimation of win, loss probabilities, and win ratio based on right-censored event data. Scandinavian Journal of Statistics, 52(1), 170–184. https://doi.org/10.1111/sjos.12734
##' 
##' @return A data.frame with an attribute \code{iid} if the argument \code{iid} was set to TRUE.
##' 

## * prodlimConditional (examples)
##' @examples
##' library(prodlim)
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
##' ## min estimator
##' \dontrun{
##' eProdlim2 <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                                 status = paste0("status",1:2,"Uncensored"),
##'                                 strata = "treatment", restriction = rep(0.75,2),
##'                                 data = d2TTE)
##' summary(eProdlim2, times = seqT)
##' }
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 2          C 0.188    912       0      0 0.888  0.0185 0.852 0.925
##' ## 3          C 0.375    836       0      0 0.809  0.0233 0.764 0.855
##' ## 4          C 0.562    755       0      0 0.726  0.0262 0.674 0.777
##' ## 5          C 0.750    688       0    688 0.688  0.0272 0.635 0.742
##' ## 6          T 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 7          T 0.188    569       0      0 0.899  0.0221 0.856 0.943
##' ## 8          T 0.375    432       0      0 0.782  0.0300 0.724 0.841
##' ## 9          T 0.562    359       0      0 0.695  0.0334 0.629 0.760
##' ## 10         T 0.750    311       0    311 0.636  0.0348 0.568 0.705
##' 
##' ## Peron (ignores conditionning)
##' \dontrun{
##' e.peron <- prodlim(Hist(eventtime2Uncensored, status2Uncensored) ~ treatment,
##'                    data = d2TTE)
##' summary(e.peron, times = seqT)
##' }
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000   1000       0      0 1.000 0.00000 1.000 1.000
##' ## 2          C 0.188    912       0      0 0.912 0.00896 0.894 0.930
##' ## 3          C 0.375    836       0      0 0.836 0.01171 0.813 0.859
##' ## 4          C 0.562    755       0      0 0.755 0.01360 0.728 0.782
##' ## 5          C 0.750    688       0      0 0.688 0.01465 0.659 0.717
##' ## 6          T 0.000   1000       0      0 1.000 0.00000 1.000 1.000
##' ## 7          T 0.188    569       0      0 0.569 0.01566 0.538 0.600
##' ## 8          T 0.375    432       0      0 0.432 0.01566 0.401 0.463
##' ## 9          T 0.562    359       0      0 0.359 0.01517 0.329 0.389
##' ## 10         T 0.750    311       0      0 0.311 0.01464 0.282 0.340
##' ## ok for C but too low for T
##' 
##' #### with censoring ####
##' 
##' ## condition on staying at risk (inefficient, incorrect se?)
##' d2TTEsurvC <- d2TTE[d2TTE$eventtime1>0.75,]
##' eNaiveC <- prodlim(Hist(eventtime2, status2) ~ treatment, data = d2TTEsurvC)
##' summary(eNaiveC, times = seqT)
##' 
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000    107       0      0 1.000  0.0000 1.000 1.000
##' ## 2          C 0.188     79       0      0 0.882  0.0322 0.819 0.945
##' ## 3          C 0.375     61       0      0 0.798  0.0420 0.716 0.880
##' ## 4          C 0.562     44       0      0 0.699  0.0507 0.600 0.799
##' ## 5          C 0.750     37       0      0 0.682  0.0524 0.579 0.785
##' ## 6          T 0.000    146       0      0 1.000  0.0000 1.000 1.000
##' ## 7          T 0.188    109       0      0 0.904  0.0253 0.855 0.954
##' ## 8          T 0.375     80       0      0 0.786  0.0378 0.711 0.860
##' ## 9          T 0.562     64       0      0 0.733  0.0419 0.651 0.816
##' ## 10         T 0.750     47       0      0 0.709  0.0440 0.622 0.795
##' 
##' ## min estimator
##' eProdlim2C <- prodlimConditional(eventtime = paste0("eventtime",1:2),
##'                              status = paste0("status",1:2),
##'                              strata = "treatment", restriction = rep(0.75,2),
##'                              data = d2TTE)
##' summary(eProdlim2C, times = seqT)
##' 
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 2          C 0.188    751       0      0 0.959  0.0384 0.883 1.034
##' ## 3          C 0.375    576       0      0 0.860  0.0493 0.763 0.957
##' ## 4          C 0.562    438       0      0 0.763  0.0558 0.654 0.873
##' ## 5          C 0.750    336       0    336 0.745  0.0592 0.629 0.861
##' ## 6          T 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 7          T 0.188    473       0      0 0.887  0.0324 0.823 0.950
##' ## 8          T 0.375    300       0      0 0.793  0.0463 0.703 0.884
##' ## 9          T 0.562    209       0      0 0.727  0.0509 0.627 0.827
##' ## 10         T 0.750    152       0    152 0.659  0.0532 0.555 0.763
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
##' 
##' ## empirical average estimator
##' d2TTEsurv.1C <- d2TTE.1C[d2TTE.1C$eventtime1Uncensored>0.75,]
##' do.call(cbind,lapply(seqT, function(iT){
##'  tapply(d2TTEsurv.1C$eventtime2Uncensored>iT,d2TTEsurv.1C$treatment,mean)
##' }))
##' ##    [,1]      [,2]      [,3]      [,4]      [,5]
##' ## C    1 0.8983051 0.8093220 0.7415254 0.6864407
##' ## T    1 0.8978873 0.7957746 0.7183099 0.6478873
##' 
##' ## min estimator
##' eProdlim2.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                                    status = paste0("status",1:2,"Uncensored"),
##'                                    strata = "treatment", restriction = rep(0.75,2),
##'                                    data = d2TTE.1C)
##' summary(eProdlim2.1C, times = seqT)
##' ##    treatment  time n.risk n.event n.lost  surv se.surv lower upper
##' ## 1          C 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 2          C 0.188    905       0      0 0.898  0.0178 0.863 0.933
##' ## 3          C 0.375    831       0      0 0.809  0.0233 0.764 0.855
##' ## 4          C 0.562    750       0      0 0.742  0.0261 0.690 0.793
##' ## 5          C 0.750    686       0    686 0.686  0.0278 0.632 0.741
##' ## 6          T 0.000   1000       0      0 1.000  0.0000 1.000 1.000
##' ## 7          T 0.188    564       0      0 0.898  0.0198 0.859 0.937
##' ## 8          T 0.375    435       0      0 0.796  0.0265 0.744 0.848
##' ## 9          T 0.562    362       0      0 0.718  0.0295 0.660 0.776
##' ## 10         T 0.750    308       0    308 0.648  0.0313 0.587 0.709
##' 
##' ## IPCW estimator
##' \dontrun{
##' eIPCW2.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2,"Uncensored"),
##'                                 status = paste0("status",1:2,"Uncensored"),
##'                                 strata = "treatment", restriction = rep(0.75,2),
##'                                 data = d2TTE.1C, times = seqT, method = "IPCW")
##' eIPCW2.1C
##' }
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
##'
##' #### with censoring ####
##' ## empirical average estimator
##' d2TTEsurvC.1C <- d2TTE.1C[d2TTE.1C$eventtime1>0.75,]
##' summary(prodlim(Hist(eventtime2, status2) ~ treatment, data = d2TTEsurvC.1C),
##'         times = seqT)
##' 
##' ## min estimator
##' eProdlim2C.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2),
##'                                    status = paste0("status",1:2),
##'                                    strata = "treatment", restriction = rep(0.75,2),
##'                                    data = d2TTE.1C, times = seqT)
##' eProdlim2C.1C
##' ##    treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.9142038 0.06255723 0.7915939 1.0368137
##' ## 4          T 0.1875 0.8998606 0.06149527 0.7793320 1.0203891
##' ## 5          C 0.3750 0.8413301 0.06041755 0.7229139 0.9597464
##' ## 6          T 0.3750 0.7929164 0.05891098 0.6774530 0.9083798
##' ## 7          C 0.5625 0.7782851 0.05833815 0.6639445 0.8926258
##' ## 8          T 0.5625 0.7331788 0.05707721 0.6213095 0.8450481
##' ## 9          C 0.7500 0.7377882 0.05705952 0.6259536 0.8496228
##' ## 10         T 0.7500 0.6800372 0.05556328 0.5711352 0.7889392
##' 
##' ## IPCW estimator
##' \dontrun{
##' eIPCW2C.1C <- prodlimConditional(eventtime = paste0("eventtime",1:2),
##'                                 status = paste0("status",1:2),
##'                                 strata = "treatment", restriction = rep(0.75,2),
##'                                 data = d2TTE.1C, times = seqT, method = "IPCW")
##' eIPCW2C.1C
##' }
##' ##    treatment   time      surv    se.surv     lower     upper
##' ## 1          C 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 2          T 0.0000 1.0000000 0.00000000 1.0000000 1.0000000
##' ## 3          C 0.1875 0.8972191 0.08997123 0.7208788 1.0735595
##' ## 4          T 0.1875 0.9251078 0.07889512 0.7704762 1.0797394
##' ## 5          C 0.3750 0.8296865 0.08661343 0.6599273 0.9994457
##' ## 6          T 0.3750 0.8214957 0.07448751 0.6755029 0.9674886
##' ## 7          C 0.5625 0.7621539 0.08309401 0.5992926 0.9250151
##' ## 8          T 0.5625 0.7400862 0.07078031 0.6013594 0.8788131
##' ## 9          C 0.7500 0.7139163 0.08046743 0.5562030 0.8716296
##' ## 10         T 0.7500 0.6734785 0.06756046 0.5410624 0.8058946
##' 

## * prodlimConditional (code)
##' @export
prodlimConditional <- function(eventtime, status, strata, data, restriction, 
                               method = "prodlim", iid = TRUE, iid.surv = "exp", type = "conditional"){

    tol <- 1e-12

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
    if(any(data[eventtime]<=0)){
        stop("Argument \'eventtimes\' should correspond to columns with strictly positive values. \n")
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
    if(any(unlist(data[status]) %in% 0:2 == FALSE)){
        stop("Argument \'status\' should only take value 0, 1, or 2. \n")
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
        restriction.endpoint <- restriction[n.endpoint]
        restriction <- restriction[1:(n.endpoint-1)]
    }else if(length(restriction)==(n.endpoint-1)){
        restriction.endpoint <- NULL
    }else{
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

    ## *** method
    method <- match.arg(method, c("IPCW","prodlim"))
    if(method == "IPCW"){
        index.censoring <- which(apply(data[status]==1, MARGIN = 1, prod)==0)
        if(length(index.censoring)>0){
            maxEventTime <- apply(data[index.censoring,eventtime]*(data[index.censoring,status]==1), MARGIN = 1, FUN = max)
            maxCensoringTime <- apply(data[index.censoring,eventtime]*(data[index.censoring,status]==0), MARGIN = 1, FUN = max)
            if(any(maxCensoringTime < maxEventTime)){
                stop("Can only apply method = \"IPCW\" when the censoring mechanism is the same for all eventtimes. \n",
                     "(events observed after censoring of another type of event for some subjects) \n")
            }
            test.notIdentical <- apply(data[index.censoring,eventtime]*(data[index.censoring,status]==0), MARGIN = 1, FUN = function(iRow){any(diff(iRow[iRow>0])!=0)})
            if(any(test.notIdentical)){
                stop("Can only apply method = \"IPCW\" when the censoring mechanism is the same for all eventtimes. \n",
                     "(censoring times differs across type of event for some subjects) \n")
            }
            
        }
    }
    
    ## *** type
    type <- match.arg(type, c("conditional","joint"))

    ## ** create a normalized dataset

    ## *** normalize names of the variables and convert strata into factor
    data2 <- cbind(data[eventtime], data[status], interaction(data[strata], sep = "..XX.."), 1:NROW(data))
    names(data2) <- c(paste0("eventtime",1:n.endpoint),
                      paste0("status",1:n.endpoint),
                      "strata", "index")
    
    lsIndex.strata <- tapply(1:NROW(data),data2$strata, FUN = identity, simplify = FALSE)
    level.strata <- names(lsIndex.strata)
    n.strata <- length(level.strata)
    
    index.Ustrata <- which(!duplicated(data2$strata))
    indexOrder.Ustrata <- index.Ustrata[match(data2[index.Ustrata,"strata"], level.strata)]

    ## *** in presence of ties between event and censoring, break the ties by putting censoring just after events
    if(method == "IPCW"){ 
        for(iE in 1:n.endpoint){ ## iE <- 1
            for(iS in 1:n.strata){ ## iS <- 1
                iEventtime <- data2[lsIndex.strata[[iS]], paste0("eventtime",iE)]
                iStatus <- data2[lsIndex.strata[[iS]], paste0("status",iE)]
                iTiestime <- intersect(unique(iEventtime[iStatus==0]), unique(iEventtime[iStatus>0]))
                if(length(iTiestime)>0){ 
                    data2[lsIndex.strata[[iS]],paste0("eventtime",iE)] <- iEventtime + tol*(iStatus==0)*(iEventtime %in% iTiestime)
                }
            }
        }
    }

    ## *** generate normalized time to events and associated min() variables
    ## normalize: X_1/\tau_1, X_2/tau_2, ...
    if(n.endpoint == 2){
        data2$eventtimeN1 <- data2$eventtime1/restriction
    }else{
        data2[paste0("eventtimeN",1:(n.endpoint-1))] <- sweep(data2[paste0("eventtime",1:(n.endpoint-1))], MARGIN = 1, FUN = "/", STATS = restriction)
    }

    ## generate min variables: initialize
    data2$eventtime.minM1 <- data2$eventtimeN1
    data2$status.minM1 <- data2$status1
    ##                         loop
    if(n.endpoint>2){
        for(iE in 2:(n.endpoint-1)){
            data2$eventtime.minM1 <- pmin(data2$eventtime.minM1, data2[[paste0("eventtimeN",iE)]])
            data2$status.minM1 <- ifelse(data2$eventtime.minM1 < data2[[paste0("eventtimeN",iE)]], data2$status.minM1, ## conditioning event occured first
                                  ifelse(data2$eventtime.minM1 > data2[[paste0("eventtimeN",iE)]], data2[[paste0("status",iE)]], ## event of interest occured first
                                  ifelse((data2$status.minM1==0) & (data2[[paste0("status",iE)]]==0), 0, pmin(data2$status.minM1,data2[[paste0("status",iE)]])))) ## same time: if both censored then censored
            ## otherwise if any event then event otherwise death
            
        }        
    }

    ## *** find jump times

    ## fit unconditional prodlim
    data.out <- cbind(data[utils::tail(eventtime,1)], data[utils::tail(status,1)],data[strata])
    if(!is.null(restriction.endpoint)){
        data.out[[utils::tail(eventtime,1)]] <- pmin(data[[utils::tail(eventtime,1)]], restriction)
        data.out[[utils::tail(status,1)]] <- data[[utils::tail(status,1)]] * (data[[utils::tail(eventtime,1)]] <= restriction)
    }
    ff.prodlim <- reformulate(termlabels = strata, response = paste0("Hist(",utils::tail(eventtime,1),",",utils::tail(status,1),"==1)"))
    out <- prodlim::prodlim(ff.prodlim, data = data.out)
    out$surv[] <- NA
    if(iid){
        out$iid.surv <- matrix(0, nrow = NROW(data2), ncol = length(out$surv))
    }

    ## extract times
    times <- out$time
    Utimes <- sort(unique(times))
    n.times <- length(times)
    
    ## *** trim jump times
    ## first time at which the 'target' event is active:
    ## eventtime(m)/times < min(eventtime(1)/restriction(m),...,eventtime(m-1)/restriction(m-1))
    ## times > eventtime(m) / min(eventtime(1)/restriction(m),...,eventtime(m-1)/restriction(m-1))
    first.times <- tapply(data2[[paste0("eventtime",n.endpoint)]]/data2$eventtime.minM1, data2$strata, min)
    Utimes.first <- Utimes[Utimes < min(first.times)]

    ## identify observations affected by each 'jump time'
    datatime2jump <- match(data[[utils::tail(eventtime,1)]],times)
    datatime2jump[is.na(datatime2jump)] <- length(times) ## when cannot be found it is because the timepoints is beyond the restriction time or censored

    ## remove plateaux

    ## *** reconstruct strata
    out.grid <- data.frame(times = out$time,
                           strata.index = do.call("c",lapply(1:length(out$size.strata), function(iS){rep(iS, length = out$size.strata[iS])})))
    out.grid$strata <- interaction(out$X, sep = "..XX..")[out.grid$strata.index]

    if(length(levels(out.grid$strata))!=length(level.strata)){
        stop("Something went wrong when retrieve the strata levels from the prodlim object. \n",
             "Mismatch in size. \n")
    }
    if(any(levels(out.grid$strata) %in% level.strata == FALSE) || any(level.strata %in% levels(out.grid$strata) == FALSE)){
        stop("Something went wrong when retrieve the strata levels from the prodlim object. \n",
             "Mismatch between the values. \n")
    }

    ## ** marginal
    ## a) Prob[X_1 > \tau_1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2 | strata] = Prob[X_1/\tau_1 > 1, X_2/\tau_2 > 1 | strata] = Prob[min(X_1/\tau_1,X_2/\tau_2) > 1 | strata]
    if(type == "conditional" || length(Utimes.first)>0){

        ## *** fit marginal model
        eKMsurv.marginal <- .prodlim2(eventtime = "eventtime.minM1", status = "status.minM1", strata = "strata",
                                      data = transform(data2, status.minM1 = status.minM1==1), time = 1, iid = iid, iid.surv = iid.surv)

        ## eKM.marginal <- prodlim::prodlim(Hist(eventtime.minM1, status.minM1) ~ strata, data = data2)
        ## eTTEM.marginal <- BuyseTTEM(eKM.marginal, treatment = "strata", iid = iid, iid.surv = iid.surv)
        ## summary(eTTEM.marginal, times = 1, iid = iid)
        
        ## *** update object at time 0 
        if(type == "conditional" && length(Utimes.first)>0){
            ## a) Prob[X_2 > 0 | X_1 > \tau_1, strata]                 = 1
            ## b) Prob[X_3 > 0 | X_1 > \tau_1, X_2 > \tau_2, strata]   = 1
            index0.out <- which(out.grid$times < min(first.times))
            
            out$surv[index0.out] <- 1
            out$se.surv[index0.out] <- 0
            out$lower[index0.out] <- 1
            out$upper[index0.out] <- 1
            out$iid.surv[,index0.out] <- 0
        }else if(type == "joint" && length(Utimes.first)>0){
            ## a) Prob[X_1 > \tau_1, X_2 > 0 | strata]                 = Prob[X_1/tau_1 > 1 | strata]
            ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > 0 | strata]   = Prob[min(X_1/tau_1, X_2/tau_2) > 1 | strata]
            index0.out <- which(out.grid$times %in% Utimes.first)
            index0.prodlim <- match(out.grid[index0.out,"strata"],eKMsurv.marginal$strata)
            
            out$surv[index0.out] <- eKMsurv.marginal[index0.prodlim,"estimate"]
            out$se.surv[index0.out] <- eKMsurv.marginal[index0.prodlim,"se"]
            out$lower[index0.out] <- eKMsurv.marginal[index0.prodlim,"lower"]
            out$upper[index0.out] <- eKMsurv.marginal[index0.prodlim,"upper"]
            out$iid.surv[,index0.out] <- attr(eKMsurv.marginal,"iid")[,index0.prodlim]
        }
    }

    ## ** joint
    ## a) Prob[X_1 > \tau_1, X_2 > t | strata]                 = Prob[min(X_1/tau_1, X_2 / t) > 1 | strata]
    ## b) Prob[X_1 > \tau_1, X_2 > \tau_2, X_3 > t | strata]   = Prob[min(X_1/tau_1, X_2/tau_2, X_3 / t) > 1 | strata]
    if(all(Utimes < min(first.times))){return(out)}

    if((method == "IPCW")){
        data2$eventtimeIPCW <- apply(data2[paste0("eventtime",1:n.endpoint)],1,max)
        data2$statusIPCW <- 1-apply(data2[paste0("status",1:n.endpoint)]>0,1,prod)
        data2$index <- 1:NROW(data2)
        
        ff.IPCW <- stats::reformulate(termlabels = strata, response = "Hist(eventtimeIPCW, statusIPCW)")
        eKM.IPCW <- prodlim(ff.IPCW, data = data2)
        eTTEM.IPCW <- BuyseTTEM(eKM.IPCW, treatment = strata[1], iid = iid, iid.surv = iid.surv)
        eKMsurv.IPCW <- summary(eTTEM.IPCW, times = pmax(times,max(restriction)), iid = iid)
    }else{
        data2$statusN.target <- data2[[paste0("status",n.endpoint)]]        
    }

    ## loop over times
    eventtime.min.save <- rep(-1, NROW(data2))
    status.min.save <- rep(-1, NROW(data2))
    Utimes.N0 <- Utimes[Utimes >= min(first.times)]

    for(iT in 1:length(Utimes.N0)){ ## iT <- 1

        ## link to output
        iIndex0.out <- which(out.grid$times==Utimes.N0[iT])
        iIndex0.prevout <- which(out.grid$times == Utimes.N0[iT-1])
        
        if(method == "prodlim"){

            ## *** create min time to event
            data2$eventtimeN.target <- data2[[paste0("eventtime",n.endpoint)]]/Utimes.N0[iT]
            data2$eventtime.min <- pmin(data2$eventtime.minM1, data2$eventtimeN.target)
            
            iIndex.changeTime <- which(data2$eventtime.min!=eventtime.min.save)
            browser()
            datatime2jump
            data2[iIndex.changeTime,"status.min"] <- with(data2[iIndex.changeTime,,drop=FALSE],
                                                          ifelse(eventtime.minM1 < eventtimeN.target, status.minM1, ## conditioning event occured first
                                                          ifelse(eventtime.minM1 > eventtimeN.target, statusN.target, ## event of interest occured first
                                                          ifelse((status.minM1==0) & (statusN.target==0), 0, pmin(status.minM1,statusN.target)))) ## same time: if both censored then censored
                                                          ##            otherwise event of interest then death
                                                          )
            

            ## *** check if something has changed since last time
            iTest.status <- all(data2$status.min == status.min.save) 
            iTest.time <- all((data2$eventtime.min<=1) == (eventtime.min.save<=1))
            iTest.strata <- identical(out.grid$strata[iIndex0.out], out.grid$strata[iIndex0.prevout])
            if(iTest.status & iTest.time & iTest.strata){
                out$surv[iIndex0.out] <- out$surv[iIndex0.prevout]
                if(iid){
                    out$iid.surv[,iIndex0.out] <- out$iid.surv[,iIndex0.prevout]
                }
                next                
            }
            eventtime.min.save <- data2$eventtime.min
            status.min.save <- data2$status.min
            
            ## *** fit KM
            iKMsurv.joint <- .prodlim2(eventtime = "eventtime.min", status = "status.min", strata = "strata",
                                       data = transform(data2, status.min = status.min==1), time = 1, iid = iid, iid.surv = iid.surv)

            ## link to output
            iIndex0.prodlim <- match(out.grid[iIndex0.out,"strata"],eKMsurv.marginal$strata)
            
            ## *** store (possibly move from joint to conditional)
            if(type == "conditional"){

                out$surv[iIndex0.out] <- (iKMsurv.joint$estimate/eKMsurv.marginal$estimate)[iIndex0.prodlim]
                
                if(iid & any(out$surv[iIndex0.out] != 1)){ ## if survival is 1 then no uncertainty
                    term1.iid <- sweep(attr(iKMsurv.joint,"iid"), MARGIN = 1, FUN = "*", STATS = (iKMsurv.joint$estimate != 1)/eKMsurv.marginal$estimate)
                    term2.iid <- - sweep(attr(eKMsurv.marginal,"iid"), MARGIN = 1, FUN = "*", STATS = (iKMsurv.joint$estimate != 1)*iKMsurv.joint$estimate/eKMsurv.marginal$estimate^2)
                    out$iid.surv[,iIndex0.out] <- (term1.iid + term2.iid)[,iIndex0.prodlim]
                    
                }
            }else if (type == "joint"){
                out$surv[iIndex0.out] <- iKMsurv.joint[iIndex0.prodlim,"estimate"]
                out$iid.surv[,iIndex0.out] <- attr(iKMsurv.joint,"iid")[,iIndex0.prodlim]
            }
            
        }else if(method == "IPCW"){
            
            ## number of events in the observed process
            Ncount <- tapply((data2$eventtime.minM1>1) * (data2[[paste0("eventtime",n.endpoint)]]>Utimes.N0[iT]), interaction(data2[strata], drop=TRUE), mean)
            
            ## assemble
            eKMsurv.joint[[iT]] <- data.frame(eKMsurv.marginal[strata],
                                              time = Utimes.N0[iT],
                                              n.risk = NA,
                                              n.event = NA,
                                              n.lost = NA,
                                              surv = Ncount/eKMsurv.IPCW$surv)
            eKMsurv.conditional[[iT]] <- data.frame(eKMsurv.marginal[strata],
                                                    time = Utimes.N0[iT],
                                                    surv = Ncount/(eKMsurv.IPCW$surv*eKMsurv.marginal$surv))
            if(iid){
                denom.iid <- sweep(attr(eKMsurv.IPCW,"iid"), MARGIN = 1, FUN = "*", STATS = eKMsurv.marginal$surv) + sweep(attr(eKMsurv.marginal,"iid"), MARGIN = 1, FUN = "*", STATS = eKMsurv.IPCW$surv)
                
                Xstrata <- model.matrix(reformulate(termlabel = strata, intercept = FALSE), data2)
                Xevent <- (data2$eventtime.minM1>1) * (data2[[paste0("eventtime",n.endpoint)]]>Utimes.N0[iT])
                Ncount.iid <- sweep(sweep(sweep(Xstrata, FUN = "*", MARGIN = 1, STATS = Xevent), FUN = "/", MARGIN = 2, STATS = colSums(Xstrata)), FUN = "-", MARGIN = 2, STATS = Ncount/NROW(data2))

                term1.iid <- sweep(Ncount.iid, MARGIN = 2, FUN = "/", STATS = eKMsurv.IPCW$surv*eKMsurv.marginal$surv)
                term2.iid <- -sweep(denom.iid, MARGIN = 2, FUN = "/", STATS = Ncount/(eKMsurv.IPCW$surv*eKMsurv.marginal$surv)^2)
            }
            
            if(iid){
                attr(eKMsurv.conditional[[iT]],"iid") <- term1.iid + term2.iid
                eKMsurv.conditional[[iT]]$se.surv <- sqrt(colSums(attr(eKMsurv.conditional[[iT]],"iid")^2))
                eKMsurv.conditional[[iT]]$lower <- eKMsurv.conditional[[iT]]$surv + stats::qnorm(0.025) * eKMsurv.conditional[[iT]]$se.surv
                eKMsurv.conditional[[iT]]$upper <- eKMsurv.conditional[[iT]]$surv + stats::qnorm(0.975) * eKMsurv.conditional[[iT]]$se.surv
            }
        }
    }

    ## ** assemble
    if(iid){
        out$se.surv <- sqrt(colSums(out$iid.surv^2))
        out$lower <- out$surv + stats::qnorm(0.025) * out$se.surv
        out$upper <- out$surv + stats::qnorm(0.975) * out$se.surv
    }

    out$hazard <- unlist(
        tapply(out$surv, out.grid$strata, function(iSurv){1 - iSurv/c(1,iSurv[-length(iSurv)])}, simplify = FALSE)
    )
    if(method == "IPCW"){
        out$IPCW <- eKMsurv.IPCW
    }
    
    ## ** export
    return(out)
}

## * .prodlim2
##' @title Kaplan-Meier estimator for a single timepoint
##' @description Fast alternative to prodlim when only a single timepoint is of interest
##'
##' @param eventtime [character] name of the column containing the event time variable.
##' @param status [character] name of the column containing the status variable (0/1).
##' @param data [data.frame] dataset.
##' @param time [numeric] timepoint of interest.
##' @param iid [logical] should the influence function, se, CI of the survival be output.
##' @param iid.surv [character] survival estimator used when evaluating the influence function  (\code{"prodlim"}, \code{"exp"})
##' 
##' @noRd
##' @keywords internal
##' @examples
##' #######################
##' #### simulate data ####
##' #######################
##' set.seed(1)
##' dd <- simBuyseTest(100)
##'
##' ######################
##' #### Kaplan Meier ####
##' ######################
##' e.prodlim <- prodlim(Hist(eventtime, status) ~ 1, data = dd)
##' summary(e.prodlim, times = 1)
##' ##   time n.risk n.event n.lost   surv se.surv  lower upper
##' ## 1    1      6       0      0 0.0822  0.0286 0.0262 0.138
##'  
##' test <- .prodlim2("eventtime","status", data = dd, time = 1, iid = TRUE, iid.surv = "exp")
##' test
##' ##   time   estimate         se     lower     upper
##' ## 1    1 0.08221101 0.02774669 0.0278285 0.1365935
##'
##' #################################
##' #### Stratified Kaplan Meier ####
##' #################################
##' eS.prodlim <- prodlim(Hist(eventtime, status) ~ treatment, data = dd)
##' summary(eS.prodlim, times = 1)
##' ## 1  treatment time n.risk n.event n.lost   surv se.surv   lower upper
##' ## 11         C    1      3       0      0 0.0940  0.0450 0.00569 0.182
##' ## 12         T    1      3       0      0 0.0713  0.0365 0.00000 0.143
##'
##' testS <- .prodlim2("eventtime","status", strata = "treatment", data = dd, time = 1, iid = TRUE, iid.surv = "exp")
##' testS
##' ##   time strata   estimate         se      lower     upper
##' ## 1    1      C 0.09396464 0.03782186 0.01983516 0.1680941
##' ## 2    1      T 0.07134397 0.03067393 0.01122418 0.1314638
##'
##' eS.TTEM <- BuyseTTEM(eS.prodlim, treatment = "treatment", iid = TRUE, iid.surv =  "exp")
##' GS <- summary(eS.TTEM, time = 1, iid = TRUE)
##' ##   treatment time n.risk n.event n.lost       surv    se.surv       lower     upper
##' ## 1         C    1      3       0      0 0.09396464 0.04212759 0.011396088 0.1765332
##' ## 2         T    1      3       0      0 0.07134397 0.03468403 0.003364515 0.1393234
##' range(attr(GS,"iid") - attr(testS,"iid"))
.prodlim2 <- function(eventtime, status, strata = NULL, data, time, iid, iid.surv){

    ## ** prepare
    ## *** only events before time are used in the KM estimator (along with the sample size)
    index.time <- which(data[[eventtime]] <= time)
    
    ## *** find strata
    if(is.null(strata)){
        n.strata <- 1
        survival <- data.frame(time = time, estimate = NA, se = NA, lower = NA, upper = NA)
    }else{
        if(is.factor(strata)){
            level.strata <- levels(data[[strata]])
        }else{
            level.strata <- sort(unique(data[[strata]]))
        }
        n.strata <- length(level.strata)
        survival <- data.frame(time = rep(time, n.strata), strata = level.strata, estimate = rep(NA, n.strata))
    }
    
    if(iid){
        attr(survival,"iid") <- matrix(0, ncol = n.strata, nrow = NROW(data))
        if(!is.null(strata)){
            colnames(attr(survival,"iid")) <- level.strata
        }
    }
    
    ## ** loop over strata
    for(iS in 1:n.strata){

        ## index
        if(is.null(strata)){
            iIndex.strata <- 1:NROW(data)
        }else{
            iIndex.strata <- which(data[[strata]]==level.strata[iS])
        }
        iN.obs <- length(iIndex.strata)

        ## *** restrict data to strata and prediction time horizon
        iIndex.strataTime <- intersect(iIndex.strata, index.time)
        iData2 <- as.data.frame(data)[iIndex.strataTime, c(eventtime,status), drop=FALSE]                        
        iN.obsTime <- NROW(iData2)
        
        ## *** normalize dataset
        names(iData2) <- c("time","status")
        iData2$index <- iIndex.strataTime ## add index of the original dataset
        iData2 <- iData2[order(iData2$time),,drop=FALSE]
        iData2$group <- rank(iData2$time) ## faster than as.numeric(as.factor(iData2$eventtime.min))    
        
        ## *** evaluate survival
        if(any(duplicated(iData2$time))){
            iDeath <- tapply(iData2$status, iData2$group, sum)
            iAtRisk <- iN.obs - c(0,cumsum(tapply(iData2$group, iData2$group, length)))[1:length(iDeath)]                            
        }else{
            iDeath <- iData2$status
            iAtRisk <- seq(from = iN.obs, to = iN.obs-iN.obsTime+1, by = -1)
        }

        ## product limit estimator
        survival$estimate[iS] <- prod(1-iDeath/iAtRisk)
        
        ## *** evaluate influence function
        ## IF_i(t) = - S(t) \int_0^t (dN_i(s) - Y_i(s)dLambda(s))/pi(s)
        ##         = - S(t) \int_0^(min(t,T_i)) (dN_i(s) - dLambda(s))/(Y(s)/n)
        ##         = - n S(t) (\Ind[T_i<=t] N_i(T_i)/Y(T_i) - \int_0^(min(t,T_i)) dN(s)/Y^2(s)
        if(iid){

            ## prepare
            iLambda <- cumsum(iDeath/iAtRisk^2)
            iSurv.for.iid <- switch(iid.surv,
                                   "prodlim" = survival$estimate[iS],
                                   "exp" = exp(-sum(iDeath/iAtRisk)))

            iCumhazard.iid <- rep(0, NROW(data))

            ## first term
            iCumhazard.iid[iData2$index] <- iData2$status/iAtRisk[iData2$group]

            ## second term
            iCumhazard.iid[iData2$index] <- iCumhazard.iid[iData2$index] - iLambda[iData2$group]
            if(iN.obsTime<iN.obs){
                iCumhazard.iid[setdiff(iIndex.strata,iIndex.strataTime)] <- - iLambda[utils::tail(iData2$group,1)]
            }
            attr(survival,"iid")[,iS] <- - iSurv.for.iid * iCumhazard.iid

        }
    }

    ## *** export
    if(iid){
        survival$se <-  sqrt(colSums(attr(survival,"iid")^2))
        survival$lower <- survival$estimate + stats::qnorm(0.025) * survival$se
        survival$upper <- survival$estimate + stats::qnorm(0.975) * survival$se
    }
    return(survival)
}


##----------------------------------------------------------------------
### prodlimConditional.R ends here
