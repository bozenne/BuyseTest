### efronlim.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 20 2025 (10:53) 
## Version: 
## Last-Updated: feb 20 2025 (13:10) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


##' @title Constrained Kaplan-Meier Estimator
##' @description Kaplan-Meier estimator, possibly stratified, constrained to 0 just after end of follow-up in each strata.
##' 
##' @param formula formula with on the right-hand side Hist(time,event) or Surv(time,event) and 1 or stratification variable(s) on the right-hand side.
##' @param data A \code{data.frame} containing the variables of present in argument \code{formula}.
##' @param ... additional arguments passed to \code{prodlim::prodlim}
##'
##' @details If in any strata there is censoring at the observed time, then the dataset is updated
##' by setting one of the censored observations to an infinitesimal later timepoint with an event instead of censoring.
##' Then the possibly stratified Kaplan-Meier estimator is run on this updated dataset.
##' 
##' @return A \code{prodlim} object.
##'
##' @examples
##' library(data.table)
##' 
##' set.seed(1)
##' dt.data <- simBuyseTest(1e2, n.strata = 2)
##' dt.data$time1 <- pmin(dt.data$eventtime, 1)
##' dt.data$status1 <- dt.data$status * (dt.data$eventtime<=1)
##'
##' ## KM
##' if(require(prodlim)){
##'    e.KM <- prodlim(Hist(time1,status1)~1, data = dt.data)
##'    plot(e.KM)
##' }
##' 
##' e.KM0 <- efronlim(Hist(time1,status1)~1, data = dt.data)
##' plot(e.KM0)
##' 
##' ## stratfied KM
##' if(require(prodlim)){
##'    e.KMS <- prodlim(Hist(time1,status1)~strata, data = dt.data)
##'    plot(e.KMS)
##' }
##' 
##' e.KMS <- efronlim(Hist(time1,status1)~strata, data = dt.data)
##' plot(e.KMS)

##' @export
efronlim <- function(formula, data, ...){

    zeroPlus <- 1e-12

    ## ** check arguments
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' must be or inherit from \"data.frame\". \n",
             "Current class: \"",paste(class(data),collapse = "\", \""),"\".\n")
    }
    
    regressors <- all.vars(stats::delete.response(terms(formula)))
    if(length(regressors)>0 && any(regressors %in% names(data) == FALSE)){
        stop("Argument \'formula\' is incompatible with argument \'data\'. \n",
             "Missing variable in data: \"",paste(setdiff(regressors,names(data)), collapse = "\", \""),"\".\n")
    }
    outcome <- setdiff(all.vars(formula),regressors)
    dataOutcome <- stats::model.frame(stats::update(formula, ".~1"), data = data)

    if(all(abs(dataOutcome$Hist[,1] - data[[outcome[1]]]) < 1e-12)){
        outcomeTime <- outcome[1]
        outcomeStatus <- outcome[2]
    }else if(all(abs(dataOutcome$Hist[,1] - data[[outcome[2]]]) < 1e-12)){
        outcomeTime <- outcome[2]
        outcomeStatus <- outcome[1]
        
    }else{
        stop("Something went wrong when identifying the time and the status variable. \n")
    }
    
    ## ** check need for efron constrain
    if(length(regressors)==0){
        ls.indexStrata <- list(1:NROW(data))
    }else{
        levelStrata <- interaction(lapply(regressors, function(iR){as.numeric(as.factor(data[[iR]]))}))
        ls.indexStrata <- tapply(1:NROW(data), INDEX = levelStrata, FUN = identity, simplify = FALSE)
    }
    indexLastStrata <- sapply(ls.indexStrata, FUN = function(iIndex){
        iIndex[max(which(data[[outcomeTime]][iIndex] == max(data[[outcomeTime]][iIndex])))]
    })
    statusLastStrata  <- data[[outcomeStatus]][indexLastStrata]
    if(any(statusLastStrata==0)){
        data.prodlim <- as.data.frame(data)
        data.prodlim[[outcomeStatus]][indexLastStrata[statusLastStrata==0]] <- 1
        data.prodlim[[outcomeTime]][indexLastStrata[statusLastStrata==0]] <- data.prodlim[[outcomeTime]][indexLastStrata[statusLastStrata==0]] + zeroPlus/10
    }else{
        data.prodlim <- data
    }

    ## ** fit prodlim
    out <- prodlim::prodlim(formula, data.prodlim, ...)
    out$efron <- any(statusLastStrata==0)
    if(out$model != "survival"){
        stop("Efron constrain can only be applied in survival setting. \n",
             "Current setting: \"",out$model,"\". \n")
    }

    ## ** export
    return(out)

}

##----------------------------------------------------------------------
### efronlim.R ends here
