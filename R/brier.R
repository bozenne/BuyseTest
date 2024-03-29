### brier.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  5 2021 (13:44) 
## Version: 
## Last-Updated: jun 27 2023 (10:04) 
##           By: Brice Ozenne
##     Update #: 192
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * brier (documentation)
#' @title Estimation of the Brier Score (EXPERIMENTAL)
#' @name brier
#' 
#' @description Estimation of the brier score, possibly after cross validation,
#' to assess the discriminant ability and calibration of a biomarker regarding a disease status.
#' 
#' @param labels [integer/character vector] the disease status (should only take two different values).
#' @param predictions [numeric vector] A vector with the same length as \code{labels} containing the biomarker values.
#' @param iid [array, optional] influence function of the prediction. For cross validation (CV) should be a 3 dimensional array (one slice per CV fold).
#' Otherwise a matrix with as many column as observations and rows as predictions.
#' @param fold [character/integer vector] If using cross validation, the index of the fold. 
#' Should have the same length as \code{labels}.
#' @param observation [integer vector] If using cross validation, the index of the corresponding observation in the original dataset.
#' Necessary to compute the standard error when using cross validation.
#' @param null [numeric, 0-1] the value against which the AUC should be compared when computing the p-value.
#' @param conf.level [numeric, 0-1] the confidence level of the confidence intervals.
#' @param transformation [logical] should a log-log transformation be used when computing the confidence intervals and the p-value.
#' 
#' @keywords models
#' 
#' @return An S3 object of class \code{BuyseTestBrier} that inherits from data.frame.

## * brier (code)
#' @export
brier <- function(labels, predictions, iid = NULL, fold = NULL, observation = NULL,
                  null = NA, conf.level = 0.95, transformation = TRUE){

    ## ** normalize user input
    if(length(unique(labels))!=2){
        stop("Argument \'labels\' must have exactly two different values. \n")
    }
    if(is.numeric(labels)){
        labels.num <- labels
        labels.factor <- as.factor(labels)
    }else if(is.factor(labels)){
        labels.num <- as.numeric(labels)-1
        labels.factor <- labels
    }else{
        stop("Argument \'labels\' must a numeric or factor vector. \n")
    }
    n.obs <- length(labels)
    if(identical(fold,0)){fold  <- NULL}
    if(!is.null(iid) && is.null(fold)){
        if(NCOL(iid) != n.obs){
            stop("Argument \'iid\' must have one column per observation. \n")
        }
    }
    if(!is.null(fold)){
        Ufold <- unique(fold)
        n.fold <- length(Ufold)
        nFold.obs <- table(fold)
        if(length(predictions)!=length(fold)){
            stop("When not NULL, argument \'fold\' must have the same length as argument \'predictions\' \n")
        }
        if(length(observation)!=length(predictions)){
            stop("When argument \'fold\' is not NULL, argument \'observation\' must have the same length as argument \'predictions\' \n")
        }
        if(!is.null(observation) && any(observation %in% 1:n.obs == FALSE)){
            stop("When not NULL, argument \'observation\' must take integer values between 1 and ",n.obs,"\n", sep = "")
        }
        if(any(sapply(tapply(observation,fold,duplicated),any))==TRUE){
            stop("Cannot quantify uncertainty when the same observation appear several times in the same fold. \n")
        }
        if(!is.null(iid)){
            if(n.fold!=dim(iid)[3]){
                stop("The third dimension of argument \'iid\' should equal the number of folds. \n")
            }
            UnFold.obs <- unique(nFold.obs)
            if(length(UnFold.obs)!=1){
                stop("The number of observations should be the same in each fold. \n")
            }
            if(UnFold.obs!=dim(iid)[2]){
                stop("The second dimension of argument \'iid\' should equal the number of observations per fold. \n")
            }
        }
    }else{
        external <- FALSE
        if(n.obs!=length(predictions)){
            stop("Argument \'labels\' and \'predictions\' must have the same length. \n")
        }
        if(identical(observation,"external")){
            observation <- NULL
            external <- TRUE
        }else if(!is.null(observation)){
            stop("Argument \'observation\' is only useful when argument \'fold\' is specified \n")
        }
        observation <- 1:n.obs
    }
    if(!is.logical(transformation)){
        stop("Argument \'transformation\' must be TRUE or FALSE \n")
    }
    if(is.na(conf.level)){
        se <- FALSE
    }else{
        se <- TRUE
    }
    
    ## ** prepare export
    if(!is.null(fold)){
        name.fold <- as.character(sort(unique(fold)))
        n.fold <- length(name.fold)
        out <- data.frame(fold = c(name.fold,"global"),
                          estimate = as.numeric(NA),
                          se = as.numeric(NA),
                          lower = as.numeric(NA),
                          upper = as.numeric(NA),
                          p.value = as.numeric(NA))
    }else{
        out <- data.frame(fold = "global",
                          estimate = as.numeric(NA),
                          se = as.numeric(NA),
                          lower = as.numeric(NA),
                          upper = as.numeric(NA),
                          p.value = as.numeric(NA))
    }

    ## ** compute brier score
    if(is.null(fold)){
        iBrier <- (predictions - labels.num)^2
        out$estimate <- mean(iBrier)
        if(se){
            iidAverage <- (iBrier-out$estimate)/(sqrt(n.obs)*sqrt(n.obs-1))
            if(is.null(iid)){
                attr(out,"iid") <- iidAverage
                ## se: stats::sd(iBrier)/sqrt(n.obs)
            }else{
                ## sqrt(crossprod(iidAverage)) - stats::sd(iBrier)/sqrt(n.obs)
                iidNuisance <-  rowMeans(.rowMultiply_cpp(iid, 2*predictions - labels.num))
                if(external){
                    browser()
                    attr(out,"iid") <- c(iidNuisance, iidAverage)
                }else{                    
                    attr(out,"iid") <- iidAverage + iidNuisance
                }
            }
            out$se <- sqrt(crossprod(attr(out,"iid")))
        }
    }

    ## ** compute brier score (CV)
    if(!is.null(fold)){
        Uobservation <- unique(sort(observation))
        n.Uobservation <- length(Uobservation)
        
        iBrier <- rep(0, length = n.obs)
        iFactor <- vector(mode = "list", length = n.obs)
        for(iObs in 1:n.obs){ ## iObs <- 1

            if(any(observation==iObs)){
                iFactor[[iObs]] <- setNames(n.Uobservation/(nFold.obs[as.character(fold[observation==iObs])]*n.fold),fold[observation==iObs])
                iBrier[iObs] <- sum((predictions[observation==iObs] - labels.num[iObs])^2*iFactor[[iObs]])
            }
        }

        out$estimate[match(name.fold,out$fold)] <- tapply((predictions-labels.num[observation])^2, fold, mean)[name.fold]
        out$estimate[out$fold=="global"] <- mean(iBrier[Uobservation])
        ## mean(iBrier[Uobservation]) - mean(out$estimate[1:10])
        if(se){
            if(is.null(iid)){
                out$se[match(name.fold,out$fold)] <- tapply((predictions-labels.num[observation])^2, fold, function(iDiff){sqrt(stats::var(iDiff)/length(iDiff))})[name.fold]
                out$se[out$fold=="global"] <- stats::sd(iBrier[Uobservation])/sqrt(n.Uobservation)
                ## out$se - mean(tapply((predictions-labels[observation])^2,fold,sd)) ## no need to be equal
                attr(out,"iid") <- rep(0, length = n.obs)
                attr(out,"iid")[Uobservation] <- (iBrier[Uobservation]-out$estimate[out$fold=="global"])/(sqrt(n.Uobservation)*sqrt(n.Uobservation-1))
                ## out$se[out$fold=="global"] - sqrt(crossprod(attr(out,"iid")))
            }else{
                iidAverage <- rep(0, length = n.obs)
                iidNuisance <- rep(0, length = n.obs)
            
                iidAverage[Uobservation] <- (iBrier[Uobservation]-out$estimate[out$fold=="global"])/(sqrt(n.Uobservation)*sqrt(n.Uobservation-1))
                ## stats::sd(iBrier[Uobservation])/sqrt(n.Uobservation) - sqrt(crossprod(iidAverage)) ## should be equal
                for(iFold in 1:n.fold){ ## iFold <- 1
                    iiFactor <- sapply(iFactor[observation[fold==name.fold[iFold]]],function(iVec){iVec[name.fold[iFold]]})
                    iStat <- 2*(predictions[fold==name.fold[iFold]] - labels.num[observation[fold==name.fold[iFold]]])
                    iidNuisance  <- iidNuisance + rowMeans(.rowMultiply_cpp(iid[,,iFold], iStat*iiFactor))

                    ## in each fold because of CV the training and test set are separate so the uncertainties are independent
                    term1 <- stats::sd((predictions[fold==name.fold[iFold]] - labels.num[observation[fold==name.fold[iFold]]])^2)
                    term2 <- sqrt(crossprod(rowMeans(.rowMultiply_cpp(iid[,,iFold], iStat)))/sum(fold==name.fold[iFold]))
                    out[out$fold==name.fold[iFold],"se"] <- term1 + term2
                }
                attr(out,"iid") <- iidAverage + iidNuisance/sqrt(n.obs)
                out$se[out$fold=="global"] <- sqrt(crossprod(attr(out,"iid")))
            }
        }
    }   

    ## ** P-value and confidence interval
    if(se){
        alpha <- 1-conf.level
        qinf <- stats::qnorm(alpha/2)
        qsup <- stats::qnorm(1-alpha/2)
    
        if(transformation){
            newse <- out$se / out$estimate
            out$lower <- as.double(out$estimate * exp(qinf * newse))
            out$upper <- as.double(out$estimate * exp(qsup * newse))

            if(!is.na(null)){
                z.stat <- (log(out$estimate) - log(null))/newse
                out$p.value <- 2*(1-stats::pnorm(abs(z.stat)))
            }
        }else{
    
            out$lower <- as.double(out[,"estimate"] + qinf * out[,"se"])
            out$upper <- as.double(out[,"estimate"] + qsup * out[,"se"])

            if(!is.na(null)){
                z.stat <- as.double((out[,"estimate"]-null)/out[,"se"])
                out$p.value <- 2*(1-stats::pnorm(abs(z.stat)))
            }
        }
    }

    ## ** Export
    class(out) <- append("BuyseTestBrier",class(out))
    attr(out, "contrast") <- levels(labels.factor)
    ## attr(out, "n.fold") <- n.fold
    return(out)
}

## * Utilitites
## ** print.BuyseTestBrier
##' @exportMethod print
print.BuyseTestBrier <- function(x, ...){
    print.data.frame(x)
}
## ** coef.BuyseTestBrier
##' @title Extract the Brier Score
##'
##' @description Extract the Brier score.
##' 
##' @param object object of class \code{BuyseTestBrier} (output of the \code{brier} function).
##' @param ... not used. For compatibility with the generic function.
##'
##' @return Estimated value for Brier score (numeric).  
##' 
##' @method coef BuyseTestBrier
##' @keywords methods
##' 
##' @export
coef.BuyseTestBrier <- function(object,...){
    object[,"estimate"]
}
## ** confint.BuyseTestBrier
##' @title Extract the Brier Score with its Confidence Interval
##'
##' @description Extract the Brier score with its Confidence Interval and possibly a p-value.
##' 
##' @param object object of class \code{BuyseTestBrier} (output of the \code{brier} function).
##' @param ... not used. For compatibility with the generic function.
##'
##' @return Estimated value for the brier score, its standard error, the lower and upper bound of the confidence interval and the p-value.
##' 
##' @method confint BuyseTestBrier
##' @keywords methods
##' 
##' @export
confint.BuyseTestBrier <- function(object,...){
    out <- object[object$fold=="global",c("estimate","se","lower","upper","p.value")]
    rownames(out) <- NULL
    return(out)
}
## ** iid.BuyseTestBrier
##' @title Extract the idd Decomposition for the Brier Score
##'
##' @description Extract the iid decompotion relative to Brier score estimate.
##' 
##' @param x object of class \code{BuyseTestBrier} (output of the \code{brier} function).
##' @param ... not used. For compatibility with the generic function.
##'
##' @return A column vector.
##' 
##' @method iid BuyseTestBrier
##' @keywords methods
##' 
##' @export
iid.BuyseTestBrier <- function(x,...){
    object <- x
    return(attr(object,"iid"))
}

##----------------------------------------------------------------------
### brier.R ends here
