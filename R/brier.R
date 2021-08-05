### brier.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  5 2021 (13:44) 
## Version: 
## Last-Updated: aug  5 2021 (19:02) 
##           By: Brice Ozenne
##     Update #: 79
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){

    warper <- function(n){
        df <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
        e.logit <- glm(Y~X1+X2, data = df, family = binomial(link="logit"))
        e.perf <- performance(e.logit, trace = FALSE, transformation = FALSE, fold.number = 0)
        e.Score <- riskRegression::Score(list(e.logit), formula = Y~1, data = df)

        auc.Score <- data.frame(model = "logit1", method = "Score", metric = "auc",
                                estimate = e.Score$AUC$score$AUC, se = e.Score$AUC$score$se, lower = e.Score$AUC$score$lower, upper = e.Score$AUC$score$upper, p.value = NA)
        brier.Score <- data.frame(model = "logit1", method = "Score", metric = "brier",
                                  estimate = e.Score$Brier$score$Brier[2], se = e.Score$Brier$score$se[2], lower = e.Score$Brier$score$lower[2], upper = e.Score$Brier$score$upper[2], p.value = NA)
        e.perf <- rbind(e.perf, auc.Score, brier.Score)
        return(e.perf)
    }

    warper(100)

    n.sim <- 100
    ls.res <- pblapply(1:n.sim, function(iSim){
        rbind(cbind(sim = iSim, warper(100)),
              cbind(sim = iSim, warper(1000)))
    })
    dt.res  <- as.data.table(do.call(rbind,ls.res))
    head(dt.res)
    dt.res[,.(se.model = mean(se), se.empirical = sd(estimate)),by = c("metric","method","model")]
    dt.res[metric == "brier" & method == "internal",estimate]
}

## * brier (code)
brier <- function(labels, predictions, iid = NULL, fold = NULL, observation = NULL,
                  null = NA, conf.level = 0.95, transformation = TRUE){

    ## ** normalize user input
    if(length(unique(labels))!=2){
        stop("Argument \'labels\' must have exactly two different values. \n")
    }
    labels.factor <- as.factor(labels)
    labels.num <- as.numeric(labels.factor)
    n.obs <- length(labels)
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
        if(n.obs!=length(predictions)){
            stop("Argument \'labels\' and \'predictions\' must have the same length. \n")
        }
        if(!is.null(observation)){
            stop("Argument \'observation\' is only useful when argument \'fold\' is specified \n")
        }
        observation <- 1:n.obs
    }
    if(!is.logical(transformation)){
        stop("Argument \'transformation\' must be TRUE or FALSE \n")
    }

    ## ** prepare export
    out <- data.frame(estimate = as.numeric(NA),
                      se = as.numeric(NA),
                      lower = as.numeric(NA),
                      upper = as.numeric(NA),
                      p.value = as.numeric(NA))

    ## ** compute brier score
    if(is.null(fold)){
        iBrier <- (predictions - labels)^2
        out$estimate <- mean(iBrier)
        if(is.null(iid)){
            out$se <- sd(iBrier)/sqrt(n.obs)
        }else{
            iidAverage <- (iBrier-out$estimate)/(sqrt(n.obs)*sqrt(n.obs-1))
            ## sqrt(crossprod(iidAverage)) - sd(iBrier)/sqrt(n.obs)
            iidNuisance <-  rowMeans(sweep(iid, FUN = "*", MARGIN = 2, STATS = 2*predictions - labels))
            attr(out,"iid") <- iidAverage + iidNuisance/sqrt(n.obs)
            out$se <- sqrt(crossprod(attr(out,"iid")))
        }
    }

    ## ** compute brier score (CV)
    if(!is.null(fold)){
        Uobservation <- unique(sort(observation))
        n.Uobservation <- length(Uobservation)
        
        iBrier <- rep(0, length = n.obs)
        iFactor <- vector(mode = "list", length = n.obs)
        for(iObs in 1:n.obs){
            if(any(observation==iObs)){
                iFactor[[iObs]] <- setNames(n.Uobservation/(nFold.obs[fold[observation==iObs]]*n.fold),fold[observation==iObs])
                iBrier[iObs] <- sum((predictions[observation==iObs] - labels[iObs])^2*iFactor[[iObs]])
            }
        }
        out$estimate <- mean(iBrier[Uobservation])
        ## out$estimate - mean(tapply((predictions-labels[observation])^2,fold,mean)) ## should be equal
        if(is.null(iid)){
            out$se <- sd(iBrier[Uobservation])/sqrt(n.Uobservation)
            ## out$se - mean(tapply((predictions-labels[observation])^2,fold,sd)) ## no need to be equal
        }else{
            iidAverage <- rep(0, length = n.obs)
            iidNuisance <- rep(0, length = n.obs)
            
            iidAverage[Uobservation] <- (iBrier[Uobservation]-out$estimate)/(sqrt(n.Uobservation)*sqrt(n.Uobservation-1))
            ## sd(iBrier[Uobservation])/sqrt(n.Uobservation) - sqrt(crossprod(iidAverage)) ## should be equal
            for(iFold in 1:n.fold){ ## iFold <- 1
                iiFactor <- sapply(iFactor[observation[fold==iFold]],function(iVec){iVec[as.character(iFold)]})
                iStat <- (2*predictions[fold==iFold] - labels[observation[fold==iFold]])*iiFactor
                iidNuisance  <- iidNuisance + rowMeans(sweep(iid[,,iFold], FUN = "*", MARGIN = 2, STATS = iStat))
            }
            attr(out,"iid") <- iidAverage + iidNuisance/sqrt(n.obs)
            out$se <- sqrt(crossprod(attr(out,"iid")))
        }
    }   

    ## ** P-value and confidence interval
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

    ## ** Export
    class(out) <- append("BuyseTestBrier",class(out))
    attr(out, "contrast") <- levels(labels.factor)
    ## attr(out, "n.fold") <- n.fold
    return(out)
}

## * Utilitites
## ** print.BuyseTestBrier
print.BuyseTestBrier <- function(x, ...){
    print.data.frame(x)
}
## ** coef.BuyseTestBrier
coef.BuyseTestBrier <- function(object,...){
    object[,"estimate"]
}
## ** confint.BuyseTestBrier
confint.BuyseTestBrier <- function(object,...){
    attr(object, "iid") <- NULL
    attr(object, "constrast") <- NULL
    ## attr(object, "n.fold") <- NULL
    return(object)
}
## ** iid.BuyseTestBrier
iid.BuyseTestBrier <- function(object,...){
    return(attr(object,"iid"))
}

##----------------------------------------------------------------------
### brier.R ends here
