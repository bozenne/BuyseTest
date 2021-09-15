### auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:29) 
## Version: 
## Last-Updated: sep 15 2021 (14:58) 
##           By: Brice Ozenne
##     Update #: 380
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - auc
#' @title Estimation of the Area Under the ROC Curve
#' @name auc
#' 
#' @description Estimation of the Area Under the ROC curve, possibly after cross validation,
#' to assess the discriminant ability of a biomarker regarding a disease status.
#' 
#' @param labels [integer/character vector] the disease status (should only take two different values).
#' @param predictions [numeric vector] A vector with the same length as \code{labels} containing the biomarker values.
#' @param fold [character/integer vector] If using cross validation, the index of the fold. 
#' Should have the same length as \code{labels}.
#' @param observation [integer vector] If using cross validation, the index of the corresponding observation in the original dataset.
#' Necessary to compute the standard error when using cross validation.
#' @param direction [character] \code{">"} lead to estimate P[Y>X],
#' \code{"<"} to estimate P[Y<X],
#' and \code{"auto"} to estimate max(P[Y>X],P[Y<X]).
#' @param add.halfNeutral [logical] should half of the neutral score be added to the favorable and unfavorable scores?
#' Useful to match the usual definition of the AUC in presence of ties.
#' @param null [numeric, 0-1] the value against which the AUC should be compared when computing the p-value.
#' @param conf.level [numeric, 0-1] the confidence level of the confidence intervals.
#' @param transformation [logical] should a log-log transformation be used when computing the confidence intervals and the p-value.
#' @param order.Hprojection [1,2] the order of the H-projection used to linear the statistic when computing the standard error.
#' 2 is involves more calculations but is more accurate in small samples. Only active when the \code{fold} argument is \code{NULL}.
#' 
#' @details The iid decomposition of the AUC is based on a first order decomposition.
#' So its squared value will not exactly match the square of the standard error estimated with a second order H-projection.
#'
#' @return A \emph{data.frame} containing for each fold the AUC value with its standard error (when computed).
#' The last line of the data.frame contains the global AUC value with its standard error.
#'
#' @references Erin LeDell, Maya Petersen, and Mark van der Laan (2015). \bold{Computationally efficient confidence intervals for cross-validated area under the ROC curve estimates}. \emph{Electron J Stat.} 9(1):1583–1607. \cr

## * Example - auc
#' @rdname auc
#' @examples
#' library(data.table)
#'
#' n <- 200
#' set.seed(10)
#' X <- rnorm(n)
#' dt <- data.table(Y = as.factor(rbinom(n, size = 1, prob = 1/(1+exp(1/2-X)))),
#'                  X = X,
#'                  fold = unlist(lapply(1:10,function(iL){rep(iL,n/10)})))
#'
#' ## compute auc
#' BuyseTest::auc(labels = dt$Y, predictions = dt$X, direction = ">")
#'
#' ## compute auc after 10-fold cross-validation
#' BuyseTest::auc(labels = dt$Y, prediction = dt$X, fold = dt$fold, observation = 1:NROW(dt))
#'

## * Code - auc
#' @export
auc <- function(labels, predictions, fold = NULL, observation = NULL,
                direction = ">", add.halfNeutral = TRUE,
                null = 0.5, conf.level = 0.95, transformation = TRUE, order.Hprojection = 2){

    ## ** Normalize user imput
    if(length(unique(labels))!=2){
        stop("Argument \'labels\' must have exactly two different values \n")
    }
    n.obs <- length(labels)    
    if(!is.null(fold)){
        if(length(fold)!=length(predictions)){
            stop("When not NULL, argument \'fold\' must have the same length as argument \'predictions\' \n")
        }
        if(length(observation)!=length(predictions)){
            stop("When argument \'fold\' is not NULL, argument \'observation\' must have the same length as argument \'predictions\' \n")
        }
        if(!is.null(observation) && any(observation %in% 1:n.obs == FALSE)){
            stop("When not NULL, argument \'observation\' must take integer values between 1 and ",n.obs,"\n", sep = "")
        }

        if(any(tapply(observation,fold, function(iObs){any(duplicated(iObs))}))){
            stop("The same observation cannot appear twice in the same fold. \n")
        }
    }else{
        if(n.obs!=length(predictions)){
            stop("Argument \'labels\' and \'predictions\' must have the same length \n")
        }
        if(!is.null(observation)){
            stop("Argument \'observation\' is only useful when argument \'fold\' is specified \n")
        }
        observation <- 1:n.obs
    }
    direction <- match.arg(direction, c(">","<","auto"), several.ok = TRUE)
    
    if(!is.logical(transformation)){
        stop("Argument \'transformation\' must be TRUE or FALSE \n")
    }
    
    df <- data.frame(Y = labels[observation],
                     X = predictions,
                     observation = observation,
                     stringsAsFactors = FALSE)
    formula0 <- Y ~ cont(X)

    if(!is.null(fold)){
        df$fold <- fold
        formula <- stats::update(formula0,.~.+fold)
        name.fold <- sort(unique(df$fold))
        n.fold <- length(name.fold)
    }else{
        df$fold <- 1
        formula <- formula0
        name.fold <- NULL
        n.fold <- 0
    }

    if(!identical(direction,"auto") && (length(direction) %in% c(1,max(n.fold,1)) == FALSE)){
        stop("Argument \'direction\' must have length 1 or the number of folds (here ",n.fold,"). \n")
    }
    
    ## ** Prepare 
    ## *** Make sure that all prediction are in the increasing means outcome direction
    direction.save <- direction
    if(direction == "auto"){
        if(sum(e.BT@count.favorable)>=sum(e.BT@count.unfavorable)){
            direction <- rep(">",max(n.fold,1))
        }else{
            direction <- rep("<",max(n.fold,1))
            df$X <- -df$X
        }
        Udirection <- as.character(NA)
    }else if(length(direction)==1){
        if(direction=="<"){
            df$X <- -df$X
        }
        direction <- rep(direction, max(n.fold,1))
        Udirection <- direction.save
    }else{
        for(iFold in 1:n.fold){
            if(direction[iFold]=="<"){
                df[df$fold==name.fold[iFold],"X"] <- -df[df$fold==name.fold[iFold],"X"]
            }
        }
        Udirection <- as.character(NA)
    }

    if(is.null(fold)){
        out <- data.frame(fold = "global",
                          direction = direction,
                          estimate = 0,
                          se = 0,
                          stringsAsFactors = FALSE)
    }else{
        out <- data.frame(fold = c(name.fold,"global"),
                          direction = c(direction,Udirection),
                          estimate = 0,
                          se = 0,
                          stringsAsFactors = FALSE)
    }
    attr(out,"iid") <- matrix(NA, nrow = n.obs, ncol = n.fold+1, dimnames = list(NULL,c(name.fold,"global")))

    ## ** Global AUC
    ## WARNING: cannot use the "global" results as if there is not the same number of pairs in all strata
    ##          it will weight differently the strata-specific AUCs
    if(is.null(fold)){
        order.save <- BuyseTest.options()$order.Hprojection
        BuyseTest.options(order.Hprojection = order.Hprojection)
    }

    e.BT <- BuyseTest(formula, method.inference = "u-statistic", data = df, trace = 0)

    if(is.null(fold)){
        BuyseTest.options(order.Hprojection = order.save)
    }

    ## store
    if(is.null(fold)){
        out[out$fold=="global","estimate"] <- as.double(coef(e.BT, statistic = "favorable", add.halfNeutral = add.halfNeutral))
        attr(out,"iid")[sort(unique(observation)),out$fold=="global"] <- getIid(e.BT, statistic = "favorable", add.halfNeutral = add.halfNeutral) ## no need for cluster argument when fold=NULL
        if(add.halfNeutral && any(abs(e.BT@iidAverage$neutral[,1])>1e-10) || any(duplicated(observation))){
            out[out$fold=="global","se"] <- as.double(sqrt(crossprod(attr(out,"iid")[,out$fold=="global"]))) ## no second order term
        }else{
            out[out$fold=="global","se"] <- as.double(confint(e.BT, statistic = "favorable")[,"se"]) ## no contribution of the neutral pairs
        }
    }else{
        attr(out,"iid")[] <- 0
    }
    
    ## ** Fold-specific AUC
    if(!is.null(fold)){
        ePOINT.BT <- coef(e.BT, statistic = "favorable", stratified = TRUE, add.halfNeutral = add.halfNeutral)[,1]
        iIID.BT <- getIid(e.BT, normalize = FALSE, statistic = "favorable", add.halfNeutral = add.halfNeutral)

        indexC <- attr(e.BT@level.treatment,"indexC")
        n.C <- length(indexC)
        indexT <- attr(e.BT@level.treatment,"indexT")
        n.T <- length(indexT)
        n.obs <- n.C + n.T

        ## *** loop over folds
        for(iFold in 1:n.fold){ ## iFold <- 1

            iIndex <- (1:n.obs)[df$fold == name.fold[iFold]]
            iIndexC <- intersect(iIndex,indexC)
            iIndexT <- intersect(iIndex,indexT)
            iObs <- observation[iIndex]
            iObsC <- observation[iIndexC]
            iObsT <- observation[iIndexT]

            out[out$fold==name.fold[iFold],"estimate"] <- as.double(ePOINT.BT[iFold])
                
            ## fold-specific iid
            ## center around the expected value and scale by the number of observations in the group of the fold
            attr(out,"iid")[iObsC,out$fold==name.fold[iFold]] <- (iIID.BT[iIndexC,1] - ePOINT.BT[iFold])/length(iIndexC)
            attr(out,"iid")[iObsT,out$fold==name.fold[iFold]] <- (iIID.BT[iIndexT,1] - ePOINT.BT[iFold])/length(iIndexT)
            out[out$fold==name.fold[iFold],"se"] <- as.double(sqrt(crossprod(attr(out,"iid")[,out$fold==name.fold[iFold]]))) ## no second order term
            ## iE.BT <- BuyseTest(formula0, method.inference = "u-statistic", data = df[df$fold==name.fold[iFold],,drop=FALSE], trace = 0)
            ## confint(iE.BT, statistic = "favorable")
            ## sqrt(crossprod(getIid(iE.BT, statistic = "favorable")))

            ## global iid
            ## center around the expected value and scale by the number of observations in the group
            attr(out,"iid")[iObsC,out$fold=="global"] <- attr(out,"iid")[iObsC,out$fold=="global"] + (iIID.BT[iIndexC,1] - ePOINT.BT[iFold])/n.C
            attr(out,"iid")[iObsT,out$fold=="global"] <- attr(out,"iid")[iObsT,out$fold=="global"] + (iIID.BT[iIndexT,1] - ePOINT.BT[iFold])/n.T
        }

        out[out$fold=="global","estimate"] <- mean(out[out$fold!="global","estimate"])
        out[out$fold=="global","se"] <- as.double(sqrt(crossprod(attr(out,"iid")[,out$fold=="global"]))) ## no second order term
    }

    ## ** P-value and confidence interval
    alpha <- 1-conf.level
    qinf <- stats::qnorm(alpha/2)
    qsup <- stats::qnorm(1-alpha/2)

    ## riskRegression:::transformCIBP(estimate = cbind(out$estimate), se = cbind(out$se), null = 1/2, conf.level =  0.95, type = "none",
    ## ci = TRUE, band = FALSE, p.value = TRUE,
    ## min.value = 0, max.value = 1)
    ## riskRegression:::transformCIBP(estimate = cbind(out$estimate), se = cbind(out$se), null = 1/2, conf.level =  0.95, type = "loglog",
    ## ci = TRUE, band = FALSE, p.value = TRUE,
    ## min.value = 0, max.value = 1)
    if(all(out$estimate==1)){
        out$lower <- 1
        out$upper <- 1
        out$p.value <- as.numeric(null==1)
    }else if(transformation){
        newse <- out$se / (- out$estimate * log(out$estimate))
        z.stat <- (log(-log(out$estimate)) - log(-log(null)))/newse
        
        out$lower <- as.double(out$estimate ^ exp(qsup * newse))
        out$upper <- as.double(out$estimate ^ exp(qinf * newse))
        out$p.value <- 2*(1-stats::pnorm(abs(z.stat)))
    }else{
        z.stat <- as.double((out[,"estimate"]-null)/out[,"se"])

        out$lower <- as.double(out[,"estimate"] + qinf * out[,"se"])
        out$upper <- as.double(out[,"estimate"] + qsup * out[,"se"])
        out$p.value <- 2*(1-stats::pnorm(abs(z.stat)))
    }

    ## ** Export
    attr(out, "n.fold") <- n.fold
    class(out) <- append("BuyseTestAuc",class(out))
    attr(out, "contrast") <- e.BT@level.treatment
    return(out)
 
}

## * Utilitites
## ** print.auc
#' @export
print.BuyseTestAuc <- function(x, ...){
    ##    if(attr(x,"n.fold")==0){
        print.data.frame(x, ...)
    ## }else{
    ##     label.upper <- paste0(attr(x,"contrast")[2],">",attr(x,"contrast")[1])
    ##     label.lower <- paste0(attr(x,"contrast")[1],">",attr(x,"contrast")[2])
    ##     x$direction <- sapply(x$direction, function(iD){
    ##         if(iD==">"){return(label.upper)}else if(iD=="<"){return(label.lower)}else{return(iD)}
    ##     })
    ##     print.data.frame(x[x$fold == "global",c("direction","estimate","se","lower","upper","p.value")], row.names = FALSE)
    ## }
}

## ** coef.auc
#' @title Extract the AUC Value
#'
#' @description Extract the AUC value.
#' 
#' @param object object of class \code{BuyseTestAUC} (output of the \code{auc} function).
#' @param ... not used. For compatibility with the generic function.
#'
#' @return Estimated value for the AUC (numeric).  
#' 
#' @method coef BuyseTestAuc
#' 
#' @export
coef.BuyseTestAuc <- function(object,...){
    object[object$fold=="global","estimate"]
}
## ** confint.auc
#' @title Extract the AUC value with its Confidence Interval
#'
#' @description Extract the AUC value with its Confidence Interval and p-value testing whether the AUC equals 0.5.
#' 
#' @param object object of class \code{BuyseTestAUC} (output of the \code{auc} function).
#' @param ... not used. For compatibility with the generic function.
#'
#' @return Estimated value for the AUC, its standard error, the lower and upper bound of the confidence interval and the p-value.
#' 
#' @method confint BuyseTestAuc
#' @export
confint.BuyseTestAuc <- function(object,...){
    out <- object[object$fold=="global",c("estimate","se","lower","upper","p.value")]
    rownames(out) <- NULL
    return(out)
}
## ** iid.auc
#' @title Extract the idd Decomposition for the AUC
#'
#' @description Extract the iid decompotion relative to AUC estimate.
#' 
#' @param object object of class \code{BuyseTestAUC} (output of the \code{auc} function).
#' @param ... not used. For compatibility with the generic function.
#'
#' @return A column vector.
#' 
#' @method iid BuyseTestAuc
#' @export
iid.BuyseTestAuc <- function(object,...){
    return(attr(object,"iid")[,"global",drop=FALSE])
}


######################################################################
### auc.R ends here
