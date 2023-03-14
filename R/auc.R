### auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:29) 
## Version: 
## Last-Updated: mar 14 2023 (13:43) 
##           By: Brice Ozenne
##     Update #: 459
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - auc
#' @title Estimation of the Area Under the ROC Curve (EXPERIMENTAL)
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
#' @param pooling [character] method used to compute the global AUC from the fold-specific AUC: either an empirical average \code{"mean"}
#' or a weighted average with weights proportional to the number of pairs of observations in each fold \code{"pairs"}.
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
#' auc(labels = dt$Y, predictions = dt$X, direction = ">")
#'
#' ## compute auc after 10-fold cross-validation
#' auc(labels = dt$Y, prediction = dt$X, fold = dt$fold, observation = 1:NROW(dt))
#'

## * Code - auc
#' @export
auc <- function(labels, predictions, fold = NULL, observation = NULL,
                direction = ">", add.halfNeutral = TRUE,
                null = 0.5, conf.level = 0.95, transformation = TRUE, order.Hprojection = 2, pooling = "mean"){

    ## ** Normalize user imput
    if(length(unique(labels))!=2){
        stop("Argument \'labels\' must have exactly two different values \n")
    }
    if(any(is.na(predictions))){
        warning("Missing values in argument \'predictions'. \n")
    }
    if(!is.null(fold)){
        if(length(fold)!=length(predictions)){
            stop("When not NULL, argument \'fold\' must have the same length as argument \'predictions\' \n")
        }
        if(length(observation)!=length(predictions)){
            stop("When argument \'fold\' is not NULL, argument \'observation\' must have the same length as argument \'predictions\' \n")
        }
        if(length(labels)!=length(predictions)){ ##  either the user provides the outcome for each observation
            n.obs <- length(labels)
            if(!is.null(observation) && any(observation %in% 1:n.obs == FALSE)){
                stop("When not NULL, argument \'observation\' must take integer values between 1 and ",n.obs,"\n", sep = "")
            }
        }else{ ## or the user provides the outcome for each prediction (i.e. several times the same value as some predictions correspond to the same obs)
            n.obs <- length(unique(observation))
            observation <- as.numeric(as.factor(observation))
        }

        if(any(tapply(observation,fold, function(iObs){any(duplicated(iObs))}))){
            stop("The same observation cannot appear twice in the same fold. \n")
        }
    }else{
        n.obs <- length(labels)    
        if(n.obs!=length(predictions)){
            stop("Argument \'labels\' and \'predictions\' must have the same length \n")
        }
        if(!is.null(observation)){
            stop("Argument \'observation\' is only useful when argument \'fold\' is specified \n")
        }
        observation <- 1:n.obs
    }
    direction <- match.arg(direction, c(">","<","auto"), several.ok = TRUE)
    pooling <- match.arg(pooling, c("pairs","mean"), several.ok = TRUE)
    
    if(!is.logical(transformation)){
        stop("Argument \'transformation\' must be TRUE or FALSE \n")
    }
    if(length(labels)==length(predictions)){
        df <- data.frame(Y = labels,
                         X = predictions,
                         observation = observation,
                         stringsAsFactors = FALSE)
    }else{
        df <- data.frame(Y = labels[observation],
                         X = predictions,
                         observation = observation,
                         stringsAsFactors = FALSE)
    }
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

    if(is.na(conf.level)){
        method.inference <- "none"
    }else{
        method.inference <- "u-statistic"
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
    if(method.inference!="none"){
        attr(out,"iid") <- matrix(NA, nrow = n.obs, ncol = n.fold+1, dimnames = list(NULL,c(name.fold,"global")))
    }else{
        out$se <- NA
    }
    
    ## ** Global AUC
    order.save <- BuyseTest.options()$order.Hprojection
    if(order.save!=order.Hprojection){
        BuyseTest.options(order.Hprojection = order.Hprojection)
        on.exit(BuyseTest.options(order.Hprojection = order.save))
    }
    e.BT <- BuyseTest(formula, method.inference = method.inference, data = df, trace = 0, add.halfNeutral = add.halfNeutral)

    ## store
    if(is.null(fold)){
        out[out$fold=="global","estimate"] <- as.double(coef(e.BT, statistic = "favorable"))
        if(method.inference!="none"){
            out[out$fold=="global","se"] <- as.double(confint(e.BT, statistic = "favorable")[,"se"])  ## may differ from iid when second order H-decomposition
            attr(out,"iid")[sort(unique(observation)),out$fold=="global"] <- getIid(e.BT, scale = TRUE, center = TRUE, statistic = "favorable") ## no need for cluster argument when fold=NULL
        }
    }else if(pooling == "mean"){ ## Here: strata have the same weigth
        ## WARNING: cannot use the "global" results as if there is not the same number of pairs in all strata
        ##          it will weight differently the strata-specific AUCs
        if(method.inference!="none"){
            attr(out,"iid")[] <- 0
        }
    }else if(pooling  == "pairs"){ ## Here: strata are weigthed according to the number of pairs
        out[out$fold=="global","estimate"] <- as.double(coef(e.BT, statistic = "favorable"))
        if(method.inference!="none"){
            out[out$fold=="global","se"] <- as.double(confint(e.BT, cluster = observation, statistic = "favorable")[,"se"])
            attr(out,"iid")[sort(unique(observation)),out$fold=="global"] <- getIid(e.BT, cluster = observation, scale = TRUE, center = TRUE, statistic = "favorable") 
            ## sqrt(as.double(crossprod(attr(out,"iid")[,out$fold=="global"])))
        }
    }
    
    ## ** Fold-specific AUC
    if(!is.null(fold)){
        if(order.Hprojection==1){
            ePOINT.BT <- coef(e.BT, statistic = "favorable", stratified = TRUE)[,1]

            normWithinStrata <- FALSE
            attr(normWithinStrata, "skipScaleCenter") <- TRUE

            out[match(name.fold,out$fold),"estimate"] <- as.double(ePOINT.BT)

            if(method.inference!="none"){
                iIID.BT <- getIid(e.BT, scale = normWithinStrata, center = normWithinStrata, statistic = "favorable")[,1]
                out[match(name.fold,out$fold),"se"] <- sqrt(as.double(tapply(iIID.BT, fold, crossprod)))
                ## iE.BT <- BuyseTest(formula0, method.inference = "u-statistic", data = df[df$fold==name.fold[2],,drop=FALSE], trace = 0, add.halfNeutral = add.halfNeutral)
                ## confint(iE.BT, statistic = "favorable")

                for(iFold in 1:n.fold){
                    attr(out,"iid")[observation[fold==name.fold[iFold]],iFold] <- iIID.BT[fold==name.fold[iFold]]
                }
            }
            
        }else{
            for(iFold in 1:n.fold){ ## iFold <- 1
                iData <- df[df$fold==name.fold[iFold],,drop=FALSE]
                iE.BT <- BuyseTest(formula0, method.inference = method.inference, data = iData,
                                   trace = 0, add.halfNeutral = add.halfNeutral)
                iConfint <- confint(iE.BT, statistic = "favorable")
                out[match(name.fold[iFold],out$fold),"estimate"] <- as.double(iConfint$estimate)
                if(method.inference!="none"){
                    out[match(name.fold[iFold],out$fold),"se"] <- as.double(iConfint$se)
                    attr(out,"iid")[iData$observation,iFold] <- getIid(iE.BT, scale = TRUE, center = TRUE, statistic = "favorable")
                }
            }
        }

        if(pooling == "mean"){ ## same weight to each fold
            out[out$fold=="global","estimate"] <- mean(out[out$fold!="global","estimate"]) 

            if(method.inference!="none"){
                attr(out,"iid")[,"global"] <- rowMeans(attr(out,"iid")[,1:n.fold,drop=FALSE])
                out[out$fold=="global","se"] <- sqrt(as.double(crossprod(attr(out,"iid")[,"global"]))) ## may not match sum(out[out$fold!="global","se"]^2)/n.fold^2 with non-independent folds
                ## also does not have 2nd order term
            }
        }

    }

    ## ** P-value and confidence interval
    if(method.inference!="none"){
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
        }else if(all(out$estimate==0)){
            out$lower <- 0
            out$upper <- 0
            out$p.value <- as.numeric(null==0)
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
    }else{
        out$lower <- NA
        out$upper <- NA
        out$p.value <- NA
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
#' @param x object of class \code{BuyseTestAUC} (output of the \code{auc} function).
#' @param ... not used. For compatibility with the generic function.
#'
#' @return A column vector.
#' 
#' @method iid BuyseTestAuc
#' @export
iid.BuyseTestAuc <- function(x,...){
    object <- x
    return(attr(object,"iid")[,"global",drop=FALSE])
}


######################################################################
### auc.R ends here
