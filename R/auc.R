### auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:29) 
## Version: 
## Last-Updated: apr  6 2020 (11:06) 
##           By: Brice Ozenne
##     Update #: 166
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
#' @param null [numeric, 0-1] the value against which the AUC should be compared when computing the p-value.
#' @param conf.level [numeric, 0-1] the confidence level of the confidence intervals.
#' @param transformation [logical] should a log-log transformation be used when computing the confidence intervals and the p-value.
#' 
#' @details Compared to other functions computing the AUC (e.g. the auc fonction of the ROCR package),
#' the AUC is defined here as P[Y>X] with a strict inequality sign (i.e. not P[Y>=X]).
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
#' @rdname auc
#' @export
auc <- function(labels, predictions, fold = NULL, observation = NULL, direction = ">",
                null = 0.5, conf.level = 0.95, transformation = FALSE){

    ## ** Normalize user imput
    if(length(unique(labels))!=2){
        stop("Argument \'labels\' must have exactly two different values \n")
    }
    n.obs <- length(labels)
    if(n.obs!=length(predictions)){
        stop("Argument \'labels\' and \'predictions\' must have the same length \n")
    }
    if(!is.null(fold)){
        if(n.obs!=length(fold)){
            stop("When not NULL, argument \'fold\' must have the same length as argument  \'labels\' \n")
        }
        if(!is.null(observation)){
            if(n.obs!=length(fold)){
                stop("When not NULL, argument \'observation\' must have the same length as argument  \'labels\' \n")
            }
            if(!is.null(observation) && any(observation %in% 1:n.obs == FALSE)){
                stop("When not NULL, argument \'observation\' must take integer values between 1 and ",n.obs,"\n", sep = "")
            }
        }
    }else{
        if(!is.null(observation)){
            stop("Argument \'observation\' is only useful when argument \'fold\' is specified \n")
        }
        observation <- 1:n.obs
    }
    direction <- match.arg(direction, c(">","<","auto","best"))
    if(!is.logical(transformation)){
        stop("Argument \'transformation\' must be TRUE or FALSE \n")
    }
    
    df <- data.frame(Y = labels,
                     X = predictions,
                     stringsAsFactors = FALSE)
    formula <- Y ~ cont(X)

    if(!is.null(fold)){
        df$fold <- fold
        formula <- stats::update(formula,.~.+fold)
    }else{
        df$fold <- 1
    }

    ## ** Perform GPC
    order.save <- BuyseTest.options()$order.Hprojection

    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(formula, method.inference = "u-statistic", data = df, trace = 0)
    BuyseTest.options(order.Hprojection = order.save)
    indexC <- attr(e.BT@level.treatment,"indexC")
    n.C <- length(indexC)
    indexT <- attr(e.BT@level.treatment,"indexT")
    n.T <- length(indexT)
    
    ## ** Extra AUC
    direction.save <- direction
    if(direction == "auto"){
        if(sum(e.BT@count.favorable)>=sum(e.BT@count.unfavorable)){
            direction <- ">"
        }else{
            direction <- "<"
        }
    }

    name.fold <- e.BT@level.strata
    n.fold <- length(name.fold)

    if(direction==">"){
        out <- data.frame(fold = c(name.fold,"global"),
                          direction = ">",
                          estimate = c(e.BT@count.favorable/e.BT@n.pairs,NA),
                          se = NA,
                          stringsAsFactors = FALSE)
        if(!is.null(observation)){
            M.iid <- do.call(cbind,tapply(observation, df$fold, function(iVec){
                iIID <- vector(mode = "numeric", length = n.obs)
                iIID[iVec] <- scale(getIid(e.BT, normalize = FALSE)[iVec,"favorable",drop=FALSE], center = TRUE, scale = FALSE)
                iIID[intersect(iVec,indexC)] <- iIID[intersect(iVec,indexC)]/n.C
                iIID[intersect(iVec,indexT)] <- iIID[intersect(iVec,indexT)]/n.T
                return(iIID)
            }))
        }
        ## colSums(M.iid^2)
    }else if(direction == "<"){
        out <- data.frame(fold = c(name.fold,"global"),
                          direction = "<",
                          estimate = c(e.BT@count.unfavorable/e.BT@n.pairs, NA),
                          se = NA,
                          stringsAsFactors = FALSE)
        if(!is.null(observation)){
            M.iid <- do.call(cbind,tapply(observation, df$fold, function(iVec){
                iIID <- vector(mode = "numeric", length = n.obs)
                iIID[iVec] <- scale(getIid(e.BT, normalize = FALSE)[iVec,"unfavorable",drop=FALSE], center = TRUE, scale = FALSE)
                iIID[intersect(iVec,indexC)] <- iIID[intersect(iVec,indexC)]/n.C
                iIID[intersect(iVec,indexT)] <- iIID[intersect(iVec,indexT)]/n.T
                return(iIID)
            }))
        }
    }else if(direction == "best"){
        out <- data.frame(fold = c(name.fold,"global"),
                          direction = "best",
                          estimate = 0,
                          se = 0,
                          stringsAsFactors = FALSE)
        if(!is.null(observation)){
            M.iid <- matrix(0, nrow = n.obs, ncol = n.fold)
        }
        for(iFold in 1:n.fold){ ## iFold <- 1
            iDirection <- c("favorable", "unfavorable")[which.max(c(e.BT@count.favorable[iFold],e.BT@count.unfavorable[iFold]))][1]
            iCount <- as.double(slot(e.BT, paste0("count.",iDirection))[iFold])
            iVec <- observation[fold == name.fold[iFold]]
            
            out$direction[iFold] <- if(iDirection=="favorable"){">"}else{"<"}
            out$estimate[iFold] <- iCount/e.BT@n.pairs[iFold]

            if(!is.null(observation)){
                M.iid[iVec,iFold] <- scale(getIid(e.BT, normalize = FALSE)[iVec,iDirection,drop=FALSE], center = TRUE, scale = FALSE)
                M.iid[intersect(iVec,indexC),iFold] <- M.iid[intersect(iVec,indexC),iFold]/n.C
                M.iid[intersect(iVec,indexT),iFold] <- M.iid[intersect(iVec,indexT),iFold]/n.T
            }
        }
    }
    out$estimate[n.fold+1] <- mean(out$estimate[1:n.fold])
    
    if(!is.null(observation)){
        if(!is.null(fold)){
            n.obsfold <- tapply(observation,fold,function(vecObs){length(unique(vecObs))})
        }else{
            n.obsfold <- n.obs
        }
        out$se[1:n.fold] <- sqrt(colSums(M.iid^2))*(n.obs/n.obsfold)
        out$se[n.fold+1] <- sqrt(sum(rowSums(M.iid)^2))
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
    if(transformation){
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
    class(out) <- append("BuyseTestAuc",class(out))
    attr(out, "contrast") <- e.BT@level.treatment
    attr(out, "n.fold") <- n.fold
    attr(out, "iid") <- M.iid
    return(out)
 
}

## * Utilitites
## ** print.auc
#' @export
print.BuyseTestAuc <- function(x, ...){
    if(NROW(x) < attr(x,"n.fold")){
        print.data.frame(x)
    }else{
        label.upper <- paste0(attr(x,"contrast")[2],">",attr(x,"contrast")[1])
        label.lower <- paste0(attr(x,"contrast")[1],">",attr(x,"contrast")[2])
        x$direction <- sapply(x$direction, function(iD){
            if(iD==">"){return(label.upper)}else if(iD=="<"){return(label.lower)}else{return(iD)}
        })
        print.data.frame(x[x$fold == "global",c("direction","estimate","se","lower","upper","p.value")], row.names = FALSE)
    }
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
    return(as.data.frame(out, stringsAsFactors = FALSE))
}


######################################################################
### auc.R ends here
