### auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:29) 
## Version: 
## Last-Updated: dec  3 2019 (17:50) 
##           By: Brice Ozenne
##     Update #: 107
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
#' @aliases auc
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
#' 
#' @details Compared to other functions computing the AUC (e.g. the auc fonction of the ROCR package),
#' the AUC is defined here as P[Y>X] with a strict inequality sign (i.e. not P[Y>=X]).
#'
#' @return A \emph{data.frame} containing for each fold the AUC value with its standard error (when computed).
#' The last line of the data.frame contains the global AUC value with its standard error.
#' 

## * Example - auc
#' @rdname auc
#' @examples
#'
#'
#' 

## * Code - auc
#' @rdname auc
#' @export
auc <- function(labels, predictions, fold = NULL, observation = NULL, direction = ">",
                null = 0.5, conf.level = 0.95){

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
    
    df <- data.frame(Y = labels,
                     X = predictions)
    formula <- Y ~ cont(X)

    if(!is.null(fold)){
        df$fold <- fold
        formula <- update(formula,.~.+fold)
    }else{
        df$fold <- 1
    }

    ## ** Perform GPC
    order.save <- BuyseTest.options()$order.Hprojection

    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(formula, method.inference = "u-statistic", data = df, trace = 0)
    BuyseTest.options(order.Hprojection = order.save)

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
                iIID[iVec] <- iid(e.BT)[iVec,"favorable",drop=FALSE]
                return(iIID)
            }))
            out$se[1:n.fold] <- sqrt(colSums(M.iid^2))
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
                iIID[iVec] <- iid(e.BT)[iVec,"unfavorable",drop=FALSE]
                return(iIID)
            }))
            out$se[1:n.fold] <- sqrt(colSums(M.iid^2))
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
            iObservation <- observation[fold == name.fold[iFold]]
            
            out$direction[iFold] <- if(iDirection=="favorable"){">"}else{"<"}
            out$estimate[iFold] <- iCount/e.BT@n.pairs[iFold]

            if(!is.null(observation)){
                M.iid[iObservation,iFold] <- iid(e.BT)[iObservation,iDirection,drop=FALSE]
                out$se[iFold] <- sqrt(sum(M.iid[,iFold]^2))
            }
        }
    }
    out$estimate[n.fold+1] <- mean(out$estimate[1:n.fold])
    if(!is.null(observation)){
        out$se[n.fold+1] <- sqrt(sum(rowSums(M.iid)^2))
    }
    
    z.stat <- as.double((out[,"estimate"]-null)/out[,"se"])
    p.value <- 2*(1-pnorm(abs(z.stat)))

    alpha <- 1-conf.level
    out <- cbind(out,
                 lower = as.double(out[,"estimate"] + qnorm(alpha/2)*out[,"se"]),
                 upper = as.double(out[,"estimate"] + qnorm(1-alpha/2)*out[,"se"]),
                 p.value = p.value)
    
    ## ** Export
    class(out) <- append("BuyseTestAuc",class(out))
    attr(out, "contrast") <- e.BT@level.treatment
    attr(out, "n.fold") <- n.fold
    if(!is.null(observation)){
        attr(out, "iid") <- M.iid
    }
    return(out)

}

######################################################################
### auc.R ends here
