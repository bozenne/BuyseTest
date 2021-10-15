### iid.logit.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  3 2021 (11:55) 
## Version: 
## Last-Updated: okt 14 2021 (19:23) 
##           By: Brice Ozenne
##     Update #: 90
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .predict.logit (documentation)
##' @title Predicted Probability with Influence Function
##' @description Compute the predicted probabilities from a logistic regression,
##' with their (robust) standard error,
##' and the corresponding influence function.
##' @name predict.logit
##' 
##' @param object a logistic model.
##' @param newdata [data.frame] dataset containing
##' @param level [character] level of the outcome for which the probability should be computed.
##' @param robust [logit] when FALSE uses the individual contribution to the modeled variance-covariance matrix as iid decomposition.
##' 
##' @keywords internal

## * .predict.logit (examples)
##' @examples
##' \dontrun{ ## will not run as .predict.logit is not exported
##' set.seed(10)
##' n <- 100
##' df <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
##' e.logit <- glm(Y~X1+X2, data = df, family = binomial(link="logit"))
##'
##' test1 <- .predict.logit(e.logit, newdata = df[1:5,])
##' test1["se",] - sqrt(diag(crossprod(attr(test1,"iid")))) ## exact
##' test1["var.se",] - diag(crossprod(attr(test1,"iid.se"))) ## exact
##' 
##' GS <- predict(e.logit, newdata = df[1:5,], se = TRUE, type =  "response")
##' test2 <- .predict.logit(e.logit, newdata = df[1:5,], robust = FALSE)
##' test2["estimate",] - GS$fit ## exact
##' test2["se",] - GS$se.fit ## exact
##' test2["se",] - sqrt(diag(crossprod(attr(test2,"iid")))) ## approximate
##' ## since it uses the robust estimator of the correlation
##' ## (and the modeled estimator of the variance)
##'
##' ## Sanity check (fully stratified model)
##' df <- data.frame(Y = rbinom(n, prob = 0.5, size = 1),
##'                  X1 = rnorm(n),
##'                  X2 = as.factor(rbinom(n, size = 1, prob = 0.5)))
##' newdata <- data.frame(X1=c(0,1),X2=as.factor(0:1))
##' 
##' e.logit <- glm(Y~X1+X2, data = df, family = binomial(link="logit"))
##' e.predlogit <- .predict.logit(e.logit, newdata = newdata)
##' cor(attr(e.predlogit,"iid"))
##' 
##' e.logitS <- glm(Y~X1*X2, data = df, family = binomial(link="logit"))
##' e.predlogitS <- .predict.logit(e.logitS, newdata = newdata)
##' cor(attr(e.predlogitS,"iid"))
##' }



## * .predict.logit (simulation study)
## warper <- function(sim,n){
##     df <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
##     e.logit <- glm(Y~X1+X2, data = df, family = binomial(link="logit"))
##     res <- .predict.logit(e.logit, newdata = data.frame(X1=2,X2=0.5))
##     out <- cbind(sim = sim, n=n,t(res))
##     return(out)
## }
## library(pbapply)
## n.sim <- 1000
## ls.res <- pblapply(1:n.sim, FUN = function(iSim){
##     rbind(warper(sim = iSim, n = 50),
##           warper(sim = iSim, n = 100),
##           warper(sim = iSim, n = 250),
##           warper(sim = iSim, n = 500),
##           warper(sim = iSim, n = 1000))
## })
## dtW.res <- as.data.table(do.call(rbind,ls.res))
## dtW.res[,empirical.se := sd(estimate), by = "n"]
## dtW.res[,empirical.var.se := var(se), by = "n"]
## gg.se <- ggplot(dtW.res, aes(x=as.factor(n))) + geom_boxplot(aes(y=se)) + geom_point(aes(y=empirical.se), shape = 2, size = 2) + geom_line(aes(y=empirical.se, group = "d"))
## gg.var.se <- ggplot(dtW.res, aes(x=as.factor(n))) + geom_boxplot(aes(y=var.se)) + geom_point(aes(y=empirical.var.se), shape = 2, size = 2) + geom_line(aes(y=empirical.var.se, group = "d"))
## gg.var.se + coord_cartesian(ylim = c(0,0.00001))

## * .predict.logit (code)
.predict.logit <- function(object, newdata, level = NULL, robust = TRUE){

    ## ** check input
    if (!inherits(object,"glm")){
        stop("Only implemented for glm objects. \n")
    }
    if (object$family$family!="binomial"){
        stop("Not implemented for other families than binomial. \n")
    }
    if (object$family$link!="logit"){
        stop("Not implemented for other link function than logit. \n")
    }
    n.newobs <- NROW(newdata)

    ## ** prepare
    ff <- stats::formula(object)
    X <- stats::model.matrix(delete.response(terms(ff)), newdata)
    beta <- coef(object)
    name.param <- names(beta)
    n.param <- length(beta)
    iid.beta <- lava::iid(object)
    n.obs <- NROW(iid.beta)

    ## ** compute predictions
    Xbeta <- as.double(X %*% beta)
    ## Xbeta - predict(object, type = "link", newdata = newdata, se = FALSE)
    pred <- 1/(1+exp(-Xbeta))
    ## pred - predict(object, type = "response", newdata = newdata, se = FALSE)

    ## ** compute variance of the predictions
    if(robust){
        Sigma.beta <- crossprod(iid.beta)
    }else{
        Sigma.beta <- stats::vcov(object)
    }

    ## var.Xbeta <- X %*% Sigma.beta %*% t(X)
    var.Xbeta <- rowSums((X %*% Sigma.beta) * X)
    ## var.pred  <- var.Xbeta * (-exp(-Xbeta)  / (1+exp(-Xbeta))^2)^2
    var.pred  <- var.Xbeta * exp(-2*Xbeta)  / (1+exp(-Xbeta))^4
    se.pred <- sqrt(var.pred)
    ## se.pred - predict(object, type = "response", newdata = newdata, se = TRUE)$se.fit
    
    ## ** compute influence function of the predictions
    if(robust == FALSE){
        seR.beta <- sqrt(colSums(iid.beta^2))
        se.beta <- sqrt(diag(Sigma.beta))
        iid.beta  <- .rowMultiply_cpp(iid.beta, se.beta/seR.beta)
    }
    iid.pred <- iid.beta %*% t(.colMultiply_cpp(X, scale = exp(-Xbeta)/(1+exp(-Xbeta))^2))
    ## colSums(iid.pred)
    ## colSums(iid.pred^2) - var.pred

    ## ** compute the influence function of the variance of the prediction
    if(robust){
        iid.Sigma.beta <- array(unlist(lapply(1:n.obs, function(i){crossprod(iid.beta[i,,drop=FALSE])-Sigma.beta/n.obs})),
                                dim = c(n.param,n.param,n.obs), dimnames = list(name.param,name.param,NULL))
    }else{
        iid.Sigma.beta <- .vcov.logit(object, indiv = TRUE, center = TRUE)
    }
    ## apply(iid.Sigma.beta,1:2,sum)
    iid.var.Xbeta <- do.call(rbind,lapply(1:n.obs, FUN = function(iObs){rowSums((X %*% iid.Sigma.beta[,,iObs]) * X)}))
    iid.var.e2XB_1eXb4 <- iid.beta %*% t(.colMultiply_cpp(X, scale = -2*exp(-2*Xbeta)/(1+exp(-Xbeta))^4 + 4*exp(-3*Xbeta)/(1+exp(-Xbeta))^5))
    ## colSums(iid.var.Xbeta)
    iid.var.pred <- .rowMultiply_cpp(iid.var.Xbeta, exp(-2*Xbeta)  / (1+exp(-Xbeta))^4) + .rowMultiply_cpp(iid.var.e2XB_1eXb4, var.Xbeta) 
    iid.se.pred <- .rowScale_cpp(iid.var.pred, 2*se.pred)
    
    ## ** set correct level
    if(!is.null(level)){
        matching.Ylevel <- table(object$data[[all.vars(formula(object))[1]]],
                                 object$y)
        all.levels <- rownames(matching.Ylevel)
        level <- match.arg(level, all.levels)

        index.level <- which(matching.Ylevel[level,]>0)
        if(length(index.level) > 1){
            stop("Unknown value for the outcome variable \n")
        }else if(index.level == 1){
            out <- 1 - out
            iid.pred <- - iid.pred
        }
    }

    
    ## ** export
    out <- rbind(estimate = pred,
                 se = se.pred,
                 var.se = colSums(iid.se.pred^2))
    attr(out,"iid") <- iid.pred
    attr(out,"iid.se") <- iid.se.pred
    return(out)
}

## * .score.logit
#' @title Score for Logistic Regressions
#' @description Compute the first derivative of the log-likelihood IPCW logistic regressions.
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total score?
#' 
#' @keywords internal
.score.logit <- function(object, indiv){
    X <- stats::model.matrix(object)
    pi <- stats::predict(object, type = "response")
    Y <- object$y
    W <- object$prior.weights
    out <- .colMultiply_cpp(X, scale = W*(Y - pi))
    colnames(out) <- colnames(X)
    if(indiv){
        return(out)
    }else{
        return(colSums(out))
    }
}

## * .information.wglm
#' @title Information for Logistic Regressions
#' @description Compute the information (i.e. opposit of the expectation of the second derivative of the log-likelihood) for logistic regressions.
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total information?
#' @param center [logical] should the individual contribution be centered around the average?
#' 
#' @keywords internal
.information.logit <- function(object, indiv, center){
    X <- stats::model.matrix(object)
    n.obs <- NROW(X)
    n.param <- NCOL(X)
    pi <- stats::predict(object, type = "response")
    W <- object$prior.weights
    tXWpi <- t(.colMultiply_cpp(X, scale = W*pi*(1-pi)))

    if(indiv){
        if(center){
            Info <- tXWpi %*% X
            out <- array(unlist(lapply(1:n.obs, function(i){tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE]-Info/n.obs})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }else{
            out <- array(unlist(lapply(1:n.obs, function(i){tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE]})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }
        return(out)
        ## apply(out, MARGIN = 1:2, FUN = sum) - tXWpi %*% X
        
    }else{
        out <- tXWpi %*% X
        rownames(out) <- colnames(out)
        return(out)
    }
}

## * .vcov.wglm
#' @title Variance-covariance matrix for Logistic Regressions
#' @description Compute the variance covariance matrix (i.e. inverse of the information) for logistic regressions.
#'
#' @param object a glm object corresponding to a logistic regression.
#' @param indiv [logical] should the individual contribution be output instead of the total variance-covariance?
#' @param center [logical] should the individual contribution be centered around the average?
#' 
#' @keywords internal
.vcov.logit <- function(object, indiv, center){
    X <- stats::model.matrix(object)
    n.obs <- NROW(X)
    n.param <- NCOL(X)
    pi <- stats::predict(object, type = "response")
    W <- object$prior.weights
    tXWpi <- t(.colMultiply_cpp(X, scale = W*pi*(1-pi)))
    Sigma <- solve(tXWpi %*% X)
    
    if(indiv){
        if(center){
            out <- array(unlist(lapply(1:n.obs, function(i){Sigma %*% tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE] %*% Sigma - Sigma/n.obs})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }else{
            out <- array(unlist(lapply(1:n.obs, function(i){Sigma %*% tXWpi[,i,drop=FALSE] %*% X[i,,drop=FALSE] %*% Sigma})),
                         dim = c(n.param,n.param,n.obs), dimnames = list(colnames(X),colnames(X),NULL))
        }
        return(out)
        ## apply(out, MARGIN = 1:2, FUN = sum) - tXWpi %*% X
        
    }else{
        out <- Sigma
        rownames(out) <- colnames(out)
        return(out)
    }
}


##----------------------------------------------------------------------
### iid.logit.R ends here
