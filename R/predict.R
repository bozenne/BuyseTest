### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2026 (14:25) 
## Version: 
## Last-Updated: jul  3 2026 (14:26) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predict.BuyseTTEM
#' @title Prediction with Time to Event Model
#' @name predict.BuyseTTEM
#' 
#' @description Evaluate the cumulative incidence function (cif) / survival in one of the treatment groups. 
#' 
#' @param object time to event model. 
#' @param times [numeric vector] time at which to evaluate the cif/survival.
#' @param treatment [character/integer] Treatment or index of the treatment group.
#' @param strata  [character/integer] Strata or index of the strata. Can be a one row data.frame where each column contain a strata variable.
#' @param newdata  [data.frame] alternative to specify the argument treatment and strata values in a data.frame format.
#' @param cause [integer] The cause relative to which the cif will be evaluated.
#' @param se [logical] Should the standard errors associated with the cif/survival be output? 
#' @param iid [logical] Should the influence function associated with the cif/survival be output? 
#' @param type [character] by default (\code{NULL}), output the cif (competing risks) or survival (no competing risks) at the timepoints specified by argument times.
#' Otherwise, when argument times is missing, output the discretisation times for the survival (\code{"jump"}),
#' or the survival at the last observation time (\code{"last"}),
#' or the last observation time (\code{"last.time"}),
#' or the change in cif or survival between discretisation times (\code{"change"}).
#' @param ... not used, for compatibility with the generic method.
#'
#' @return a list containing the survival (element \code{survival}) or the cumulative incidence function (element \code{cif}),
#' and possible standard errors (element \code{.se}) and influence function (element \code{.iid}).
#' 
#' @keywords methods
#' 
#' @export
predict.BuyseTTEM <- function(object, times, treatment, strata, newdata = NULL, cause = 1, type = NULL,
                              se = NULL, iid = FALSE, ...){

    ## ** normalize input
    dots <- list(...)
    strata.var <- attr(object$peron$strata,"original")
    
    ## *** treatment
    if(missing(treatment) & missing(strata) & !is.null(newdata)){
        if(dots$level.chaos==0){
            object.original <- object
            class(object.original) <- setdiff(class(object.original), "BuyseTTEM")
            out <- do.call(predict, c(list(object.original, times = times, newdata = newdata, cause = cause, type = type),dots))
            return(out)
        }
        treatment.var <- object$peron$treatment
        if(any(c(treatment.var,strata.var) %in% names(newdata) == FALSE)){
            stop("Argument \'newdata\' must contain variable(s) \"",paste(setdiff(c(treatment.var,strata.var),names(newdata)), collapse = "\", \""),"\" when argument \'treatment\' is missing. \n")
        }
        ls.pred <- lapply(1:NROW(newdata), function(iD){ ## iD <- 1
            iData <- newdata[iD,,drop=FALSE]
            rownames(iData) <- NULL
            if(is.null(strata.var)){
                iPred <- predict(object, times = times, treatment = iData[[treatment.var]], cause = cause, type = NULL, se = se, iid = FALSE)
            }else{
                iPred <- predict(object, times = times, treatment = iData[[treatment.var]], strata = iData[strata.var], cause = cause, type = NULL, se = se, iid = FALSE)
            }
            return(cbind(iData,as.data.frame(iPred)))                
        })
        df.pred <- do.call(rbind,ls.pred)
        return(df.pred)
    }else if(missing(treatment)){
        stop("Argument \'treatment\' is missing. \n",
             "Should take value \"",paste(object$peron$level.treatment, collapse = "\" or \""),"\". \n")
    }else if(!is.null(newdata)){
        message("Argument \'newdata\ is ignored when arguemnt \'treatment\' is specified. \n")
    }
    

    if(!is.numeric(treatment)){
        treatment <- which(treatment == object$peron$level.treatment)
    }

    ## *** strata
    if(missing(strata)){
        if(object$peron$n.strata == 1){
            strata <- 1
        }else{
            stop("Argument \'strata\' argument is missing: needed when the time to event model contains covariates. \n",
                 "Should be  data.frame or take value \"",paste(object$peron$level.strata, collapse = "\" or \""),"\". \n")
        }
    }else{
        if(inherits(strata,"data.frame")){

            if(NROW(strata)!=1){
                stop("Argument \'strata\' should have a single row when being a data.frame. \n")
            }
            if(any(strata.var %in% names(strata) == FALSE)){
                stop("Missing column(s) \"",paste(setdiff(attr(object$peron$strata.var,"original"), names(strata)), collapse = "\", \""),"\" in argument \'strata\'. \n")
            }
            
            if(length(strata.var)==1){
                strata <- strata[[attr(object$peron$strata.var,"original")]]
            }else{
                strata <- interaction(as.data.frame(strata)[strata.var], sep = ".")
            }
            
            if(length(strata.var)==1 && strata.var == "..strata.."){
                ## do nothing as it should already index the strata 
            }else if(any(strata %in% object$peron$level.strata == FALSE)){
                stop("Unknown strata defined by argument \'strata\'. \n",
                     "Should be one of ",paste(object$peron$level.strata, collapse = ", ")," \n")
            }else{
                strata <- which(object$peron$level.strata==strata)
            }
            
        }else{

            if(length(strata)!=1){
                stop("Argument \'strata\' should either be a data.frame specifying the strata variables \n",
                     "or have length 1, e.g. be an integer indexing the strata. \n")
            }
            if(is.numeric(strata)){
                if(strata %in% 1:object$peron$n.strata == FALSE){
                    stop("When numeric argument \'strata\' should be one of ",paste(1:object$peron$n.strata, collapse = ", ")," \n")
                }
            }else if(is.character(strata) || is.factor(strata)){
                if(strata %in% object$peron$level.strata == FALSE){
                    stop("When character/factor argument \'strata\' should be one of ",paste(object$peron$level.strata, collapse = ", ")," \n")
                }
                strata <- which(object$peron$level.strata==strata)
            }else{
                stop("Unknown type for argument \'strata\': should be numeric/character/factor of length 1 or a data.frame. \n")
            }
        }
    }

    ## *** type
    if(!is.null(type)){
        type <- match.arg(type, c("jump","last","change","last.time"))
        if(!missing(times)){
            stop("Argument \'type\' cannot be \"",type,"\" unless argument \'times\' is missing. \n")
        }
        if(iid){
            stop("Argument \'iid\' cannot be TRUE unless argument \'type\' is NULL. \n")
        }
        if(!is.null(se) && se){
            stop("Argument \'se\' cannot be TRUE unless argument \'type\' is NULL. \n")
        }
        se <- FALSE
    }else if(is.null(se)){
        se <- any(lengths(object$peron$iid.cif[[1]][[1]])>0)
    }
    
    ## ** output last cif estimate
    type.CR <- ifelse(object$peron$n.CR==1,"survival","competing.risks")
    if(!is.null(type) && type == "last.time"){        
        return(object$peron$last.time[strata,treatment])
    }else if(!is.null(type) && type == "last"){
        if(type.CR=="survival"){
            return(1-object$peron$last.estimate[strata,treatment,cause])
        }else{
            return(object$peron$last.estimate[strata,treatment,cause])
        }
    }else if(!is.null(type) && type == "jump"){
        return(object$peron$jumpSurvHaz[[strata]][[treatment]]$time.jump)
    }

    ## ** evaluate the cif at the requested times
    if(is.null(type)){
        table.peron <- object$peron$cif[[strata]][[treatment]][[cause]]
        index.table <- pmax(1,prodlim::sindex(jump.time = table.peron$time, eval.time = times)) ## pmin since 1 is taking care of negative times
        if(iid || se){
            tableIID.peron <- lava::iid(object, strata = strata, treatment = treatment, cause = cause)
        }
        
    }else if(type == "change"){
        if(type.CR=="survival"){ ## cif  = 1-survival
            table.peron <- object$peron$jumpSurvHaz[[strata]][[treatment]]
            table.peron$indexBefore.jump <- c(0,table.peron$index.jump[-NROW(table.peron)])
            table.peron$cif <- (1-table.peron$survival[table.peron$index.jump]) - (1-c(1,table.peron$survival)[table.peron$indexBefore.jump+1])
        }else{ ## competing risk case: can only be prodlim object
            table.peron <- object$peron$cif[[strata]][[treatment]][[cause]][!is.na(object$peron$cif[[strata]][[treatment]][[cause]]$obs),c("obs","index","time","cif"),drop=FALSE]
            names(table.peron)[1:3] <- c("index.obs", "index.jump", "time.jump")
            table.peron$indexBefore.jump <- c(0,table.peron$index.jump[-NROW(table.peron)])
            table.peron$cif <- table.peron$cif[table.peron$index.jump] - c(0,table.peron$cif)[table.peron$indexBefore.jump+1]
        }
        index.table <- 1:NROW(table.peron)        
    }

    ## ** export
    if(is.null(type)){
        out <- list(time = table.peron[index.table,"time"],
                    index = table.peron[index.table,"index"])
    }else if(type == "change"){
        out <- list(time = table.peron[index.table,"time.jump"],
                    indexBefore = table.peron[index.table,"indexBefore.jump"],
                    indexAfter = table.peron[index.table,"index.jump"])
    }
    if(type.CR=="survival"){
        if(is.null(type)){
            out$survival <- 1-table.peron[index.table,"cif"]
        }else if(type == "change"){
            out$survival <- -table.peron[index.table,"cif"]
        }
        if(se){ ## only when type is NULL
            out$survival.se <- sqrt(colSums(tableIID.peron[,out$index+1,drop=FALSE]^2))
        }
        if(iid){ ## only when type is NULL
            out$survival.iid <- tableIID.peron[,out$index+1,drop=FALSE]
        }
        
    }else{
        out$cif <- table.peron[index.table,"cif"]
        if(se){## only when type is NULL
            out$cif.se <- sqrt(colSums(tableIID.peron[,out$index+1,drop=FALSE]^2))
        }
        if(iid){ ## only when type is NULL
            out$cif.iid <- tableIID.peron[,out$index+1,drop=FALSE]
        }        
    }
    return(out)
}


## * .predict.logit (documentation)
##' @title Predicted Probability with Influence Function
##' @description Compute the predicted probabilities from a logistic regression,
##' with their (robust) standard error,
##' and the corresponding influence function.
##' @noRd
##' 
##' @param object a logistic model.
##' @param newdata [data.frame] dataset containing
##' @param level [character] level of the outcome for which the probability should be computed.
##' @param robust [logit] when FALSE uses the individual contribution to the modeled variance-covariance matrix as iid decomposition.
##' 

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
    if (object$family$family %in% c("binomial","quasibinomial") == FALSE){
        stop("Not implemented for other families than binomial or quasibinomial. \n")
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


##----------------------------------------------------------------------
### predict.R ends here
