### performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  3 2021 (11:17) 
## Version: 
## Last-Updated: aug  5 2021 (18:02) 
##           By: Brice Ozenne
##     Update #: 108
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @examples
##'if(FALSE){
##' ## Simulate data
##' set.seed(10)
##' n <- 100
##' df <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
##'
##' ## fit logistic model
##' e.logit <- glm(Y~X1+X2, data = df, family = binomial(link="logit"))
##'
##' performance(e.logit, fold.number = 10)
##' performance(e.logit, fold.number = 20)
##' auc(labels = df$Y, predictions = fitted(e.logit))
##' }

## * performance
performance <- function(object, data = NULL, newdata = NA, fold.size = 1/10, fold.number = 0,
                 name.response = NULL,
                 null = c(brier = NA, AUC = 0.5), conf.level = 0.95, transformation = TRUE,
                 trace = TRUE){

    ## ** normalize user input
    if(is.null(data)){
        data <- stats::model.frame(object)
    }
    if(is.null(name.response)){
        name.response <- all.vars(formula(object))[1]
    }
    if(inherits(object,"ranger")){
        object <- list(ranger1 = object)
    }else  if(inherits(object,"glm")){
        if(object$family$family!="binomial"){
            stop("Cannot handle glm with family other than binomial. \n")
        }
        if(object$family$link!="logit"){
            stop("Cannot handle glm with link other than logit. \n")
        }
        object <- list(logit1 = object)
    }else if(!is.list(object)){
        stop("Argument \'object\' should be a \'glm\' object, a \'ranger\' object, or a list containing \'glm\' or \'ranger\' objects. \n")
    }else{
        possible.names <- lapply(object, function(iO){
            if(inherits(iO,"ranger")){
                return("ranger")
            }else if(inherits(object,"glm")){
                if(iO$family$family!="binomial"){
                    stop("Cannot handle glm with family other than binomial. \n")
                }
                if(iO$family$link!="logit"){
                    stop("Cannot handle glm with link other than logit. \n")
                }
                return("logit")
            }else{
                stop("Argument \'object\' should be a \'glm\' object, a \'ranger\' object, or a list containing \'glm\' or \'ranger\' objects. \n")
            }
        })
        if(is.null(names(object))){
            names(object) <- possible.names
        }
    }
    if("brier" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called brier. \n",
             "(corresponding to the null hypothesis for the brier score, possibly NA)")
    }
    if("AUC" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called AUC. \n",
             "(corresponding to the null hypothesis for the AUC, possibly NA)")
    }
    nobs.object <- unname(unique(sapply(object,stats::nobs)))
    if(length(nobs.object)!=1){
        stop("The training set seems to differ in size between models. \n")
    }

    if(fold.number!=0 && !is.na(fold.number) && !is.null(fold.number)){
        if(!inherits(data,"data.frame") && !inherits(data,"matrix")){
            stop("Argument \'data\' must be specified when using cross-validation. \n")
        }
        data <- as.data.frame(data)
        if(fold.size<0){
            stop("Argument \'fold.size\' must be strictly positive \n")
        }
        if(fold.size<1){
            fold.size <- ceiling(nobs.object * fold.size)
        }
    }else if((!is.na(data) && !identical(data,FALSE))){
        data <- as.data.frame(data)
    }
    names.object <- names(object)
    n.object <- length(object)


    out <- NULL
    ## ** internal performance
    if(!is.na(data) && !identical(data,FALSE)){
        if(trace){cat("- assess internal performance: ")}
        nData.obs <- NROW(data)

        ## *** get predictions
        internal.predictions <- matrix(NA, nrow = nData.obs, ncol = n.object,
                                       dimnames = list(NULL, names.object))
        internal.se.predictions <- matrix(NA, nrow = nData.obs, ncol = n.object,
                                          dimnames = list(NULL, names.object))
        internal.iid <- setNames(vector(mode = "list", length = n.object), names.object)
        
        for(iO in 1:n.object){
            if(trace){cat(".")}

            if(inherits(object[[iO]],"ranger")){
                iPred <- predict(object[[iO]], data = data, type = "response")
                browser()
            }else if(inherits(object[[iO]],"glm")){
                iPred <- .predict.logit(object[[iO]], newdata = data)
                internal.predictions[,iO] <- iPred["estimate",]
                internal.se.predictions[,iO] <- iPred["se",]
                internal.iid[[iO]] <- attr(iPred,"iid")
            }
        }
        if(trace){cat(" ")}

        ## *** assess performance
        internal.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        internal.iid.auc <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        internal.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        internal.iid.brier <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iAUC <- auc(labels = data[[name.response]], predictions = internal.predictions[,iO],
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            internal.auc[iO,] <- confint(iAUC)
            internal.iid.auc[,iO] <- iid(iAUC)

            iBrier <- brier(labels = data[[name.response]], predictions = internal.predictions[,iO], iid = internal.iid[[iO]],
                           null = null["Brier"], conf.level = conf.level, transformation = transformation)
            internal.brier[iO,] <- confint(iBrier)
            internal.iid.brier[,iO] <- iid(iBrier)
        }

        ## *** export
        out <- rbind(out,
                     cbind(model = rownames(internal.auc), method = "internal", metric = "auc", internal.auc),
                     cbind(model = rownames(internal.brier), method = "internal", metric = "brier", internal.brier)
                     )
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["internal"]] <- internal.iid.auc
        attr(out,"iid")$brier[["internal"]] <- internal.iid.brier
        if(trace){cat("\n  done. \n")}
    }

    ## ** external performance
    if(!is.null(newdata) && !is.na(newdata) && !identical(newdata,FALSE)){
        if(trace){cat("- assess external performance: ")}
        nNewdata.obs <- NROW(newdata)

        ## *** get predictions
        external.predictions <- matrix(NA, nrow = nNewdata.obs, ncol = n.object,
                                       dimnames = list(NULL, names.object))
        external.se.predictions <- matrix(NA, nrow = nNewdata.obs, ncol = n.object,
                                          dimnames = list(NULL, names.object))
        external.iid <- setNames(vector(mode = "list", length = n.object), names.object)
        
        for(iO in 1:n.object){
            if(trace){cat(".")}

            if(inherits(object[[iO]],"ranger")){
                iPred <- predict(object[[iO]], data = newdata, type = "response")
                browser()
            }else if(inherits(object[[iO]],"glm")){
                iPred <- .predict.logit(object[[iO]], newdata = newdata)
                external.predictions[,iO] <- iPred["estimate",]
                external.se.predictions[,iO] <- iPred["se",]
                external.iid[[iO]] <- attr(iPred,"iid")
            }
        }
        if(trace){cat(" ")}
        
        ## *** assess performance
        external.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        external.iid.auc <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        external.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        external.iid.brier <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iAUC <- auc(labels = newdata[[name.response]], predictions = external.predictions[,iO],
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            external.auc[iO,] <- confint(iAUC)
            external.iid.auc[,iO] <- iid(iAUC)

            iBrier <- brier(labels = newdata[[name.response]], predictions = external.predictions[,iO], iid = external.iid[[iO]],
                            null = null["Brier"], conf.level = conf.level, transformation = transformation)
            external.brier[iO,] <- confint(iBrier)
            external.iid.brier[,iO] <- iid(iBrier)
        }
        
        ## *** export
        out <- rbind(out,
                     cbind(model = rownames(external.auc), method = "external", metric = "auc", external.auc),
                     cbind(model = rownames(external.brier), method = "external", metric = "brier", external.brier)
                     )
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["external"]] <- external.iid.auc
        attr(out,"iid")$brier[["external"]] <- external.iid.brier
        if(trace){cat("\n  done. \n")}
    }

    ## ** CV performance
    if(fold.number>0){
        if(trace){cat("- assess performance using cross-validation: \n")}
        nData.obs <- NROW(data)

        ## *** identify folds
        fold.test <- do.call(cbind,lapply(1:fold.number, function(iFold){
            sample(1:nData.obs, replace = FALSE, size = fold.size)
        }))
        fold.train <- do.call(cbind,lapply(1:fold.number, function(iFold){
            setdiff(1:nData.obs,fold.test[,iFold])
        }))

        ## *** get predictions
        cv.indexing <- array(NA, dim = c(fold.size, 2, fold.number),
                             dimnames = list(NULL, c("observation","fold"), NULL))
        cv.predictions <- array(NA, dim = c(fold.size, n.object, fold.number),
                                dimnames = list(NULL, names.object, NULL))
        cv.se.predictions <- array(NA, dim = c(fold.size, n.object, fold.number),
                                   dimnames = list(NULL, names.object, NULL))
        cv.iid <- setNames(lapply(1:n.object, function(iFold){
            array(0, dim = c(nData.obs, fold.size, fold.number))
        }), names.object)

        if(trace){
            pb <- utils::txtProgressBar(max = fold.number, style = 3)
        }
        for(iFold in 1:fold.number){ ## iFold <- 1
            if(trace){
                utils::setTxtProgressBar(pb, iFold)
            }
            iDataTrain <- data[fold.train[,iFold],,drop=FALSE]
            iDataTest <- data[fold.test[,iFold],,drop=FALSE]
            cv.indexing[,"observation",iFold] <- fold.test[,iFold]
            cv.indexing[,"fold",iFold] <- iFold
            
            for(iO in 1:n.object){ ## iO <- 1
                iObject <- update(object[[iO]], data = iDataTrain)
                
                if(inherits(iObject,"ranger")){
                    iPred <- predict(iObject, data = iDataTest, type = "response")
                    browser()
                }else if(inherits(iObject,"glm")){
                    iPred <- .predict.logit(iObject, newdata = iDataTest)
                    cv.predictions[,iO,iFold] <- iPred["estimate",]
                    cv.se.predictions[,iO,iFold] <- iPred["se",]
                    cv.iid[[iO]][fold.train[,iFold],,iFold] <- attr(iPred,"iid")
                }
            }
        }
        if(trace){close(pb)}
        
        ## *** compute auc
        cv.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        cv.iid.auc <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        cv.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        cv.iid.brier <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iObs <- as.double(cv.indexing[,"observation",])
            iFold <- as.double(cv.indexing[,"fold",])
            iPred <- as.double(cv.predictions[,iO,])

            iAUC <- auc(labels = data[[name.response]], predictions = iPred, fold = iFold, observation = iObs,
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            cv.auc[iO,] <- confint(iAUC)
            cv.iid.auc[,iO] <- iid(iAUC)
            ## sqrt(crossprod(rowMeans(attr(iAUC,"iid")))) - confint(iAUC)["se"]
            
            iBrier <- brier(labels = data[[name.response]], predictions = iPred, fold = iFold, observation = iObs, iid = cv.iid[[iO]],
                           null = null["Brier"], conf.level = conf.level, transformation = transformation)
            cv.brier[iO,] <- confint(iBrier)
            cv.iid.brier[,iO] <- iid(iBrier)
        }

        ## *** export
        out <- rbind(out,
                     cbind(model = rownames(cv.auc), method = "cv", metric = "auc", cv.auc),
                     cbind(model = rownames(cv.brier), method = "cv", metric = "brier", cv.brier)
                     )
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["cv"]] <- cv.iid.auc
        attr(out,"iid")$brier[["cv"]] <- cv.iid.brier
        if(trace){cat("\n  done. \n")}
        
    }


    ## ** export
    rownames(out) <- NULL
    return(out)
}

##----------------------------------------------------------------------
### performance.R ends here
