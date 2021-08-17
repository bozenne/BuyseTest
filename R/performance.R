### performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  3 2021 (11:17) 
## Version: 
## Last-Updated: aug 17 2021 (19:06) 
##           By: Brice Ozenne
##     Update #: 224
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Assess Performance of a Classifier
##' @description Assess the performance in term of AUC and brier score of one or several binary classifiers.
##' Currently limited to logistic regressions and random forest.
##'
##' @param object a \code{glm} or \code{range} object, or a list of such object.
##' @param data [data.frame] the training data.
##' @param newdata [data.frame] an external data used to assess the performance.
##' @param fold.size [double, >0] Either the size of the test dataset (when >1) or the fraction of the dataset (when <1) to be used for testing when using cross-validation.
##' @param fold.number [integer] When strictly positive, the number of fold used in the cross-validation. If 0 then no cross validation is performed.
##' @param name.response [character] The name of the response variable (i.e. the one containing the categories).
##' @param null [numeric vector of length 2] the right-hand side of the null hypothesis relative to each metric.
##' @param conf.level [numeric] confidence level for the confidence intervals.
##' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
##' Otherwise they are computed without any transformation.
##' @param trace [logical] Should the execution of the function be traced.
##' 
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
##' @export
performance <- function(object, data = NULL, newdata = NA, fold.size = 1/10, fold.number = 0,
                 name.response = NULL,
                 null = c(brier = NA, AUC = 0.5), conf.level = 0.95, transformation = TRUE,
                 trace = TRUE){

    ## ** normalize user input
    ## convert object into a list of models with names
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
        possible.names <- lapply(1:length(object), function(iO){
            if(inherits(object[[iO]],"ranger")){
                return(paste0("ranger",iO))
            }else if(inherits(object[[iO]],"glm")){
                if(object[[iO]]$family$family!="binomial"){
                    stop("Cannot handle glm with family other than binomial. \n")
                }
                if(object[[iO]]$family$link!="logit"){
                    stop("Cannot handle glm with link other than logit. \n")
                }
                return(paste0("logit",iO))
            }else{
                stop("Argument \'object\' should be a \'glm\' object, a \'ranger\' object, or a list containing \'glm\' or \'ranger\' objects. \n")
            }
        })
        if(is.null(names(object))){
            names(object) <- possible.names
        }
    }
    ## extract full training data
    if(is.null(data)){
        ls.data <- lapply(object, stats::model.frame)
        if(length(unique(sapply(object, NROW)))>1){
            stop("All models should be fitted using the same dataset. \n")
        }
        ls.names <- lapply(ls.data,names)
        index.max <- which.max(sapply(ls.names, length))[1]
        max.names <- ls.names[[index.max]]
        test.max <- sapply(ls.names, function(iNames){all(iNames %in% max.names)})
        data <- ls.data[[index.max]]
        if(all(test.max)==FALSE){
            for(iObj in 1:n.object){
                add.cols <- ls.names[[iObj]][ls.names[[iObj]] %in% colnames(data)]
                if(length(add.cols)>1){
                    data <- cbind(data,ls.data[[iObj]][,add.cols])
                }
            }
        }
    }
    data <- as.data.frame(data)

    ## extract name of the outcome variable
    if(is.null(name.response)){
        name.response <- all.vars(formula(object[[1]]))[1]
    }
    ref.response <- sort(data[[name.response]], decreasing = TRUE)[1]

    if(is.character(data[[name.response]])){
        data[[name.response]] <- as.factor(data[[name.response]])
    }
    if(is.factor(data[[name.response]])){
        if(length(levels(data[[name.response]]))==2){
            data[[name.response]] <- as.numeric(data[[name.response]])-1
        }else stop("The column corresponding to the argument \'name.response\' should take exactly two different values. \n",
                   "Unique values found: \"",paste0(levels(data[[name.response]]), collapse = "\" \""),"\".\n")
    }else if(any(data[[name.response]] %in% 0:1 == FALSE)){
        stop("The column corresponding to the argument \'name.response\' should correspond to a binary variable. \n",
             "Unique values found: \"",paste0(levels(data[[name.response]]), collapse = "\" \""),"\".\n")
    }
    
    ## null hypothesis
    if("brier" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called brier. \n",
             "(corresponding to the null hypothesis for the brier score, possibly NA)")
    }
    if("AUC" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called AUC. \n",
             "(corresponding to the null hypothesis for the AUC, possibly NA)")
    }
    ## extract sample size of the training set
    nobs.object <- unname(unique(sapply(object,function(iO){ ## iO <- object$RF
        if(inherits(iO,"glm")){
            return(stats::nobs(iO))
        }else if(inherits(iO,"ranger")){
            return(iO$num.samples)
        }
    })))
    if(length(nobs.object)!=1){
        stop("The training set seems to differ in size between models. \n")
    }

    ## arguments for cross validation
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
        if(fold.size<2){
            stop("Argument \'fold.size\' must be greater or equal to 2 (at least one sample in each category) \n")
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
                iPred <- predict(object[[iO]], data = data, type = "se")
                internal.predictions[,iO] <- iPred$predictions
                internal.se.predictions[,iO] <- iPred$se
                internal.iid[[iO]] <- matrix(0, nrow = nData.obs, ncol = nData.obs)
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
        if(n.object>1){
            internal.auc <- cbind(internal.auc, p.value_comp = NA)
            internal.brier <- cbind(internal.brier, p.value_comp = NA)
        }
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iAUC <- auc(labels = data[[name.response]], predictions = internal.predictions[,iO],
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            internal.auc[iO,] <- confint(iAUC)
            internal.iid.auc[,iO] <- iid(iAUC)
            if(iO>1){
                iStat <- (internal.auc[iO,"estimate"] - internal.auc[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.auc[,iO]-internal.iid.auc[,iO-1]))
                internal.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            iBrier <- brier(labels = data[[name.response]], predictions = internal.predictions[,iO], iid = internal.iid[[iO]],
                            null = null["Brier"], conf.level = conf.level, transformation = transformation)
            internal.brier[iO,] <- confint(iBrier)
            internal.iid.brier[,iO] <- iid(iBrier)
            if(iO>1){
                iStat <- (internal.brier[iO,"estimate"] - internal.brier[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.brier[,iO]-internal.iid.brier[,iO-1]))
                internal.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            
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
                iPred <- predict(object[[iO]], data = newdata, type = "se")
                external.predictions[,iO] <- iPred$predictions
                external.se.predictions[,iO] <- iPred$se
                external.iid[[iO]] <- matrix(0, nrow = nData.obs, ncol = nNewdata.obs)
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
        external.iid.auc <- matrix(NA, nrow = nNewdata.obs, ncol = n.object, dimnames = list(NULL,names.object))
        external.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 5, dimnames = list(names.object,c("estimate","se","lower","upper","p.value"))))
        external.iid.brier <- matrix(NA, nrow = nobs.object + nNewdata.obs, ncol = n.object, dimnames = list(NULL,names.object))
        if(n.object>1){
            external.auc <- cbind(external.auc, p.value_comp = NA)
            external.brier <- cbind(external.brier, p.value_comp = NA)
        }
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iBrier <- brier(labels = newdata[[name.response]], predictions = external.predictions[,iO], iid = external.iid[[iO]], observation = "external",
                            null = null["Brier"], conf.level = conf.level, transformation = transformation)

            external.brier[iO,] <- confint(iBrier)
            external.iid.brier[,iO] <- iid(iBrier)
            if(iO>1){
                iStat <- (external.brier[iO,"estimate"] - external.brier[iO-1,"estimate"]) / sqrt(crossprod(external.iid.brier[,iO]-external.iid.brier[,iO-1]))
                external.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            iAUC <- auc(labels = newdata[[name.response]], predictions = external.predictions[,iO],
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            external.auc[iO,] <- confint(iAUC)
            external.iid.auc[,iO] <- iid(iAUC)
            if(iO>1){
                iStat <- (external.auc[iO,"estimate"] - external.auc[iO-1,"estimate"]) / sqrt(crossprod(external.iid.auc[,iO]-external.iid.auc[,iO-1]))
                external.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

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
        index.response1 <- which(data[[name.response]]==ref.response)
        index.response0 <- which(data[[name.response]]!=ref.response)
        fold.test <- do.call(cbind,lapply(1:fold.number, function(iFold){
            iSample <- c(sample(index.response0, size = 1), sample(index.response1, size = 1))
            if(fold.size>2){
                c(iSample,sample(setdiff(1:nData.obs,iSample), replace = FALSE, size = fold.size-2))
            }else{
                return(iSample)
            }
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
                iObject <- stats::update(object[[iO]], data = iDataTrain)
                if(inherits(iObject,"ranger")){
                    iPred <- predict(iObject, data = iDataTest, type = "se")
                    cv.predictions[,iO,iFold] <- iPred$predictions
                    cv.se.predictions[,iO,iFold] <- iPred$se
                    cv.iid[[iO]][fold.train[,iFold],,iFold] <- matrix(0, nrow = nData.obs-fold.size, ncol = fold.size)
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
        if(n.object>1){
            cv.auc <- cbind(cv.auc, p.value_comp = NA)
            cv.brier <- cbind(cv.brier, p.value_comp = NA)
        }
        ls.auc <- setNames(vector(mode = "list", length = n.object), names.object)
        ls.brier <- setNames(vector(mode = "list", length = n.object), names.object)

        for(iO in 1:n.object){ ## iO <- 1
            if(trace){cat("*")}
            iObs <- as.double(cv.indexing[,"observation",])
            iFold <- as.double(cv.indexing[,"fold",])
            iPred <- as.double(cv.predictions[,iO,])
            ls.auc[[iO]] <- auc(labels = data[[name.response]], predictions = iPred, fold = iFold, observation = iObs,
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            cv.auc[iO,] <- confint(ls.auc[[iO]])
            cv.iid.auc[,iO] <- iid(ls.auc[[iO]])
            if(iO>1){
                iStat <- (cv.auc[iO,"estimate"] - cv.auc[iO-1,"estimate"]) / sqrt(crossprod(cv.iid.auc[,iO]-cv.iid.auc[,iO-1]))
                cv.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            ## sqrt(crossprod(rowMeans(attr(ls.auc[[iO]],"iid")))) - confint(ls.auc[[iO]])["se"]
            ls.brier[[iO]] <- brier(labels = data[[name.response]], predictions = iPred, fold = iFold, observation = iObs, iid = cv.iid[[iO]],
                           null = null["Brier"], conf.level = conf.level, transformation = transformation)
            cv.brier[iO,] <- confint(ls.brier[[iO]])
            cv.iid.brier[,iO] <- iid(ls.brier[[iO]])
            if(iO>1){
                iStat <- (cv.brier[iO,"estimate"] - cv.brier[iO-1,"estimate"]) / sqrt(crossprod(cv.iid.brier[,iO]-cv.iid.brier[,iO-1]))
                cv.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }
        }

        ## *** export
        out <- rbind(out,
                     cbind(model = rownames(cv.auc), method = "cv", metric = "auc", cv.auc),
                     cbind(model = rownames(cv.brier), method = "cv", metric = "brier", cv.brier)
                     )
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["cv"]] <- cv.iid.auc
        attr(out,"iid")$brier[["cv"]] <- cv.iid.brier
        attr(out,"auc") <- ls.auc
        attr(out,"brier") <- ls.brier
        if(trace){cat("\n  done. \n")}
        
    }


    ## ** export
    rownames(out) <- NULL
    return(out)
}

##----------------------------------------------------------------------
### performance.R ends here
