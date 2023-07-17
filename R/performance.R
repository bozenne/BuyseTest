### performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  3 2021 (11:17) 
## Version: 
## Last-Updated: jul  4 2023 (18:46) 
##           By: Brice Ozenne
##     Update #: 1191
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
##' @param fold.size [double, >0] either the size of the test dataset (when >1) or the fraction of the dataset (when <1) to be used for testing when using cross-validation.
##' @param fold.repetition [integer] when strictly positive, the number of folds used in the cross-validation. If 0 then no cross validation is performed.
##' @param fold.balance [logical] should the outcome distribution in the folds of the cross-validation be similar to the one of the original dataset?
##' @param impute [character] in presence of missing value in the regressors of the training dataset, should a complete case analysis be performed (\code{"none"})
##' or should the median/mean (\code{"median"}/\code{"mean"}) value be imputed. For categorical variables, the most frequent value is imputed.
##' @param individual.fit [logical] if \code{TRUE} the predictive model is refit for each individual using only the predictors with non missing values.
##' @param name.response [character] the name of the response variable (i.e. the one containing the categories).
##' @param null [numeric vector of length 2] the right-hand side of the null hypothesis relative to each metric.
##' @param se [logical] should the uncertainty about AUC/brier be computed?
##' When \code{TRUE} adapt the method of LeDell et al. (2015) to repeated cross-validation for the AUC and the brier score.
##' @param conf.level [numeric] confidence level for the confidence intervals.
##' @param auc.type [character] should the auc be computed approximating the predicted probability by a dirac (\code{"classical"}, usual AUC formula)
##' or approximating the predicted probability by a normal distribution.
##' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
##' Otherwise they are computed without any transformation.
##' @param trace [logical] Should the execution of the function be traced.
##' @param simplify [logical] should the number of fold and the size of the fold used for the cross validation be removed from the output?
##' @param seed [integer, >0] Random number generator (RNG) state used when starting data spliting.
##' If \code{NULL} no state is set.
##'
##' @references LeDell E, Petersen M, van der Laan M. Computationally efficient confidence intervals for cross-validated area under the ROC curve estimates. Electron J Stat. 2015;9(1):1583-1607. doi:10.1214/15-EJS1035 
##' 
##' @return An S3 object of class \code{performance}.
##' @keywords model
##' 
##' @examples
##' ## Simulate data
##' set.seed(10)
##' n <- 100
##' df.train <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
##' df.test <- data.frame(Y = rbinom(n, prob = 0.5, size = 1), X1 = rnorm(n), X2 = rnorm(n))
##'
##' ## fit logistic model
##' e.null <- glm(Y~1, data = df.train, family = binomial(link="logit"))
##' e.logit1 <- glm(Y~X1, data = df.train, family = binomial(link="logit"))
##' e.logit2 <- glm(Y~X1+X2, data = df.train, family = binomial(link="logit"))
##'
##' ## assess performance on the training set (biased)
##' ## and external dataset
##' performance(e.logit1, newdata = df.test)
##' e.perf <- performance(list(null = e.null, p1 = e.logit1, p2 = e.logit2),
##'                       newdata = df.test)
##' e.perf
##' summary(e.perf, order.model = c("null","p2","p1"))
##' 
##' ## assess performance using cross validation
##' \dontrun{
##' set.seed(10)
##' performance(e.logit1, fold.repetition = 10, se = FALSE)
##' set.seed(10)
##' performance(list(null = e.null, prop = e.logit1), fold.repetition = 10)
##' performance(e.logit1, fold.repetition = c(50,20,10))
##' }

## * performance
##' @export
performance <- function(object, data = NULL, newdata = NA, individual.fit = FALSE, impute = "none", name.response = NULL,
                        fold.size = 1/10, fold.repetition = 0, fold.balance = FALSE,
                        null = c(brier = NA, AUC = 0.5), conf.level = 0.95, se = TRUE, transformation = TRUE, auc.type = "classical",
                        simplify = TRUE, trace = TRUE, seed = NULL){

    out <- list(call = match.call(),
                response = list(),
                performance = NULL,
                prediction = list(),
                iid.auc = list(),
                iid.brier = list())
    if(se==FALSE){
        out$iid.auc <- NULL
        out$iid.brier <- NULL
    }
    
    ## ** fix randomness
    stop.after.init <- FALSE
    if(identical(seed,"only")){
        stop.after.init <- TRUE
    } else if(!is.null(seed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)
    }

    ## ** normalize user input
    init.args <- .performance_args(object = object,
                                   data = data,
                                   newdata = newdata,
                                   individual.fit = individual.fit,
                                   impute = impute,
                                   name.response = name.response,
                                   fold.size = fold.size,
                                   fold.repetition = fold.repetition,
                                   fold.balance = fold.balance,
                                   null = null,
                                   conf.level = conf.level,
                                   se = se,
                                   transformation = transformation,
                                   auc.type = auc.type,
                                   simplify = simplify,
                                   trace = trace,
                                   seed = seed)

    for(iL in init.args$args){
        assign(x = iL, init.args[[iL]])
    }
    internal <- init.args$internal
    fold.allnumber <- init.args$fold.allnumber
    fold.group <- init.args$fold.group
    fold.test <- init.args$fold.test
    save.data <- init.args$save.data
    names.object <- names(object)
    n.object <- length(object)

    ## ** initialize data
    init.data <- .performance_init(object = object,
                                   data = data,
                                   newdata = newdata,
                                   individual.fit = individual.fit,
                                   name.response = name.response,
                                   fold.size = fold.size,
                                   fold.repetition = fold.repetition,
                                   fold.balance = fold.balance,
                                   fold.group = fold.group,
                                   fold.test = fold.test,
                                   internal = internal,
                                   trace = trace)

    data <- init.data$data
    newdata <- init.data$newdata
    ref.response <- init.data$ref.response
    ref.response.num <- init.data$ref.response.num
    nobs.object <- init.data$nobs.object
    data.missingPattern <- init.data$data.missingPattern  ## only data
    newdata.missingPattern <- init.data$newdata.missingPattern ## data, data+newdata, or newdata acorrding to argument internal
    external <- !is.null(newdata)
    nobs.newdata <- NROW(newdata)
    fold.size <- init.data$fold.size
    fold.test <- init.data$fold.test
    ## any(sapply(fold.test, function(x){any(duplicated(x))}))
    ## sapply(fold.test, function(x){min(x, na.rm=TRUE)})
    ## sapply(fold.test, function(x){max(x, na.rm=TRUE)})
    ## sapply(fold.test, function(x){sum(!is.na(x))})

    if(stop.after.init){
        return(fold.test)
    }
    
    ## ** predictions
    if(trace){
        cat("    Assessment of the predictive performance of ",n.object," model", if(n.object>1){"s"},"\n\n",sep="")
        cat("- Prediction: ")
    }

    ## *** internal/external
    if(internal||external){
        
        if(trace){
            txt <- NULL
            if(internal){txt <- c(txt,"internal")}
            if(external){txt <- c(txt,"external")}
            cat(paste(txt, collapse = " and "), sep = "")
        }
        if(internal){
            data.test <- rbind(data, newdata[,colnames(data),drop=FALSE])
        }else{
            data.test <- newdata
        }
        perf.intext <- .performance_runfold(object, names.object = names.object, n.object = n.object,
                                            data.train = data, n.train = nobs.object, impute = impute,
                                            data.test = data.test, n.test = NROW(data.test), missingPattern.test = newdata.missingPattern,
                                            individual.fit = individual.fit, auc.type = auc.type, se = se, trace = trace-1)

        if(internal){
            internal.predictions <- perf.intext$estimate[1:NROW(data),,drop=FALSE]
            if(auc.type == "probabilistic"){
                internal.se.predictions <- perf.intext$se[1:NROW(data),,drop=FALSE]
            }else{
                internal.se.predictions <- NULL
            }
            if(se>1){
                internal.iid.predictions <- perf.intext$iid[1:NROW(data),,,drop=FALSE]
            }else{
                internal.iid.predictions <- NULL
            }
        }
        if(external){
            if(internal){
                external.predictions <- perf.intext$estimate[(NROW(data)+1):(NROW(data)+NROW(newdata)),,drop=FALSE]
                if(auc.type == "probabilistic"){
                    external.se.predictions <- perf.intext$se[(NROW(data)+1):(NROW(data)+NROW(newdata)),,drop=FALSE]
                }else{
                    external.se.predictions <- NULL
                }
                if(se>1){
                    external.iid.predictions <- perf.intext$iid[(NROW(data)+1):(NROW(data)+NROW(newdata)),,,drop=FALSE]
                }else{
                    external.iid.predictions <- NULL
                }
            }else{
                external.predictions <- perf.intext$estimate
                if(auc.type == "probabilistic"){
                    external.se.predictions <- perf.intext$se
                }else{
                    external.se.predictions <- NULL
                }
                if(se>1){
                    external.iid.predictions <- perf.intext$se
                }else{
                    external.iid.predictions <- NULL
                }
            }
        }
        

        if(trace){
            cat("\n")
        }
    }
    
    ## *** cross validation
    if(any(fold.repetition>0)){

        if(trace){
            if(internal||external){
                space <- rep(" ",14)
            }else{
                space <- ""
            }
            if(fold.group>1){                
                cat(space,fold.group," folds cross-validation repeated ",fold.repetition," times \n",sep="")
            }else{
                cat(space," 1 fold cross-validation with ",fold.size," samples repeated ",fold.repetition," times \n",sep="")
            }
        }

        ## *** get predictions
        cv.indexing <- array(NA, dim = c(sum(fold.size), 3, fold.repetition),
                             dimnames = list(NULL, c("observation","repetition","fold"), NULL))
        cv.predictions <- array(NA, dim = c(sum(fold.size), n.object, fold.repetition),
                                dimnames = list(NULL, names.object, NULL))
        if(auc.type == "probabilistic"){
            cv.se.predictions <- array(NA, dim = c(sum(fold.size), n.object, fold.repetition),
                                       dimnames = list(NULL, names.object, NULL))
        }
        if(se>1){
            cv.iid.predictions <- setNames(lapply(1:n.object, function(iFold){
                array(0, dim = c(nobs.object, sum(fold.size), fold.repetition))
            }), names.object)
        }
        if(trace){
            pb <- utils::txtProgressBar(max = fold.repetition, style = 3)
        }
        for(iRepeat in 1:fold.repetition){ ## iRepeat <- 1
            if(trace){
                utils::setTxtProgressBar(pb, iRepeat)
            }
            cv.indexing[,"repetition",iRepeat] <- iRepeat
            
            for(iFold in 1:fold.group){ ## iFold <- 1
                indexData.iFoldTest <- na.omit(fold.test[[iRepeat]][iFold,])
                indexData.iFoldTrain <- setdiff(1:nobs.object,indexData.iFoldTest)
                
                indexStore.iFoldTest <- (1+sum(c(0,fold.size)[1:iFold])):sum(fold.size[1:iFold])
                cv.indexing[indexStore.iFoldTest,"observation",iRepeat] <- indexData.iFoldTest
                cv.indexing[indexStore.iFoldTest,"fold",iRepeat] <- iFold

                if(!is.null(data.missingPattern)){
                    iData.missingPattern <- lapply(data.missingPattern, function(iModel){ ## iModel <- newdata.missingPattern[[1]]
                        iOut <- droplevels(iModel[indexData.iFoldTest])
                        attr(iOut, "index") <- lapply(attr(iModel, "index"), FUN = function(iVec){which(indexData.iFoldTest %in% iVec)})
                        attr(iOut, "index")[sapply(attr(iOut, "index"), length)==0] <- NULL
                        attr(iOut, "formula") <- attr(iModel, "formula")[names(attr(iOut, "index"))]
                        return(iOut)
                    })
                }else{
                    iData.missingPattern <- NULL
                }
                iPerf.fold <- .performance_runfold(object, names.object = names.object, n.object = n.object,
                                                   data.train = data[indexData.iFoldTrain,,drop=FALSE], n.train = length(indexData.iFoldTrain), impute = impute,
                                                   data.test = data[indexData.iFoldTest,,drop=FALSE], n.test = length(indexData.iFoldTest), missingPattern.test = iData.missingPattern,
                                                   individual.fit = individual.fit, auc.type = auc.type, se = se, trace = FALSE)
                cv.predictions[indexStore.iFoldTest,,iRepeat] <- iPerf.fold$estimate
                if(auc.type == "probabilistic"){
                    cv.se.predictions[indexStore.iFoldTest,,iRepeat] <- iPerf.fold$se
                }
                if(se>1){
                    for(iO in 1:n.object){
                        cv.iid.predictions[indexData.iFoldTrain,indexStore.iFoldTest,iRepeat] <- iPerf.fold[[iO]]
                    }
                }
            
            }
        }
    }
    if(trace){cat("\n")}

    ## ** Performance
    if(trace){cat("- Performance:")}
    ls.auc <- list()
    ls.brier <- list()

    ## *** internal
    if(internal){
        if(trace){cat(" internal")}

        internal.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 7, dimnames = list(names.object,c("model","estimate","se","lower","upper","p.value","p.value_comp"))))
        internal.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 7, dimnames = list(names.object,c("model","estimate","se","lower","upper","p.value","p.value_comp"))))
        if(se>0){
            internal.iid.auc <- matrix(NA, nrow = NROW(data), ncol = n.object, dimnames = list(NULL,names.object))
            internal.iid.brier <- matrix(NA, nrow = NROW(data), ncol = n.object, dimnames = list(NULL,names.object))
        }
        ls.auc$internal <- setNames(vector(mode = "list", length = n.object), names.object)
        ls.brier$internal <- setNames(vector(mode = "list", length = n.object), names.object)

        for(iO in 1:n.object){ ## iO <- 1
            
            if(any(is.na(internal.predictions[,iO]))){
                internal.auc[iO,"model"] <- names.object[iO]
                internal.brier[iO,"model"] <- names.object[iO]
                next
            }

            ## AUC
            ls.auc$internal[[iO]] <- auc(labels = data$XXresponseXX, predictions = internal.predictions[,iO],
                                         add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            internal.auc[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO], confint(ls.auc$internal[[iO]]))
            if(se>0){
                internal.iid.auc[,iO] <- iid(ls.auc$internal[[iO]])
                if(iO>1){
                    iStat <- (internal.auc[iO,"estimate"] - internal.auc[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.auc[,iO]-internal.iid.auc[,iO-1]))
                    internal.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                }
            }

            ## Brier score
            ls.brier$internal[[iO]] <- brier(labels = data$XXresponseXX, predictions = internal.predictions[,iO], iid = internal.iid.predictions[[iO]],
                                             null = null["Brier"], conf.level = conf.level, transformation = transformation)
            internal.brier[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(ls.brier$internal[[iO]]))
            if(se>0){
                internal.iid.brier[,iO] <- iid(ls.brier$internal[[iO]])
                if(iO>1){
                    iStat <- (internal.brier[iO,"estimate"] - internal.brier[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.brier[,iO]-internal.iid.brier[,iO-1]))
                    internal.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                }
            }
        }
        ## export
        out$performance <- rbind(out$performance,
                                 cbind(method = "internal", metric = "auc", internal.auc),
                                 cbind(method = "internal", metric = "brier", internal.brier)
                                 )
        out$response[["internal"]] <- data$XXresponseXX
        out$prediction[["internal"]] <- internal.predictions
        if(se>0){
            out$iid.auc[["internal"]] <- internal.iid.auc
            out$iid.brier[["internal"]] <- internal.iid.brier
        }
        if(trace){cat("(done)")}
    }

    ## *** external
    if(external){
        if(trace){cat(" external")}
    
        external.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 7, dimnames = list(NULL,c("model","estimate","se","lower","upper","p.value","p.value_comp"))))
        external.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 7, dimnames = list(NULL,c("model","estimate","se","lower","upper","p.value","p.value_comp"))))
        if(se>0){
            external.iid.auc <- matrix(NA, nrow = nobs.newdata, ncol = n.object, dimnames = list(NULL,names.object))
            external.iid.brier <- matrix(NA, nrow = NROW(data) + nobs.newdata, ncol = n.object, dimnames = list(NULL,names.object))
        }
        ls.auc$external <- setNames(vector(mode = "list", length = n.object), names.object)
        ls.brier$external <- setNames(vector(mode = "list", length = n.object), names.object)

        for(iO in 1:n.object){
            if(any(is.na(internal.predictions[,iO]))){
                external.auc[iO,"model"] <- names.object[iO]
                external.brier[iO,"model"] <- names.object[iO]
                next
            }

            ## AUC
            ls.auc$external[[iO]] <- auc(labels = newdata$XXresponseXX, predictions = external.predictions[,iO],
                                         add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            external.auc[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(ls.auc$external[[iO]]))
            if(se>0){
                external.iid.auc[,iO] <- iid(ls.auc$external[[iO]])
                if(iO>1){
                    iStat <- (external.auc[iO,"estimate"] - external.auc[iO-1,"estimate"]) / sqrt(crossprod(external.iid.auc[,iO]-external.iid.auc[,iO-1]))
                    external.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                }
            }

            ## Brier score
            ls.brier$external[[iO]] <- brier(labels = newdata$XXresponseXX, predictions = external.predictions[,iO], iid = external.iid.predictions[[iO]], observation = "external",
                                             null = NA, conf.level = conf.level, transformation = transformation)
            external.brier[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(ls.brier$external[[iO]]))
            if(se>0){
                external.iid.brier[,iO] <- iid(ls.brier$external[[iO]])
                if(iO>1){
                    iStat <- (external.brier[iO,"estimate"] - external.brier[iO-1,"estimate"]) / sqrt(crossprod(external.iid.brier[,iO]-external.iid.brier[,iO-1]))
                    external.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                }
            }            
        }
        
        ## export
        out$performance <- rbind(out$performance,
                                 cbind(method = "external", metric = "auc", external.auc),
                                 cbind(method = "external", metric = "brier", external.brier)
                                 )
        out$response[["external"]] <- newdata$XXresponseXX
        out$prediction[["external"]] <- external.predictions
        if(se>0){
            out$iid.auc[["external"]] <- external.iid.auc
            out$iid.brier[["external"]] <- external.iid.brier
        }
        if(trace){cat("(done)")}
    }

    ## *** cross-validation
    if(fold.repetition>0){
        if(trace){cat(" CV")}

        n.number <- length(fold.allnumber)
        name.col <- c("model","estimate","se","lower","upper","p.value","p.value_comp","foldCV.number","foldCV.size")
        cv.auc <- data.frame(matrix(NA, nrow = n.object*n.number, ncol = length(name.col), dimnames = list(NULL,name.col)))
        cv.brier <- data.frame(matrix(NA, nrow = n.object*n.number, ncol = length(name.col), dimnames = list(NULL,name.col)))
        if(se>0){
            cv.iid.auc <- matrix(NA, nrow = NROW(data), ncol = n.object, dimnames = list(NULL,names.object))
            cv.iid.brier <- matrix(NA, nrow = NROW(data), ncol = n.object, dimnames = list(NULL,names.object))
        }
        ls.auc$cv <- setNames(vector(mode = "list", length = n.object), names.object)
        ls.brier$cv <- setNames(vector(mode = "list", length = n.object), names.object)

        for(iO in 1:n.object){ ## iO <- 2
            for(iNumber in 1:n.number){ ## iNumber <- 1

                ## prepare
                iObs <- as.double(cv.indexing[,"observation",1:fold.allnumber[iNumber]])
                iRepeat <- as.double(cv.indexing[,"repetition",1:fold.allnumber[iNumber]])
                iPred <- as.double(cv.predictions[,iO,1:fold.allnumber[iNumber]])
                if(se>1){
                    iCV.iid <- cv.iid.predictions[[iO]][,,1:fold.allnumber[iNumber],drop=FALSE]
                }else{
                    iCV.iid <- NULL
                }
                iCurrent <- (iO-1)*n.number+iNumber
                iPrevious <- (iO-2)*n.number+iNumber
                ## remove folds with missing values
                iRepeat.NA <- unique(iRepeat[is.na(iPred)])
                if(length(iRepeat.NA)>0){
                    iPred <- iPred[iRepeat %in% iRepeat.NA == FALSE]
                    iObs <- iObs[iRepeat %in% iRepeat.NA == FALSE]
                    if(se>1){
                        iCV.iid <- iCV.iid[,,-iRepeat.NA,drop=FALSE]
                    }
                    iRepeat <- iRepeat[iRepeat %in% iRepeat.NA == FALSE]
                    if(length(iRepeat)==0){
                        cv.auc[iCurrent,"model"] <- names.object[iO]
                        cv.brier[iCurrent,"model"] <- names.object[iO]
                        next
                    }
                }

                ## AUC
                iAUC <- auc(labels = data$XXresponseXX, predictions = iPred, fold = iRepeat, observation = iObs,
                            add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
                cv.auc[iCurrent,setdiff(name.col,"p.value_comp")] <- cbind(model = names.object[iO], confint(iAUC), foldCV.number = fold.allnumber[iNumber], foldCV.size = sum(fold.size))
                if(iNumber==1){
                    ls.auc$cv[[iO]] <- iAUC
                    if(se>0){
                        cv.iid.auc[,iO] <- iid(iAUC)
                        if(iO>1){
                            iStat <- (cv.auc[iCurrent,"estimate"] - cv.auc[iPrevious,"estimate"]) / sqrt(crossprod(cv.iid.auc[,iO]-cv.iid.auc[,iO-1]))
                            cv.auc[iCurrent,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                        }
                    }
                }
                
                ## sqrt(crossprod(rowMeans(attr(ls.auc[[iO]],"iid")))) - confint(ls.auc[[iO]])["se"]
                ## Brier score
                iBrier <- brier(labels = data$XXresponseXX, predictions = iPred, fold = iRepeat, observation = iObs, iid = iCV.iid,
                                null = null["Brier"], conf.level = conf.level, transformation = transformation)
                cv.brier[iCurrent,setdiff(name.col,"p.value_comp")] <- cbind(model = names.object[iO], confint(iBrier), foldCV.number = fold.allnumber[iNumber], foldCV.size = sum(fold.size))
                if(iNumber==1){
                    ls.brier$cv[[iO]] <- iBrier
                    if(se>0){
                        cv.iid.brier[,iO] <- iid(iBrier)
                        if(iO>1){
                            iStat <- (cv.brier[iCurrent,"estimate"] - cv.brier[iPrevious,"estimate"]) / sqrt(crossprod(cv.iid.brier[,iO]-cv.iid.brier[,iO-1]))
                            cv.brier[iCurrent,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                        }
                    }
                }
            }
        }
        ## export
        if(!is.null(out$performance)){
            out$performance$foldCV.size <- NA
            out$performance$foldCV.number <- NA
        }
        
        out$performance <- rbind(out$performance,
                                 cbind(method = "cv", metric = "auc", cv.auc),
                                 cbind(method = "cv", metric = "brier", cv.brier)
                                 )

        if(simplify){
            out$performance$foldCV.size <- NULL
            if(n.number==1){
                out$performance$foldCV.number <- NULL
            }
        }

        out$response[["cv"]] <- data$XXresponseXX
        out$prediction[["cv"]] <- cv.predictions
        attr(out$prediction[["cv"]],"index") <- cv.indexing
        if(se>0){
            out$iid.auc[["cv"]] <- cv.iid.auc
            out$iid.brier[["cv"]] <- cv.iid.brier
        }
        if(trace){cat("(done)")}
    }

    out$auc <- ls.auc
    out$brier <- ls.brier
    if(trace){cat("\n")}
    
    ## ** export
    if(save.data){out$data <- data}
    out$args <- list(individual.fit = individual.fit,
                     impute = impute,
                     name.response = name.response,
                     fold.size = fold.size,
                     fold.repetition = fold.repetition,
                     fold.balance = fold.balance,
                     null = null,
                     conf.level = conf.level,
                     transformation = transformation,
                     auc.type = auc.type,
                     seed = seed)
    rownames(out$performance) <- NULL
    class(out) <- append("performance",class(out))
    return(out)
}


## * .performance_args
## normalize user input
.performance_args <- function(object,
                              data,
                              newdata,
                              individual.fit,
                              impute,
                              name.response,
                              fold.size,
                              fold.repetition,
                              fold.balance,
                              null,
                              conf.level,
                              se,
                              transformation,
                              auc.type,
                              simplify,
                              trace,
                              seed){
    
    ## ** object argument
    ## convert object into a list of models with names
    if(inherits(object,"ranger")){
        if(any(names(object$call)[-1]=="")){
            stop("All arguments must be named when calling ranger. \n")
        }
        if("probability" %in% names(object$call) == FALSE || object$call$probability == FALSE){
            stop("Argument \'probability\' must be set to TRUE when calling ranger. \n")
        }
        object <- list(ranger1 = object)
    }else  if(inherits(object,"randomForest")){
        object <- list(randomForest1 = object)
    }else if(inherits(object,"glm")){
        if(object$family$family %in% c("binomial","quasibinomial") == FALSE){
            stop("Cannot handle glm with family other than binomial or quasibinomial. \n")
        }
        if(object$family$link!="logit"){
            stop("Cannot handle glm with link other than logit. \n")
        }
        object <- list(logit1 = object)
    }else if(inherits(object,"miss.glm")){
        object <- list(logit1 = object)
    } else if(!is.list(object)){
        stop("Argument \'object\' should be a \'glm\' object, a \'miss.glm\' object, or a \'ranger\' object. \n")
    }else{
        possible.names <- lapply(1:length(object), function(iO){
            if(inherits(object[[iO]],"ranger")){
                if(any(names(object[[iO]]$call)[-1]=="")){
                    stop("All arguments must be named when calling ranger. \n")
                }
                if("probability" %in% names(object[[iO]]$call) == FALSE || object[[iO]]$call$probability == FALSE){
                    stop("Argument \'probability\' must be set to TRUE when calling ranger. \n")
                }
                return(paste0("ranger",iO))
            }else if(inherits(object[[iO]],"randomForest")){
                return(paste0("randomForest",iO))
            }else if(inherits(object[[iO]],"glm")){
                if(object[[iO]]$family$family %in% c("binomial","quasibinomial") == FALSE){
                    stop("Cannot handle glm with family other than binomial or quasibinomial. \n")
                }
                if(object[[iO]]$family$link!="logit"){
                    stop("Cannot handle glm with link other than logit. \n")
                }
                return(paste0("logit",iO))
            }else if(inherits(object[[iO]],"miss.glm")){
                ## nothing
            }else{
                stop("Argument \'object\' should be a list containing \'glm\', \'miss.glm\' or \'ranger\' objects. \n")
            }
        })
        if(is.null(names(object))){
            names(object) <- possible.names
        }
    }
    
    ## *** impute argument
    impute <- match.arg(impute, c("none","median","mean"))
    
    ## *** name.response argument
    ## extract name of the outcome variable
    if(is.null(name.response)){
        name.response <- all.vars(formula(object[[1]]))[1]
    }

    ## *** data argument
    if(identical(data,FALSE) || identical(data,NA)){
        internal <- FALSE
        data <- NULL
    }else if(identical(attr(data,"internal"), FALSE)){
        internal <- FALSE
    }else{
        internal <- TRUE
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
            for(iObj in 1:length(object)){
                add.cols <- ls.names[[iObj]][ls.names[[iObj]] %in% colnames(data)]
                if(length(add.cols)>1){
                    data <- cbind(data,ls.data[[iObj]][,add.cols])
                }
            }
        }
        save.data <- TRUE
    }else{
        save.data <- FALSE
    }
    data <- as.data.frame(data)
    if(name.response %in% names(data) == FALSE){
        stop("Could not find the variable ",name.response," in argument \'data\'. \n")
    }
    if(any(is.na(data[[name.response]]))){
        stop("Cannot handle missing data in the outcome variable (argument data). \n")
    }

    ## *** newdata argument
    if(identical(newdata,NA) || identical(newdata,FALSE)){
        newdata <- NULL
    }
    if(!is.null(newdata)){
        if(name.response %in% names(newdata) == FALSE){
            stop("Could not find the variable ",name.response," in argument \'newdata\'. \n")
        }
        if(any(is.na(newdata[[name.response]]))){
            stop("Cannot handle missing data in the outcome variable (argument newdata). \n")
        }
        if(any(names(data) %in% names(newdata) == FALSE)){
            stop("Argument \'newdata\' should contain column(s) \"",paste(names(data)[names(data) %in% names(newdata) == FALSE], collapse = "\" \""),"\" \n")
        }
    }


    ## *** null argument
    ## null hypothesis
    if("brier" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called brier. \n",
             "(corresponding to the null hypothesis for the brier score, possibly NA)")
    }
    if("AUC" %in% names(null) == FALSE){
        stop("Argument \'null\' should be a vector with an element called AUC. \n",
             "(corresponding to the null hypothesis for the AUC, possibly NA)")
    }

    ## *** auc.type argument
    auc.type <- match.arg(auc.type,c("classical","probabilistic"))
    if(auc.type %in% "probabilistic"){
        stop("Probabilistic AUC not yet implemented. \n",
             "Consider setting the argument \'auc.type\' to \"classical\". \n")
    }

    ## *** conf.level argument
    if(se==FALSE){
        conf.level <- NA
    }else if(is.na(conf.level)){
        se <- FALSE
    }
    

    ## *** fold.repetition arguments
    if(is.data.frame(fold.repetition)){
        if(any(names(fold.repetition) %in% c("observation","fold","repetition") == FALSE)){
            stop("When a data.frame, argument \'fold.repetition\' should contain columns \"observation\", \"fold\", \"repetition\". \n")
        }
        df.fold <- fold.repetition
        fold.repetition <- length(unique(df.fold$repetition))
        fold.size <- max(tapply(df.fold$repetition,interaction(df.fold$repetition,df.fold$fold), length))
        fold.group <- length(unique(df.fold$fold))
        
        if(any(df.fold$repetition %in% 1:fold.repetition==FALSE)){
            stop("Incorrect values for argument \'fold.repetition\' in column \"repetition\". \n",
                 "Should contain values between 1 and ",fold.repetition,". \n")
        }
        if(any(df.fold$fold %in% 1:fold.group==FALSE)){
            stop("Incorrect values for argument \'fold.repetition\' in column \"fold\". \n",
                 "Should contain values between 1 and ",fold.group,". \n")
        }
        if(any(df.fold$observation %in% 1:NROW(data)==FALSE)){
            stop("Incorrect values for argument \'fold.repetition\' in column \"observation\". \n",
                 "Should contain values between 1 and ",NROW(data),". \n")
        }
        size.tempo <- unique(tapply(df.fold$repetition,df.fold$repetition, length))
        if(length(size.tempo)>1){
            stop("Incorrect structure for argument \'fold.repetition\'. \n",
                 "Should contain the same number of lines for each repetitions. \n")
        }
        
        fold.test <- lapply(1:fold.repetition, function(iRep){
            matrix(NA, nrow = fold.group, ncol = fold.size)
        })
        for(iRep in 1:fold.repetition){
            iIndex <- df.fold[df.fold$repetition==iRep,c("observation","fold","repetition")]
            for(iFold in 1:fold.group){
                fold.test[[iRep]][iFold,1:sum(iIndex$fold == iFold)] <- iIndex[iIndex$fold == iFold,"observation"]
            }
        }
    }else{
        fold.test <- NULL
    }
    if(is.na(fold.repetition) || is.null(fold.repetition)){
        fold.repetition <- 0
    }
    if(any(fold.repetition<0)){
        stop("Argument \'fold.repetition\' must be positive \n")
    }
    if(any(fold.repetition %% 1 > 0)){
        stop("Argument \'fold.repetition\' must be an integer \n")
    }

    ## *** fold.size arguments
    if(fold.size<=0){
        stop("Argument \'fold.size\' must be strictly positive \n")
    }

    ## *** fold.repetition arguments
    if(any(fold.repetition>0)){
        fold.allnumber <- fold.repetition
        if(length(fold.repetition)>1){ ## enable to run once performance and get the results with various number of repetitions
            fold.repetition <- max(fold.allnumber)
        }

        nobs.object <- NROW(data)
        if(fold.size<1){
            fold.group <- ceiling(1/fold.size)
            fold.size <- rep(ceiling(nobs.object / fold.group),fold.group)
            if(sum(fold.size)>nobs.object){
                n.extra <- sum(fold.size)-nobs.object
                fold.size[(fold.group-n.extra+1):fold.group] <- fold.size[(fold.group-n.extra+1):fold.group]-1
            }
        }else{
            fold.group <- 1
        }
    }else{
        fold.allnumber <- 0
        fold.group <- NA
    }

    ## ** export
    args <- names(match.call()[-1])
    n.args <- length(args)
    out <- stats::setNames(vector(mode = "list", length = n.args), args)
    for(iL in 1:n.args){
        iE <- environment()[[args[iL]]]
        if(!is.null(iE)){
            out[[args[iL]]] <- iE
        }
    }
    out$args <- args
    out$internal <- internal
    out$fold.allnumber <- fold.allnumber
    out$fold.group <- fold.group
    out$fold.test <- fold.test
    out$save.data <- save.data
    return(out)
}

## * .performance_init
## initialize data and missing data patterns
.performance_init <- function(object, data, newdata,
                              individual.fit, name.response,
                              fold.size, fold.repetition, fold.balance, fold.group, fold.test, internal, trace){

    ## ** extract reference level of the outcome variable
    ref.response <- sort(data[[name.response]], decreasing = TRUE)[1]

    ## ** normalize the outcome variable in data
    if("XXresponseXX" %in% names(data)){
        stop("No column called \"XXresponseXX\" should exist in argument \'data\'. \n",
             "This name is used internally.\n")
    }

    if(is.character(data[[name.response]])){
        data$XXresponseXX <- as.factor(data[[name.response]])
    }else if(is.logical(data[[name.response]])){
        data$XXresponseXX <- as.numeric(data[[name.response]])
    }else{
        data$XXresponseXX <- data[[name.response]]
    }
    if(is.factor(data$XXresponseXX)){
        if(length(levels(data$XXresponseXX))==2){
            data$XXresponseXX <- as.numeric(data$XXresponseXX)-1
            ref.response.num <- as.numeric(ref.response)-1
        }else stop("In argument \'data\', the column corresponding to the argument \'name.response\' should take exactly two different values. \n",
                   "Unique values found: \"",paste0(levels(data$XXresponseXX), collapse = "\" \""),"\".\n")
    }else if(any(data$XXresponseXX %in% 0:1 == FALSE)){
        stop("In argument \'data\', the column corresponding to the argument \'name.response\' should correspond to a binary variable. \n",
             "Unique values found: \"",paste0(levels(data$XXresponseXX), collapse = "\" \""),"\".\n")
    }else{
        ref.response.num <- ref.response
    }

    ## ** normalize the outcome variable in newdata
    if(!is.null(newdata)){
        newdata <- as.data.frame(newdata)
        if("XXresponseXX" %in% names(newdata)){
            stop("No column called \"XXresponseXX\" should exist in argument \'newdata\'. \n",
                 "This name is used internally.\n")
        }
        if(is.character(newdata[[name.response]])){
            newdata$XXresponseXX <- as.factor(newdata[[name.response]])
        }else if(is.logical(newdata[[name.response]])){
            newdata$XXresponseXX <- as.numeric(newdata[[name.response]])
        }else{
            newdata$XXresponseXX <- newdata[[name.response]]
        }
        if(is.factor(newdata$XXresponseXX)){
            if(length(levels(newdata$XXresponseXX))==2){
                newdata$XXresponseXX <- as.numeric(newdata$XXresponseXX)-1
            }else stop("In argument \'newdata\', the column corresponding to the argument \'name.response\' should take exactly two different values. \n",
                       "Unique values found: \"",paste0(levels(newdata$XXresponseXX), collapse = "\" \""),"\".\n")
        }else if(any(newdata$XXresponseXX %in% 0:1 == FALSE)){
            stop("In argument \'newdata\', the column corresponding to the argument \'name.response\' should correspond to a binary variable. \n",
                 "Unique values found: \"",paste0(levels(newdata$XXresponseXX), collapse = "\" \""),"\".\n")
        }
    }


    ## ** extract sample size of the training set
    nobs.object <- unname(unique(sapply(object,function(iO){ ## iO <- object$RF
        if(inherits(iO,"glm")){
            return(stats::nobs(iO))
        }else if(inherits(iO,"ranger")){
            return(iO$num.samples)
        }else if(inherits(iO,"randomForest")){
            return(NROW(iO$votes))
        }
    })))
    if(length(nobs.object)!=1){
        if(trace){
            message("The training set seems to differ in size between models: ",paste0(nobs.object, collapse = ", "),". \n")
        }
        nobs.object <- max(nobs.object)
    }

    ## ** define missing data patterns
    if(individual.fit){  ## and formula/index of observations for each missing data pattern
        object.formula <- lapply(object,function(iO){
            ff <- try(formula(iO), silent = TRUE)
            if(inherits(ff,"try-error")){
                return(eval(iO$call[[2]]))
            }else{
                return(ff)
            }
        })
        object.xvar <- lapply(object.formula, function(iF){all.vars(stats::delete.response(stats::terms(iF)))})

        ## data
        nobs.object <- NROW(data)
        object.iformula <- setNames(vector(mode = "list", length = length(object)), names(object))
        data.missingPattern <- setNames(vector(mode = "list", length = length(object)), names(object))

        for(iO in 1:length(object)){ ## iO <- 1
            iTest.na <- is.na(data[,object.xvar[[iO]],drop=FALSE])
            ## fields::image.plot(iTest.na)
            object.iformula[[iO]] <- apply(iTest.na,1,function(iRow){
                if(any(iRow)){
                    iNewFF <- as.formula(paste(".~.-",paste(colnames(iTest.na)[iRow],collapse="-")))
                    return(stats::update(object.formula[[iO]],iNewFF))
                }else{
                    return(object.formula[[iO]])
                }
            })
            data.missingPattern[[iO]] <- interaction(as.data.frame(1*iTest.na),drop=TRUE)
            attr(data.missingPattern[[iO]],"index") <- tapply(1:length(data.missingPattern[[iO]]),data.missingPattern[[iO]],list)
            attr(data.missingPattern[[iO]],"formula") <- lapply(attr(data.missingPattern[[iO]],"index"),function(iVec){object.iformula[[iO]][[iVec[1]]]})
        }

        ## newdata
        if(!is.null(newdata)){
            newdata.missingPattern <- setNames(vector(mode = "list", length = length(object)), names(object))
            if(internal){
                newdata2 <- rbind(data,newdata)
            }else{
                newdata2 <- newdata
            }

            for(iO in 1:length(object)){ ## iO <- 1
                iTest.na <- is.na(newdata2[,object.xvar[[iO]],drop=FALSE])
                newdata.missingPattern[[iO]] <- interaction(as.data.frame(1*iTest.na),drop=TRUE)
                attr(newdata.missingPattern[[iO]],"index") <- tapply(1:length(newdata.missingPattern[[iO]]),newdata.missingPattern[[iO]],list)
                attr(newdata.missingPattern[[iO]],"formula") <- lapply(attr(newdata.missingPattern[[iO]],"index"),function(iObs){ ## iObs <- 3
                    iVar.rm <- colnames(iTest.na)[iTest.na[iObs[1],]]
                    if(length(iVar.rm)>0){
                        return(stats::update(object.formula[[iO]],as.formula(paste(".~.-",paste(iVar.rm,collapse="-")))))
                    }else{
                        return(object.formula[[iO]])
                    }
                })
            }
        }else if(internal){
            newdata.missingPattern <- data.missingPattern
        }else{
            newdata.missingPattern <- NULL
        }

    }else{
        data.missingPattern <- NULL
        newdata.missingPattern <- NULL        
    }

    ## ** prepare folds for cross validation
    if(!is.null(fold.test)){
        ## nothing to do        
    }else if(fold.repetition>0){
        index.response1 <- which(data$XXresponseXX==ref.response.num)
        index.response0 <- which(data$XXresponseXX!=ref.response.num)
        prevalence <- length(index.response1)/length(index.response0)

        if(fold.balance){ ## position of sampled observations in each fold (such that prevalence is preserved)
            index.sample0 <- .balanceFold(n.obs = length(index.response0), n.fold = fold.group)
            index.sample1 <- .balanceFold(n.obs = length(index.response1), n.fold = fold.group)
            fold.size <- sapply(index.sample0,length)+sapply(index.sample1,length)
        }
        fold.test <- lapply(1:fold.repetition, function(iRepeat){ ## iRepeat <- 1
            iM <- matrix(NA, nrow = fold.group, ncol = max(fold.size))
            if(fold.balance){
                iSample0 <- sample(index.response0, replace = FALSE)
                iSample1 <- sample(index.response1, replace = FALSE)
                for(iFold in 1:fold.group){ ## iFold <- 1
                    iM[iFold,1:(length(index.sample0[[iFold]])+length(index.sample1[[iFold]]))] <- c(iSample0[index.sample0[[iFold]]], iSample1[index.sample1[[iFold]]])
                }
            }else{
                ## make sure there is at least one 0 and one 1 in the test dataset
                iSample0 <- sample(index.response0, replace = FALSE, size = fold.group)
                iSample1 <- sample(index.response1, replace = FALSE, size = fold.group)
                ## sample remaining observations
                iSample01 <- sample(setdiff(1:NROW(data),c(iSample0,iSample1)), replace = FALSE, size = sum(fold.size)-2*fold.group)
                ## collect into a matrix
                iM[,1] <- iSample0
                iM[,2] <- iSample1
                if(any(fold.size>2)){
                    iM[,-(1:2)] <- c(iSample01,rep(NA,length(iM)-sum(fold.size)))
                }
            }
            return(iM)
        })

    }else{
        fold.size <- 0
        fold.test <- NULL
    }

    ## ** export
    out <- list(data = data,
                newdata = newdata,
                ref.response = ref.response,
                ref.response.num = ref.response.num,
                nobs.object = nobs.object,
                data.missingPattern = data.missingPattern,
                newdata.missingPattern = newdata.missingPattern,
                fold.size = fold.size,
                fold.test = fold.test)
}


## * .performance_predict
##' @description Compute prediction with uncertainty for various type of models
##' @param object model from which prediction should be evaluated.
##' @param n.obs [integer] number of observations in the training set.
##' @param newdata [data.frame] test set.
##' @param se [logical] should the uncertainty (i.e. standard error, influence function)
##' associated to the predictions be extracted when possible?
##' @noRd
.performance_predict <- function(object, n.obs, newdata, se){

    ## ** prepare output
    out <- list(estimate = NULL)
    out$se <- rep(NA, ncol = NROW(newdata))
    out$iid <- matrix(0, nrow = n.obs, ncol = NROW(newdata))

    ## ** predictions
    if(inherits(object,"ranger")){
        if(se>0){
            iPred <- predict(object, data = newdata, type = "se")
            out$estimate <- iPred$predictions[,which(object$forest$class.values==1)]
            out$se <- iPred$se
        }else{
            iPred <- predict(object, data = newdata, type = "response")
            out$estimate <- iPred$predictions[,which(object$forest$class.values==1)]
        }
    }else if(inherits(object,"randomForest")){
        out$estimate <- stats::predict(object, newdata = newdata, type = "prob")[,2]
    }else if(inherits(object,"glm")){
        if(object$converged==FALSE){return(NULL)}
        if(se>0){
            iPred <- .predict.logit(object, newdata = newdata)
            out$estimate <- iPred["estimate",]
            out$se <- iPred["se",]
            out$iid <- attr(iPred,"iid")
        }else{
            out$estimate <- stats::predict(object, newdata = newdata, type = "response")
        }
    }else if(inherits(object,"miss.glm")){
        Xb <- stats::model.matrix(stats::formula(object), newdata) %*% stats::coef(object)
        out$estimate <- as.double(1/(1+exp(-Xb)))
    }

    ## ** export
    return(out)

}

## * .performance_runfold
.performance_runfold <- function(object, names.object, n.object,
                                 data.train, n.train, impute,
                                 data.test, n.test, missingPattern.test,
                                 individual.fit, auc.type, se, trace){


    ## ** prepare output
    out <- list(estimate = matrix(NA, nrow = n.test, ncol = n.object,
                                  dimnames = list(NULL, names.object)))
    if(auc.type == "probabilistic"){
        out$se <- matrix(NA, nrow = n.test, ncol = n.object,
                         dimnames = list(NULL, names.object))
    }
    if(se>1){
        out$iid <- setNames(vector(mode = "list", length = n.object), names.object)
    }

    ## ** imputation
    if(any(is.na(data.train)) && impute != "none"){
        data.train0 <- data.train
        col.impute <- names(which(colSums(is.na(data.train))>0))
        for(iCol in col.impute){ ## iCol <- "X1"
            if(is.numeric(data.train[[iCol]])){
                data.train[is.na(data.train[[iCol]]),iCol] <- do.call(impute, args = list(data.train[[iCol]], na.rm = TRUE))
            }else{
                data.train[is.na(data.train[[iCol]]),iCol] <- names(sort(table(data.train[[iCol]], useNA = "no"), decreasing = TRUE))[1]
            }
        }
    }


    ## ** evaluate predictions over models
    for(iO in 1:n.object){ ## iO <- 3
        
        if(individual.fit){
            ## *** Missing data
            ## Refit model for each observation depending on the available predictors
            ## Median imputation for the predictors

            iIndex.pattern <- attr(missingPattern.test[[iO]],"index")
            iFormula.pattern <- attr(missingPattern.test[[iO]],"formula")
            iPred.pattern <- attr(missingPattern.test[[iO]],"prediction") ## previously fit model
            if(se>1){
                out$iid[[iO]] <- matrix(0, nrow = n.train, ncol = n.test)
            }

            ## For each set of non-missing predictors
            for(iPattern in 1:length(iIndex.pattern)){ ## iPattern <- 1
                iName.pattern <- names(iIndex.pattern)[iPattern]

                ## update model and compute predictions
                iData.test <- data.test[iIndex.pattern[[iPattern]],all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                    
                if(inherits(object[[iO]],"miss.glm")){
                    iData.train <- data.train0[,all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                    iIndex.train <- 1:NROW(iData.train)
                }else{
                    iData.train <- data.train[,all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                    iIndex.train <- which(rowSums(is.na(iData.train)) == 0)
                }
                iObject <- stats::update(object[[iO]], formula = iFormula.pattern[[iPattern]], data = iData.train[iIndex.train,,drop=FALSE])
                iPred <- .performance_predict(iObject, n.obs = length(iIndex.train), newdata = iData.test, se = se>1)
                
                ## store results
                if(!is.null(iPred)){ ## convergence check
                    out$estimate[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$estimate)
                    if(auc.type == "probabilistic"){
                        out$se[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$se)
                    }
                    if(se>1){
                        out$iid[[iO]][iIndex.pattern[[iPattern]],iIndex.train] <- iPred$iid
                    }
                }
            }
        }else{
            ## *** Complete case
            if(inherits(object[[iO]],"miss.glm")){
                iData.train <- data.train0                
            }else{
                iData.train <- data.train
            }
            if(!is.null(data.train)){
                iPred <- .performance_predict(stats::update(object[[iO]], data = data.train),
                                              n.obs = n.train, newdata = data.test, se = se>1)
            }else{
                iPred <- .performance_predict(object[[iO]],
                                              n.obs = n.train, newdata = data.test, se = se>1)
            }
            if(!is.null(iPred)){ ## convergence check
                out$estimate[,iO] <- as.double(iPred$estimate)
                if(auc.type == "probabilistic"){
                    out$se[,iO] <- as.double(iPred$se)
                }
                if(se>1){
                    out$iid[[iO]] <- iPred$iid
                }
            }
        }
    }
    if(trace){cat(" ")}


    ## ** export
    return(out)
}

## * .balanceFold
##' @description Generate indexes relative to each fold.
##' @param n.obs number of observations.
##' @param n.fold number of folds.
##' @noRd
.balanceFold <- function(n.obs, n.fold){
    if(n.obs < n.fold){
        return(c(lapply(1:(n.fold-n.obs), function(i){NULL}),as.list(1:n.obs)))
    }

    ## number of observations per fold
    nobs.fold <- rep(floor(n.obs/n.fold),n.fold)
    nobs.fold <- nobs.fold + c(rep(1,n.obs-sum(nobs.fold)),rep(0,n.fold-(n.obs-sum(nobs.fold)))) 
    ## list of indexes
    out <- mapply(x = cumsum(c(0,nobs.fold[-n.fold]-1)+1), y = cumsum(nobs.fold), function(x,y){list(x:y)})

    return(out)
}
##----------------------------------------------------------------------
### performance.R ends here
