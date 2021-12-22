### performance.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  3 2021 (11:17) 
## Version: 
## Last-Updated: Dec 21 2021 (16:49) 
##           By: Brice Ozenne
##     Update #: 518
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Assess Performance of a Classifier (EXPERIMENTAL)
##' @description Assess the performance in term of AUC and brier score of one or several binary classifiers.
##' Currently limited to logistic regressions and random forest.
##'
##' @param object a \code{glm} or \code{range} object, or a list of such object.
##' @param data [data.frame] the training data.
##' @param newdata [data.frame] an external data used to assess the performance.
##' @param fold.size [double, >0] Either the size of the test dataset (when >1) or the fraction of the dataset (when <1) to be used for testing when using cross-validation.
##' @param fold.number [integer] When strictly positive, the number of fold used in the cross-validation. If 0 then no cross validation is performed.
##' @param individual.fit [logical] If \code{TRUE} the predictive model is refit for each individual using only the predictors with non missing values.
##' @param name.response [character] The name of the response variable (i.e. the one containing the categories).
##' @param null [numeric vector of length 2] the right-hand side of the null hypothesis relative to each metric.
##' @param conf.level [numeric] confidence level for the confidence intervals.
##' @param auc.type [character] should the auc be computed approximating the predicted probability by a dirac (\code{"classical"}, usual AUC formula)
##' or approximating the predicted probability by a normal distribution.
##' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
##' Otherwise they are computed without any transformation.
##' @param trace [logical] Should the execution of the function be traced.
##' @param simplify [logical] Should the number of fold and the size of the fold used for the cross validation be removed from the output?
##'
##' @details WARNING: this function is still in development. In particular standard errors, confidence intervals, and p-values should not be trusted.
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
##' e.logit <- glm(Y~X1+X2, data = df.train, family = binomial(link="logit"))
##'
##' ## assess performance on the training set (biased)
##' ## and external dataset
##' BuyseTest::performance(e.logit, newdata = df.test)
##' BuyseTest::performance(list(null = e.null, prop = e.logit), newdata = df.test)
##' 
##' ## assess performance using cross validation
##' set.seed(10)
##' BuyseTest::performance(e.logit, fold.number = 10)
##' set.seed(10)
##' \dontrun{
##' BuyseTest::performance(list(null = e.null, prop = e.logit), fold.number = 10)
##' BuyseTest::performance(e.logit, fold.number = c(50,20,10))
##' }

## * performance
##' @export
performance <- function(object, data = NULL, newdata = NA, fold.size = 1/10, fold.number = 0,
                        individual.fit = FALSE, name.response = NULL,
                        null = c(brier = NA, AUC = 0.5), conf.level = 0.95, transformation = TRUE, auc.type = "classical",
                        simplify = TRUE, trace = TRUE){

    ## ** normalize user input
    ## convert object into a list of models with names
    if(inherits(object,"ranger")){
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
    }else if(!is.list(object)){
        stop("Argument \'object\' should be a \'glm\' object, a \'ranger\' object, or a list containing \'glm\' or \'ranger\' objects. \n")
    }else{
        possible.names <- lapply(1:length(object), function(iO){
            if(inherits(object[[iO]],"ranger")){
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
            }else{
                stop("Argument \'object\' should be a \'glm\' object, a \'ranger\' object, or a list containing \'glm\' or \'ranger\' objects. \n")
            }
        })
        if(is.null(names(object))){
            names(object) <- possible.names
        }
    }

    if(individual.fit){
        for(iO in 1:length(object)){
            test.formula <- try(stats::formula(object[[iO]]),silent = TRUE)
            if(inherits(test.formula,"try-error")){
                stop("Could not extract the formula from one of the model in argument \'object\'. \n",
                     "Consider explicitely naming the argument \'formula\' when fitting the model. \n",
                     "Error when calling stats::formula: ",test.formula[1],"\n")
            }
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
            for(iObj in 1:length(object)){
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
    ## missing data in outcome
    if(any(is.na(data[[name.response]]))){
        stop("Cannot handle missing data in the outcome variable. \n")
    }
    ## extract reference level of the outcome variable
    ref.response <- sort(data[[name.response]], decreasing = TRUE)[1]

    ## normalize the outcome variable in data
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

    ## normalize the outcome variable in newdata
    if(!is.null(newdata) && !identical(newdata,NA) && !identical(newdata,FALSE)){
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

            for(iO in 1:length(object)){ ## iO <- 1
                iTest.na <- is.na(newdata[,object.xvar[[iO]],drop=FALSE])
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
                attr(newdata.missingPattern[[iO]],"prediction") <- list()
                 
            }
        }

    }else{
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
            stop("The training set seems to differ in size between models. \n")
        }
    }

    ## arguments for cross validation
    fold.allnumber <- fold.number
    if(length(fold.number)>1){
        fold.number <- max(fold.allnumber)
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
            fold.group <- ceiling(1/fold.size)
            fold.size <- rep(ceiling(nobs.object / fold.group),fold.group)
            if(sum(fold.size)>nobs.object){
                n.extra <- sum(fold.size)-nobs.object
                fold.size[(fold.group-n.extra+1):fold.group] <- fold.size[(fold.group-n.extra+1):fold.group]-1
            }
        }else{
            fold.group <- 1
        }
        
        if(any(table(data$XXresponseXX)<fold.group)){
            stop("Argument \'fold.size\' too large as there are less than ",fold.group," observations with outcome ",names(table(data$XXresponseXX))[which.min(table(data$XXresponseXX))],"\n")
        }
        if(any(fold.size<2)){
            stop("Argument \'fold.size\' must be greater or equal to 2 (at least one sample in each category) \n")
        }
    }else if(!identical(data,NA) && !identical(data,FALSE)){
        data <- as.data.frame(data)
    }
    auc.type <- match.arg(auc.type,c("classical","probabilistic"))
    if(auc.type %in% "probabilistic"){
        stop("Probabilistic AUC not yet implemented. \n",
             "Consider setting the argument \'auc.type\' to \"classical\". \n")
    }
    
    names.object <- names(object)
    n.object <- length(object)
    out <- NULL
    
    ## ** internal performance
    if(!identical(data,NA) && !identical(data,FALSE)){
        if(trace){cat("- assess internal performance: ")}

        ## *** get predictions
        internal.predictions <- matrix(NA, nrow = nobs.object, ncol = n.object,
                                       dimnames = list(NULL, names.object))
        if(auc.type == "probabilistic"){
            internal.se.predictions <- matrix(NA, nrow = nobs.object, ncol = n.object,
                                              dimnames = list(NULL, names.object))
        }
        internal.iid <- setNames(vector(mode = "list", length = n.object), names.object)
        
        for(iO in 1:n.object){ ## iO <- 3
            if(trace){cat(".")}
            if(individual.fit){
                internal.iid[[iO]] <- matrix(0, nrow = nobs.object, ncol = nobs.object)
                iIndex.pattern <- attr(data.missingPattern[[iO]],"index")
                iFormula.pattern <- attr(data.missingPattern[[iO]],"formula")
                
                for(iPattern in 1:length(iIndex.pattern)){ ## iPattern <- 1
                    iName.pattern <- names(iIndex.pattern)[iPattern]
                    iData <- data[,all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                    iIndex <- which(rowSums(is.na(iData)) == 0)
                    iObject <- stats::update(object[[iO]], formula = iFormula.pattern[[iPattern]], data = iData[iIndex,,drop=FALSE])
                    iPred <- .performance_predict(iObject, n.obs = length(iIndex), newdata = data[iIndex.pattern[[iPattern]],,drop=FALSE], auc.type = auc.type)

                    if(!is.null(iPred)){ ## convergence check
                        internal.predictions[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$estimate)
                        if(auc.type == "probabilistic"){
                            internal.se.predictions[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$se)
                        }
                        internal.iid[[iO]][iIndex.pattern[[iPattern]],iIndex] <- iPred$iid
                    }

                    if(!is.null(newdata) && iName.pattern %in% levels(newdata.missingPattern[[iO]])){ ## re-use the object for prediction on other dataset
                        iNewdata <- newdata[attr(newdata.missingPattern[[iO]],"index")[[iName.pattern]],,drop=FALSE]
                        attr(newdata.missingPattern[[iO]],"prediction")[[iName.pattern]] <- .performance_predict(iObject, n.obs = NROW(iData), newdata = iNewdata, auc.type = auc.type)
                    }
                    
                }
            }else{
                iPred <- .performance_predict(object[[iO]], n.obs = nobs.object, newdata = data, auc.type = auc.type)
                if(!is.null(iPred)){ ## convergence check
                    internal.predictions[,iO] <- as.double(iPred$estimate)
                    if(auc.type == "probabilistic"){
                        internal.se.predictions[,iO] <- as.double(iPred$se)
                    }
                    internal.iid[[iO]] <- iPred$iid
                }
            }
        }
        if(trace){cat(" ")}

        ## *** assess performance
        internal.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 6, dimnames = list(names.object,c("model","estimate","se","lower","upper","p.value"))))
        internal.iid.auc <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        internal.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 6, dimnames = list(names.object,c("model","estimate","se","lower","upper","p.value"))))
        internal.iid.brier <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        if(n.object>1){
            internal.auc <- cbind(internal.auc, p.value_comp = NA)
            internal.brier <- cbind(internal.brier, p.value_comp = NA)
        }
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iAUC <- auc(labels = data$XXresponseXX, predictions = internal.predictions[,iO],
                        add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            internal.auc[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO], confint(iAUC))
            internal.iid.auc[,iO] <- iid(iAUC)
            if(iO>1){
                iStat <- (internal.auc[iO,"estimate"] - internal.auc[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.auc[,iO]-internal.iid.auc[,iO-1]))
                internal.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            iBrier <- brier(labels = data$XXresponseXX, predictions = internal.predictions[,iO], iid = internal.iid[[iO]],
                            null = null["Brier"], conf.level = conf.level, transformation = transformation)
            internal.brier[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(iBrier))
            internal.iid.brier[,iO] <- iid(iBrier)
            if(iO>1){
                iStat <- (internal.brier[iO,"estimate"] - internal.brier[iO-1,"estimate"]) / sqrt(crossprod(internal.iid.brier[,iO]-internal.iid.brier[,iO-1]))
                internal.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            
        }
        
        ## *** export
        out <- rbind(out,
                     cbind(method = "internal", metric = "auc", internal.auc),
                     cbind(method = "internal", metric = "brier", internal.brier)
                     )
        if(is.null(attr(out,"response"))){attr(out,"response") <- list()}
        attr(out,"response")[["internal"]] <- data$XXresponseXX
        if(is.null(attr(out,"prediction"))){attr(out,"prediction") <- list()}
        attr(out,"prediction")[["internal"]] <- internal.predictions
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["internal"]] <- internal.iid.auc
        attr(out,"iid")$brier[["internal"]] <- internal.iid.brier
        if(trace){cat(" done. \n")}
    }

    ## ** external performance
    if(!is.null(newdata) && !identical(newdata,NA) && !identical(newdata,FALSE)){
        if(trace){cat("- assess external performance: ")}
        nNewdata.obs <- NROW(newdata)

        ## *** get predictions
        external.predictions <- matrix(NA, nrow = nNewdata.obs, ncol = n.object,
                                       dimnames = list(NULL, names.object))
        if(auc.type == "probabilistic"){
            external.se.predictions <- matrix(NA, nrow = nNewdata.obs, ncol = n.object,
                                              dimnames = list(NULL, names.object))
        }
        external.iid <- setNames(vector(mode = "list", length = n.object), names.object)
        
        for(iO in 1:n.object){
            if(trace){cat(".")}

            if(individual.fit){
                external.iid[[iO]] <- matrix(0, nrow = nobs.object, ncol = nNewdata.obs)
                iIndex.pattern <- attr(newdata.missingPattern[[iO]],"index")
                iFormula.pattern <- attr(newdata.missingPattern[[iO]],"formula")
                iPred.pattern <- attr(newdata.missingPattern[[iO]],"prediction")

                for(iPattern in 1:length(iIndex.pattern)){ ## iPattern <- 1
                    iName.pattern <- names(iIndex.pattern)[iPattern]
                    if(is.null(iPred.pattern[[iName.pattern]])){ ## new pattern
                        iData <- data[,all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                        iIndex2.pattern <- which(rowSums(is.na(iData))==0)
                        iObject <- stats::update(object[[iO]], formula = iFormula.pattern[[iPattern]], data = iData[iIndex2.pattern,,drop=FALSE])

                        iNewdata <- newdata[iIndex.pattern[[iPattern]],all.vars(iFormula.pattern[[iPattern]]),drop=FALSE]
                        iPred <- .performance_predict(iObject, n.obs = length(iIndex2.pattern), newdata = iNewdata, auc.type = auc.type)
                    }else{ ## already seen pattern
                        iIndex2.pattern <- attr(data.missingPattern[[iO]],"index")[[iName.pattern]]
                        iPred <- iPred.pattern[[iName.pattern]]
                    }
                    if(!is.null(iPred)){ ## convergence check
                        external.predictions[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$estimate)
                        if(auc.type == "probabilistic"){
                            external.se.predictions[iIndex.pattern[[iPattern]],iO] <- as.double(iPred$se)
                        }
                        external.iid[[iO]][iIndex.pattern[[iPattern]],iIndex2.pattern] <- iPred$iid
                    }
                }
            }else{
                iPred <- .performance_predict(object[[iO]], n.obs = nobs.object, newdata = newdata, auc.type = auc.type)
                if(!is.null(iPred)){ ## convergence check
                    external.predictions[,iO] <- as.double(iPred$estimate)
                    if(auc.type == "probabilistic"){
                        external.se.predictions[,iO] <- as.double(iPred$se)
                    }
                    external.iid[[iO]] <- iPred$iid
                }
            }
            
        }
        if(trace){cat(" ")}

        ## *** assess performance
        external.auc <- data.frame(matrix(NA, nrow = n.object, ncol = 6, dimnames = list(NULL,c("model","estimate","se","lower","upper","p.value"))))
        external.iid.auc <- matrix(NA, nrow = nNewdata.obs, ncol = n.object, dimnames = list(NULL,names.object))
        external.brier <- data.frame(matrix(NA, nrow = n.object, ncol = 6, dimnames = list(NULL,c("model","estimate","se","lower","upper","p.value"))))
        external.iid.brier <- matrix(NA, nrow = nobs.object + nNewdata.obs, ncol = n.object, dimnames = list(NULL,names.object))
        if(n.object>1){
            external.auc <- cbind(external.auc, p.value_comp = NA)
            external.brier <- cbind(external.brier, p.value_comp = NA)
        }
        for(iO in 1:n.object){
            if(trace){cat("*")}
            iBrier <- brier(labels = newdata$XXresponseXX, predictions = external.predictions[,iO], iid = external.iid[[iO]], observation = "external",
                            null = NA, conf.level = conf.level, transformation = transformation)
            external.brier[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(iBrier))
            external.iid.brier[,iO] <- iid(iBrier)
            if(iO>1){
                iStat <- (external.brier[iO,"estimate"] - external.brier[iO-1,"estimate"]) / sqrt(crossprod(external.iid.brier[,iO]-external.iid.brier[,iO-1]))
                external.brier[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

            iAUC <- auc(labels = newdata$XXresponseXX, predictions = external.predictions[,iO],
                                   add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
            external.auc[iO,c("model","estimate","se","lower","upper","p.value")] <- cbind(model = names.object[iO],confint(iAUC))
            external.iid.auc[,iO] <- iid(iAUC)
            if(iO>1){
                iStat <- (external.auc[iO,"estimate"] - external.auc[iO-1,"estimate"]) / sqrt(crossprod(external.iid.auc[,iO]-external.iid.auc[,iO-1]))
                external.auc[iO,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
            }

        }
        
        ## *** export
        out <- rbind(out,
                     cbind(method = "external", metric = "auc", external.auc),
                     cbind(method = "external", metric = "brier", external.brier)
                     )
        if(is.null(attr(out,"response"))){attr(out,"response") <- list()}
        attr(out,"response")[["external"]] <- newdata$XXresponseXX
        if(is.null(attr(out,"prediction"))){attr(out,"prediction") <- list()}
        attr(out,"prediction")[["external"]] <- external.predictions
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["external"]] <- external.iid.auc
        attr(out,"iid")$brier[["external"]] <- external.iid.brier
        if(trace){cat(" done. \n")}
    }

    ## ** CV performance
    if(fold.number>0){
        if(trace){
            if(fold.group>1){                
                cat("- assess performance using ",fold.group," folds cross-validation repeated ",fold.number," times: \n",sep="")
            }else{
                cat("- assess performance using 1 fold cross-validation with ",fold.size," samples repeated ",fold.number," times: \n",sep="")
            }
        }
        out$foldCV.number <- as.numeric(NA)
        out$foldCV.size <- as.numeric(NA)
            
        ## *** identify folds
        index.response1 <- which(data$XXresponseXX==ref.response.num)
        index.response0 <- which(data$XXresponseXX!=ref.response.num)
        fold.test <- lapply(1:fold.number, function(iFold){
            ## make sure there is at least one 0 and one 1 in the test dataset
            iSample0 <- sample(index.response0, replace = FALSE, size = fold.group)
            iSample1 <- sample(index.response1, replace = FALSE, size = fold.group)
            ## sample remaining observations
            iSample01 <- sample(setdiff(1:nobs.object,c(iSample0,iSample1)), replace = FALSE, size = sum(fold.size)-2*fold.group)
            ## collect into a matrix
            iM <- matrix(NA, nrow = fold.group, ncol = max(fold.size))
            iM[,1] <- iSample0
            iM[,2] <- iSample1
            if(any(fold.size>2)){
                iM[,-(1:2)] <- c(iSample01,rep(NA,length(iM)-sum(fold.size)))
            }
            return(iM)
        })

        ## *** get predictions
        cv.indexing <- array(NA, dim = c(sum(fold.size), 2, fold.number),
                             dimnames = list(NULL, c("observation","fold"), NULL))
        cv.predictions <- array(NA, dim = c(sum(fold.size), n.object, fold.number),
                                dimnames = list(NULL, names.object, NULL))
        if(auc.type == "probabilistic"){
            cv.se.predictions <- array(NA, dim = c(sum(fold.size), n.object, fold.number),
                                       dimnames = list(NULL, names.object, NULL))
        }
        cv.iid <- setNames(lapply(1:n.object, function(iFold){
            array(0, dim = c(nobs.object, sum(fold.size), fold.number))
        }), names.object)

        if(trace){
            pb <- utils::txtProgressBar(max = fold.number, style = 3)
        }
        for(iFold in 1:fold.number){ ## iFold <- 1
            if(trace){
                utils::setTxtProgressBar(pb, iFold)
            }
            cv.indexing[,"fold",iFold] <- iFold
            
            for(iSubFold in 1:fold.group){ ## iFold <- 1
                iFoldTest <- na.omit(fold.test[[iFold]][iSubFold,])
                iFoldTrain <- setdiff(1:nobs.object,iFoldTest)

                iDataTrain <- data[iFoldTrain,,drop=FALSE]
                iDataTest <- data[iFoldTest,,drop=FALSE]

                index.iFoldTest <- (1+sum(c(0,fold.size)[1:iSubFold])):sum(fold.size[1:iSubFold])
                cv.indexing[index.iFoldTest,"observation",iFold] <- iFoldTest
                
                for(iO in 1:n.object){ ## iO <- 1

                    if(individual.fit){
                        for(iObs in 1:fold.size){ ## iObs <- 3
                            iObject <- stats::update(object[[iO]], formula = object.iformula[[iO]][[iFoldTest[iObs]]], data = iDataTrain)
                            iPerf <- .performance_predict(iObject, n.obs = nobs.object-fold.size[iSubFold], newdata = iDataTest[iObs,], auc.type = auc.type)
                            if(!is.null(iPerf)){ ## convergence check
                                cv.predictions[index.iFoldTest[iObs],iO,iFold] <- as.double(iPerf$estimate)
                                if(auc.type == "probabilistic"){
                                    cv.se.predictions[index.iFoldTest[iObs],iO,iFold] <- as.double(iPerf$se)
                                }
                                cv.iid[[iO]][iFoldTrain[,iFold][setdiff(1:NROW(iDataTrain),iObject$na.action)],index.iFoldTest[iObs],iFold] <- iPerf$iid
                            }
                        }
                    }else{
                        iObject <- stats::update(object[[iO]], data = iDataTrain)
                        iPred <- .performance_predict(iObject, n.obs = nobs.object-fold.size[iSubFold], newdata = iDataTest, auc.type = auc.type)
                        if(!is.null(iPred)){ ## convergence check
                            cv.predictions[index.iFoldTest,iO,iFold] <- as.double(iPred$estimate)
                            if(auc.type == "probabilistic"){
                                cv.se.predictions[index.iFoldTest,iO,iFold] <- as.double(iPred$se)
                            }
                            cv.iid[[iO]][iFoldTrain,index.iFoldTest,iFold] <- iPred$iid
                        }
                    }
                }
            }
        }
        if(trace){close(pb)}

        ## *** assess predictions
        n.number <- length(fold.allnumber)
        name.col <- c("model","estimate","se","lower","upper","p.value","foldCV.number","foldCV.size")
        cv.auc <- data.frame(matrix(NA, nrow = n.object*n.number, ncol = length(name.col), dimnames = list(NULL,name.col)))
        cv.iid.auc <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        cv.brier <- data.frame(matrix(NA, nrow = n.object*n.number, ncol = length(name.col), dimnames = list(NULL,name.col)))
        cv.iid.brier <- matrix(NA, nrow = nobs.object, ncol = n.object, dimnames = list(NULL,names.object))
        if(n.object>1){
            cv.auc <- cbind(cv.auc, p.value_comp = NA)
            cv.brier <- cbind(cv.brier, p.value_comp = NA)
        }
        ls.auc <- setNames(vector(mode = "list", length = n.object), names.object)
        ls.brier <- setNames(vector(mode = "list", length = n.object), names.object)
        
        for(iO in 1:n.object){ ## iO <- 2
            if(trace){cat("*")}
            for(iNumber in 1:n.number){ ## iNumber <- 1

                ## prepare
                iObs <- as.double(cv.indexing[,"observation",1:fold.allnumber[iNumber]])
                iFold <- as.double(cv.indexing[,"fold",1:fold.allnumber[iNumber]])
                iPred <- as.double(cv.predictions[,iO,1:fold.allnumber[iNumber]])
                iCV.iid <- cv.iid[[iO]][,,1:fold.allnumber[iNumber],drop=FALSE]
                iCurrent <- (iO-1)*n.number+iNumber
                iPrevious <- (iO-2)*n.number+iNumber

                ## remove folds with missing values
                iFold.NA <- unique(iFold[is.na(iPred)])
                if(length(iFold.NA)>0){
                    iPred <- iPred[iFold %in% iFold.NA == FALSE]
                    iObs <- iObs[iFold %in% iFold.NA == FALSE]
                    iCV.iid <- iCV.iid[,,-iFold.NA,drop=FALSE]
                    iFold <- iFold[iFold %in% iFold.NA == FALSE]
                }
                
                ## auc
                iAUC <- auc(labels = data$XXresponseXX, predictions = iPred, fold = iFold, observation = iObs,
                            add.halfNeutral = TRUE, null = null["AUC"], conf.level = conf.level, transformation = transformation)
                cv.auc[iCurrent,name.col] <- cbind(model = names.object[iO], confint(iAUC), foldCV.number = fold.allnumber[iNumber], foldCV.size = sum(fold.size))
                if(iNumber==1){
                    ls.auc[[iO]] <- iAUC
                    cv.iid.auc[,iO] <- iid(iAUC)
                    if(iO>1){
                        iStat <- (cv.auc[iCurrent,"estimate"] - cv.auc[iPrevious,"estimate"]) / sqrt(crossprod(cv.iid.auc[,iO]-cv.iid.auc[,iO-1]))
                        cv.auc[iCurrent,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                    }
                }
                
                ## sqrt(crossprod(rowMeans(attr(ls.auc[[iO]],"iid")))) - confint(ls.auc[[iO]])["se"]
                ## brier
                iBrier <- brier(labels = data$XXresponseXX, predictions = iPred, fold = iFold, observation = iObs, iid = iCV.iid,
                                null = null["Brier"], conf.level = conf.level, transformation = transformation)
                cv.brier[iCurrent,name.col] <- cbind(model = names.object[iO], confint(iBrier), foldCV.number = fold.allnumber[iNumber], foldCV.size = sum(fold.size))
                if(iNumber==1){
                    ls.brier[[iO]] <- iBrier
                    cv.iid.brier[,iO] <- iid(iBrier)
                    if(iO>1){
                        iStat <- (cv.brier[iCurrent,"estimate"] - cv.brier[iPrevious,"estimate"]) / sqrt(crossprod(cv.iid.brier[,iO]-cv.iid.brier[,iO-1]))
                        cv.brier[iCurrent,"p.value_comp"] <- 2*(1-stats::pnorm(abs(iStat)))
                    }
                }
            }
        }

        ## *** export
        out <- rbind(out,
                     cbind(method = "cv", metric = "auc", cv.auc),
                     cbind(method = "cv", metric = "brier", cv.brier)
                     )
        if(simplify>0){
            out$foldCV.size <- NULL
            if(n.number==1){
                out$foldCV.number <- NULL
            }
        }
        if(is.null(attr(out,"response"))){attr(out,"response") <- list()}
        attr(out,"response")[["cv"]] <- data$XXresponseXX
        if(is.null(attr(out,"prediction"))){attr(out,"prediction") <- list()}
        attr(out,"prediction")[["cv"]] <- cv.predictions
        attr(attr(out,"prediction")[["cv"]],"index") <- cv.indexing
        if(is.null(attr(out,"iid"))){attr(out,"iid") <- list(auc = NULL, brier = NULL)}
        attr(out,"iid")$auc[["cv"]] <- cv.iid.auc
        attr(out,"iid")$brier[["cv"]] <- cv.iid.brier
        attr(out,"auc") <- ls.auc
        attr(out,"brier") <- ls.brier
        if(trace){cat(" done. \n")}
        
    }


    ## ** export
    rownames(out) <- NULL
    class(out) <- append("performance",class(out))
    return(out)
}

## * .performance_predict
.performance_predict <- function(object, n.obs, newdata, auc.type){

    out <- list(estimate = NULL,
                se = NULL,
                iid = matrix(0, nrow = n.obs, ncol = NROW(newdata)))
    
    if(inherits(object,"ranger")){
        if(auc.type == "classical"){
            iPred <- predict(object, data = newdata, type = "response")
            if(is.matrix(object$predictions)){
                out$estimate <- iPred$predictions[,which(object$forest$class.values==1)]
            }else{
                warning("0/1 predictions are used instead of probabilities. \n")
                out$estimate <- iPred$predictions
            }
        }else if(auc.type == "probabilistic"){
            iPred <- predict(object, data = newdata, type = "se")
            out$estimate <- iPred$predictions[,which(object$forest$class.values==1)]
            out$se <- iPred$se
        }
    }else if(inherits(object,"randomForest")){
        out$estimate <- predict(object, newdata = newdata, type = "prob")[,2]
    }else if(inherits(object,"glm")){
        if(object$converged==FALSE){return(NULL)}
        iPred <- .predict.logit(object, newdata = newdata)
        out$estimate <- iPred["estimate",]
        if(auc.type == "probabilistic"){
            out$se <- iPred["se",]
        }
        out$iid <- attr(iPred,"iid")
    }

    ##
    return(out)

}

##----------------------------------------------------------------------
### performance.R ends here
