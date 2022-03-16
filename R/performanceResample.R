### performanceResample.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  3 2022 (12:01) 
## Version: 
## Last-Updated: mar 16 2022 (13:57) 
##           By: Brice Ozenne
##     Update #: 61
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * performanceResample (documentation)
##' @title Uncertainty About Performance of a Classifier (EXPERIMENTAL)
##' @description Use resampling to quantify uncertainties about the performance of one or several binary classifiers evaluated via cross-validation.
##'
##' @param object a \code{glm} or \code{range} object, or a list of such object.
##' @param data [data.frame] the training data.
##' @param name.response [character] The name of the response variable (i.e. the one containing the categories).
##' @param type.resampling [character] Should non-parametric bootstrap (\code{"bootstrap"}) or permutation of the outcome (\code{"permutation"}) be used.
##' @param n.resampling [integer,>0] Nnumber of bootstrap samples or permutations.
##' @param fold.number [integer,>0] Nnumber of folds used in the cross-validation. Should be strictly positive.
##' @param conf.level [numeric, 0-1] confidence level for the confidence intervals.
##' @param cpus [integer, >0] the number of CPU to use. If strictly greater than 1, resampling is perform in parallel. 
##' @param seed [integer, >0] seed used to ensure reproducibility.
##' @param trace [logical] Should the execution of the function be traced.
##' @param filename [character] Prefix for the files containing each result.
##' @param ... arguments passed to \code{\link{performance}}.
##' 

## * performanceResample (code)
##' @export
performanceResample <- function(object, data = NULL, name.response = NULL,
                                type.resampling, n.resampling = 1000, fold.number = 100, conf.level = 0.95,
                                cpus = 1, seed = NULL, trace = TRUE, filename = NULL, ...){

    ## ** fix randomness
    if(!is.null(seed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)
    }

    ## ** normalize arguments
    if(fold.number<=0){
        stop("Argument \'fold.number\' should be strictly positive. \n")
    }
    type.resampling <- match.arg(type.resampling, c("permutation", "bootstrap"))

    if(is.null(data)){            
        initData <- FALSE
    }else{
        initData <- data
        attr(initData,"internal") <- FALSE
    }
    initPerf <- performance(object, data = initData, name.response = name.response,
                            fold.number = fold.number, conf.level = NA, trace = FALSE, seed = NULL, ...)
    if(!is.null(filename)){
        if(!is.null(seed)){
            filename <- paste0(filename,"-seed",seed)
        }
        saveRDS(initPerf, file = paste0(filename,".rds"))
    }
                

    if(is.null(data)){
        data <- attr(initPerf,"data")
    }
    if(is.null(name.response)){
        name.response <- attr(initPerf,"name.response")
    }
    data[["XXresponseXX"]] <- NULL

    ## ** single run function
    if(type.resampling=="permutation"){
        dataResample <- as.data.frame(data)
        attr(dataResample,"internal") <- FALSE ## only do CV

        warperResampling <- function(i){
            dataResample[[name.response]] <- sample(data[[name.response]])
            iPerf <- try(suppressWarnings(performance(object, data = dataResample, name.response = name.response, fold.number = fold.number, trace = trace-1, conf.level = NA, seed = NULL, ...)),
                         silent = FALSE)
            if(inherits(iPerf, "try-error")){
                return(NULL)
            }else{
                if(!is.null(filename)){
                    saveRDS(cbind(sample = i, iPerf[,c("metric","model","estimate")]), file = paste0(filename,"-",type.resampling,i,".rds"))
                }
                return(cbind(sample = i, iPerf[,c("metric","model","estimate")]))
            }
        }
    }else if(type.resampling=="bootstrap"){
        warperResampling <- function(i){
            dataResample <- data[sample(NROW(data), size = NROW(data), replace = TRUE),,drop=FALSE]
            attr(dataResample,"internal") <- FALSE ## only do CV
            iPerf <- try(suppressWarnings(performance(object, data = dataResample, name.response = name.response, fold.number = fold.number, trace = trace-1, conf.level = NA, seed = NULL, ...)),
                         silent = FALSE)
            if(inherits(iPerf, "try-error")){
                return(NULL)
            }else{
                if(!is.null(filename)){
                    saveRDS(cbind(sample = i, iPerf[,c("metric","model","estimate")]), file = paste0(filename,"-",type.resampling,i,".rds"))
                }
                return(cbind(sample = i, iPerf[,c("metric","model","estimate")]))
            }
        }
    }
    ## warperResampling(5)
    
    ## serial calculations
    if(cpus==1){
        ## ** method to loop
        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }

        ## ** loop
        ls.resampling <- do.call(method.loop,
                                 args = list(X = 1:n.resampling,
                                             FUN = warperResampling)
                                 )

        if(!is.null(seed)){rm(.Random.seed, envir=.GlobalEnv)} # restaure original seed

    }else{ ## parallel calculations
        ## define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.resampling, style = 3)          
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        ## link to foreach
        doParallel::registerDoParallel(cl)

        ## export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(library(BuyseTest, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE))
        })
        toExport <- NULL
        iB <- NULL ## [:forCRANcheck:] foreach        
        ls.resampling <- foreach::`%dopar%`(
                                      foreach::foreach(iB=1:n.resampling,
                                                       .export = toExport,
                                                       .packages = "data.table"),                                            
                                      {                                           
                                          if(trace>0){utils::setTxtProgressBar(pb, iB)}

                                           return(warperResampling(iB))
                      
                                       })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }

    ## ** gather results
    alpha <- 1-conf.level
    dt.resampling <- data.table::as.data.table(do.call(rbind, ls.resampling))
    if(type.resampling=="permutation"){
        data.table::setnames(dt.resampling, old = "estimate", new ="estimate.perm")
        dt.initPerf <- as.data.table(as.list(initPerf))
        dt.merge <- dt.resampling[dt.initPerf, on = c("metric","model")]
        out <- rbind(dt.merge[dt.merge$metric == "auc", list(estimate = mean(.SD$estimate), se = NA, lower = NA, upper = NA,
                                                             p.value = mean(abs(.SD$estimate)<=abs(.SD$estimate.perm)),
                                                             p.value_comp = NA), by= c("metric","model")],
                     dt.merge[dt.merge$metric == "brier", list(estimate = mean(.SD$estimate), se = NA, lower = NA, upper = NA,
                                                               p.value = mean(.SD$estimate>.SD$estimate.perm),
                                                               p.value_comp = NA), by= c("metric","model")])
    }else if(type.resampling=="bootstrap"){
        vec.estimate <- initPerf[,"estimate"]
        M.resampling <- as.matrix(dcast(dt.resampling, value.var = "estimate", formula = sample~metric+model))[,-1,drop=FALSE]
        out.metric <- sapply(strsplit(colnames(M.resampling), split = "_", fixed = TRUE),"[[",1)
        out.model <- as.character(mapply(x = paste0("^",out.metric,"_"), y = colnames(M.resampling), FUN = function(x,y){gsub(pattern = x, replacement = "", x = y)}))

        out <- data.frame(metric = out.metric, model = out.model,
                          confint_percentileBootstrap(Delta = vec.estimate,
                                                      Delta.resampling = M.resampling,
                                                      null = c(auc = 0.5, brier = NA)[out.metric], alternative = "two.sided", alpha = 1-conf.level,
                                                      endpoint = colnames(M.resampling), backtransform.delta = function(x){x}))


        vecauc.estimate <- vec.estimate[out.metric=="auc"]
        vecbrier.estimate <- vec.estimate[out.metric=="brier"]
        Mauc.resampling <- M.resampling[,out.metric=="auc",drop=FALSE]
        Mbrier.resampling <- M.resampling[,out.metric=="brier",drop=FALSE]
        if(length(vecauc.estimate)>1){
                
            deltaout <- confint_percentileBootstrap(Delta = c(0,vecauc.estimate[-length(vecauc.estimate)] - vecauc.estimate[-1],
                                                              0,vecbrier.estimate[-length(vecbrier.estimate)] - vecbrier.estimate[-1]),
                                                    Delta.resampling = cbind(0,Mauc.resampling[,-ncol(Mauc.resampling)] - Mauc.resampling[,-1],
                                                                             0,Mbrier.resampling[,-ncol(Mbrier.resampling)] - Mbrier.resampling[,-1]),
                                                    null = c(NA, rep(0, length(vecauc.estimate)-1), NA, rep(0, length(vecbrier.estimate)-1)), alternative = "two.sided", alpha = 1-conf.level,
                                                    endpoint = colnames(M.resampling), backtransform.delta = function(x){x})
            out$p.value_comp <- deltaout[,"p.value"]
        }else{
            out$p.value_comp <- NA
        }

        names(out)[names(out) == "lower.ci"] <- "lower"
        names(out)[names(out) == "upper.ci"] <- "upper"
        out$null <- NULL
    }
    attr(out, "resample") <- dt.resampling
    attr(out, "original") <- initPerf


    ## ** export
    class(out) <- append("performanceResample",class(out))
    return(out)
}



##----------------------------------------------------------------------
### performanceResample.R ends here
