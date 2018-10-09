### evalBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: okt  9 2018 (08:24) 
##           By: Brice Ozenne
##     Update #: 133
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## * Documentation - evalBuyseTest
#' @name evalBuyseTest
#' @title Performing simulation studies with BuyseTest
#' @aliases evalBuyseTest
#' 
#' @description Performs a simulation studies for several sample sizes.
#' Returns estimates, standard errors, confidence intervals and p.values.
#'
#' @param sim [function] take two arguments:
#' the sample size in the control group (\code{n.C}) and the sample size in the treatment group (\code{n.C})
#' and generate datasets. The datasets must be data.table objects.
#' @param sample.size the various sample sizes at which the simulation should be perform.
#' Disregarded if any of the arguments \code{sample.sizeC} or \code{sample.sizeT} are specified.
#' @param sample.sizeC the various sample sizes in the control group.
#' @param sample.sizeT the various sample sizes in the treatment group.
#' @param n.rep the number of simulations.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the alternative hypothesis.
#' Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. 
#' Default value read from \code{BuyseTest.options()}.
#' @param seed [integer, >0] the seed to consider for the simulation study.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param trace [integer] should the execution of the function be traced?
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}.
#' @param ... parameters from \code{BuyseTest}.
#' 

## * evalBuyseTest (examples)
##' @rdname evalBuyseTest
##' @examples
##'
##' simFCT <- function(n.C, n.T){
##'     out <- data.table(Y=rnorm(n.C+n.T),
##'                       T=c(rep(1,n.C),rep(0,n.T))
##'                      )
##' return(out)
##' }
##'
##' evalBuyseTest(sim = simFCT, sample.size = c(100), n.rep = 2,
##'               formula = T ~ cont(Y), method.inference = "asymptotic", trace = 4)
##' 

## * evalBuyseTest (code)
##' @rdname evalBuyseTest
##' @export
evalBuyseTest <- function(sim, sample.size, sample.sizeC = NULL, sample.sizeT = NULL, n.rep, cpus = 1,                          
                          alternative = NULL, seed = 10, conf.level = NULL, transformation = NULL, trace = 1, ...){

    name.call <- names(match.call())
    option <- BuyseTest.options()
    alpha <- 1 - conf.level
    if(is.null(transformation)){
        transformation <- option$transformation
    }
    if(is.null(conf.level)){
        conf.level <- option$conf.level
    }
    if(is.null(alternative)){
        alternative <- option$alternative
    }
    
    if("keep.pairScore" %in% name.call){
        stop("\'keep.pairScore\' is not an argument of evalBuyseTest \n")
    }
    if(is.null(sample.sizeC) && is.null(sample.sizeT)){
        sample.sizeC <- sample.size
        sample.sizeT <- sample.size
    }

    validInteger(sample.sizeC,
                 valid.length = NULL,
                 min = 1, refuse.duplicates = TRUE,
                 method = "BuyseTest")
    validInteger(sample.sizeT,
                 valid.length = NULL,
                 min = 1, refuse.duplicates = TRUE,
                 method = "BuyseTest")
    if(length(sample.sizeT)!=length(sample.sizeC)){
        stop("Arguments \'sample.sizeT\ and \'sample.sizeC\' must have the same length \n")
    }
    n.sample.size <- length(sample.sizeT)
    sample.sizeCmax <- sample.sizeC[n.sample.size]
    sample.sizeTmax <- sample.sizeT[n.sample.size]

    dt.tempo <- sim(n.C = sample.sizeC[1], n.T = sample.sizeT[1])
    if(!data.table::is.data.table(dt.tempo)){
        stop("The function defined by the argument \'sim\' must return a data.table object.\n")
    }
    
    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(cpus = cpus, option = option, name.call = name.call, alternative = alternative,
                              data = NULL, keep.pairScore = TRUE, ...)

    if(any(outArgs$operator!=">0")){
        stop("Cannot use argument \'operator\' with evalBuyseTest \n")
    }
    if(!is.null(outArgs$strata)){
        stop("Cannot use argument \'strata\' with evalBuyseTest \n")
    }
    if(outArgs$method.inference %in% c("none","asymptotic") == FALSE){
        stop("Argument \'method.inference\' must be \"none\" or \"asymptotic\" \n")
    }

    cpus <- outArgs$cpus
    outArgs$cpus <- 1
    outArgs$trace <- 0

    ## ** initialization data
    outArgs$level.treatment <- levels(as.factor(dt.tempo[[outArgs$treatment]]))
    outArgs$n.strata <- 1
    outArgs$level.strata <- "1"
    outArgs$allstrata <- NULL

    ## ** create weights matrix for survival endpoints
    if(outArgs$D.TTE>0 && outArgs$D>1){
        ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
        outArgs[c("Wscheme","index.survivalM1","threshold.TTEM1")] <- buildWscheme(endpoint = outArgs$endpoint,
                                                                                   D = outArgs$D,
                                                                                   D.TTE = outArgs$D.TTE,
                                                                                   type = outArgs$type,
                                                                                   threshold = outArgs$threshold)
    }else{ #  factice arguments. Will be sent to the C++ arguments to fill the argument but not used by the function.
        outArgs$Wscheme <- matrix(nrow=0,ncol=0)
        outArgs$index.survivalM1 <- numeric(0)
        outArgs$threshold.TTEM1 <- numeric(0)
    }

    ## ** Display
    if (trace > 1) {
        cat("         Simulation study with BuyseTest \n\n")

        if(trace > 2){
            do.call(printGeneral, args = outArgs)
            if(outArgs$method.inference!="none"){
                do.call(printInference, args = outArgs)
            }
        }
        
        cat("Simulation\n",
            "   - sample size: ",paste(sample.size, collapse = " "),"\n",
            "   - repetitions: ",n.rep,"\n",
            "   - cpus       : ",cpus,"\n",
            sep = "")
        cat(" \n")

     
        
    }
    ## ** define environment
    name.copy <- c("initializeSurvival_Peron", "initializeData", ".BuyseTest", "sim",
                   "outArgs", "sample.sizeTmax", "sample.sizeCmax", "n.sample.size",
                   "sample.size", "sample.sizeC", "sample.sizeT", "n.rep", "alternative", "seed")
    envirBT <- new.env()
    for(iObject in name.copy){ ## iObject <- name.copy[1]
        envirBT[[iObject]] <- eval(parse(text = iObject))
    }
    names(envirBT)

    ## ** warper
    warper <- function(i, envir){
        iOut <- matrix(NA, nrow = n.sample.size, ncol = 13,
                       dimnames = list(NULL, c("simulation","n.T","n.C",
                                               "netChance","netChance.se","netChance.lower","netChance.upper","netChance.p.value",
                                               "winRatio","winRatio.se","winRatio.lower","winRatio.upper","winRatio.p.value")))

        ## *** Simulate data
        envir$outArgs$data <- envir$sim(n.T = sample.sizeTmax, n.C = sample.sizeCmax)
        envir$outArgs$n.obs <- NROW(envir$outArgs$data)
        envir$outArgs$n.obsStrata <- envir$outArgs$n.obs

        ## *** Point estimate
        outPoint <- .BuyseTest(envir = envir,
                               keep.pairScore = TRUE,
                               method.inference = "none")

        ## *** put results into a data.table
        envir$indexT <- which(envir$outArgs$data[[envir$outArgs$treatment]]==1)
        envir$indexC <- which(envir$outArgs$data[[envir$outArgs$treatment]]==0)
        ## outPoint$tableScore[[1]]

        tablePairScore <- pairScore2dt(outPoint$tableScore,
                                       level.treatment = envir$outArgs$level.treatment,
                                       level.strata = envir$outArgs$level.strata,
                                       n.strata = envir$outArgs$n.strata,
                                       endpoint = envir$outArgs$endpoint,
                                       threshold = envir$outArgs$threshold,
                                       indexT = envir$indexT,
                                       indexC = envir$indexC)

        for(iSample in 1:n.sample.size){ ## iSample <- 1
            iOut[iSample,"simulation"] <- i
            iOut[iSample,"n.C"] <- envir$sample.sizeC[iSample]
            iOut[iSample,"n.T"] <- envir$sample.sizeT[iSample]

            tableSample <- lapply(tablePairScore, function(iEndpoint){ ## iEndpoint <- tablePairScore[[1]]
                index <- intersect(which(iEndpoint$indexWithinStrata.C <= envir$sample.sizeC[iSample]),
                                   which(iEndpoint$indexWithinStrata.T <= envir$sample.sizeT[iSample]))
                return(iEndpoint[index])
            })

            ## *** Point estimate
            MresSample <- do.call(rbind,lapply(tableSample, function(iEndpoint){
                return(c("npairs" = NROW(iEndpoint),
                         "favorable" = sum(iEndpoint$favorable),
                         "unfavorable" = sum(iEndpoint$unfavorable)))
            }))
            iOut[iSample,"netChance"] <- (sum(MresSample[,"favorable"]) - sum(MresSample[,"unfavorable"]))/MresSample[1,"npairs"]
            iOut[iSample,"winRatio"] <- sum(MresSample[,"favorable"]) / sum(MresSample[,"unfavorable"])

            ## *** Inference 
            if(outArgs$method.inference %in% c("asymptotic")){
                outCovariance <- inferenceUstatistic(tableSample,
                                                     count.favorable = matrix(MresSample[,"favorable"], nrow = 1),
                                                     count.unfavorable = matrix(MresSample[,"unfavorable"], nrow = 1),
                                                     n.pairs = envir$sample.sizeC[iSample]*envir$sample.sizeT[iSample],
                                                     n.C = envir$sample.sizeC[iSample],
                                                     n.T = envir$sample.sizeT[iSample],
                                                     n.strata = envir$outArgs$n.strata,
                                                     n.endpoint = length(envir$outArgs$endpoint),
                                                     endpoint = envir$outArgs$endpoint)

                outCI <- confint_Ustatistic(Delta = iOut[iSample,"netChance"], covariance = outCovariance$Sigma, statistic = "netChance",
                                            null = 0, alternative = alternative, alpha = alpha,
                                            endpoint = envir$outArgs$endpoint, transformation = transformation)
                iOut[iSample,"netChance.se"] <- outCI[1,"se"]
                iOut[iSample,"netChance.lower"] <- outCI[1,"lower.ci"]
                iOut[iSample,"netChance.upper"] <- outCI[1,"upper.ci"]
                iOut[iSample,"netChance.p.value"] <- outCI[1,"p.value"]

            }

        }
        
        return(iOut)
    }
    

    ## ** simulation study
    if (cpus == 1) { ## *** sequential permutation test
           
        if (!is.null(seed)) {set.seed(seed)} # set the seed

        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.simulation <- do.call(method.loop,
                                 args = list(X = 1:n.rep,
                                             FUN = function(X){
                                                 return(warper(i = X, envir = envirBT))                                                  
                                             })
                                 )
        dt.out <- as.data.table(do.call(rbind, ls.simulation))
    }else { ## *** parallel permutation test

        ## define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.rep, style = 3)          
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        ## link to foreach
        doParallel::registerDoParallel(cl)

        ## export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(library(BuyseTest, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE))
        })
        ## export functions
        toExport <- c(".BuyseTest","initializeSurvival_Peron","pairScore2dt","inferenceUstatistic","confint_Ustatistic")

        i <- NULL ## [:forCRANcheck:] foreach        
        ls.permutation <- foreach::`%dopar%`(
                                       foreach::foreach(i=1:n.rep,
                                                        .export = toExport),                                            
                                       {                                           
                                           if(trace>0){utils::setTxtProgressBar(pb, i)}
                                           return(warper(i = i, envir = envirBT))
                      
                                       })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }

    ## ** export
    return(dt.out)
}

######################################################################
### evalBuyseTest.R ends here
