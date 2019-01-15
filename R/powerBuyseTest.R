### powerBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: jan 15 2019 (10:16) 
##           By: Brice Ozenne
##     Update #: 337
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## * Documentation - powerBuyseTest
#' @name powerBuyseTest
#' @title Performing simulation studies with BuyseTest
#' @aliases powerBuyseTest
#' 
#' @description Performs a simulation studies for several sample sizes.
#' Returns estimates, standard errors, confidence intervals and p.values.
#'
#' @param sim [function] take two arguments:
#' the sample size in the control group (\code{n.C}) and the sample size in the treatment group (\code{n.C})
#' and generate datasets. The datasets must be data.table objects.
#' @param sample.size [integer vector, >0] the various sample sizes at which the simulation should be perform.
#' Disregarded if any of the arguments \code{sample.sizeC} or \code{sample.sizeT} are specified.
#' @param sample.sizeC [integer vector, >0] the various sample sizes in the control group.
#' @param sample.sizeT [integer vector, >0] the various sample sizes in the treatment group.
#' @param n.rep [integer, >0] the number of simulations.
#' @param null [numeric vector] the null hypothesis to be tested for the net benefit (first element) and the win ratio (second element).
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
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the asymptotic variance.
#' @param ... parameters from \code{BuyseTest}.
#' 

## * powerBuyseTest (examples)
##' @rdname powerBuyseTest
##' @examples
##'
##' ## using simBuyseTest
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(100), n.rep = 2,
##'               formula = Treatment ~ tte(eventtime, censoring = status),
##'               method.inference = "asymptotic", trace = 4)
##'
##' ## using user defined simulation function
##' simFCT <- function(n.C, n.T){
##'     out <- data.table(Y=rnorm(n.C+n.T),
##'                       T=c(rep(1,n.C),rep(0,n.T))
##'                      )
##' return(out)
##' }
##'
##' powerBuyseTest(sim = simFCT, sample.size = c(100), n.rep = 2,
##'               formula = T ~ cont(Y), method.inference = "asymptotic", trace = 4)
##' 

## * powerBuyseTest (code)
##' @rdname powerBuyseTest
##' @export
powerBuyseTest <- function(sim, sample.size, sample.sizeC = NULL, sample.sizeT = NULL, n.rep, null = c(0,1), cpus = 1,                          
                           alternative = NULL, seed = 10, conf.level = NULL, order.Hprojection = NULL, transformation = NULL, trace = 1,
                           ...){

    call <- match.call()$sim

    ## ** normalize and check arguments
    name.call <- names(match.call())
    option <- BuyseTest.options()
    if(is.null(transformation)){
        transformation <- option$transformation
    }
    if(is.null(order.Hprojection)){
        order.Hprojection <- option$order.Hprojection
    }
    if(is.null(conf.level)){
        conf.level <- option$conf.level
    }
    if(is.null(alternative)){
        alternative <- option$alternative
    }
    alpha <- 1 - conf.level
    
    if("keep.pairScore" %in% name.call){
        stop("\'keep.pairScore\' is not an argument of powerBuyseTest \n")
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
    validNumeric(null,
                 valid.length = 2,
                 method = "BuyseTest")
    names(null) <- c("netBenefit","winRatio")
        
    n.sample.size <- length(sample.sizeT)
    grid.inference <- expand.grid(order = order.Hprojection,
                                  transformation = transformation)
    n.inference <- NROW(grid.inference)
    sample.sizeCmax <- sample.sizeC[n.sample.size]
    sample.sizeTmax <- sample.sizeT[n.sample.size]

    dt.tempo <- sim(n.C = sample.sizeC[1], n.T = sample.sizeT[1])
    if(!data.table::is.data.table(dt.tempo)){
        stop("The function defined by the argument \'sim\' must return a data.table object.\n")
    }
    
    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(cpus = cpus, option = option, name.call = name.call, alternative = alternative,
                              data = NULL, model.tte = NULL, keep.pairScore = TRUE, ...)
    if(outArgs$method.tte==1 && n.sample.size >1){
        stop("Peron correction not compatible with powerBuyseTest for more than one sample size \n")
    }
    if(any(outArgs$operator!=">0")){
        stop("Cannot use argument \'operator\' with powerBuyseTest \n")
    }
    if(!is.null(outArgs$strata)){
        stop("Cannot use argument \'strata\' with powerBuyseTest \n")
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
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    out.name <- c("Wscheme","endpoint.UTTE","index.UTTE","D.UTTE","reanalyzed","outSurv")
    outArgs[out.name] <- buildWscheme(method.tte = outArgs$method.tte,
                                      endpoint = outArgs$endpoint,
                                      D = outArgs$D,
                                      D.TTE = outArgs$D.TTE,
                                      n.strata = outArgs$n.strata,
                                      type = outArgs$type,
                                      threshold = outArgs$threshold)

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
    name.copy <- c("call", "sim", "option",
                   "outArgs", "sample.sizeTmax", "sample.sizeCmax", "n.sample.size",
                   "sample.size", "sample.sizeC", "sample.sizeT", "n.rep", "alternative", "seed")
    envirBT <- new.env()
    for(iObject in name.copy){ ## iObject <- name.copy[1]
        envirBT[[iObject]] <- eval(parse(text = iObject))
    }

    ## ** warper
    warper <- function(i, envir){
        iOut <- matrix(NA, nrow = n.sample.size * n.inference, ncol = 14,
                       dimnames = list(NULL, c("simulation","n.T","n.C","method.inference",
                                               "netBenefit","netBenefit.se","netBenefit.lower","netBenefit.upper","netBenefit.p.value",
                                               "winRatio","winRatio.se","winRatio.lower","winRatio.upper","winRatio.p.value")))
        iOut <- as.data.frame(iOut)
        iOut[,"simulation"] <- i

        ## *** Initialize data
        out.name <- c("data","M.endpoint","M.censoring",
                      "index.C","index.T","index.strata",
                      "index.endpoint","index.censoring","level.treatment","level.strata",
                      "n.strata","n.obs","n.obsStrata","cumn.obsStrata")
        envir$outArgs[out.name] <- initializeData(data = do.call(eval(envir$call), args = list(n.T = sample.sizeTmax, n.C = sample.sizeCmax)),
                                                  type = envir$outArgs$type,
                                                  endpoint = envir$outArgs$endpoint,
                                                  method.tte = envir$outArgs$method.tte,
                                                  censoring = envir$outArgs$censoring,
                                                  operator = envir$outArgs$operator,
                                                  strata = envir$outArgs$strata,
                                                  treatment = envir$outArgs$treatment,
                                                  copy = FALSE)

        ## *** Point estimate
        outPoint <- .BuyseTest(envir = envir,
                               keep.pairScore = TRUE,
                               method.inference = "none")

        ## *** put results into a data.table
        tablePairScore <- pairScore2dt(outPoint$tableScore,
                                       level.treatment = envir$outArgs$level.treatment,
                                       level.strata = envir$outArgs$level.strata,
                                       n.strata = envir$outArgs$n.strata,
                                       endpoint = envir$outArgs$endpoint,
                                       threshold = envir$outArgs$threshold)

        for(iSample in 1:n.sample.size){ ## iSample <- 1
            iIndex.store <- seq(1 + (iSample-1)*n.inference, iSample*n.inference)
            iOut[iIndex.store,"n.C"] <- envir$sample.sizeC[iSample]
            iOut[iIndex.store,"n.T"] <- envir$sample.sizeT[iSample]
            iPairs <- envir$sample.sizeC[iSample]*envir$sample.sizeT[iSample]
            
            tableSample <- lapply(tablePairScore, function(iEndpoint){ ## iEndpoint <- tablePairScore[[1]]
                ## restrict pairs 
                index <- intersect(which(iEndpoint$indexWithinStrata.C <= envir$sample.sizeC[iSample]),
                                   which(iEndpoint$indexWithinStrata.T <= envir$sample.sizeT[iSample]))
                iEndpoint <- iEndpoint[index]

                ## update index in dataset
                old.position <- sort(c(unique(iEndpoint$index.T),unique(iEndpoint$index.C)))
                new.position <- rank(old.position)
                old2new <- rep(NA, max(old.position))
                old2new[old.position] <- new.position
                iEndpoint[, c("index.T") := old2new[.SD$index.T]]
                iEndpoint[, c("index.C") := old2new[.SD$index.C]]

                ## update correction (no strata)
                if(envir$outArgs$correction.uninf>0 && iSample < n.sample.size && sum(iEndpoint$uninf)>0){
## new weighting                    
                        mfactor <- sum(iEndpoint$favorable + iEndpoint$unfavorable + iEndpoint$neutral + iEndpoint$uninf) / sum(iEndpoint$favorable + iEndpoint$unfavorable + iEndpoint$neutral)
                        iEndpoint[, c("favorableC") :=.SD$favorable * .SD$weight * mfactor]
                        iEndpoint[, c("unfavorableC") :=.SD$unfavorable * .SD$weight * mfactor]
                        iEndpoint[, c("neutralC") :=.SD$neutral * .SD$weight * mfactor]
                    
                }
                
                ## export
                return(iEndpoint)
            })

            ## *** Point estimate
            MresSample <- do.call(rbind,lapply(tableSample, function(iEndpoint){
                return(c("npairs" = NROW(iEndpoint),
                         "favorable" = sum(iEndpoint$favorableC),
                         "unfavorable" = sum(iEndpoint$unfavorableC)))
            }))
            
            iOut[iIndex.store,"netBenefit"] <- (sum(MresSample[,"favorable"]) - sum(MresSample[,"unfavorable"]))/MresSample[1,"npairs"]
            iOut[iIndex.store,"winRatio"] <- sum(MresSample[,"favorable"]) / sum(MresSample[,"unfavorable"])

            ## *** Inference 
            if(outArgs$method.inference %in% c("asymptotic")){
                ## warning: only work if no strata, otherwise n.pairs/count.favorable/count.unfavorable needs to be sum over strata
                ## see BuyseTest.R
                outCovariance <- inferenceUstatistic(tableSample,
                                                     order = max(order.Hprojection), ## if order 1 and 2 are requested by the user then feed order 2
                                                     count.favorable = matrix(MresSample[,"favorable"], nrow = 1),
                                                     count.unfavorable = matrix(MresSample[,"unfavorable"], nrow = 1),
                                                     n.pairs = iPairs,
                                                     n.C = envir$sample.sizeC[iSample],
                                                     n.T = envir$sample.sizeT[iSample],
                                                     level.strata = envir$outArgs$level.strata,
                                                     n.strata = envir$outArgs$n.strata,
                                                     n.endpoint = length(envir$outArgs$endpoint),
                                                     endpoint = envir$outArgs$endpoint)
                
                for(iStatistic in c("netBenefit","winRatio")){
                    for(iInference in 1:n.inference){ ## iInference <- 1
                        ## cat("\n",iStatistic," ",iInference," ",iSample,"\n")
                            
                        iTransformation <- grid.inference[iInference,"transformation"]
                        iOrder <- grid.inference[iInference,"order"]
                        if(iOrder == max(order.Hprojection)){
                            iCovariance <- outCovariance$Sigma                         
                        }else{
                            ## recompute the covariance matrix to remove second order term 
                            iCovariance <- .iid2cov(A.iid = outCovariance$iid1, A2.iid = NULL,
                                                    order = iOrder, endpoint = envir$outArgs$endpoint, n.endpoint = length(envir$outArgs$endpoint))
                        }
                        outCI <- confint_Ustatistic(Delta = iOut[iIndex.store[iInference],iStatistic],
                                                    pc.favorable = MresSample[,"favorable"]/MresSample[,"npairs"],
                                                    pc.unfavorable =  MresSample[,"unfavorable"]/MresSample[,"npairs"],
                                                    covariance = iCovariance, statistic = iStatistic,
                                                    alternative = alternative, alpha = alpha, null = null[iStatistic],
                                                    endpoint = envir$outArgs$endpoint, transformation = iTransformation,
                                                    continuity.correction = envir$option$continuity.correction,
                                                    n.pairs = iPairs)

                        iOut[iIndex.store[iInference],"method.inference"] <- paste0("order=",iOrder," - transformation=",iTransformation)
                        iOut[iIndex.store[iInference],paste0(iStatistic,".se")] <- outCI[1,"se"]
                        iOut[iIndex.store[iInference],paste0(iStatistic,".lower")] <- outCI[1,"lower.ci"]
                        iOut[iIndex.store[iInference],paste0(iStatistic,".upper")] <- outCI[1,"upper.ci"]
                        iOut[iIndex.store[iInference],paste0(iStatistic,".p.value")] <- outCI[1,"p.value"]
                    }
                }
            

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
    }else { ## *** parallel permutation test
        n.block <- max(cpus,round(sqrt(n.rep)))
        rep.perBlock0 <- max(1,floor(n.rep/n.block))
        rep.perBlock <- c(rep(rep.perBlock0, n.block-1), n.rep - (n.block-1)*rep.perBlock0)
        cumsum.rep.perBlock <- c(0,cumsum(rep.perBlock))

        ## define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.block, style = 3)          
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        ## link to foreach
        doParallel::registerDoParallel(cl)

        ## seed
        if (!is.null(seed)) {
            set.seed(seed)
            seqSeed <- sample.int(1e3, size = cpus)
            parallel::clusterApply(cl, seqSeed, function(x){
                set.seed(x)
            })
        }         

        ## export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(library(BuyseTest, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE))
        })
        ## export functions
        toExport <- c(".BuyseTest","initializeData","pairScore2dt","inferenceUstatistic","confint_Ustatistic", ".iid2cov", "validNumeric")
        
        i <- NULL ## [:forCRANcheck:] foreach
        ls.simulation <- foreach::`%dopar%`(
                                      foreach::foreach(i=1:n.block,
                                                       .export = toExport),                                            
                                      {                                           
                                          if(trace>0){utils::setTxtProgressBar(pb, i)}
                                          ls.out <- list()
                                          for(j in 1:rep.perBlock[i]){
                                            ls.out[[j]] <- warper(i = j + cumsum.rep.perBlock[i], envir = envirBT)
                                          }
                                          return(do.call(rbind, ls.out))
                      
                                       })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }
    dt.out <- as.data.table(do.call(rbind, ls.simulation))

    ## ** export
    attr(outArgs$method.inference,"continuity.correction") <- option$continuity.correction
        
    BuyseSim.object <- BuyseSim(
        alternative = alternative,      
        method.inference = outArgs$method.inference,
        conf.level = conf.level,
        null = null,
        n.rep = n.rep,      
        results = dt.out
    )

    return(BuyseSim.object)
}

######################################################################
### powerBuyseTest.R ends here
