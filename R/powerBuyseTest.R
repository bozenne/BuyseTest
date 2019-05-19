### powerBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: maj 19 2019 (16:50) 
##           By: Brice Ozenne
##     Update #: 468
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
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#' @param ... parameters from \code{BuyseTest}.
#' 

## * powerBuyseTest (examples)
##' @rdname powerBuyseTest
##' @examples
##' library(data.table)
##' 
##' ## using simBuyseTest
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(100), n.rep = 2,
##'                formula = Treatment ~ bin(toxicity),
##'                method.inference = "u-statistic", trace = 4)
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
##'               formula = T ~ cont(Y), method.inference = "u-statistic", trace = 4)
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
                 min = 1, refuse.duplicates = FALSE, ## accept duplicates for checking the software
                 method = "BuyseTest")
    validInteger(sample.sizeT,
                 valid.length = NULL,
                 min = 1, refuse.duplicates = FALSE, ## accept duplicates for checking the software
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
                                  index.sample.size = 1:n.sample.size,
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
    if(outArgs$scoring.rule==1 && n.sample.size > 1){
        stop("Peron correction not compatible with powerBuyseTest for more than one sample size\n")
    }
    if(any(outArgs$operator!=">0")){
        stop("Cannot use argument \'operator\' with powerBuyseTest \n")
    }
    if(!is.null(outArgs$strata)){
        stop("Cannot use argument \'strata\' with powerBuyseTest \n")
    }
    if(outArgs$method.inference %in% c("none","u-statistic") == FALSE){
        stop("Argument \'method.inference\' must be \"none\" or \"u-statistic\" \n")
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
    outArgs[out.name] <- buildWscheme(scoring.rule = outArgs$scoring.rule,
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
        iOut <- matrix(NA, nrow = n.inference, ncol = 14,
                       dimnames = list(NULL, c("simulation","n.T","n.C","method.inference",
                                               "netBenefit","netBenefit.se","netBenefit.lower","netBenefit.upper","netBenefit.p.value",
                                               "winRatio","winRatio.se","winRatio.lower","winRatio.upper","winRatio.p.value")))
        iOut <- as.data.frame(iOut)
        iOut[,"simulation"] <- i

        ## *** Initialize data
        out.name <- c("data","M.endpoint","M.censoring",
                      "index.C","index.T","index.strata",
                      "index.endpoint","index.censoring","level.treatment","level.strata", "method.score",
                      "n.strata","n.obs","n.obsStrata","cumn.obsStrata")

        envir$outArgs[out.name] <- initializeData(data = do.call(eval(envir$call), args = list(n.T = sample.sizeTmax, n.C = sample.sizeCmax)),
                                                  type = envir$outArgs$type,
                                                  endpoint = envir$outArgs$endpoint,
                                                  scoring.rule = envir$outArgs$scoring.rule,
                                                  censoring = envir$outArgs$censoring,
                                                  operator = envir$outArgs$operator,
                                                  strata = envir$outArgs$strata,
                                                  treatment = envir$outArgs$treatment,
                                                  copy = FALSE)

        ## *** Point estimate
        outPoint <- .BuyseTest(envir = envir,
                               method.inference = "none",
                               iid = TRUE,
                               pointEstimation = TRUE)

        ## *** generate BT results based on all sample sizes
        ls.BT <- .createSubBT(outPoint, order = order.Hprojection,
                              type = envir$outArgs$type,
                              endpoint = envir$outArgs$endpoint,
                              level.treatment = envir$outArgs$level.treatment,
                              scoring.rule = envir$outArgs$scoring.rule,
                              method.inference = envir$outArgs$method.inference,
                              hierarchical = envir$outArgs$hierarchical,
                              correction.uninf = envir$outArgs$correction.uninf,
                              threshold = envir$outArgs$threshold,
                              weight = envir$outArgs$weight,
                              sample.sizeC = envir$sample.sizeC, sample.sizeT = envir$sample.sizeT, n.sample.size = n.sample.size)

        for(iInference in 1:n.inference){ ## iInference <- 1

            ##
            iTransformation <- grid.inference[iInference,"transformation"]
            iOrder <- grid.inference[iInference,"order"]
            iSample.size <- grid.inference[iInference,"index.sample.size"]

            ##
            iOut[iInference,"n.T"] <- sample.sizeT[iSample.size]
            iOut[iInference,"n.C"] <- sample.sizeC[iSample.size]
            iOut[iInference,"method.inference"] <- paste0("order=",iOrder," - transformation=",iTransformation)

            ##
            BT.tempo <- ls.BT[[iSample.size]]
            ## ls.BT[[iSample.size]]@covariance
            if(length(order.Hprojection)==2 && iOrder==1){
                BT.tempo@covariance <- attr(BT.tempo@covariance,"first.order")
            }
                
            for(iStatistic in c("netBenefit","winRatio")){ ## iStatistic <- "winRatio"

                outCI <- suppressWarnings(confint(BT.tempo, transformation = iTransformation, statistic = iStatistic))

                iOut[iInference,paste0(iStatistic)] <- outCI[NROW(outCI),"estimate"]
                iOut[iInference,paste0(iStatistic,".se")] <- outCI[NROW(outCI),"se"]
                iOut[iInference,paste0(iStatistic,".lower")] <- outCI[NROW(outCI),"lower.ci"]
                iOut[iInference,paste0(iStatistic,".upper")] <- outCI[NROW(outCI),"upper.ci"]
                iOut[iInference,paste0(iStatistic,".p.value")] <- outCI[NROW(outCI),"p.value"]
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
        if(!is.null(seed)){rm(.Random.seed, envir=.GlobalEnv)} # restaure original seed
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

## * .createSubBT
.createSubBT <- function(object, order,
                         type, endpoint, level.treatment, scoring.rule, method.inference, hierarchical, correction.uninf,
                         threshold, weight,
                         sample.sizeT, sample.sizeC, n.sample.size){

    n.endpoint <- length(endpoint)
    out <- vector(mode = "list", length = n.sample.size)

    ## ** put pairwise comparisons into a data.table
    table <- pairScore2dt(object$tableScore,
                          level.treatment = level.treatment,
                          level.strata = "1",
                          n.strata = 1,
                          endpoint = endpoint,
                          threshold = threshold)

    ## ** create BT for each sample size
    for(iSample in 1:n.sample.size){ ## iSample <- 1

        if(iSample!=n.sample.size || 2 %in% order){
            iCount_favorable <- matrix(NA, nrow = 1, ncol = n.endpoint)
            iCount_unfavorable <- matrix(NA, nrow = 1, ncol = n.endpoint)
            iCount_neutral <- matrix(NA, nrow = 1, ncol = n.endpoint)
            iCount_uninf <- matrix(NA, nrow = 1, ncol = n.endpoint)

            iLS.table <- vector(mode = "list", length = n.endpoint)
            
            ## *** subset of index
            iN.C <- sample.sizeC[iSample]
            iN.T <- sample.sizeT[iSample]
            iN_pairs <- iN.C * iN.T

            ## *** update scores
            for(iEndpoint in 1:length(endpoint)){ ## iEndpoint <- 1

                ## restrict pairs 
                index <- intersect(which(table[[iEndpoint]]$indexWithinStrata.C <= iN.C),
                                   which(table[[iEndpoint]]$indexWithinStrata.T <= iN.T))
                iTable <- table[[iEndpoint]][index]

                ## update index in dataset
                old.position <- sort(c(unique(iTable$index.T),unique(iTable$index.C)))
                new.position <- rank(old.position)
                old2new <- rep(NA, max(old.position))
                old2new[old.position] <- new.position
                iTable[, c("index.T") := old2new[.SD$index.T]]
                iTable[, c("index.C") := old2new[.SD$index.C]]

                ## update correction (no strata)
                if(correction.uninf>0 && iSample < n.sample.size && sum(iTable$uninf)>0){
                    ## new weighting
                    if(correction.uninf == 1){
                        mfactorFavorable <- sum(iTable$favorable) / sum(iTable$favorable + iTable$unfavorable + iTable$neutral)
                        mfactorUnfavorable <- sum(iTable$unfavorable) / sum(iTable$favorable + iTable$unfavorable + iTable$neutral)
                        mfactorNeutral <- sum(iTable$neutral) / sum(iTable$favorable + iTable$unfavorable + iTable$neutral)
                        iTable[, c("favorableC") := (.SD$favorable + .SD$uninf * mfactorFavorable) * .SD$weight]
                        iTable[, c("unfavorableC") := (.SD$unfavorable + .SD$uninf * mfactorUnfavorable) * .SD$weight]
                        iTable[, c("neutralC") := (.SD$neutral + .SD$uninf * mfactorNeutral) * .SD$weight]
                    }else if(correction.uninf == 2){
                        mfactor <- sum(iTable$favorable + iTable$unfavorable + iTable$neutral + iTable$uninf) / sum(iTable$favorable + iTable$unfavorable + iTable$neutral)
                        iTable[, c("favorableC") :=.SD$favorable * .SD$weight * mfactor]
                        iTable[, c("unfavorableC") :=.SD$unfavorable * .SD$weight * mfactor]
                        iTable[, c("neutralC") :=.SD$neutral * .SD$weight * mfactor]
                    }
                }

                ## Point estimates
                iCount_favorable[1,iEndpoint] <- sum(iTable$favorableC)
                iCount_unfavorable[1,iEndpoint] <- sum(iTable$unfavorableC)
                iCount_neutral[1,iEndpoint] <- sum(iTable$neutralC)
                iCount_uninf[1,iEndpoint] <- sum(iTable$uninfC)

                iLS.table[[iEndpoint]] <- iTable
            }

            ## *** new point estimate
            idelta_netBenefit <- iCount_favorable/iN_pairs-iCount_unfavorable/iN_pairs
            idelta_winRatio <- iCount_favorable/iCount_unfavorable

            iDelta_netBenefit <- cumsum(iCount_favorable[1,]*weight)/iN_pairs-cumsum(iCount_unfavorable[1,]*weight)/iN_pairs
            iDelta_winRatio <- cumsum(iCount_favorable[1,]*weight)/cumsum(iCount_unfavorable[1,]*weight)

            ## *** compute variance            
            iSigma <- inferenceUstatistic(iLS.table, order = order,
                                           weight = weight, count.favorable = iCount_favorable, count.unfavorable = iCount_unfavorable,
                                           n.pairs = iN_pairs, n.C = iN.C, n.T = iN.T,
                                           level.strata = "1", n.strata = 1, n.endpoint = n.endpoint, endpoint = endpoint)$Sigma

        }else{
            iN_pairs <- object$n_pairs

            iCount_favorable <- object$count_favorable
            iCount_unfavorable <- object$count_unfavorable
            iCount_neutral <- object$count_neutral
            iCount_uninf <- object$count_uninf

            idelta_netBenefit <- object$delta_netBenefit
            idelta_winRatio <- object$delta_winRatio
            iDelta_netBenefit <- object$Delta_netBenefit
            iDelta_winRatio <- object$Delta_winRatio

            iSigma <- object$Mvar
        }
            
        ## *** Create object
        out[[iSample]] <- BuyseRes(
            count.favorable = iCount_favorable,      
            count.unfavorable = iCount_unfavorable,
            count.neutral = iCount_neutral,    
            count.uninf = iCount_uninf,
            n.pairs = iN_pairs,
            delta.netBenefit = idelta_netBenefit,
            delta.winRatio = idelta_winRatio,
            Delta.netBenefit = iDelta_netBenefit,
            Delta.winRatio = iDelta_winRatio,
            type = type,
            endpoint = endpoint,
            level.treatment = level.treatment,
            scoring.rule = switch(as.character(scoring.rule),
                                  "0" = "Gehan",
                                  "1" = "Peron"),
            hierarchical = hierarchical,
            correction.uninf = correction.uninf,
            method.inference = method.inference,
            strata = NULL,
            level.strata = "1",
            threshold = threshold,
            n.resampling = as.double(NA),
            deltaResampling.netBenefit = array(dim=c(0,0,0)),
            deltaResampling.winRatio = array(dim=c(0,0,0)),
            DeltaResampling.netBenefit = matrix(NA, nrow = 0, ncol = 0),
            DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
            covarianceResampling = array(NA, dim = c(0,0,0)),
            covariance = iSigma,
            weight = weight,
            iid_favorable = NULL,
            iid_unfavorable = NULL,
            iidNuisance_favorable = NULL,
            iidNuisance_unfavorable = NULL,
            tablePairScore = list(),
            tableSurvival = list()
        )
            
    }

    return(out)
    
}
2
######################################################################
### powerBuyseTest.R ends here
