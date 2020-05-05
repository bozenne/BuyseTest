### powerBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: maj  5 2020 (10:37) 
##           By: Brice Ozenne
##     Update #: 841
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
#' @param null [numeric vector] For each statistic of interest, the null hypothesis to be tested.
#' The vector should be named with the names of the statistics.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' Default value read from \code{BuyseTest.options()}.
#' @param seed [integer, >0] the seed to consider for the simulation study.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param trace [integer] should the execution of the function be traced?
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}.
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#' Default value read from \code{BuyseTest.options()}.
#' @param ... other arguments (e.g. \code{scoring.rule}, \code{method.inference}) to be passed to \code{initializeArgs}.
#' @author Brice Ozenne

## * powerBuyseTest (examples)
##' @rdname powerBuyseTest
##' @examples
##' library(data.table)
##' 
##' #### Using simBuyseTest ####
##' ## only point estimate
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(10, 50, 100), n.rep = 10,
##'                formula = treatment ~ bin(toxicity), seed = 10,
##'                method.inference = "none", trace = 4)
##'
##' ## point estimate with rejection rate
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(10, 50, 100), n.rep = 10,
##'                formula = treatment ~ bin(toxicity), seed = 10,
##'                method.inference = "u-statistic", trace = 4)
##'
##' #### Using user defined simulation function ####
##' ## Example of power calculation for Wilcoxon test
##' simFCT <- function(n.C, n.T){
##'     out <- rbind(cbind(Y=stats::rt(n.C, df = 5), group=0),
##'                  cbind(Y=stats::rt(n.T, df = 5), group=1) + 1)
##'     return(data.table::as.data.table(out))
##' }
##'
##' \dontshow{
##' powerW <- powerBuyseTest(sim = simFCT, sample.size = c(5, 10,20,30,50,100),
##'                          n.rep = 10, formula = group ~ cont(Y))
##' summary(powerW)
##' }
##' \dontrun{
##' powerW <- powerBuyseTest(sim = simFCT, sample.size = c(5, 10,20,30,50,100),
##'                          n.rep = 1000, formula = group ~ cont(Y), cpus = "all")
##' summary(powerW)
##' } 
##' 

## * powerBuyseTest (code)
##' @export
powerBuyseTest <- function(sim,
                           sample.size,
                           sample.sizeC = NULL,
                           sample.sizeT = NULL,
                           n.rep,
                           null = c("netBenefit" = 0),
                           cpus = 1,                          
                           seed = NULL,
                           conf.level = NULL,
                           alternative = NULL,
                           order.Hprojection = NULL,
                           transformation = NULL,
                           trace = 1,
                           ...){

    call <- match.call()$sim

    ## ** normalize and check arguments
    name.call <- names(match.call())
    option <- BuyseTest.options()
    if(is.null(conf.level)){
        conf.level <- option$conf.level
    }
    if(is.null(alternative)){
        alternative <- option$alternative
    }
    if(is.null(order.Hprojection)){
        order.Hprojection <- option$order.Hprojection
    }
    if(is.null(transformation)){
        transformation <- option$transformation
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
    statistic <- names(null)
    validCharacter(statistic,
                   name1 = "names(null)",
                   valid.length = 1:4,
                   valid.values = c("favorable","unfavorable","netBenefit","winRatio"),
                   refuse.NULL = TRUE,
                   refuse.duplicates = TRUE,
                   method = "BuyseTest")
    
    n.sample.size <- length(sample.sizeT)
    sample.sizeCmax <- sample.sizeC[n.sample.size]
    sample.sizeTmax <- sample.sizeT[n.sample.size]

    dt.tempo <- sim(n.C = sample.sizeC[1], n.T = sample.sizeT[1])
    if(!data.table::is.data.table(dt.tempo)){
        stop("The function defined by the argument \'sim\' must return a data.table object.\n")
    }
    
    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(cpus = cpus, option = option, name.call = name.call, 
                              data = NULL, model.tte = NULL, ...)

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = c(outArgs[setdiff(names(outArgs),"data")], list(data = dt.tempo)))        
        if(any(outArgs$operator!=">0")){
            stop("Cannot use argument \'operator\' with powerBuyseTest \n")
        }
        if(!is.null(outArgs$strata)){
            stop("Cannot use argument \'strata\' with powerBuyseTest \n")
        }
        if(outArgs$method.inference %in% c("none","u-statistic") == FALSE){
            stop("Argument \'method.inference\' must be \"none\" or \"u-statistic\" \n")
        }
    }

    cpus <- outArgs$cpus
    outArgs$cpus <- 1
    outArgs$trace <- 0
    
    ## ** initialization data
    outArgs$level.treatment <- levels(as.factor(dt.tempo[[outArgs$treatment]]))
    outArgs$n.strata <- 1
    outArgs$level.strata <- "1"
    outArgs$allstrata <- NULL
    
    ## ** Display
    if (trace > 1) {
        cat("         Simulation study with BuyseTest \n\n")

        if(trace > 2){
            argsInit <- setdiff(names(as.list(args(initializeData))), c("","copy","data"))
            M.status <- do.call(initializeData, args = c(outArgs[argsInit], list(copy = FALSE, data = dt.tempo)))$M.status
            do.call(printGeneral, args = c(outArgs, list(M.status = M.status)))
            if(outArgs$method.inference!="none"){
                do.call(printInference, args = outArgs)
            }
        }
        if(!missing(sample.size) && !is.null(sample.size)){
            text.sample.size <- paste0("   - sample size: ",paste(sample.size, collapse = " "),"\n")
        }else{
            text.sample.size <- paste0("   - sample size: ",paste(sample.sizeC, collapse = " ")," (control)\n",
                                       "                : ",paste(sample.sizeT, collapse = " ")," (treatment)\n")
        }
        cat("Simulation\n",
            "   - repetitions: ",n.rep,"\n",
            "   - cpus       : ",cpus,"\n",
            sep = "")
        cat(" \n")
    }
    
    ## ** define environment
    envirBT <- new.env()
    ## envirBT[[deparse(call)]] <- sim
    name.copy <- c("sim", "option",
                   "outArgs", "sample.sizeTmax", "sample.sizeCmax", "n.sample.size",
                   "sample.sizeC", "sample.sizeT", "n.rep", "seed",
                    "statistic", "null", "conf.level", "alternative", "transformation", "order.Hprojection", 
                   ".powerBuyseTest", ".createSubBT")
    for(iObject in name.copy){ ## iObject <- name.copy[2]
        envirBT[[iObject]] <- eval(parse(text = iObject))
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
                                                 return(.powerBuyseTest(i = X,
                                                                        envir = envirBT,
                                                                        statistic = statistic,
                                                                        null = null,
                                                                        conf.level = conf.level,
                                                                        alternative = alternative,
                                                                        transformation = transformation,
                                                                        order.Hprojection = order.Hprojection))                                                  
                                             })
                                 )
        if(!is.null(seed)){rm(.Random.seed, envir=.GlobalEnv)} # restaure original seed
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
        toExport <- c(".BuyseTest",
                      ".powerBuyseTest",
                      ".createSubBT",
                      "wsumPairScore",
                      "S4BuyseTest",
                      "initializeData",
                      "calcSample",
                      "calcPeron",
                      "pairScore2dt",
                      "inferenceUstatistic",
                      "confint_Ustatistic",
                      ".iid2cov",
                      "validNumeric")

        ## try sim
        test <- try(parallel::clusterCall(cl, fun = function(x){
            sim(n.T = sample.sizeTmax, n.C = sample.sizeCmax)
        }), silent = TRUE)
        if(inherits(test,"try-error")){
            stop(paste0("Could not run argument \'sim\' when using multiple CPUs. \n Consider trying first to run powerBuyseTest with cpus=1. \n If it runs, make sure that \'sim\' does not depend on any variable in the global environment or package without explicit mention of the namespace. \n",test))
        }

        ## run simul
        i <- NULL ## [:forCRANcheck:] foreach
        ls.simulation <- foreach::`%dopar%`(
                                      foreach::foreach(i=1:n.rep, .export = toExport), {                                           
                                          if(trace>0){utils::setTxtProgressBar(pb, i)}
                                          .powerBuyseTest(i = i,
                                                          envir = envirBT,
                                                          statistic = statistic,
                                                          null = null,
                                                          conf.level = conf.level,
                                                          alternative = alternative,
                                                          transformation = transformation,
                                                          order.Hprojection = order.Hprojection)
                                      })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }
    dt.out <- data.table::as.data.table(do.call(rbind, ls.simulation))

    ## ** export
    BuysePower.object <- S4BuysePower(
        alternative = alternative,      
        method.inference = outArgs$method.inference,
        conf.level = conf.level,
        endpoint =  outArgs$endpoint,
        threshold =  outArgs$threshold,
        type =  outArgs$type,
        null = null,
        n.rep = n.rep,      
        results = dt.out
    )

    return(BuysePower.object)
}

## * .powerBuyseTest
.powerBuyseTest <- function(i, envir, statistic, null, conf.level, alternative, transformation, order.Hprojection){

    out <- NULL
    allBT <- vector(mode = "list", length = envir$n.sample.size)
    
    n.endpoint <- length(envir$outArgs$endpoint)
    n.statistic <- length(statistic)
    rerun <- (envir$n.sample.size>1) * (1 + (envir$outArgs$scoring.rule>0)) ## 0 no, 1 yes but only via tablePairScore, 2 yes normal one

    ## ** Initialize data
    data <- envir$sim(n.T = envir$sample.sizeTmax, n.C = envir$sample.sizeCmax)

    out.name <- c("data","M.endpoint","M.status",
                  "index.C","index.T","index.strata",
                  "level.treatment","level.strata", "method.score",
                  "n.strata","n.obs","n.obsStrata","n.obsStrataResampling","cumn.obsStrataResampling","skeletonPeron",
                  "scoring.rule", "iidNuisance", "nUTTE.analyzedPeron_M1", "endpoint.UTTE", "status.UTTE", "D.UTTE","index.UTTE","keep.pairScore")


    envir$outArgs[out.name] <- initializeData(data = data,
                                              type = envir$outArgs$type,
                                              endpoint = envir$outArgs$endpoint,
                                              Uendpoint = envir$outArgs$Uendpoint,
                                              D = envir$outArgs$D,
                                              scoring.rule = envir$outArgs$scoring.rule,
                                              status = envir$outArgs$status,
                                              Ustatus = envir$outArgs$Ustatus,
                                              method.inference = envir$outArgs$method.inference,
                                              operator = envir$outArgs$operator,
                                              censoring = envir$outArgs$censoring,
                                              strata = envir$outArgs$strata,
                                              treatment = envir$outArgs$treatment,
                                              hierarchical = envir$outArgs$hierarchical,
                                              copy = FALSE,
                                              keep.pairScore = envir$outArgs$keep.pairScore,
                                              endpoint.TTE = envir$outArgs$endpoint.TTE,
                                              status.TTE = envir$outArgs$status.TTE,
                                              iidNuisance = envir$outArgs$iidNuisance)

    ## ** Point estimate for the largest sample size
    ## largest sample size
    if(rerun==1){envir$outArgs$keep.pairScore <- TRUE}
    outPoint <- .BuyseTest(envir = envir,
                           method.inference = "none",
                           iid = envir$outArgs$iid,
                           pointEstimation = TRUE)
    keep.args <- c("index.C", "index.T", "type","endpoint","level.strata","level.treatment","scoring.rule","hierarchical","neutral.as.uninf",
                   "correction.uninf","method.inference","method.score","strata","threshold","weight","n.resampling")
    allBT[[envir$n.sample.size]] <- do.call("S4BuyseTest", args = c(outPoint, envir$outArgs[keep.args]))

    ## ** Loop over other sample sizes
    if(rerun>0){

        for(iSize in 1:(envir$n.sample.size-1)){

            if(rerun==1){ ## Gehan's scoring rule
                outPoint <- .createSubBT(allBT[[envir$n.sample.size]], 
                                         sample.sizeC = envir$sample.sizeC[iSize], sample.sizeT = envir$sample.sizeT[iSize])
            }else if(rerun==2){ ## Peron's scoring rule or correction
                envir$outArgs[out.name] <- initializeData(data = rbind(data[envir$outArgs$index.C[1:envir$sample.sizeC[iSize]]],
                                                                       data[envir$outArgs$index.T[1:envir$sample.sizeT[iSize]]]),
                                                          type = envir$outArgs$type,
                                                          endpoint = envir$outArgs$endpoint,
                                                          Uendpoint = envir$outArgs$Uendpoint,
                                                          D = envir$outArgs$D,
                                                          scoring.rule = envir$outArgs$scoring.rule,
                                                          status = envir$outArgs$status,
                                                          Ustatus = envir$outArgs$Ustatus,
                                                          method.inference = envir$outArgs$method.inference,
                                                          operator = envir$outArgs$operator,
                                                          censoring = envir$outArgs$censoring,
                                                          strata = envir$outArgs$strata,
                                                          treatment = envir$outArgs$treatment,
                                                          hierarchical = envir$outArgs$hierarchical,
                                                          copy = FALSE,
                                                          keep.pairScore = envir$outArgs$keep.pairScore,
                                                          endpoint.TTE = envir$outArgs$endpoint.TTE,
                                                          status.TTE = envir$outArgs$status.TTE,
                                                          iidNuisance = envir$outArgs$iidNuisance)

                outPoint <- .BuyseTest(envir = envir,
                                       iid = envir$outArgs$iid,
                                       method.inference = "none",
                                       pointEstimation = TRUE)
            }
            allBT[[iSize]] <- do.call("S4BuyseTest", args = c(outPoint, envir$outArgs[keep.args]))
        }
    }

    ## ** Inference
    for(iSize in 1:envir$n.sample.size){
        for(iStatistic in statistic){
            for(iTransformation in transformation){
                for(iOrder.Hprojection in order.Hprojection){
                    iCI <- confint(allBT[[iSize]],
                                   statistic = iStatistic,
                                   null  = null[iStatistic],
                                   conf.level  = conf.level,
                                   alternative = alternative,
                                   order.Hprojection = iOrder.Hprojection,
                                   transformation = iTransformation)

                    out <- rbind(out,
                                 data.table::data.table(n.T = envir$sample.sizeC[[iSize]],
                                                        n.C = envir$sample.sizeT[[iSize]],
                                                        endpoint = rownames(iCI),
                                                        statistic = iStatistic,
                                                        transformation = iTransformation,
                                                        order.Hprojection = iOrder.Hprojection,
                                                        iCI,
                                                        stringsAsFactors = FALSE)
                                 )
                }
            }
        }
    }
    ## ** Export
    rownames(out) <- NULL
    return(out)
}

## * .createSubBT
.createSubBT <- function(object, sample.sizeT, sample.sizeC){

    ## ** extract information
    endpoint <- object@endpoint
    n.endpoint <- length(endpoint)
    n.strata <- length(object@level.strata)
    weight <- object@weight
    method.inference <- object@method.inference
    correction.uninf <- object@correction.uninf
    tableScore <- object@tablePairScore
    order.Hprojection <- attr(method.inference,"hprojection")
    
    ## ** prepare export
    out <- list(count_favorable = matrix(NA, nrow = 1, ncol = n.endpoint),
                count_unfavorable = matrix(NA, nrow = 1, ncol = n.endpoint),
                count_neutral = matrix(NA, nrow = 1, ncol = n.endpoint),
                count_uninf = matrix(NA, nrow = 1, ncol = n.endpoint),
                delta = array(NA, dim = c(n.strata,n.endpoint,4)),
                Delta = NULL,
                n_pairs = sample.sizeC * sample.sizeT,
                iidAverage_favorable = NULL,
                iidAverage_unfavorable = NULL,
                iidNuisance_favorable = NULL,
                iidNuisance_unfavorable = NULL,
                covariance = NULL,
                tableScore = vector(mode = "list", length = n.endpoint)
                )
            
    ## ** update scores
    for(iEndpoint in 1:n.endpoint){ ## iEndpoint <- 1
        ## cat("endpoint ",iEndpoint,"\n")

        ## restrict pairs 
        index <- intersect(which(tableScore[[iEndpoint]]$indexWithinStrata.C <= sample.sizeC),
                           which(tableScore[[iEndpoint]]$indexWithinStrata.T <= sample.sizeT))
        iTable <- tableScore[[iEndpoint]][index]

        ## update index in dataset
        old.position <- sort(c(unique(iTable$index.T),unique(iTable$index.C)))
        new.position <- rank(old.position)
        old2new <- rep(NA, max(old.position))
        old2new[old.position] <- new.position
        iTable[, c("index.T") := old2new[.SD$index.T]]
        iTable[, c("index.C") := old2new[.SD$index.C]]

        ## update correction (no strata)
        if(correction.uninf>0){
            ## new weighting
            if(iEndpoint>1){
                iTable[, c("weight") := weightM1[match(iTable$index.pair,as.numeric(names(weightM1)))]]
            }
            
            if(correction.uninf == 1){
                mfactorFavorable <- sum(iTable$favorable * iTable$weight) / sum((iTable$favorable + iTable$unfavorable + iTable$neutral) * iTable$weight)
                mfactorUnfavorable <- sum(iTable$unfavorable * iTable$weight) / sum((iTable$favorable + iTable$unfavorable + iTable$neutral) * iTable$weight)
                mfactorNeutral <- sum(iTable$neutral * iTable$weight) / sum((iTable$favorable + iTable$unfavorable + iTable$neutral) * iTable$weight)
                iTable[, c("favorableC") := (.SD$favorable + .SD$uninf * mfactorFavorable) * .SD$weight]
                iTable[, c("unfavorableC") := (.SD$unfavorable + .SD$uninf * mfactorUnfavorable) * .SD$weight]
                iTable[, c("neutralC") := (.SD$neutral + .SD$uninf * mfactorNeutral) * .SD$weight]
            }else if(correction.uninf == 2){                
                mfactor <- sum(iTable$favorable + iTable$unfavorable + iTable$neutral + iTable$uninf) / sum(iTable$favorable + iTable$unfavorable + iTable$neutral)
                iTable[, c("favorableC") :=.SD$favorable * mfactor * .SD$weight]
                iTable[, c("unfavorableC") :=.SD$unfavorable * mfactor * .SD$weight]
                iTable[, c("neutralC") :=.SD$neutral * mfactor * .SD$weight]
            }
            weightM1 <-  setNames(iTable$neutralC,iTable$index.pair)
        }
        ## new counts
        out$count_favorable[1,iEndpoint] <- sum(iTable$favorableC)
        out$count_unfavorable[1,iEndpoint] <- sum(iTable$unfavorableC)
        out$count_neutral[1,iEndpoint] <- sum(iTable$neutralC)
        out$count_uninf[1,iEndpoint] <- sum(iTable$uninfC)

        out$tableScore[[iEndpoint]] <- iTable
    }

    ## new point estimate
    out$delta[,,1] <- out$count_favorable/out$n_pairs
    out$delta[,,2] <- out$count_unfavorable/out$n_pairs
    out$delta[,,3] <- out$delta[,,1]-out$delta[,,2]
    out$delta[,,4] <- out$delta[,,1]/out$delta[,,2]
    
    out$Delta <- cbind(cumsum(out$count_favorable[1,]*weight)/out$n_pairs,
                       cumsum(out$count_unfavorable[1,]*weight)/out$n_pairs,
                       cumsum(out$count_favorable[1,] - out$count_unfavorable[1,])*weight/out$n_pairs,
                       cumsum(out$count_favorable[1,])/cumsum(out$count_unfavorable[1,]))

    ## compute variance
    if(method.inference == "u-statistic"){
        iInference <- inferenceUstatistic(out$tableScore, order = order.Hprojection,
                                          weight = weight, count.favorable = out$count_favorable, count.unfavorable = out$count_unfavorable,
                                          n.pairs = out$n_pairs, n.C = sample.sizeC, n.T =sample.sizeT, 
                                          level.strata = "1", n.strata = 1, n.endpoint = n.endpoint, endpoint = endpoint)

        ##  extract variance
        out$iidAverage_favorable <- iInference$iidAverage_favorable
        out$iidAverage_unfavorable <- iInference$iidAverage_unfavorable
        out$covariance <- iInference$covariance
    }

    ## ** export
    return(out)
    
}

######################################################################
### powerBuyseTest.R ends here
