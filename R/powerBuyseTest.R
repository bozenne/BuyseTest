### powerBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: mar 20 2023 (14:08) 
##           By: Brice Ozenne
##     Update #: 1076
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
#' Returns estimates, their standard deviation, the average estimated standard error, and the rejection rate.
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
#' @param cpus [integer, >0] the number of CPU to use. Default value is 1.
#' @param seed [integer, >0] the seed to consider for the simulation study.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric, 0-1] type 1 error level.
#' Default value read from \code{BuyseTest.options()}.
#' @param power [numeric, 0-1] type 2 error level used to determine the sample size. Only relevant when \code{sample.size}, \code{sample.sizeC}, and \code{sample.sizeT} are not given. See details.
#' @param max.sample.size [interger, 0-1] sample size used to approximate the sample size achieving the requested type 1 and type 2 error (see details).
#' Can have length 2 to indicate the sample in each group when the groups have unequal sample size.
#' @param trace [integer] should the execution of the function be traced?
#' @param transformation [logical] should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}.
#' @param order.Hprojection [integer 1,2] the order of the H-project to be used to compute the variance of the net benefit/win ratio.
#' Default value read from \code{BuyseTest.options()}.
#' @param ... other arguments (e.g. \code{scoring.rule}, \code{method.inference}) to be passed to \code{initializeArgs}.
#'
#' @details \bold{Sample size calculation}: to approximate the sample size achieving the requested type 1 (\eqn{\alpha}) and type 2 error (\eqn{\beta}),
#' GPC are applied on a large sample (as defined by the argument \code{max.sample.size}): \eqn{N^*=m^*+n^*} where \eqn{m^*} is the sample size in the control group and \eqn{n^*} is the sample size in the active group.
#' Then the effect (\eqn{\delta}) and the asymptotic variance of the estimator (\eqn{\sigma^2}) are estimated. The total sample size is then deduced as (two-sided case):
#' \deqn{\hat{N} = \hat{\sigma}^2\frac{(u_{1-\alpha/2}+u_{1-\beta})^2}{\hat{\delta}^2}} from which the group specific sample sizes are deduced: \eqn{\hat{m}=\hat{N}\frac{m^*}{N^*}} and \eqn{\hat{n}=\hat{N}\frac{n^*}{N^*}}. Here \(u_x\) denotes the x-quantile of the normal distribution.
#' Since this is an approximation, a simulation study is performed with this sample size to provide the estimated power. It may not be exactly the requested power but should provide a reasonnable guess which can be refined with further simulation studies. Note that the maximumal sample size should be very high for the estimation to be accurate (ideally 10000 or more).
#' 
#' @author Brice Ozenne

## * powerBuyseTest (examples)
##' @rdname powerBuyseTest
##' @examples
##' library(data.table)
##' 
##' #### Using simBuyseTest ####
##' ## only point estimate
##' \dontrun{
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(10, 25, 50, 75, 100), 
##'                formula = treatment ~ bin(toxicity), seed = 10, n.rep = 1000,
##'                method.inference = "none", trace = 2, keep.pairScore = FALSE)
##' }
##' 
##' ## point estimate with rejection rate
##' \dontshow{
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(10, 50, 100), 
##'                formula = treatment ~ bin(toxicity), seed = 10, n.rep = 10,
##'                method.inference = "u-statistic", trace = 4)
##' }
##' \dontrun{
##' powerBuyseTest(sim = simBuyseTest, sample.size = c(10, 50, 100), 
##'                formula = treatment ~ bin(toxicity), seed = 10, n.rep = 1000,
##'                method.inference = "u-statistic", trace = 4)
##' }
##'
##' #### Using user defined simulation function ####
##' ## power calculation for Wilcoxon test
##' simFCT <- function(n.C, n.T){
##'     out <- rbind(cbind(Y=stats::rt(n.C, df = 5), group=0),
##'                  cbind(Y=stats::rt(n.T, df = 5), group=1) + 1)
##'     return(data.table::as.data.table(out))
##' }
##' simFCT2 <- function(n.C, n.T){
##'     out <- rbind(cbind(Y=stats::rt(n.C, df = 5), group=0),
##'                  cbind(Y=stats::rt(n.T, df = 5), group=1) + 0.25)
##'     return(data.table::as.data.table(out))
##' }
##'
##' \dontshow{
##' powerW <- powerBuyseTest(sim = simFCT, sample.size = c(5, 10,20,30,50,100),
##'                          n.rep = 10, formula = group ~ cont(Y))
##' summary(powerW)
##' }
##' \dontrun{
##' powerW <- powerBuyseTest(sim = simFCT, sample.size = c(5,10,20,30,50,100),
##'                          n.rep = 1000, formula = group ~ cont(Y), cpus = "all")
##' summary(powerW)
##' } 
##'
##' ## sample size needed to reach (approximately) a power
##' ## based on summary statistics obtained on a large sample 
##' \dontrun{
##' sampleW <- powerBuyseTest(sim = simFCT, power = 0.8, max.sample.size = 10000,
##'                          n.rep = 1000, formula = group ~ cont(Y), cpus = "all")
##' summary(sampleW) ## not very accurate but gives an order of magnitude
##' 
##' sampleW2 <- powerBuyseTest(sim = simFCT2,
##'                            power = 0.8, max.sample.size = 10000,
##'                            n.rep = 1000, formula = group ~ cont(Y), cpus = "all")
##' summary(sampleW2) ## more accurate with larger samples
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
                           power = NULL,
                           max.sample.size = 5000,
                           alternative = NULL,
                           order.Hprojection = NULL,
                           transformation = NULL,
                           trace = 1,
                           ...){

    call <- match.call()
    
    ## ** normalize and check arguments
    name.call <- names(call)
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
    
    if(is.null(sample.sizeC) && is.null(sample.sizeT)){
        if((missing(sample.size) || is.null(sample.size)) && !is.null(power)){
            if(length(max.sample.size)==1){
                max.sample.size <- c(max.sample.size,max.sample.size)
            }
            if (trace > 1) {
                cat("         Determination of the sample using a large sample (m=",max.sample.size[1],", n=",max.sample.size[2],")  \n\n",sep="")
            }
            e.BTmax <- BuyseTest(..., data = sim(n.C = max.sample.size[1], n.T = max.sample.size[2]), trace = 0)
            if(e.BTmax@method.inference == "u-statistic"){
                DeltaMax <- coef(e.BTmax, statistic = names(null)) - null
                IidMax <- getIid(e.BTmax, statistic = names(null), scale = FALSE)
                ratio <- c(max.sample.size[1]/sum(max.sample.size),max.sample.size[2]/sum(max.sample.size))
                sigma2Max <- mean(IidMax[attr(e.BTmax@level.treatment,"indexC"),]^2)/ratio[1] + mean(IidMax[attr(e.BTmax@level.treatment,"indexT"),]^2)/ratio[2]
                if(alternative=="two.sided"){
                    n.approx <- sigma2Max*(stats::qnorm(1-alpha/2) + stats::qnorm(power))^2/DeltaMax^2
                }else if(alternative=="less"){
                    if(DeltaMax<0){
                        n.approx <- sigma2Max*(stats::qnorm(1-alpha) + stats::qnorm(power))^2/DeltaMax^2
                    }else{
                        message("No power: positive effect detected. \n")
                        return(invisible(DeltaMax))
                    }
                }else if(alternative=="greater"){
                    if(DeltaMax>0){
                        n.approx <- sigma2Max*(stats::qnorm(1-alpha) + stats::qnorm(power))^2/DeltaMax^2
                    }else{
                        message("No power: positive effect detected. \n")
                        return(invisible(DeltaMax))
                    }
                }
                sample.sizeC <- ceiling(n.approx*ratio[1])
                sample.sizeT <- ceiling(n.approx*ratio[2])

                ## (mean(IidMax[attr(e.BTmax@level.treatment,"indexC"),]^2) + mean(IidMax[attr(e.BTmax@level.treatment,"indexT"),]^2))*(stats::qnorm(1-alpha/2)+stats::qnorm(power))^2/DeltaMax^2
                if (trace > 1) {
                    cat("   - estimated effect (variance): ",as.double(DeltaMax)," (",sigma2Max,")\n",sep="")
                    cat("   - estimated sample size      : (m=",sample.sizeC,", n=",sample.sizeT,")\n\n",sep="")
                }
            
            }
        }else{
            sample.sizeC <- sample.size
            sample.sizeT <- sample.size
        }
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
    
    dt.tempo <- sim(n.C = sample.sizeC[1], n.T = sample.sizeT[1])
    if(!data.table::is.data.table(dt.tempo)){
        stop("The function defined by the argument \'sim\' must return a data.table object.\n")
    }

    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(cpus = cpus, option = option, name.call = name.call, 
                              data = NULL, model.tte = NULL, ...)
    outArgs$call <- setNames(as.list(call),names(call))
    ## if((outArgs$keep.pairScore == FALSE) && ("keep.pairScore" %in% names(call) == FALSE) && (outArgs$method.inference %in% c("none","u-statistic")) && (outArgs$correction.uninf==0) && all(outArgs$scoring.rule<=0) ){
    ##     outArgps$keep.pairScore <- TRUE
    ## }

    ## ** test arguments
    if(option$check){
        index.pb <- which(outArgs$status[outArgs$type=="tte"] == "..NA..")
        if(length(index.pb)>0){
            if(any(attr(outArgs$censoring,"original") %in% c("left","right") == FALSE)){
                stop("BuyseTest: wrong specification of \'status\'. \n",
                     "\'status\' must indicate a variable in data for TTE endpoints. \n",
                     "\'censoring\' is used to indicate whether there is left or right censoring. \n",
                     "Consider changing \'censoring =\' into \'status =\' when in the argument \'formula\' \n")
            }else{        
                stop("BuyseTest: wrong specification of \'status\'. \n",
                     "\'status\' must indicate a variable in data for TTE endpoints. \n",
                     "TTE endoints: ",paste(outArgs$endpoint[outArgs$type=="tte"],collapse=" "),"\n",
                     "proposed \'status\' for these endoints: ",paste(outArgs$status[outArgs$type=="tte"],collapse=" "),"\n")
            }
        }
        ## outTest <- do.call(testArgs, args = outArgs)        
        ##     outTest <- do.call(testArgs, args = c(outArgs[setdiff(names(outArgs),"data")], list(data = dt.tempo)))        
        ##     if(!is.null(outArgs$strata)){
        ##         stop("Cannot use argument \'strata\' with powerBuyseTest \n")
        ##     }
        ##     if(outArgs$method.inference %in% c("none","u-statistic") == FALSE){
        ##         stop("Argument \'method.inference\' must be \"none\" or \"u-statistic\" \n")
        ##     }
    }

    cpus <- outArgs$cpus
    outArgs$cpus <- 1
    outArgs$trace <- 0

    ## ** initialization sample size and data
    n.sample.size <- length(sample.sizeT)
    sample.sizeCmax <- sample.sizeC[n.sample.size]
    sample.sizeTmax <- sample.sizeT[n.sample.size]

    outArgs$level.treatment <- levels(as.factor(dt.tempo[[outArgs$treatment]]))
    ## outArgs$n.strata <- 1
    ## outArgs$level.strata <- "1"
    ## outArgs$allstrata <- NULL

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
                   ".powerBuyseTest")
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
                      "wsumPairScore",
                      "S4BuyseTest",
                      "initializeData",
                      "calcSample",
                      "calcPeron",
                      "pairScore2dt",
                      "confint_Ustatistic",
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
        restriction =  outArgs$restriction,
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
    
    scoring.rule <- envir$outArgs$scoring.rule
    iidNuisance <- envir$outArgs$iidNuisance
    n.endpoint <- length(envir$outArgs$endpoint)
    n.statistic <- length(statistic)
    rerun <- (envir$n.sample.size>1)

    ## when creating S4 object
    keep.args <- c("index.C", "index.T", "index.strata", "type","endpoint","level.strata","level.treatment","scoring.rule","hierarchical","neutral.as.uninf","add.halfNeutral",
                   "correction.uninf","method.inference","method.score","strata","threshold","restriction","weightObs","weightEndpoint","pool.strata","n.resampling","call")

    ## ** Simulate data
    data <- envir$sim(n.T = envir$sample.sizeTmax, n.C = envir$sample.sizeCmax)

    iInit <- initializeData(data = data,
                            type = envir$outArgs$type,
                            endpoint = envir$outArgs$endpoint,
                            Uendpoint = envir$outArgs$Uendpoint,
                            D = envir$outArgs$D,
                            scoring.rule = envir$outArgs$scoring.rule,
                            status = envir$outArgs$status,
                            Ustatus = envir$outArgs$Ustatus,
                            method.inference = envir$outArgs$method.inference,
                            censoring = envir$outArgs$censoring,
                            strata = envir$outArgs$strata,
                            treatment = envir$outArgs$treatment,
                            hierarchical = envir$outArgs$hierarchical,
                            copy = FALSE,
                            keep.pairScore = envir$outArgs$keep.pairScore,
                            endpoint.TTE = envir$outArgs$endpoint.TTE,
                            status.TTE = envir$outArgs$status.TTE,
                            iidNuisance = envir$outArgs$iidNuisance)
    out.name <- names(iInit)
    envir$outArgs[out.name] <- iInit

    ## save for subsetting the data set with other sample sizes
    index.C <- envir$outArgs$index.C
    index.T <- envir$outArgs$index.T

    ## ** Point estimate for the largest sample size
    if(envir$outArgs$method.inference %in% c("none","u-statistic")){
        ## compute estimate and possibly uncertainty
        outPoint <- .BuyseTest(envir = envir,
                               method.inference = envir$outArgs$method.inference,
                               iid = envir$outArgs$iid,
                               pointEstimation = TRUE)

        ## create S4 object
        allBT[[envir$n.sample.size]] <- do.call("S4BuyseTest", args = c(outPoint, envir$outArgs[keep.args]))
    }else{
        data[["..strata.."]]  <- NULL
        data[["..rowIndex.."]]  <- NULL
        data[["..NA.."]]  <- NULL
        allBT[[envir$n.sample.size]] <- BuyseTest(data = data, scoring.rule = envir$outArgs$scoring.rule, pool.strata = envir$outArgs$pool.strata, correction.uninf = envir$outArgs$correction.uninf, 
                                                  model.tte = envir$outArgs$model.tte, method.inference = envir$outArgs$method.inference, n.resampling = envir$outArgs$n.resampling, 
                                                  strata.resampling = envir$outArgs$strata.resampling, hierarchical = envir$outArgs$hierarchical, weightEndpoint = envir$outArgs$weightEndpoint, 
                                                  neutral.as.uninf = envir$outArgs$neutral.as.uninf, add.halfNeutral = envir$outArgs$add.halfNeutral, 
                                                  trace = FALSE, treatment = envir$outArgs$treatment, endpoint = envir$outArgs$endpoint, 
                                                  type = envir$outArgs$type, threshold = envir$outArgs$threshold, restriction = envir$outArgs$restriction, status = envir$outArgs$status, operator = envir$outArgs$operator, 
                                                  censoring = envir$outArgs$censoring, strata = envir$outArgs$strata)
    }
    
    ## ** Loop over other sample sizes
    if(rerun>0){
        ## test.bebu <- envir$outArgs$keep.pairScore && (envir$outArgs$method.inference %in% c("none","u-statistic")) && all(scoring.rule <= 4) && (envir$outArgs$correction.uninf == 0)

        for(iSize in 1:(envir$n.sample.size-1)){

            ## if(test.bebu){ REMOVED AS IT IS SLOWER TO KEEP pairScore THAN RE-RUN THE C++ CODE

            ##     outCov2 <- inferenceUstatisticBebu(tablePairScore = allBT[[envir$n.sample.size]]@tablePairScore,
            ##                                        subset.C = 1:envir$sample.sizeC[iSize],
            ##                                        subset.T = 1:envir$sample.sizeT[iSize],
            ##                                        order = envir$outArgs$order.Hprojection,
            ##                                        weightEndpoint = envir$outArgs$weightEndpoint,
            ##                                        n.pairs = envir$sample.sizeC[iSize]*envir$sample.sizeT[iSize],
            ##                                        n.C = envir$sample.sizeC[iSize],
            ##                                        n.T = envir$sample.sizeT[iSize],
            ##                                        level.strata = envir$outArgs$level.strata,
            ##                                        n.strata = envir$outArgs$n.strata,
            ##                                        n.endpoint = n.endpoint,
            ##                                        endpoint = envir$outArgs$endpoint)
                
            ##     outPoint2 <- list(count_favorable = outCov2$count_favorable,
            ##                       count_unfavorable = outCov2$count_unfavorable,
            ##                       count_neutral = outCov2$count_neutral,
            ##                       count_uninf = outCov2$count_uninf, ## outPoint$count_uninf
            ##                       delta = outCov2$delta,
            ##                       Delta = outCov2$Delta, ## outPoint$Delta
            ##                       n_pairs = outCov2$n.pairs,
            ##                       iidAverage_favorable = matrix(nrow = 0, ncol = 0),
            ##                       iidAverage_unfavorable = matrix(nrow = 0, ncol = 0),
            ##                       iidAverage_neutral = matrix(nrow = 0, ncol = 0),
            ##                       iidNuisance_favorable = matrix(nrow = 0, ncol = 0),
            ##                       iidNuisance_unfavorable = matrix(nrow = 0, ncol = 0),
            ##                       iidNuisance_neutral = matrix(nrow = 0, ncol = 0),
            ##                       covariance = outCov2$Sigma,
            ##                       tableScore = list()                                  
            ##                       )

            ##     allBT[[iSize]] <- do.call("S4BuyseTest", args = c(outPoint2, envir$outArgs[keep.args]))
            ## }else{
                iData <- rbind(data[index.C[1:envir$sample.sizeC[iSize]]],
                               data[index.T[1:envir$sample.sizeT[iSize]]])

                if(envir$outArgs$method.inference %in% c("none","u-statistic")){
                    envir$outArgs[out.name] <- initializeData(data = iData,
                                                              type = envir$outArgs$type,
                                                              endpoint = envir$outArgs$endpoint,
                                                              Uendpoint = envir$outArgs$Uendpoint,
                                                              D = envir$outArgs$D,
                                                              scoring.rule = scoring.rule,
                                                              status = envir$outArgs$status,
                                                              Ustatus = envir$outArgs$Ustatus,
                                                              method.inference = envir$outArgs$method.inference,
                                                              censoring = envir$outArgs$censoring,
                                                              strata = envir$outArgs$strata,
                                                              treatment = envir$outArgs$treatment,
                                                              hierarchical = envir$outArgs$hierarchical,
                                                              copy = FALSE,
                                                              keep.pairScore = envir$outArgs$keep.pairScore,
                                                              endpoint.TTE = envir$outArgs$endpoint.TTE,
                                                              status.TTE = envir$outArgs$status.TTE,
                                                              iidNuisance = iidNuisance)

                    outPoint <- .BuyseTest(envir = envir,
                                           iid = envir$outArgs$iid,
                                           method.inference = envir$outArgs$method.inference,
                                           pointEstimation = TRUE)

                    allBT[[iSize]] <- do.call("S4BuyseTest", args = c(outPoint, envir$outArgs[keep.args]))
                }else{
                    iData[["..strata.."]]  <- NULL
                    iData[["..rowIndex.."]]  <- NULL
                    iData[["..NA.."]]  <- NULL
                    allBT[[iSize]] <- BuyseTest(data = iData,
                                                scoring.rule = envir$outArgs$scoring.rule,
                                                pool.strata = envir$outArgs$pool.strata,
                                                correction.uninf = envir$outArgs$correction.uninf, 
                                                model.tte = envir$outArgs$model.tte,
                                                method.inference = envir$outArgs$method.inference,
                                                n.resampling = envir$outArgs$n.resampling, 
                                                strata.resampling = envir$outArgs$strata.resampling,
                                                hierarchical = envir$outArgs$hierarchical,
                                                weightEndpoint = envir$outArgs$weightEndpoint, 
                                                neutral.as.uninf = envir$outArgs$neutral.as.uninf, 
                                                add.halfNeutral = envir$outArgs$add.halfNeutral, 
                                                trace = FALSE,
                                                treatment = envir$outArgs$treatment,
                                                endpoint = envir$outArgs$endpoint, 
                                                type = envir$outArgs$type,
                                                threshold = envir$outArgs$threshold,
                                                restriction = envir$outArgs$restriction,
                                                status = envir$outArgs$status,
                                                operator = envir$outArgs$operator, 
                                                censoring = envir$outArgs$censoring,
                                                strata = envir$outArgs$strata)
                }
            ## }
        }
    }

    ## ** Inference
    for(iSize in 1:envir$n.sample.size){
        for(iStatistic in statistic){
            for(iTransformation in transformation){
                for(iOrder.Hprojection in order.Hprojection){
                    iCI <- suppressMessages(confint(allBT[[iSize]],
                                                    statistic = iStatistic,
                                                    null  = null[iStatistic],
                                                    conf.level  = conf.level,
                                                    alternative = alternative,
                                                    order.Hprojection = iOrder.Hprojection,
                                                    transformation = iTransformation))

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


######################################################################
### powerBuyseTest.R ends here
