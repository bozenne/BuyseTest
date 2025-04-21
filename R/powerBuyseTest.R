### powerBuyseTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2018 (12:57) 
## Version: 
## Last-Updated: Apr 21 2025 (11:47) 
##           By: Brice Ozenne
##     Update #: 1333
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
#' Can also be use for power calculation or to approximate the sample size needed to reach a specific power.
#'
#' @param formula [formula] a symbolic description of the GPC model,
#' typically \code{treatment ~ type1(endpoint1) + type2(endpoint2, threshold2) + strata}.
#' See Details in the documentation of the \code{\link{BuyseTest}} function, section "Specification of the GPC model".
#' @param sim [function] take two arguments:
#' the sample size in the control group (\code{n.C}) and the sample size in the treatment group (\code{n.C})
#' and generate datasets. The datasets must be data.frame objects or inherits from data.frame.
#' @param sample.size [integer vector or matrix, >0] the group specific sample sizes relative to which the simulations should be perform.
#' When a vector, the same sample size is used for each group. Alternatively can be a matrix with two columns, one for each group (respectively T and C).
#' @param n.rep [integer, >0] the number of simulations.
#' When specifying the power instead of the sample size, should be a vector of length 2 where the second element indicates the number of simulations used to identify the sample size.
#' @param null [numeric vector] For each statistic of interest, the null hypothesis to be tested.
#' The vector should be named with the names of the statistics.
#' @param cpus [integer, >0] the number of CPU to use. Default value is 1.
#' @param export.cpus [character vector] name of the variables to export to each cluster.
#' @param seed [integer, >0] Random number generator (RNG) state used when starting the simulation study.
#' If \code{NULL} no state is set.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric, 0-1] type 1 error level.
#' Default value read from \code{BuyseTest.options()}.
#' @param power [numeric, 0-1] type 2 error level used to determine the sample size. Only relevant when \code{sample.size} is not given. See details.
#' @param max.sample.size [interger, 0-1] sample size used to approximate the sample size achieving the requested type 1 and type 2 error (see details).
#' Can have length 2 to indicate the sample in each group (respectively T and C) when the groups have unequal sample size.
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
#' \deqn{\hat{N} = \hat{\sigma}^2\frac{(u_{1-\alpha/2}+u_{1-\beta})^2}{\hat{\delta}^2}} from which the group specific sample sizes are deduced: \eqn{\hat{m}=\hat{N}\frac{m^*}{N^*}} and \eqn{\hat{n}=\hat{N}\frac{n^*}{N^*}}. Here \eqn{u_x} denotes the x-quantile of the normal distribution. \cr
#' This approximation can be improved by increasing the sample size (argument \code{max.sample.size}) and/or by performing it multiple times based on a different dataset and average estimated sample size per group (second element of argument \code{n.rep}). \cr
#' To evaluate the approximation, a simulation study is then performed with the estimated sample size. It will not exactly match the requested power but should provide a reasonnable guess which can be refined with further simulation studies. The larger the sample size (and/or number of CPUs) the more accurate the approximation.
#'
#' \bold{seed}: the seed is used to generate one seed per simulation. These simulation seeds are the same whether one or several CPUs are used.
#' 
#' @return An S4 object of class  \code{\linkS4class{S4BuysePower}}.
#' @keywords htest
#' 
#' @author Brice Ozenne

## * powerBuyseTest (examples)
##' @rdname powerBuyseTest
##' @examples
##' library(data.table)
##' 
##' #### Using simBuyseTest ####
##' ## save time by not generating TTE outcomes
##' simBuyseTest2 <- function(...){simBuyseTest(..., argsCont = NULL, argsTTE = NULL)}
##' 
##' ## only point estimate
##' \dontrun{
##' pBT <- powerBuyseTest(sim = simBuyseTest2, sample.size = c(10, 25, 50, 75, 100), 
##'                   formula = treatment ~ bin(toxicity), seed = 10, n.rep = 1000,
##'                   method.inference = "none", keep.pairScore = FALSE, cpus = 5)
##' summary(pBT)
##' model.tables(pBT)
##' }
##' 
##' ## point estimate with rejection rate
##' \dontshow{
##' powerBuyseTest(sim = simBuyseTest2, sample.size = c(10, 50, 100), 
##'                formula = treatment ~ bin(toxicity), seed = 10, n.rep = 10,
##'                method.inference = "u-statistic", trace = 4)
##' }
##' \dontrun{
##' powerBuyseTest(sim = simBuyseTest2, sample.size = c(10, 50, 100), 
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
##' sampleW <- powerBuyseTest(sim = simFCT, power = 0.8, formula = group ~ cont(Y), 
##'                          n.rep = c(1000,10), max.sample.size = 2000, cpus = 5,
##'                          seed = 10)
##' nobs(sampleW)
##' summary(sampleW) ## not very accurate but gives an order of magnitude
##' 
##' sampleW2 <- powerBuyseTest(sim = simFCT2, power = 0.8, formula = group ~ cont(Y), 
##'                          n.rep = c(1000,10), max.sample.size = 2000, cpus = 5,
##'                          seed = 10)
##' summary(sampleW2) ## more accurate when the sample size needed is not too small
##' }
##'

## * powerBuyseTest (code)
##' @export
powerBuyseTest <- function(formula,
                           sim,
                           sample.size,
                           n.rep = c(1000,10),
                           null = c("netBenefit" = 0),
                           cpus = 1,
                           export.cpus = NULL,
                           seed = NULL,
                           conf.level = NULL,
                           power = NULL,
                           max.sample.size = 2000,
                           alternative = NULL,
                           order.Hprojection = NULL,
                           transformation = NULL,
                           trace = 1,
                           ...){

    mycall <- match.call()

    ## ** normalize and check arguments
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
    outArgs <- initializeArgs(formula = formula, cpus = cpus, option = option, call = mycall, 
                              data = NULL, model.tte = NULL, ...)
    outArgs$call <- setNames(as.list(mycall),names(mycall))

    ## power
    if(!is.null(power) && (!missing(sample.size) && !is.null(sample.size))){
        warning("Argument power is disregarded when arguments \'sample.size\' is specified. \n")
        power <- NULL
    }else if(is.null(power) && (missing(sample.size) || is.null(sample.size))){
        stop("Argument \'sample.size\' or argument \'power\' must be specified. \n")
    }

    ## sample size
    if(is.null(power)){
        if(!is.numeric(sample.size) || (!is.vector(sample.size) && !is.matrix(sample.size))){
            stop("Argument \'sample.size\' must be a vector of integers or a matrix of integers with two-columns (one for each group). \n")
        }
        if(any(sample.size<=0) || any(sample.size %% 1 != 0)){
            stop("Argument \'sample.size\' must only contain strictly positive integers. \n")
        }
        if(is.matrix(sample.size)){
            if(NCOL(sample.size)!=2){
                stop("When a matrix, argument \'sample.size\' must have two-columns (one for each group). \n")
            }
            if(is.null(colnames(sample.size))){
                sample.sizeT <- sample.size[,1]
                sample.sizeC <- sample.size[,2]
            }else{
                if(any(c("C","T") %in% colnames(sample.size) == FALSE)){
                    stop("When a matrix, argument \'sample.size\' must have column names \"C\" and \"T\" \n")
                }
                sample.sizeT <- sample.size[,"T"]
                sample.sizeC <- sample.size[,"C"]
            }
        }else{
            sample.sizeT <- sample.size
            sample.sizeC <- sample.size
        }

    }else{
        if(length(n.rep)==1){
            rep2rep <- function(x){sapply(x, function(iX){ceiling(log10(iX) + 3*pmax(0,log10(iX)-1) + pmax(0,log10(iX)-2) + 5*pmax(0,log10(iX)-3))})}
            ## rep2rep(10^(1:5))
            ## [1]  1  5 10 20 30

            n.rep <- 10000
            n.rep <- c(n.rep, ceiling(log10(n.rep) + 3*pmax(0,log10(n.rep)-1) + pmax(0,log10(n.rep)-2) + 5*pmax(0,log10(n.rep)-3)))

        }
        if(!is.vector(max.sample.size)){
            stop("Argument \'max.sample.size\' must be a vector. \n")
        }
        if(length(max.sample.size) %in% 1:2 == FALSE){
            stop("Argument \'max.sample.size\' must have length 1 or 2. \n")
        }
        if(length(max.sample.size) == 2){
            if(is.null(names(max.sample.size))){
                names(max.sample.size) <- c("T","C")
            }else if(any(c("C","T") %in% names(max.sample.size) == FALSE)){
                stop("When a vector of length 2, argument \'max.sample.size\' must have names T and C. \n")
            }
        }
    }
    
    ## statistic
    statistic <- names(null)
    validCharacter(statistic,
                   name1 = "names(null)",
                   valid.length = 1:4,
                   valid.values = c("favorable","unfavorable","netBenefit","winRatio"),
                   refuse.NULL = TRUE,
                   refuse.duplicates = TRUE,
                   method = "BuyseTest")

    ## sim
    dt.tempo <- sim(n.C = 10, n.T = 10)
    if(!inherits(dt.tempo, "data.frame")){
        stop("The function defined by the argument \'sim\' must return a data.frame or an object that inherits from data.frame.\n")
    }

    ## cluster
    if(identical(cpus,"all")){
        cpus <- parallel::detectCores()
    }else if(cpus>1){
        validInteger(cpus,
                     valid.length = 1,
                     min = 1,
                     max = parallel::detectCores(),
                     method = "powerBuyseTest")
    }

    ## seed
    if (!is.null(seed)) {
        tol.seed <- 10^(floor(log10(.Machine$integer.max))-1)
        if(n.rep[1]>tol.seed){
            stop("Cannot set a seed per simulation when considering more than ",tol.seed," similations. \n")
        }
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)
        seqSeed <- sample.int(tol.seed, max(n.rep),  replace = FALSE)        
    }else{
        seqSeed <- NULL
    }

    ## trace
    if (trace > 0) {
        requireNamespace("pbapply")
        method.loop <- pbapply::pblapply
    }else{
        method.loop <- lapply
    }
    
    ## ** initialize cluster
    if(cpus>1){
        cl <- parallel::makeCluster(cpus)
        ## link to foreach
        doSNOW::registerDoSNOW(cl)
        ## export
        if(!is.null(export.cpus)){
            parallel::clusterExport(cl, export.cpus)
        }
        ## seed
        if (!is.null(seed)) {
            parallel::clusterExport(cl, varlist = "seqSeed", envir = environment())
        }         
        ## export all BuyseTest functions (except C++)
        BT.fct <- ls(getNamespace("BuyseTest"))
        toExport <- c(".BuyseTest", ".powerBuyseTest", "wsumPairScore", "S4BuyseTest", "initializeData", "calcSample", "calcPeron", "pairScore2dt", "confint_Ustatistic",
                      grep("^valid",BT.fct,value = TRUE))

        ## parallel::clusterExport(cl, varlist = c(BT.fct[!grepl("_cpp$",BT.fct)], ".powerBuyseTest", ".BuyseTest"))
        ## [:forCRANcheck:] foreach
        iB <- NULL        
    }

    ## ** initialize sample size
    if(!is.null(power)){
            
        if(length(max.sample.size)==1){
            max.sample.size <- c(T = max.sample.size, C = max.sample.size)
        }
            
        if (trace > 1) {
            cat("         Determination of the sample using a large sample (T=",max.sample.size[1],", C=",max.sample.size[2],")  \n\n",sep="")
        }

        if (cpus == 1) {            
            ls.BTmax <- do.call(method.loop,
                                args = list(X = 1:n.rep[2],
                                            FUN = function(X){
                                                if(!is.null(seed)){set.seed(seqSeed[X])}
                                                iOut <- BuyseTest(formula, data = sim(n.T = max.sample.size["T"], n.C = max.sample.size["C"]), trace = 0, ...)
                                                return(iOut)
                                            })
                                )
        }else if(cpus > 1){

            ## split into a 100 jobs
            split.resampling <- parallel::splitIndices(nx = n.rep[2], ncl = min(max(100,10*cpus), n.rep[2]))
            nsplit.resampling <- length(split.resampling)

            ## define progress bar
            if(trace>0){
                pb <- utils::txtProgressBar(max = nsplit.resampling, style = 3)          
                progress <- function(n){utils::setTxtProgressBar(pb, n)}
                opts <- list(progress = progress)
            }else{
                opts <- list()
            }

            ## export functions
            ls2.BTmax <- foreach::`%dopar%`(
                                      foreach::foreach(iB=1:nsplit.resampling,
                                                       .export = toExport,
                                                       .packages = c("data.table","BuyseTest","lava"),
                                                       .options.snow = opts), {
                                         iOut <- lapply(split.resampling[[iB]], function(iSplit){
                                             if(!is.null(seed)){set.seed(seqSeed[iSplit])}
                                             iBT <- BuyseTest(formula, data = sim(n.T = max.sample.size["T"], n.C = max.sample.size["C"]), trace = 0, ...)
                                             return(iBT)                                             
                                         })
                                         return(iOut)                                         
                                     })
            if(trace>0){close(pb)}
            if(n.rep[1]<=0){parallel::stopCluster(cl)}

            ## collect
            ls.BTmax <- do.call("c",ls2.BTmax)
        }

        if(ls.BTmax[[1]]@method.inference == "u statistic"){
            DeltaMax <- sapply(ls.BTmax, function(iBT){utils::tail(coef(iBT, statistic = names(null)),1) - null})
            IidMax <- do.call(cbind,lapply(ls.BTmax, FUN = getIid, statistic = names(null), scale = FALSE))
                
            ratio <- c(T = as.double(max.sample.size["T"]/sum(max.sample.size)),
                       C = as.double(max.sample.size["C"]/sum(max.sample.size)))
            indexT <- attr(ls.BTmax[[1]]@level.treatment,"indexT")
            indexC <- attr(ls.BTmax[[1]]@level.treatment,"indexC")
                
            sigma2Max <- colMeans(IidMax[indexC,,drop=FALSE]^2)/ratio["C"] + colMeans(IidMax[indexT,,drop=FALSE]^2)/ratio["T"]
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
            sample.sizeC <- ceiling(mean(n.approx*ratio["C"]))
            attr(sample.sizeC,"sample") <- unname(n.approx*ratio["C"])
            sample.sizeT <- ceiling(mean(n.approx*ratio["T"]))
            attr(sample.sizeT,"sample") <- unname(n.approx*ratio["C"])

            ## (mean(IidMax[attr(e.BTmax@level.treatment,"indexC"),]^2) + mean(IidMax[attr(e.BTmax@level.treatment,"indexT"),]^2))*(stats::qnorm(1-alpha/2)+stats::qnorm(power))^2/DeltaMax^2
            if (trace > 1) {
                if(n.rep[2]==1){
                    cat("   - estimated effect (variance): ",unname(DeltaMax)," (",sigma2Max,")\n",sep="")
                    cat("   - estimated sample size      : m=",sample.sizeC,", n=",sample.sizeT,"\n\n",sep="")
                }else{
                    cat("   - average estimated effect (average asymptotic variance)    : ",unname(mean(DeltaMax))," (",mean(sigma2Max),")\n",sep="")
                    cat("   - average estimated sample size [min;max]: m=",
                        sample.sizeC," [",ceiling(min(n.approx*ratio["C"])),";",ceiling(max(n.approx*ratio["C"])),"], n=",
                        sample.sizeT," [",ceiling(min(n.approx*ratio["T"])),";",ceiling(max(n.approx*ratio["T"])),"]\n\n",sep="")
                }
            }
            
        }else{
            stop("Can only determine the sample size when argument \'method.inference\' equals \"u statistic\". \n")
        }
    }
    
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
        ##     if(outArgs$method.inference %in% c("none","u statistic") == FALSE){
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
            resInitData <- do.call(initializeData, args = c(outArgs[argsInit], list(copy = FALSE, data = dt.tempo)))
            do.call(printGeneral, args = c(outArgs, list(M.status = resInitData$M.status, paired = resInitData$paired)))
            if(outArgs$method.inference!="none"){
                do.call(printInference, args = c(outArgs, list(paired = resInitData$paired)))
            }
        }
        if(!missing(sample.size) && !is.null(sample.size)){
            text.sample.size <- paste0("   - sample size: ",paste(sample.size, collapse = " "),"\n")
        }else{
            text.sample.size <- paste0("   - sample size: ",paste(sample.sizeC, collapse = " ")," (control)\n",
                                       "                : ",paste(sample.sizeT, collapse = " ")," (treatment)\n")
        }
        cat("Simulation\n",
            "   - repetitions: ",n.rep[1],"\n",
            "   - cpus       : ",cpus,"\n",
            sep = "")
        cat(" \n")
    }
    
    ## ** define environment
    envirBT <- new.env()
    ## envirBT[[deparse(call)]] <- sim
    name.copy <- c("sim", "option",
                   "outArgs", "sample.sizeTmax", "sample.sizeCmax", "n.sample.size",
                   "sample.sizeC", "sample.sizeT", "n.rep", 
                   "statistic", "null", "conf.level", "alternative", "transformation", "order.Hprojection",
                   ".BuyseTest",".powerBuyseTest")
    for(iObject in name.copy){ ## iObject <- name.copy[2]
        envirBT[[iObject]] <- eval(parse(text = iObject))
    }

    ## ** simulation study
    if (cpus == 1) { ## *** sequential simulation

        ls.simulation <- do.call(method.loop,
                                 args = list(X = 1:n.rep[1],
                                             FUN = function(X){
                                                 if(!is.null(seed)){set.seed(seqSeed[X])}
                                                 iOut <- .powerBuyseTest(i = X,
                                                                         envir = envirBT,
                                                                         statistic = statistic,
                                                                         null = null,
                                                                         conf.level = conf.level,
                                                                         alternative = alternative,
                                                                         transformation = transformation,
                                                                         order.Hprojection = order.Hprojection)
                                                 if(!is.null(seed)){
                                                     return(cbind(iOut, seed = seqSeed[X]))
                                                 }else{
                                                     return(iOut)
                                                 }
                                             })
                                 )
    }else { ## *** parallel simulation

        ## split into a 100 jobs
        split.resampling <- parallel::splitIndices(nx = n.rep[1], ncl = min(max(100,10*cpus), n.rep[1]))
        nsplit.resampling <- length(split.resampling)

        ## define progress bar
        if(trace>0){
            pb <- utils::txtProgressBar(max = nsplit.resampling, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }


        ## try sim
        test <- try(foreach::`%dopar%`(
                                 foreach::foreach(iB=1:cpus,
                                                  .export = toExport,
                                                  .packages = c("data.table","BuyseTest","lava"),
                                                  .options.snow = opts), {
                                                      sim(n.T = sample.sizeTmax, n.C = sample.sizeCmax)
                                                  }),
                    silent = TRUE)
        if(inherits(test,"try-error")){
            stop(paste0("Could not run argument \'sim\' when using multiple CPUs. \n Consider trying first to run powerBuyseTest with cpus=1. \n If it runs, make sure that \'sim\' does not depend on any variable in the global environment or package without explicit mention of the namespace. \n",test))
        }

        ## run simul
        ls2.simulation <- foreach::`%dopar%`(
                                       foreach::foreach(iB=1:nsplit.resampling,
                                                        .export = toExport,
                                                        .packages = c("data.table","BuyseTest","lava"),
                                                        .options.snow = opts), {

                                                            iOut <- lapply(split.resampling[[iB]], function(iSplit){
                                                                if(!is.null(seed)){set.seed(seqSeed[iSplit])}
                                                                iBT <- .powerBuyseTest(i = iSplit,
                                                                                       envir = envirBT,
                                                                                       statistic = statistic,
                                                                                       null = null,
                                                                                       conf.level = conf.level,
                                                                                       alternative = alternative,
                                                                                       transformation = transformation,
                                                                                       order.Hprojection = order.Hprojection)
                                                                if(!is.null(seed)){
                                                                    return(cbind(iBT,seed = seqSeed[iSplit]))
                                                                }else{
                                                                    return(iBT)
                                                                }
                                                            })
                                                        })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}

        ## collect
        ls.simulation <- do.call("c",ls2.simulation)
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
        max.sample.size = max.sample.size,
        power = power,
        n.rep = n.rep,      
        results = dt.out,
        sample.sizeT =  sample.sizeT,
        sample.sizeC =  sample.sizeC,
        seed = seqSeed
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
                   "correction.uninf","method.inference","method.score","paired","strata","grid.strata","threshold","restriction","weightObs","weightEndpoint","pool.strata","n.resampling","call")

    ## ** Simulate data
    data <- data.table::as.data.table(envir$sim(n.T = envir$sample.sizeTmax, n.C = envir$sample.sizeCmax))

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
                            pool.strata = envir$outArgs$pool.strata,
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
    if(envir$outArgs$method.inference %in% c("none","u statistic")){
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

                if(envir$outArgs$method.inference %in% c("none","u statistic")){
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
                                                              pool.strata = envir$outArgs$pool.strata,
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
                    iTransform <- "none"
                    if(!is.null(attr(iCI,"nametransform"))){
                        iTransform <- attr(iCI,"nametransform")
                    }else{
                        iTransform <- "none"
                    }
                    
                    out <- rbind(out,
                                 cbind(data.table::data.table(n.T = envir$sample.sizeT[[iSize]],
                                                              n.C = envir$sample.sizeC[[iSize]],
                                                              endpoint = rownames(iCI),
                                                              statistic = iStatistic,
                                                              transformation = iTransform,
                                                              order.Hprojection = iOrder.Hprojection,
                                                              stringsAsFactors = FALSE),
                                       iCI)
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
