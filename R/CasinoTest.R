### CasinoTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2023 (15:15) 
## Version: 
## Last-Updated: jul 23 2025 (17:09) 
##           By: Brice Ozenne
##     Update #: 184
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * CasinoTest (documentation)
##' @title Multi-group GPC (EXPERIMENTAL)
##' @description Perform Generalized Pairwise Comparisons (GPC) for two or more groups.
##' Can handle one or several binary, continuous and time-to-event endpoints.
##' 
##' @param formula [formula] a symbolic description of the GPC model, see the \code{BuyseTest} function
##' @param data [data.frame] dataset.
##' @param type [character] Type of estimator: can be \code{"unweighted"} or \code{"weighted"}.
##' @param statistic [character] the statistic summarizing the pairwise comparison: \code{"netBenefit"} (Net Treatment Benefit) or \code{"favorable"} (Mann-Withney parameter).
##' Default value read from \code{BuyseTest.options()}. 
##' @param method.inference [character] method used to compute confidence intervals and p-values.
##' Can be \code{"none"}, \code{"u-statistic"}, or \code{"rank"}. 
##' Default value read from \code{BuyseTest.options()}.
##' @param method.multcomp [character] method used to adjust for multiple comparisons.
##' Can be any element of ‘p.adjust.methods’ (e.g. "holm"), "maxT-integration", or "maxT-simulation".
##' @param conf.level [numeric] confidence level for the confidence intervals.
##' Default value read from \code{BuyseTest.options()}.
##' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
##' Default value read from \code{BuyseTest.options()}.
##' @param transformation [logical]  should the CI be computed on the inverse hyperbolic tangent scale / log scale for the net benefit / win ratio and backtransformed.
##' Otherwise they are computed without any transformation.
##' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
##' @param seed [integer, >0] Random number generator (RNG) state used when adjusting for multiple comparisons.
##' If \code{NULL} no state is set.
##' @param ... arguments passed to \code{\link{BuyseTest}}
##' 
##' @details Require to have installed the package riskRegression and BuyseTest
##'
##' Setting argument \code{method.inference} to \code{"rank"} uses a U-statistic approach with a small sample correction to match the variance estimator derived in Result 4.16 page 228 of Brunner (2018).
##' 
##' @return An S3 object of class \code{CasinoTest} that inherits from data.frame.
##' @references Edgar Brunner, Arne C Bathke, and Frank Konietschke (2018). \bold{Rank and pseudo-rank procedures for independent observations in factorial designs}. Springer.
##' @keywords models
##' 

## * CasinoTest (example)
##' @examples
##' library(data.table)
##' library(BuyseTest)
##'
##' #### simulate data ####
##' set.seed(11)
##' n <- 10
##' dt <- rbind(data.table(score = rnorm(n), group = "A"),
##'             data.table(score = rnorm(2*n), group = "B"),
##'             data.table(score = rnorm(3*n), group = "C"))
##' dt$index <- 1:NROW(dt)
##'
##' #### estimation ####
##' score.casino <- dt$score
##'
##' ## naive casino (by hand)
##' M.score <- outer(dt[group=="A",score],score.casino,function(x,y){(x>y)+0.5*(x==y)})
##' mean(M.score) ## 0.4766667
##'
##' ## naive casino (via BuyseTest)
##' CasinoTest(group ~ cont(score), data = dt, type = "weighted", statistic = "favorable")
##' CasinoTest(group ~ cont(score), data = dt, type = "weighted", statistic = "netBenefit")
##' 
##' ## harmonic casino (by hand)
##' hweight <- unlist(tapply(dt$group, dt$group, function(x){rep(1/length(x),length(x))}))
##' M.scoreW <- sweep(M.score, MARGIN = 2, FUN = "*", STATS = NROW(dt)*hweight/3)
##' mean(M.scoreW) ## 0.3680556
##' 
##' ## harmonic casino (via BuyseTest)
##' CasinoTest(group ~ cont(score), data = dt, type = "unweighted", statistic = "favorable")
##'
##' #### Relative liver weights data (Brunner 2018, table 4.1, page 183) ####
##' liverW <- rbind(
##'   data.frame(value = c(3.78, 3.40, 3.29, 3.14, 3.55, 3.76, 3.23, 3.31),
##'              group = "Placebo"),
##'   data.frame(value = c(3.46,3.98,3.09,3.49,3.31,3.73,3.23),
##'              group = "Dose 1"),
##'   data.frame(value = c(3.71, 3.36, 3.38, 3.64, 3.41, 3.29, 3.61, 3.87),
##'              group = "Dose 2"),
##'   data.frame(value = c(3.86,3.80,4.14,3.62,3.95,4.12,4.54),
##'              group = "Dose 3"),
##'   data.frame(value = c(4.14,4.11,3.89,4.21,4.81,3.91,4.19, 5.05),
##'              group = "Dose 4")
##' )
##' liverW$valueU <- liverW$value + (1:NROW(liverW))/1e6
##'
##' ## same as table 4.1, page 183 in Brunner et al (2018)
##' CasinoTest(group ~ cont(value), data = liverW, type = "weighted", statistic = "favorable")
##' CasinoTest(group ~ cont(valueU), data = liverW, type = "unweighted", statistic = "favorable")

## * CasinoTest (code)
##' @export
CasinoTest <- function(formula, data, statistic = NULL, method.inference = NULL, type = "unweighted", 
                       conf.level = NULL, transformation = NULL, alternative = NULL, method.multcomp = "none",
                       seed = NA, ...){

    requireNamespace("riskRegression") ## for confidence bands

    ## ** normalize arguments
    option <- BuyseTest.options()

    ## *** statistic
    if(is.null(statistic)){
        statistic <- option$statistic
    }
    statistic <- match.arg(statistic, c("favorable","netBenefit"))
    
    ## *** method.inference
    if(is.null(method.inference)){
        method.inference <- option$method.inference
    }
    method.inference <- match.arg(gsub("-"," ",tolower(method.inference), fixed = TRUE), c("none","u statistic","rank"))
    if(method.inference=="rank"){
        method.inference <- "u statistic"
        ssc <- TRUE
    }else{
        ssc <- FALSE
    }
    ## *** conf.level
    if(is.null(conf.level)){
        conf.level <- option$conf.level
    }
    ## *** transformation
    if(is.null(transformation)){
        transformation <- option$transformation
    }
    ## *** alternative
    if(is.null(alternative)){
        alternative <- option$alternative
    }
    ## *** type
    if(tolower(type)=="un-weighted"){
        type <- "unweighted"
    }
    type <- match.arg(type, c("weighted","unweighted"))
    if("XXindexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXindexXX\" as this name is used internally by the CasinoTest function. \n")
    }
    if("XXweightXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXweightXX\" as this name is used internally by the CasinoTest function. \n")
    }

    ## *** data
    data <- as.data.frame(data)
    data$XXindexXX <- 1:NROW(data)
    
    ## ** read formula
    details.formula <- initializeFormula(formula, hierarchical = TRUE, envir = environment())
    name.treatment <- details.formula$treatment
    name.endpoint <- details.formula$endpoint
    n.endpoint <- length(name.endpoint)
    if(!is.factor(data[[name.treatment]])){
        data[[name.treatment]] <- as.factor(data[[name.treatment]])
    }else{
        data[[name.treatment]] <- droplevels(data[[name.treatment]])
    }
    level.treatment <- levels(data[[name.treatment]])
    n.treatment <- length(level.treatment)
    elevel.treatement <- paste(level.treatment,collapse=".")
    n.obs <- NROW(data)
    ## prepare normalization
    n.group <- table(data[[name.treatment]])
    norm.groupvar <- 1/n.group
    norm.group <- norm.groupvar[data[[name.treatment]]]
    
    ## ** re-organize data (split by treatment and duplicate with new treatment level)
    ls.data <- by(data,data[[name.treatment]],function(x){x})
    ls.data2 <- lapply(ls.data, function(x){
        x[[name.treatment]] <- elevel.treatement
        return(x)
    })
    data2 <- do.call(rbind,ls.data2)

    ## ** prepare grid and storage
    ## grid
    grid <- .unorderedPairs(level.treatment)
    n.grid <- NCOL(grid)
    grid.BT <- vector(length = n.grid, mode = "list")

    ## storage
    M.estimate <- array(NA, dim = c(n.treatment, n.treatment, n.endpoint),
                        dimnames = list(level.treatment, level.treatment, name.endpoint))
    M.null <- array(NA, dim = c(n.treatment, n.treatment, n.endpoint),
                    dimnames = list(level.treatment, level.treatment, name.endpoint))
    if(method.inference!="none"){
        M.iid <- array(0, dim = c(n.obs,n.endpoint, n.treatment, n.treatment),
                       dimnames = list(NULL, name.endpoint, level.treatment, level.treatment))
    }else{
        M.iid <- NULL
    }

    ## ** Generalized Pairwise Comparisons
    for(iGrid in 1:n.grid){ ## iGrid <- 2
        iTreat1 <- grid[1,iGrid]
        iTreat2 <- grid[2,iGrid]
        iData <- rbind(ls.data[[iTreat1]],ls.data2[[iTreat2]])
        iData[[name.treatment]] <- droplevels(stats::relevel(iData[[name.treatment]], iTreat1))

        ## GPC
        grid.BT[[iGrid]] <- BuyseTest(formula, data = iData, method.inference = method.inference, ..., add.halfNeutral = (statistic=="favorable"), trace = FALSE)
        if(attr(grid.BT[[1]]@scoring.rule,"test.match")){
            stop("CasinoTest does not currently handle matching. \n")
        }
        if(attr(grid.BT[[1]]@weightStrata,"type") %in%  c("standarisation","standarization")){
            stop("CasinoTest does not currently handle standarisation. \n")
        }

        ## store estimate
        iInference <- confint(grid.BT[[iGrid]], statistic = statistic)
        M.estimate[iTreat1,iTreat2,] <- iInference$estimate
        M.null[iTreat1,iTreat2,] <- iInference$null
        if(iTreat1!=iTreat2){
            if(statistic=="favorable"){
                M.estimate[iTreat2,iTreat1,] <- coef(grid.BT[[iGrid]], statistic = "unfavorable")
            }else if(statistic=="netBenefit"){
                M.estimate[iTreat2,iTreat1,] <- -iInference$estimate
            }
            M.null[iTreat2,iTreat1,] <- iInference$null
        }

        ## store iid
        if(method.inference!="none"){
            if(iTreat1==iTreat2){
                iIndex <- unique(sort(iData$XXindexXX))
                M.iid[iIndex,,iTreat1,iTreat1] <- getIid(grid.BT[[iGrid]], statistic = statistic, scale = FALSE, center = TRUE, cluster = iData$XXindexXX)/n.group[iTreat1]
            }else{
                iIndex <- iData$XXindexXX
                M.iid[iIndex,,iTreat1,iTreat2] <- getIid(grid.BT[[iGrid]], statistic = statistic, scale = TRUE, center = TRUE)
                if(statistic=="favorable"){
                    M.iid[iIndex,,iTreat2,iTreat1] <- getIid(grid.BT[[iGrid]], statistic = "unfavorable", scale = TRUE, center = TRUE)
                }else if(statistic=="netBenefit"){
                    M.iid[iIndex,,iTreat2,iTreat1] <- -getIid(grid.BT[[iGrid]], statistic = "netBenefit", scale = TRUE, center = TRUE)
                }
            }
        }
    }

    ## ** collect results averaged over all treatments
    grid.treatmentXendpoint <- expand.grid(level.treatment,name.endpoint)
    out.estimate <- as.data.frame(matrix(NA, nrow = n.treatment*n.endpoint, ncol = 8,
                                         dimnames = list(interaction(grid.treatmentXendpoint, sep = ": "), c("endpoint",name.treatment,"estimate","se","lower.ci","upper.ci","null","p.value"))))
    out.estimate[[name.treatment]] <- grid.treatmentXendpoint[,1]
    out.estimate$endpoint <- grid.treatmentXendpoint[,2]

    if(method.inference!="none"){
        out.iid <- array(NA, dim = c(n.obs, n.treatment, n.endpoint),
                         dimnames = list(NULL, level.treatment, name.endpoint))
    }else{
        out.iid <- NULL
    }

    if(type=="weighted"){
        weight.GPC <- n.group/n.obs
    }else if(type=="unweighted"){
        weight.GPC <- rep(1/n.treatment, n.treatment)
    }
    for(iE in 1:n.endpoint){ ## iE <- 1
        out.estimate[out.estimate$endpoint == name.endpoint[iE],"estimate"] <- colSums(.colMultiply_cpp(M.estimate[,,iE], weight.GPC))
        out.estimate[out.estimate$endpoint == name.endpoint[iE],"null"] <- colSums(.colMultiply_cpp(M.null[,,iE], weight.GPC))

        if(method.inference!="none"){
            for(iT in 1:n.treatment){ ## iT <- 1
                out.iid[,iT,iE] <- rowSums(.rowMultiply_cpp(M.iid[,iE,,iT], weight.GPC))
                if(ssc){
                    out.iid[,iT,iE] <- out.iid[,iT,iE]*sqrt(n.group/(n.group-1))[data[[name.treatment]]]
                }
                out.estimate[out.estimate$endpoint == name.endpoint[iE] & out.estimate[[name.treatment]] == level.treatment[iT],"se"] <- sqrt(sum(out.iid[,iT,iE]^2))
            }
        }
    }

    ## ** statistical inference
    if(method.inference!="none"){
        if(transformation){
            type.transformation <- switch(statistic,
                                          favorable = "atanh2",
                                          netBenefit = "atanh")
        }else{
            type.transformation <- "none"
        }
        e.Band <- riskRegression::transformCIBP(estimate = matrix(out.estimate$estimate, nrow = n.treatment, ncol = n.endpoint, byrow = FALSE,
                                                                  dimnames = list(level.treatment, name.endpoint)),
                                                se = matrix(out.estimate$se, nrow = n.treatment, ncol = n.endpoint, byrow = FALSE,
                                                            dimnames = list(level.treatment, name.endpoint)),
                                                iid = out.iid,
                                                null = matrix(out.estimate$null, nrow = n.treatment, ncol = n.endpoint, byrow = FALSE,
                                                              dimnames = list(level.treatment, name.endpoint)),
                                                conf.level = conf.level,
                                                alternative = alternative,
                                                ci = TRUE, type = type.transformation, min.value = 0, max.value = 1,
                                                p.value = TRUE, band = method.multcomp!="none", 
                                                method.band = method.multcomp,
                                                seed = seed)

        out.estimate$lower.ci <- as.vector(e.Band$lower)
        out.estimate$upper.ci <- as.vector(e.Band$upper)
        out.estimate$p.value <- as.vector(e.Band$p.value)
        if(method.multcomp!="none"){
            out.estimate$lower.band <- as.vector(e.Band$lowerBand)
            out.estimate$upper.band <- as.vector(e.Band$upperBand)
            out.estimate$adj.p.value <- as.vector(e.Band$adj.p.value)
        }
        attr(out.estimate,"iid") <- out.iid
    }

    ## ** export
    class(out.estimate) <- append("CasinoTest",class(out.estimate))
    return(out.estimate)
}

## * .unorderedPairs (from LMMstar package)
.unorderedPairs <- function (x, distinct = FALSE){
    n.x <- length(x)
    out <- do.call(cbind, lapply(1:n.x, function(iK) {
        rbind(x[iK], x[iK:n.x])
    }))
    if (distinct) {
        return(out[, out[1, ] != out[2, ], drop = FALSE])
    }
    else {
        return(out)
    }
}


## alternative implementation 
##     if(method.inference == "u statistic"){
##         if(iTreat1==iTreat2){
##             iIndex <- unique(sort(iData$XXindexXX))
##             M.iid[iIndex,iTreat1,iTreat1] <- getIid(grid.BT[[iGrid]], statistic = "favorable", scale = FALSE, center = FALSE, cluster = iData$XXindexXX)
##         }else{
##             iIndex <- iData$XXindexXX
##             M.iid[iIndex,iTreat1,iTreat2] <- getIid(grid.BT[[iGrid]], statistic = "favorable", scale = FALSE, center = FALSE)
##             M.iid[iIndex,iTreat2,iTreat1] <- getIid(grid.BT[[iGrid]], statistic = "unfavorable", scale = FALSE, center = FALSE)
##         }
##     }
## }

## ## ** collect results averaged over all treatments
## out.estimate <- as.data.frame(matrix(NA, nrow = n.treatment, ncol = 6,
##                                      dimnames = list(level.treatment, c("estimate","se","lower.ci","upper.ci","null","p.value"))))
## out.iid <- matrix(NA, nrow = n.obs, ncol = n.treatment,
##                   dimnames = list(NULL, level.treatment))

## if(type=="weighted"){
##     weight.GPC <- n.group/n.obs
## }else if(type=="unweighted"){
##     weight.GPC <- rep(1/n.treatment, n.treatment)
## }
## out.estimate$estimate <- colSums(.colMultiply_cpp(M.estimate, weight.GPC))
## out.estimate$null <- colSums(.colMultiply_cpp(M.null, weight.GPC))

## for(iT in 1:n.treatment){ ## iT <- 1

##     iExpectation <- rowSums(.rowMultiply_cpp(M.iid[,,iT], weight.GPC))
##     iCenter <- tapply(iExpectation, data[[name.treatment]], mean)
##     if(ssc){
##         ## one n.group is used to evaluate the variance of the H-decomposition, i.e. becomes n.group-1
##         ## one n.group is from the H-decomposition (1/n \sum H)
##         out.iid[,iT] <- (iExpectation-iCenter[data[[name.treatment]]])/sqrt(n.group[data[[name.treatment]]]*(n.group[data[[name.treatment]]]-1))
##     }else{
##         out.iid[,iT] <- (iExpectation-iCenter[data[[name.treatment]]])/n.group[data[[name.treatment]]]
##     }

## }

##----------------------------------------------------------------------
### CasinoTest.R ends here
