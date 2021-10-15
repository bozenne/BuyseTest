### multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  4 2021 (16:17) 
## Version: 
## Last-Updated: okt 15 2021 (11:40) 
##           By: Brice Ozenne
##     Update #: 143
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * BuyseMultComp (documentation)
##' @title Adjustment for Multiple Comparisons
##' @description Adjustment p-values and confidence estimated via GPC for multiple comparisons.
##' @name BuyseMultComp
##'
##' @param object A BuyseTest object or a list of BuyseTest objects. All objects should contain the same endpoints.
##' @param cluster [character] name of the variable identifying the observations in the dataset used by each BuyseTest model.
##' Only relevant when using a list of BuyseTest objects to correctly combine the influence functions.
##' If NULL, then it is assumed that the BuyseTest objects correspond to different groups of individuals.
##' @param linfct [numeric matrix] a contrast matrix of size the number of endpoints times the number of BuyseTest models.
##' @param rhs [numeric vector] the values for which the test statistic should be tested against. Should have the same number of rows as \code{linfct}.
##' @param endpoint [character or numeric vector] the endpoint(s) to be considered.
##' @param statistic [character] the statistic summarizing the pairwise comparison:
##' \code{"netBenefit"} displays the net benefit, as described in Buyse (2010) and Peron et al. (2016)),
##' \code{"winRatio"} displays the win ratio, as described in Wang et al. (2016),
##' \code{"favorable"} displays the proportion in favor of the treatment (also called Mann-Whitney parameter), as described in Fay et al. (2018).
##' \code{"unfavorable"} displays the proportion in favor of the control.
##' Default value read from \code{BuyseTest.options()}.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param ... argument passsed to the function \code{transformCIBP} of the riskRegression package.
#'
#' @details Simulateneous confidence intervals and adjusted p-values are computed using a single-step max-test approach via the function \code{transformCIBP} of the riskRegression package.
#'
#' @examples
#' ## simulate data
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 3)
#' df.data$id <- paste0("Id",1:NROW(df.data))
#'
#' ## adjustment over endpoints
#' ff1 <- treatment ~ TTE(eventtime, status = status, threshold = 0.1)
#' ff2 <- update(ff1, .~. + cont(score, threshold = 1))
#' BT2 <- BuyseTest(ff2, data= df.data, trace = FALSE)
#'
#' ## (require riskRegression >= 2021.10.04 to match)
#' confint(BT2) ## not adjusted
#' BuyseMultComp(BT2, endpoint = 1:2) ## adjusted
#' 
#' ## adjustment strata
#' BT1a <- BuyseTest(treatment ~ TTE(eventtime, status = status, threshold = 0.1),
#'                  data= df.data[strata == "a",], trace = FALSE)
#' BT1b <- BuyseTest(treatment ~ TTE(eventtime, status = status, threshold = 0.1),
#'                  data= df.data[strata == "b",], trace = FALSE)
#' BT1c <- BuyseTest(treatment ~ TTE(eventtime, status = status, threshold = 0.1),
#'                  data= df.data[strata == "c",], trace = FALSE)
#' rbind(a = confint(BT1a), b = confint(BT1b), c = confint(BT1c)) ## not adjusted
#' BuyseMultComp(list(a = BT1a, b = BT1b, c = BT1c))
#' BuyseMultComp(list(a = BT1a, b = BT1b, c = BT1c), cluster = "id")

## * BuyseMultComp (code)
##' @rdname BuyseMultComp
##' @export
BuyseMultComp <- function(object, cluster = NULL, linfct = NULL, rhs = NULL, endpoint = NULL, statistic = NULL, 
                          conf.level = NULL, alternative = NULL, transformation = NULL, ...){

    ## ** normalize arguments
    option <- BuyseTest.options()

    ## object
    if(inherits(object,"S4BuyseTest")){
        test.list <- FALSE
    }else if(all(sapply(object,inherits,"S4BuyseTest"))){
        n.object <- length(object)
        test.list <- TRUE
        if(is.null(object)){
            names(object) <- paste0("test",1:n.object)
        }
        name.object <- names(object)
    }else{
        stop("Incorrect \'object\': should be a BuyseTest object or a list of BuyseTest objects. \n")
    }

    ## statistic
    if(is.null(statistic)){
        statistic <- tolower(option$statistic)
    }else{
        if(length(statistic)>1){
            stop("Argument \'statistic\' must have length 1. \n")
        }
        statistic <- match.arg(statistic, c("netbenefit","winratio","favorable","unfavorable"), several.ok = FALSE)
    }

    ## endpoint
    if(test.list){
        ls.valid.endpoint <- unique(lapply(object,function(iBT){paste0(iBT@endpoint,"_",iBT@threshold)}))
        if(length(ls.valid.endpoint)==1){
            valid.endpoint <- ls.valid.endpoint[[1]]
        }else{
            stop("Cannot handle BuyseTest objects with different endpoints. \n")
        }
    }else{
        valid.endpoint <- paste0(object@endpoint,"_",object@threshold)
    }
    if(is.null(endpoint)){
        endpoint <- valid.endpoint[length(valid.endpoint)]
    }else if(is.numeric(endpoint)){
        validInteger(endpoint,
                     name1 = "endpoint",
                     min = 1, max = length(valid.endpoint),
                     valid.length = NULL,
                     method = "iid[BuyseTest]")
        endpoint <- valid.endpoint[endpoint]
    }else{
        validCharacter(endpoint,
                       valid.length = 1:length(valid.endpoint),
                       valid.values = valid.endpoint,
                       refuse.NULL = FALSE)
    }
    n.endpoint <- length(endpoint)

    ## conf.level
    if(is.null(conf.level)){
        conf.level <- option$conf.level
    }
    
    ## altenative
    if(is.null(alternative)){
        alternative <- option$alternative
    }

    ## transformation
    if(is.null(transformation)){
        transformation <- option$transformation
    }
    ## ** extract iid and coefficients
    if(test.list){
        if(n.endpoint>1){
            ls.iName <- lapply(name.object,paste0,paste0(": ",endpoint))
            iName <- unlist(ls.iName)
        }else{
            ls.iName <- as.list(name.object)
            iName <- name.object
        }
        vec.beta <- stats::setNames(unlist(lapply(object, coef, endpoint = endpoint, statistic = statistic)), iName)
        vec.se <- stats::setNames(unlist(lapply(object, function(iO){confint(iO, endpoint = endpoint, statistic = statistic)$se})), iName)
        ls.iid <- lapply(1:length(object), function(iO){ ## iO <- 1
            iIID <- do.call(cbind,getIid(object[[iO]], endpoint = endpoint, statistic = statistic))
            colnames(iIID) <- ls.iName[[iO]]
            return(iIID)
        })
        if(is.null(cluster)){
            seqn.id <- sapply(ls.iid,NROW)
            cumseqn.id <- cumsum(seqn.id)
            n.id <- cumseqn.id[n.object]
            M.iid <- matrix(0, nrow = n.id, ncol = n.object*n.endpoint, dimnames = list(NULL, iName))
            for(iObject in 1:length(name.object)){
                iStart <- c(1,cumseqn.id+1)[iObject]
                iStop <- cumseqn.id[iObject]
                M.iid[iStart:iStop,ls.iName[[iObject]]] <- ls.iid[[iObject]]
            }
        }else{
            ## try to retrieve data
            cluster.var <- cluster
            cluster <- try(lapply(object, function(iO){as.character(eval(iO@call$data)[[cluster.var]]) }), silent = TRUE)
            if(inherits(cluster, "try-error")){
                stop("Could not retrieve the column \"cluster\" from the evaluation of the call of the BuyseTest objects. \n")
            }
            ## find unique clusters
            Ucluster <- unique(unlist(cluster))
            n.id <- length(Ucluster)
            ## store iid according to the clusters
            M.iid <- matrix(0, nrow = n.id, ncol = n.object*n.endpoint, dimnames = list(NULL, iName))
            for(iObject in 1:length(name.object)){
                M.iid[match(cluster[[iObject]],Ucluster),ls.iName[[iObject]]] <- ls.iid[[iObject]]
            }
        }
    }else{
        iName <- endpoint

        vec.beta <- stats::setNames(coef(object, endpoint = endpoint, statistic = statistic), iName)
        vec.se <- stats::setNames(confint(object, endpoint = endpoint, statistic = statistic)$se, iName)
        M.iid <- do.call(cbind,getIid(object, endpoint = endpoint, statistic = statistic))
        colnames(M.iid) <- iName
    }

    n.beta <- length(vec.beta)
    vec.name <- names(vec.beta)

    ## ** create linfct matrix
    if(is.null(linfct)){
        linfct <- diag(1, nrow = n.beta, ncol = n.beta)
        dimnames(linfct) = list(vec.name,vec.name)
    }else{
        if(NCOL(linfct) != n.beta){
            stop("Incorrect argument \'linfct\': should have ",n.beta," columns. \n")
        }
        if(is.null(colnames(linfct))){
            colnames(linfct) <- names(vec.beta)
        }else{
            if(any(vec.name %in% colnames(linfct) == FALSE)){
                stop("Missing column \"",paste0(vec.name[vec.name%in%colnames(linfct)==FALSE],collapse="\" \""),"\" in argument \'linfct\'. \n")
            }
            linfct <- linfct[,vec.name,drop=FALSE]
        }
    }

    ## ** create rhs vector
    if(is.null(rhs)){
        rhs <- rep(switch(statistic,
                          "netbenefit" = 0,
                          "winratio" = 1,
                          "favorable" = 0.5,
                          "unfavorable" = 0.5), n.beta)
    }else{
        if(length(rhs) != NROW(linfct)){
            stop("Length of argument \'rhs\' does not match the number of row of argument \'contrast\' (",length(rhs)," vs. ",NROW(linfct),"). \n")
        }
    }
    ## ** perform adjustment for multiple comparisons
    requireNamespace("riskRegression")
    A.iid <- array(NA, dim = c(NROW(M.iid), NCOL(M.iid),1))
    A.iid[,,1] <- M.iid

    min.value <- switch(statistic,
                        "netBenefit" = -1,
                        "winRatio" = 0,
                        "favorable" = 0,
                        "unfavorable" = 0)
    max.value <- switch(statistic,
                        "netBenefit" = 1,
                        "winRatio" = Inf,
                        "favorable" = 1,
                        "unfavorable" = 1)

    if(statistic %in% c("none","netBenefit","winRatio") || inherits(try(riskRegression::transformCIBP(estimate = 0.5, se = cbind(0.1), type = "atanh2", seed = NA, band = FALSE, alternative = "two.sided", p.value = TRUE, ci = TRUE, conf.level = 0.95, min.value = -Inf, max.value = Inf, null = 0.5),silent=TRUE),"try-error")){
        type <- switch(statistic,
                       "netbenefit" = "atanh",
                       "winratio" = "log",
                       "favorable" = "cloglog",## note: not the same transformation as confint
                       "unfavorable" = "cloglog", ## note: not the same transformation as confint
                       "none" = "none") 
    }else{
        type <- switch(statistic,
                       "netbenefit" = "atanh",
                       "winratio" = "log",
                       "favorable" = "atanh2",
                       "unfavorable" = "atanh2",
                       "none" = "none") 
    }

    dots <- list(...)
    if("seed" %in% names(dots) == FALSE){
        dots$seed <- NA
    }
    if("method.band" %in% names(dots) == FALSE){
        dots$method.band <- "maxT-integration"
    }

    iBand <- do.call(riskRegression::transformCIBP,
                     args = c(list(estimate = rbind(vec.beta),
                                   se = rbind(vec.se),
                                   iid = A.iid,
                                   null = rhs,
                                   conf.level = conf.level,
                                   alternative = alternative,
                                   ci = TRUE, type = type, min.value = min.value, max.value = max.value,
                                   band = TRUE, p.value = TRUE),
                              dots))

    ## ** export
    out <- data.frame(estimate = as.double(vec.beta), se = as.double(vec.se),
                      lower.ci = as.double(iBand$lower), upper.ci = as.double(iBand$upper), null = rhs, p.value = as.double(iBand$p.value),
                      lower.band = as.double(iBand$lowerBand), upper.band = as.double(iBand$upperBand), adj.p.value = as.double(iBand$adj.p.value))
    if(test.list && n.endpoint > 1){
        out <- cbind(model = sapply(name.object, rep, length(endpoint)), endpoint = rep(endpoint, length(object)), out)
        rownames(out) <- NULL
    }else{
        rownames(out) <- iName
    }
    attr(out,"iid") <- A.iid[,,1]
    attr(out,"linfct") <- linfct
    attr(out,"quantileBand") <- iBand$quantile
    attr(out,"quantileBand") <- iBand$quantile
    return(out)
}
        

##----------------------------------------------------------------------
### multcomp.R ends here
