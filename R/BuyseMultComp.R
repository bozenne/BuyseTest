### multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  4 2021 (16:17) 
## Version: 
## Last-Updated: jul 17 2023 (17:29) 
##           By: Brice Ozenne
##     Update #: 294
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
##' @description Adjust p-values and confidence intervals estimated via GPC for multiple comparisons.
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
#' @param cumulative [logical] should the summary statistic be cumulated over endpoints?
#' Otherwise display the contribution of each endpoint.
#' @param conf.level [numeric] confidence level for the confidence intervals.
#' Default value read from \code{BuyseTest.options()}.
#' @param band [logical] Should confidence intervals and p-values adjusted for multiple comparisons be computed.
#' @param global [logical] Should global test (intersection of all null hypotheses) be made? 
#' @param alternative [character] the type of alternative hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param transformation [logical]  should the CI be computed on the logit scale / log scale for the net benefit / win ratio and backtransformed.
#' Otherwise they are computed without any transformation.
#' Default value read from \code{BuyseTest.options()}. Not relevant when using permutations or percentile bootstrap.
#' @param ... argument passsed to the function \code{transformCIBP} of the riskRegression package.
#'
#' @details Simulateneous confidence intervals and adjusted p-values are computed using a single-step max-test approach via the function \code{transformCIBP} of the riskRegression package.
#' This corresponds to the single-step Dunnett described in Dmitrienko et al (2013) in table 2 and section 7.
#'
#' @return An S3 object of class \code{BuyseMultComp}.
#' 
#' @keywords htest
#' 
#' @references Dmitrienko, A. and D'Agostino, R., Sr (2013), Traditional multiplicity adjustment methods in clinical trials. Statist. Med., 32: 5172-5218. https://doi.org/10.1002/sim.5990
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 3)
#'
#' #### adjustment for all univariate analyses ####
#' ff1 <- treatment ~ TTE(eventtime, status = status, threshold = 0.1)
#' ff2 <- update(ff1, .~. + cont(score, threshold = 1))
#' BT2 <- BuyseTest(ff2, data= df.data, trace = FALSE)
#'
#' ## (require riskRegression >= 2021.10.04 to match)
#' confint(BT2, cumulative = FALSE) ## not adjusted
#' confintAdj <- BuyseMultComp(BT2, cumulative = FALSE, endpoint = 1:2) ## adjusted
#' confintAdj
#' if(require(lava)){
#' cor(lava::iid(confintAdj)) ## correlation between test-statistic
#' }
#' 
#' #### 2- adjustment for multi-arm trial ####
#' ## case where we have more than two treatment groups
#' ## here strata will represent the treatment groups
#' df.data$strata <- as.character(df.data$strata)
#' df.data$id <- paste0("Id",1:NROW(df.data)) ## define id variable
#' 
#' BT1ba <- BuyseTest(strata ~ TTE(eventtime, status = status, threshold = 1),
#'                    data= df.data[strata %in% c("a","b"),], trace = FALSE)
#' BT1ca <- BuyseTest(strata ~ TTE(eventtime, status = status, threshold = 0.1),
#'                    data= df.data[strata %in% c("a","c"),], trace = FALSE)
#' BT1cb <- BuyseTest(strata ~ TTE(eventtime, status = status, threshold = 0.1),
#'                    data= df.data[strata %in% c("b","c"),], trace = FALSE)
#' rbind("b-a" = confint(BT1ba),
#'       "c-a" = confint(BT1ca),
#'       "c-b" = confint(BT1cb)) ## not adjusted
#' confintAdj <- BuyseMultComp(list("b-a" = BT1ba, "c-a" = BT1ca, "c-b" = BT1cb),
#'                             cluster = "id", global = TRUE)
#' confintAdj
#' if(require(lava)){
#' cor(lava::iid(confintAdj))
#' }

## * BuyseMultComp (code)
##' @rdname BuyseMultComp
##' @export
BuyseMultComp <- function(object, cluster = NULL, linfct = NULL, rhs = NULL, endpoint = NULL, statistic = NULL, cumulative = TRUE,
                          conf.level = NULL, band = TRUE, global = FALSE, alternative = NULL, transformation = NULL, ...){

    ## ** normalize arguments
    option <- BuyseTest.options()
    call <- match.call()
    ## object
    if(inherits(object,"S4BuyseTest")){
        test.list <- FALSE
        if(any(object@weightObs!=1)){
            stop("Cannot not currently handle weighted observations. \n") 
        }
        name.object <- NULL
    }else if(all(sapply(object,inherits,"S4BuyseTest"))){
        n.object <- length(object)
        test.list <- TRUE
        if(is.null(object)){
            names(object) <- paste0("test",1:n.object)
        }
        if(any(sapply(object, function(iO){any(iO@weightObs!=1)}))){
            stop("Cannot not currently handle weighted observations. \n") 
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
        valid.endpoint <- lapply(object,function(iBT){names(iBT@endpoint)})
        if(is.null(endpoint)){
            endpoint <- unlist(lapply(valid.endpoint, function(iE){iE[length(iE)]}))
        }else if(is.vector(endpoint)){
            if(length(endpoint)==1){
                endpoint <- rep(endpoint, n.object)
            }else if(length(endpoint)!=n.object){
                stop("Argument \'endpoint\' misspecified. \n",
                     "Should have length 1 or length ",n.object," (i.e. the number of objects in the list). \n")
            }
            if(is.numeric(endpoint)){
                for(iE in 1:length(endpoint)){
                    if(endpoint[iE]<=0 || endpoint[iE]>=length(valid.endpoint[[iE]])){
                        stop("The ",iE," element in argument \'endpoint\' should be between 1 and ",length(valid.endpoint[[iE]]),". \n")
                    }
                }
            }else{
                for(iE in 1:length(endpoint)){
                    if(endpoint[iE] %in% valid.endpoint[[iE]] == FALSE){
                        stop("The ",iE," element in argument \'endpoint\' should one \"",paste(valid.endpoint[[iE]],collapse = "\" \""),"\". \n")
                    }
                }
            }
        }else{
            stop("Argument \'endpoint\' should be a vector. \n")
        }
    }else{
        valid.endpoint <- names(object@endpoint)
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
    }
    

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

    ## type of transformation
    if(transformation){
        type <- switch(statistic,
                       "netbenefit" = "atanh",
                       "winratio" = "log",
                       "favorable" = "atanh2",
                       "unfavorable" = "atanh2",
                       "none" = "none") 
    }else{
        type <- "none"
    }

    ## weights

    ## ** extract iid and coefficients
    if(test.list){
        if(is.null(name.object)){
            iName <- endpoint
        }else{
            if(any(duplicated(name.object))){
                iName <- paste0(name.object,": ",endpoint)
            }else{
                iName <- name.object
            }
        }
        if(any(duplicated(iName))){
            stop("Duplicated names for the estimates: provide unique name to each element of the list containing the S4BuyseTest objects. \n")
        }
        ls.beta <- lapply(1:n.object, function(iO){coef(object[[iO]], endpoint = endpoint[iO], statistic = statistic, cumulative = cumulative, strata = "global")})
        vec.beta <- stats::setNames(unlist(ls.beta), iName)

        ls.iid <- lapply(1:n.object, function(iO){ ## iO <- 1
            iIID <- getIid(object[[iO]], endpoint = endpoint[iO], statistic = statistic, cumulative = cumulative, strata = "global", simplify = FALSE)[["global"]]
            colnames(iIID) <- iName[iO]
            return(iIID)
        })
        if(is.null(cluster)){
            ## seqn.id <- sapply(ls.iid,NROW)
            ## cumseqn.id <- cumsum(seqn.id)
            ## n.id <- cumseqn.id[n.object]
            ## M.iid <- matrix(0, nrow = n.id, ncol = n.object*n.endpoint, dimnames = list(NULL, iName))
            ## for(iObject in 1:length(name.object)){
            ##     iStart <- c(1,cumseqn.id+1)[iObject]
            ##     iStop <- cumseqn.id[iObject]
            ##     M.iid[iStart:iStop,iName[iObject]] <- ls.iid[[iObject]]
            ## }
            stop("The argument \'cluster\' must be specified to identify the common individuals across the BuyseTest object. \n")
        }else{
            ## retrieve data
            if(is.numeric(cluster)){
                if(length(unique(sapply(ls.iid,length)))>1){
                    stop("Argument \'cluster\' cannot be numeric when the BuyseTest are performed on different number of observations. \n",
                         "Number of observations per BuyseTest: ",paste(sapply(ls.iid,length), collapse = ", "),".\n")
                }
                if(any(length(cluster) != length(ls.iid)[[1]])){
                    stop("Incorrect length for argument \'cluster\': when numeric it should have length the number of observations ",length(ls.iid[[1]]),".\n",
                         "Current length: ",length(cluster),".\n")
                }
                if(any(cluster %in% 1:length(ls.iid[[1]]) == FALSE)){
                    stop("Incorrect values for argument \'cluster\': when numeric it should takes integer values between 1 and ",length(ls.iid[[1]]),".\n",
                         "Example of incorrect value: ",cluster[cluster %in% 1:length(ls.iid[[1]]) == FALSE][1],".\n")
                }
                cluster <- lapply(1:n.object, function(iO){cluster})
            }else if(is.list(cluster)){
                if(any(sapply(cluster,is.numeric) == FALSE)){
                    stop("When a list, argument \'cluster\' should contain integers indexing the clusters")
                }
                if(any(sapply(cluster,length) != sapply(ls.iid,length))){
                    stop("Incorrect argument \'cluster\': when a list, the length of each element should match the number of observations the BuyseTest.\n",
                         "Length argument vs. BuyseTest: ",paste(paste(sapply(cluster,length), sapply(ls.iid,length), sep = " vs. "), collapse = ", "),".\n")
                }
                if(any(sapply(cluster,max) < 1) || any(unlist(lapply(cluster, `%%`, 1))!=0)){
                    stop("Incorrect argument \'cluster\': when a list, element should contain integers indexing the clusters.\n")
                }
            }else{
                cluster.var <- cluster
                cluster <- lapply(object, function(iO){try(as.character(iO@call$data[[cluster.var]]), silent = TRUE) })
                if(any(sapply(cluster, inherits, "try-error"))){
                    indexPb <- which(sapply(cluster, inherits, "try-error"))
                    stop("Could not retrieve the cluster column \"",cluster.var,"\" from the evaluation of the call of the BuyseTest objects. \n",
                         "Problematic object(s): ",paste(iName[indexPb], collapse = ","),"\n")
                }
                if(any(sapply(cluster, length)==0)){
                    indexPb <- which(sapply(cluster, length)==0)
                    stop("Could not retrieve the cluster column \"",cluster.var,"\" from the evaluation of the call of the BuyseTest objects. \n",
                         "Problematic object(s): ",paste(iName[indexPb], collapse = ","),"\n")
                }
            }
            ## find unique clusters
            Ucluster <- unique(unlist(cluster))
            cluster.factor <- lapply(cluster, factor, level = Ucluster)
            n.id <- length(Ucluster)
            ## store iid according to the clusters
            M.iid <- matrix(0, nrow = n.id, ncol = n.object, dimnames = list(Ucluster, iName))
            for(iObject in 1:n.object){ ## iObject <- 1
                M.iid[,iName[iObject]] <- tapply(ls.iid[[iObject]],cluster.factor[[iObject]],sum, default = 0)
            }
        }
    }else{
        iName <- endpoint
        vec.beta <- stats::setNames(coef(object, endpoint = endpoint, statistic = statistic, cumulative = cumulative, strata = "global"), iName)
        M.iid <- getIid(object, endpoint = endpoint, statistic = statistic, cumulative = cumulative, strata = "global")
        colnames(M.iid) <- iName
    }

    n.beta <- length(vec.beta)
    vec.name <- names(vec.beta)

    ## ** create and apply linfct matrix
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
                stop("Missing column \"",paste0(vec.name[vec.name %in% colnames(linfct)==FALSE],collapse="\" \""),"\" in argument \'linfct\'. \n")
            }
            linfct <- linfct[,vec.name,drop=FALSE]
        }
    }
    vec.Cbeta <- t(linfct %*% vec.beta)
    M.Ciid <- M.iid  %*% t(linfct)
    n.C <- NROW(linfct)

    vec.Cse <- rbind(sqrt(diag(crossprod(M.Ciid))))
    A.Ciid <- array(NA, dim = c(NROW(M.Ciid), NCOL(M.Ciid),1))
    A.Ciid[,,1] <- M.Ciid

    ## ** create rhs vector
    if(is.null(rhs)){
        rhs <- rep(switch(statistic,
                          "netbenefit" = 0,
                          "winratio" = 1,
                          "favorable" = 0.5,
                          "unfavorable" = 0.5), n.C)
    }else{
        if(length(rhs) != n.C){
            stop("Length of argument \'rhs\' does not match the number of row of argument \'contrast\' (",length(rhs)," vs. ",NROW(linfct),"). \n")
        }
    }
    

    ## ** perform global test (single multivariate Wald test)
    if(global){
        if(transformation){
            vec.CbetaTrans <- riskRegression_transformT(estimate = vec.Cbeta, se = 1 , null = rhs, type = type, alternative = alternative)
            M.cSigmaTrans <- crossprod(riskRegression_transformIID(estimate = vec.Cbeta, iid = A.Ciid, type = type)[,,1])
            dimnames(M.cSigmaTrans) <- list(colnames(M.iid),colnames(M.iid))
        }else{
            vec.CbetaTrans <- vec.Cbeta
            M.CsigmaTrans <- crossprod(M.Ciid)
        }
        iStat <- try(as.double(vec.CbetaTrans %*% solve(M.cSigmaTrans) %*% t(vec.CbetaTrans), silent = TRUE))
        if(inherits(iStat,"try-error")){
            out.chisq <- iStat
        }else{
            out.chisq <- data.frame(statistic = iStat, df = n.C, p.value = 1 - stats::pchisq(iStat, df = n.C))
        }
    }else{
        out.chisq <- data.frame(statistic = NA, df = NA, p.value = NA)
    }

    ## ** perform adjustment for multiple comparisons (several univariate Wald test)
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

    
    dots <- list(...)
    if("seed" %in% names(dots) == FALSE){
        dots$seed <- NA
    }
    if("method.band" %in% names(dots) == FALSE){
        dots$method.band <- "maxT-integration"
    }
    requireNamespace("riskRegression")

    iBand <- do.call(riskRegression::transformCIBP,
                     args = c(list(estimate = vec.Cbeta,
                                   se = vec.Cse,
                                   iid = A.Ciid,
                                   null = rhs,
                                   conf.level = conf.level,
                                   alternative = alternative,
                                   ci = TRUE, type = type, min.value = min.value, max.value = max.value,
                                   band = band, p.value = TRUE),
                              dots))

    ## ** export
    out <- list(table.uni = NULL,
                table.multi = out.chisq,
                iid = M.Ciid,
                linfct = linfct,
                quantileBand = iBand$quantile
                )
    if(band){
        out$table.uni <- data.frame(estimate = as.double(vec.Cbeta), se = as.double(vec.Cse),
                                    lower.ci = as.double(iBand$lower), upper.ci = as.double(iBand$upper), null = rhs, p.value = as.double(iBand$p.value),
                                    lower.band = as.double(iBand$lowerBand), upper.band = as.double(iBand$upperBand), adj.p.value = as.double(iBand$adj.p.value))
    }else{
        out$table.uni <- data.frame(estimate = as.double(vec.Cbeta), se = as.double(vec.Cse),
                                    lower.ci = as.double(iBand$lower), upper.ci = as.double(iBand$upper), null = rhs, p.value = as.double(iBand$p.value),
                                    lower.band = NA, upper.band = NA, adj.p.value = NA)
    }
    rownames(out$table.uni) <- iName
    class(out) <- append("BuyseMultComp",class(out))
    return(out)
}
        
## * as.data.frame.BuyseMultComp
##' @method as.data.frame BuyseMultComp
##' @export
as.data.frame.BuyseMultComp <- function(x, row.names = NULL, optional = FALSE, ...){
    return(as.data.frame(x$table.uni, row.names = row.names, optional = optional, ...))
}

## * as.data.table.BuyseMultComp
##' @method as.data.table BuyseMultComp
##' @export
as.data.table.BuyseMultComp <- function(x, keep.rownames = NULL, ...){
    return(as.data.table(x$table.uni, keep.rownames = keep.rownames, ...))
}

## * coef.BuyseMultComp
##' @method coef BuyseMultComp
##' @export
coef.BuyseMultComp <- function(object, ...){
    out <- stats::setNames(object$table.uni$estimate,rownames(object$table.uni))
    return(out) 
}

## * confint.BuyseMultComp
##' @method confint BuyseMultComp
##' @export
confint.BuyseMultComp <- function(object, parm, level = 0.95, ...){
    out <- object$table.uni[,c("estimate","lower.band","upper.band","adj.p.value")]
    return(out) 
}

## * iid.BuyseMultComp
##' @method iid BuyseMultComp
##' @export
iid.BuyseMultComp <- function(x, keep.rownames = NULL, ...){
    return(x$iid)
}

## * model.tables.BuyseMultComp
##' @method model.tables BuyseMultComp
##' @export
model.tables.BuyseMultComp <- function(x, type = "univariate", ...){
    type <- match.arg(type, c("univariate","multivariate"))
    if(type == "univariate"){
        return(x$table.uni)
    }else if(type == "multivariate"){
        return(x$table.multi)
    }
}

## * print.BuyseMultComp
##' @method print BuyseMultComp
##' @export
print.BuyseMultComp <- function(x, ...){
    dots <- list(...)
    
    if(!all(is.na(x$table.multi))){
        cat("  - Multivariate test: p.value = ",x$table.multi[,"p.value"]," (df = ",x$table.multi[,"df"],")\n",sep="")
    }
    cat("  - Univariate tests:\n",sep="")
    if(any(dots$cols %in% names(x$table.uni) == FALSE) ){
        stop("Incorrect argument \'cols\'. \n",
             "Valid values: \"",paste(names(x$table.uni), collapse = "\" \""),"\".\n")
    }
    if(!is.null(dots$cols)){
        print(x$table.uni[,dots$cols,drop=FALSE])
    }else{
        print(x$table.uni)
    }
    
    return(invisible(NULL))
}



##----------------------------------------------------------------------
### multcomp.R ends here
