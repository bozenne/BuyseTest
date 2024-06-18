## * inferenceResampling
inferenceResampling <- function(envir){

    cpus <- envir$outArgs$cpus
    D <- envir$outArgs$D
    endpoint <- envir$outArgs$endpoint
    iid <- envir$outArgs$iid
    level.strata <- envir$outArgs$level.strata
    method.inference <- envir$outArgs$method.inference

    n.resampling <- envir$outArgs$n.resampling
    n.strata <- envir$outArgs$n.strata
    seed <- envir$outArgs$seed

    if (!is.null(seed)) {
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        tol.seed <- attr(seed,"max")
        set.seed(seed)
        seqSeed <- sample.int(tol.seed, n.resampling,  replace = FALSE)        
    }

    trace <- envir$outArgs$trace

    ## re-order dataset according to the strata used when resampling
    if(!is.na(attr(method.inference,"resampling-strata"))){
        envir$outArgs$data[,c("..rowIndex..") := 1:.N]
        data.table::setkeyv(envir$outArgs$data, cols = attr(method.inference,"resampling-strata"))

        envir$outArgs$M.endpoint <- envir$outArgs$M.endpoint[envir$outArgs$data[["..rowIndex.."]],,drop=FALSE]
        envir$outArgs$M.status <- envir$outArgs$M.status[envir$outArgs$data[["..rowIndex.."]],,drop=FALSE]
        envir$outArgs$index.C <- which(envir$outArgs$data[[envir$outArgs$treatment]] == 0)
        envir$outArgs$index.T <- which(envir$outArgs$data[[envir$outArgs$treatment]] == 1)
        envir$outArgs$index.strata <- tapply(1:NROW(envir$outArgs$data), envir$outArgs$data[["..strata.."]], list)
        envir$outArgs$data[,c("..rowIndex..") := NULL,]
    }
    
    ## ** computation
    if (cpus == 1) { ## *** sequential resampling test

        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.resampling <- do.call(method.loop,
                                 args = list(X = 1:n.resampling,
                                             FUN = function(iB){
                                                 if(!is.null(seed)){set.seed(seqSeed[iB])}
                                                 iOut <- .BuyseTest(envir = envir,
                                                                    iid = iid,
                                                                    method.inference = method.inference,
                                                                    pointEstimation = FALSE
                                                                    )
                                                 if(!is.null(seed)){
                                                     return(c(iOut,list(seed = seqSeed[iB])))
                                                 }else{
                                                     return(iOut)
                                                 }
                                             })
                                 )
    }else { ## *** parallel resampling test

        ## split into a 100 jobs
        split.resampling <- parallel::splitIndices(nx = n.resampling, ncl = min(max(100,10*cpus), n.resampling))
        nsplit.resampling <- length(split.resampling)

        ## define cluster
        cl <- parallel::makeCluster(cpus)
        if(trace>0){
            pb <- utils::txtProgressBar(max = nsplit.resampling, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }
        ## link to foreach
        doSNOW::registerDoSNOW(cl)

        ## seed
        if (!is.null(seed)) {
            parallel::clusterExport(cl, varlist = "seqSeed", envir = environment())
        }         

        ## export functions
        toExport <- c(".BuyseTest","calcPeron","calcSample")
        iB <- NULL ## [:forCRANcheck:] foreach

        ls2.resampling <- foreach::`%dopar%`(
                                       foreach::foreach(iB=1:nsplit.resampling,
                                                        .export = toExport,
                                                        .packages = c("data.table","BuyseTest"),
                                                        .options.snow = opts), {
                                           iOut <- lapply(split.resampling[[iB]], function(iSplit){
                                               if(!is.null(seed)){set.seed(seqSeed[iSplit])}
                                               iBT <- .BuyseTest(envir = envir,
                                                                 iid = iid,
                                                                 method.inference = method.inference,
                                                                 pointEstimation = FALSE)
                                               if(!is.null(seed)){
                                                   return(c(iBT,list(seed = seqSeed[iSplit])))
                                               }else{
                                                   return(iBT)
                                               }
                                           })
                                           return(iOut)
                                       })        

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}

        ## collect
        ls.resampling <- do.call("c",ls2.resampling)
    }

    ## ** post treatment
    test.resampling <- which(unlist(lapply(ls.resampling,is.null)) == FALSE)
    if(length(test.resampling) != n.resampling){
        n.failure <- n.resampling - length(test.resampling) 
        warning("The resampling procedure failed for ",n.failure," samples (",round(100*n.failure/n.resampling,2),"%)")
    }
    
    dim.delta <- c(n.resampling, n.strata, D, 6)
    dimnames.delta <- list(as.character(1:n.resampling), level.strata, endpoint, c("favorable","unfavorable","neutral","uninf","netBenefit","winRatio"))
    out <- list(deltaResampling = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                DeltaResampling = array(NA, dim = dim.delta[c(1,3,4)], dimnames = dimnames.delta[c(1,3,4)]),
                weightStrataResampling = matrix(NA, nrow = n.resampling, ncol = n.strata,
                                                dimnames = list(NULL, level.strata))
                )
    if(!is.null(seed)){
        out$seed <- rep(NA, n.resampling)
    }
    if(iid){
        out$covarianceResampling = array(NA, dim = c(n.resampling, D, 5),
                                         dimnames = list(as.character(1:n.resampling), endpoint, c("favorable", "unfavorable", "covariance", "netBenefit", "winRatio")))
    }else{
        out$covarianceResampling <- array(NA, dim = c(0,0,0))
    }
    
    for(iR in test.resampling){ ## iR <- 1
        out$deltaResampling[iR,,,] <- ls.resampling[[iR]]$delta
        out$DeltaResampling[iR,,] <- ls.resampling[[iR]]$Delta
        out$weightStrataResampling[iR,] <- ls.resampling[[iR]]$weightStrata
        if(!is.null(seed)){
            out$seed[iR] <- seqSeed[iR]
        }
        if(iid){
            out$covarianceResampling[iR,,] <- ls.resampling[[iR]]$covariance
        }

    }

    ## ** export
    return(out)
}


## * inference U-statistic (Bebu et al 2015)
##' @description Implement the computation of the asymptotic variance as described in Bebu et al (2015)
##' Should give results equivalent to inferenceUstatistic
##' NOTE: arguments subset.C and subset.T were used for BuysePower to re-compute statistics on a subset of the data
##'       but this happens to be slower than just re-running the test so is not used
##' @noRd
##' @references Large sample inference for a win ratio analysis of a composite outcome based on prioritized components Biostatistics (2015), pp. 1–10 doi:10.1093/biostatistics/kxv032
inferenceUstatisticBebu <- function(tablePairScore, subset.C = NULL, subset.T = NULL,
                                    order, weightEndpoint, 
                                    n.pairs, n.C, n.T, level.strata, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    out <- list()
    
    ## ** extract informations
    n.endpoint <- length(endpoint)
    ntot.pairs <- sum(n.pairs)
        
    ## ** merge tables
    ls.table <- wsumPairScore(tablePairScore, weightEndpoint = weightEndpoint,
                              subset.C = subset.C, subset.T = subset.T)

    out$n.pairs <- NROW(ls.table[[1]])
    count.favorable <- do.call(cbind,lapply(ls.table, function(iTable){sum(iTable$favorable)}))
    count.unfavorable <- do.call(cbind,lapply(ls.table, function(iTable){sum(iTable$unfavorable)}))

    p.favorable <- as.double(count.favorable/sum(n.pairs))
    p.unfavorable <- as.double(count.unfavorable/sum(n.pairs))

    ## out$Delta <- unname(cbind(p.favorable,
    ##                           p.unfavorable,
    ##                           p.favorable - p.unfavorable,
    ##                           count.favorable / count.unfavorable))
    
    ## out$delta <- array(NA, dim = c(n.strata, n.endpoint, 4))
    ## out$count_favorable <- matrix(NA, nrow = n.strata, ncol = n.endpoint)
    ## out$count_unfavorable <- matrix(NA, nrow = n.strata, ncol = n.endpoint)
    ## out$count_neutral <- matrix(NA, nrow = n.strata, ncol = n.endpoint)
    ## out$count_uninf <- matrix(NA, nrow = n.strata, ncol = n.endpoint)

    ## ** compute variance component over strata
    M.cov <- matrix(0, nrow = n.endpoint, ncol = 3,
                    dimnames = list(endpoint, c("favorable","unfavorable","covariance")))
    
    strataSum <- matrix(NA, nrow = 6, ncol = n.endpoint,
                        dimnames = list(c("favorableT","favorableC","unfavorableT","unfavorableC","mixedC","mixedT"),
                                        endpoint))

    for(iStrata in 1:n.strata){ ## iStrata <- 1

        iN.strata <- n.pairs[iStrata]
        iLS.table <- lapply(ls.table, function(iT){iT[iT$strata == level.strata[iStrata]]})
        iDT.nCT <- iLS.table[[1]][,.(n.C = length(unique(.SD$indexWithinStrata.C)),n.T = length(unique(.SD$indexWithinStrata.T)))]
        iN.C <- iDT.nCT$n.C
        iN.T <- iDT.nCT$n.T

        iP.favorable <- unlist(lapply(iLS.table, function(iT){sum(iT$favorable)}))/n.pairs[iStrata]
        iP.unfavorable <- unlist(lapply(iLS.table, function(iT){sum(iT$unfavorable)}))/n.pairs[iStrata]

        for(iE in 1:n.endpoint){ ## iE <- 1

            iTable <- iLS.table[[iE]]
            
            index2originalOrder.C <- iTable[!duplicated(iTable$index.C),
                                            stats::setNames(.SD$index.C,.SD$indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(iTable$index.T),
                                            stats::setNames(.SD$index.T,.SD$indexWithinStrata.T)]
            ## *** Hajek projection
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.T"]
            sumPair.T[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.T[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
            iN.setT <- NROW(sumPair.T)
                
            ## \E[X_i>=Y_j+\tau|Y_j] and \E[X_i+\tau<=Y_j|Y_j]
            sumPair.C <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.C"]
            sumPair.C[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.C[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
            iN.setC <- NROW(sumPair.C)
            
            ## *** variance 
            ## P[1(X_i,Y_j)1(X_i,Y_k)] = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) sum_k neq j 1(X_i,Y_j)
            ##                         = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) ( sum_k 1(X_i,Y_k) - 1(X_i,Y_j) )
            ## here we compute sum_k 1(X_i,Y_k) and m-1
            iTable[, c("sumFavorable.T") := sumPair.T$favorable[.SD$indexWithinStrata.T]]
            iTable[, c("sumUnfavorable.T") := sumPair.T$unfavorable[.SD$indexWithinStrata.T]]
            
            iTable[, c("sumFavorable.C") := sumPair.C$favorable[.SD$indexWithinStrata.C]]
            iTable[, c("sumUnfavorable.C") := sumPair.C$unfavorable[.SD$indexWithinStrata.C]]

            
            if(iN.setT > 0){
                ## E[ 1(X_i>Y_j) 1(X_i>Y_k) ]
                strataSum["favorableT",iE] <- iTable[,sum(.SD$sumFavorable.T * .SD$favorable)] / (iN.strata*iN.setT)
                ## E[ 1(X_i<Y_j) 1(X_i<Y_k) ]
                strataSum["unfavorableT",iE] <- iTable[,sum(.SD$sumUnfavorable.T * .SD$unfavorable)] / (iN.strata*iN.setT)
                ## E[ 1(X_i>Y_j) 1(X_i<Y_k) ]
                strataSum["mixedT",iE] <- iTable[,sum(.SD$sumUnfavorable.T * .SD$favorable)] / (iN.strata*iN.setT)
            }
            if(iN.setC > 0){
                ## E[ 1(X_i>Y_j) 1(X_k>Y_j) ]
                strataSum["favorableC",iE] <- iTable[,sum(.SD$sumFavorable.C * .SD$favorable)] / (iN.strata*iN.setC)
                ## E[ 1(X_i<Y_j) 1(X_k<Y_j) ]
                strataSum["unfavorableC",iE] <- iTable[,sum(.SD$sumUnfavorable.C * .SD$unfavorable)] / (iN.strata*iN.setC)
                ## E[ 1(X_i>Y_j) 1(X_k<Y_j) ]
                strataSum["mixedC",iE] <- iTable[,sum(.SD$sumUnfavorable.C * .SD$favorable)] / (iN.strata*iN.setC)
            }

        }

        ## *** first order terms: compute xi
        ## P[X1>Y1 & X1>Y1'] - P[X1>Y1]^2
        xi_10_11 <- strataSum["favorableT",] - iP.favorable^2
        ## P[X1>Y1 & X1'>Y1] - P[X1>Y1]^2
        xi_01_11 <- strataSum["favorableC",] - iP.favorable^2

        ## P[X1<Y1 & X1<Y1'] - P[X1<Y1]^2
        xi_10_22 <- strataSum["unfavorableT",] - iP.unfavorable^2
        ## P[X1<Y1 & X1'<Y1] - P[X1<Y1]^2
        xi_01_22 <- strataSum["unfavorableC",] - iP.unfavorable^2
    
        ## P[X1>Y1 & X1<Y1'] - P[X1>Y1]*P[X1<Y1]
        xi_10_12 <- strataSum["mixedC",] - iP.favorable * iP.unfavorable
        ## P[X1>Y1 & X1'<Y1] - P[X1>Y1]*P[X1<Y1]
        xi_01_12 <- strataSum["mixedT",] - iP.favorable * iP.unfavorable

        ## *** second order terms
        if(order == 2){
            Mfav <- do.call(cbind,lapply(iLS.table,"[[","favorable"))
            Munfav <- do.call(cbind,lapply(iLS.table,"[[","unfavorable"))
            
            varUfav <- colMeans(Mfav^2) - colMeans(Mfav)^2 ##instead of p1.favorable*(1-p1.favorable)
            varUunfav <- colMeans(Munfav^2) - colMeans(Munfav)^2 ##instead of p1.unfavorable*(1-p1.unfavorable)
            covUfavunfav <- colMeans(Mfav*Munfav) - colMeans(Mfav)*colMeans(Munfav) ##instead of -p1.favorable*p1.unfavorable
            
            H2.favorable <- (varUfav - xi_10_11 - xi_01_11)/(iN.strata)
            H2.unfavorable <- (varUunfav - xi_10_22 - xi_01_22)/(iN.strata)
            H2.covariance <- (covUfavunfav - xi_10_12 - xi_01_12)/(iN.strata)
        }else{
            H2.favorable <- 0
            H2.unfavorable <- 0
            H2.covariance <- 0
        }
        
        ## ** compute sigma
        ## NO STRATA:
        ## N.TC = N.T+N.C
        ## Sigma = N.TC/N.C SigmaC + N.TT/N.T SigmaT
        ## asymptotic variance i.e. sqrt(N.TC)(Uhat - U) \sim N(0,Sigma)
        ## scaled asymptotic variance i.e. (Uhat - U) \sim N(0,Sigma/N.TC) = N(0,1/N.C SigmaC + 1/N.T SigmaT)
        ##
        ## STRATA:
        ## same but adding a factor n.strata / N.TC to accound for pooling
        M.cov <- M.cov + (iN.strata/ntot.pairs)^2 * cbind(favorable =  xi_10_11 / iN.C + xi_01_11 / iN.T + H2.favorable,
                                                          unfavorable = xi_10_22 / iN.C + xi_01_22 / iN.T + H2.unfavorable,
                                                          covariance = xi_10_12 / iN.C + xi_01_12 / iN.T + H2.covariance)
    }

    ## ** export
    out$Sigma <- cbind(M.cov,
                       "netBenefit" = M.cov[,"favorable"] + M.cov[,"unfavorable"] - 2 * M.cov[,"covariance"],
                       "winRatio" =  M.cov[,"favorable"]/p.unfavorable^2 + M.cov[,"unfavorable"]*p.favorable^2/p.unfavorable^4 - 2 * M.cov[,"covariance"]*p.favorable/p.unfavorable^3
                       )
    return(out)
}

## * inference permutation (Anderson and Verbeeck 2023)
##' @description Implement the computation of the variance of the permutation distribution as described in Anderson and Verbeeck (2023)
##' exactVarPermutation, the core of this function, is essentially a copy of the code provided by the authors (Anderson and Verbeeck).
##' @noRd
##' @references William N Anderson and Johan Verbeeck. Exact permutation and bootstrap distribution of generalized pairwise comparisons statistics (2023), Mathematics.
inferenceVarPermutation <- function(data, treatment, level.treatment, weightStrata, ...){

    ## ** extend data
    n.obs <- NROW(data)
    if("XXgroupXX" %in% names(data)){
        stop("BuyseTest: Argument \'data\' must not contain a column \"XXgroupXX\". \n")
    }
    if("XXindexXX" %in% names(data)){
        stop("BuyseTest: Argument \'data\' must not contain a column \"XXindexXX\". \n")
    }
    data.ext <- rbind(cbind(XXgroupXX = "A", XXindexXX = 1:n.obs, data), 
                      cbind(XXgroupXX = "B", XXindexXX = 1:n.obs, data))

    ## ** run GPC on all pairs
    eBT.all <- BuyseTest(data = data.ext, treatment = "XXgroupXX", ,
                         method.inference = "none", keep.pairScore = TRUE,
                         trace = FALSE, ...)
    endpoint <- eBT.all@endpoint
    strata <- eBT.all@level.strata
    if(any("Peron" %in% eBT.all@scoring.rule)){
        warning("BuyseTest: the current implementation of the exact variance of the permutation distribution will not provide type 1 error control when using the Peron scoring rule. \n")
    }else if(any(eBT.all@correction.uninf>0)){
        warning("BuyseTest: the current implementation of the exact variance of the permutation distribution will not provide type 1 error control when using a correction for uninformative pairs. \n")
    }

    ## ** extract all pairwise scores
    eBTscore.all <- getPairScore(eBT.all, endpoint = endpoint, cumulative = TRUE, unlist = FALSE)        

    ## ** create score matrix
    U.favorable <- lapply(endpoint, function(iE){matrix(0, nrow = n.obs, ncol = n.obs)})
    U.unfavorable <- lapply(endpoint, function(iE){matrix(0, nrow = n.obs, ncol = n.obs)})
    U.score <- lapply(endpoint, function(iE){matrix(0, nrow = n.obs, ncol = n.obs)})
        
    for(iEndpoint in endpoint){
        eBTscore.all[[iEndpoint]]$pos.A <- data.ext$XXindexXX[eBTscore.all[[iEndpoint]]$index.A] ## location in the original dataset
        eBTscore.all[[iEndpoint]]$pos.B <- data.ext$XXindexXX[eBTscore.all[[iEndpoint]]$index.B] ## location in the original dataset

        U.favorable[[iEndpoint]][(eBTscore.all[[iEndpoint]]$pos.B-1)*n.obs + eBTscore.all[[iEndpoint]]$pos.A] <- eBTscore.all[[iEndpoint]]$favorable
        U.unfavorable[[iEndpoint]][(eBTscore.all[[iEndpoint]]$pos.B-1)*n.obs + eBTscore.all[[iEndpoint]]$pos.A] <- eBTscore.all[[iEndpoint]]$unfavorable
        U.score[[iEndpoint]] <- U.favorable[[iEndpoint]] - U.unfavorable[[iEndpoint]]
    }

    ## ** evaluate permutation variance
    if(length(strata)==1){
        ls.permvar <- lapply(endpoint, function(iEndpoint){
            exactVarPermutation(U.score[[iEndpoint]], treatment = data[[treatment]], level.treatment = level.treatment)
        })
    }else{
        index.strata <- lapply(attr(strata,"index"), intersect, 1:n.obs) ## to retrieve original size since the dataset was duplicated
        ls.permvar <- lapply(endpoint, function(iEndpoint){ ## iEndpoint <- "score"
            iM.permvar <- do.call(rbind,lapply(index.strata, function(iIndex){ ## iIndex <- index.strata[[1]]
                exactVarPermutation(U.score[[iEndpoint]][iIndex,iIndex], treatment = data[iIndex,treatment], level.treatment = level.treatment)
            }))
            rownames(iM.permvar) <- strata
            iOut <- colSums(sweep(iM.permvar, MARGIN = 1, FUN = "*", STATS = weightStrata^2))
            attr(iOut,"strata") <- iM.permvar
            return(iOut)
        })
    }

    ## ** reshape
    out <- do.call(rbind, ls.permvar)
    rownames(out) <- endpoint
    if(length(strata)>1){
        attr(out,"strata") <- lapply(ls.permvar,attr,"strata")
    }
    return(out)
}

## * wsumPairScore
##' @description cumulate over endpoint the scores
##' @noRd
wsumPairScore <- function(pairScore, weightEndpoint, subset.C, subset.T){

    keep.col <- c("strata","index.C","index.T","index.pair","indexWithinStrata.C", "indexWithinStrata.T","favorableC","unfavorableC","neutralC","uninfC")
    old.col <- c("favorableC","unfavorableC","neutralC","uninfC")
    new.col <- c("favorable","unfavorable","neutral","uninf")

    n.endpoint <- length(pairScore)
    out <- vector(mode = "list", length = n.endpoint)
    if(!is.null(subset.C)){
        subset.numC <- sort(unique(pairScore[[1]]$index.C))[subset.C]
    }
    if(!is.null(subset.T)){
        subset.numT <- sort(unique(pairScore[[1]]$index.T))[subset.T]
    }
    ## indexPair <- stats::setNames(1:NROW(pairScore[[1]]),pairScore[[1]]$index.pair)
    for(iE in 1:n.endpoint){ ## iE <- 2
        iTable <- data.table::copy(pairScore[[iE]][,.SD,.SDcols = keep.col])
        if(!is.null(subset.C)){
            iTable <- iTable[iTable$index.C %in% subset.numC]
        }
        if(!is.null(subset.T)){
            iTable <- iTable[iTable$index.T %in% subset.numT]
        }
        data.table::setnames(iTable, old = old.col, new = new.col)
        iTable[,c("favorable") := .SD$favorable * weightEndpoint[iE]]
        iTable[,c("unfavorable") := .SD$unfavorable * weightEndpoint[iE]]
        iTable[,c("neutral") := .SD$neutral * weightEndpoint[iE]]
        iTable[,c("uninf") := .SD$uninf * weightEndpoint[iE]]
        
        if(iE==1){
            out[[iE]] <- iTable
        }else{
            out[[iE]] <- data.table::copy(out[[iE-1]])
            indexMatch <- match(iTable$index.pair,out[[iE]]$index.pair) ## indexMatch - iTable$index.pair
            out[[iE]][indexMatch, c("favorable") := .SD$favorable + iTable$favorable]
            out[[iE]][indexMatch, c("unfavorable") := .SD$unfavorable + iTable$unfavorable]
            out[[iE]][indexMatch, c("neutral") := .SD$neutral + iTable$neutral]
            out[[iE]][indexMatch, c("uninf") := .SD$uninf + iTable$uninf]
        }
    }
    return(out)
}

## * exactVarPermutation
## Function adapted from William Anderson and Johan Verbeeck
exactVarPermutation <- function(winmatrix, treatment, level.treatment){
    if (any(winmatrix + t(winmatrix)  != 0)) stop("Win matrix not skew")
    m <- sum(treatment == level.treatment[1])
    n <- sum(treatment == level.treatment[2])
    N <- m + n   
  
    ## compute various vertex values
    ## computations are O(N^2), because the matrix multiplications are by component 
    ## sum of values of edges pointing in to the vertex -- real analog of indegree
    invalues <- colSums(winmatrix*as.integer(winmatrix > 0))  
    insqvalues <- colSums(winmatrix^2*as.integer(winmatrix > 0))  
    ## sum of values of edges pointing out of the vertex -- real analog of outdegree 
    outvalues <- rowSums(winmatrix*as.integer(winmatrix > 0))  
    ## sum of squares of values of edges pointing out of the vertex
    outsqvalues <- rowSums(winmatrix^2*as.integer(winmatrix > 0))
    edgesum <- sum(invalues) ## could also use outvalues
    ## count the cases at each vertex v -- case definitions in manuscript
    ## the case matrix has 5 rows (1 for each case), and N columns
    ## case 6 is not represented in this matrix, because case 6 situations do not belong to a specific vertex
    ## the matrix is perhaps useful in understanding the algorithm, but only the rowSums are actually used
    casematrix <- rbind(insqvalues,
                        invalues^2 - insqvalues,
                        outvalues^2 - outsqvalues,
                        invalues*outvalues,
                        invalues*outvalues)
    ## add the cases from all the vertices
    casecounts <- rowSums(casematrix)
    casecounts[6] <- edgesum^2 - sum(casecounts) # now we have case 6
    expected =  rep(edgesum*m*n/(N*(N - 1)), 2)
    names(expected) <- c("Test Wins", "Control Wins")
    TTerms <- c(m*n/(N*(N - 1)), m*n*(m - 1)/(N*(N - 1)*(N - 2)),
                m*n*(n - 1)/(N*(N - 1)*(N - 2)), 0, 0, 
                m*n*(m - 1)*(n - 1)/(N*(N - 1)*(N - 2)*(N - 3)))
    CTerms <- c(m*n/(N*(N - 1)), m*n*(n - 1)/(N*(N - 1)*(N - 2)),
                m*n*(m - 1)/(N*(N - 1)*(N - 2)), 0, 0,
                m*n*(m - 1)*(n - 1)/(N*(N - 1)*(N - 2)*(N - 3)))
    TCTerms <- c(0, 0, 0, m*n*(n - 1)/(N*(N - 1)*(N - 2)),
                 m*n*(m - 1)/(N*(N - 1)*(N - 2)),
                 m*n*(m - 1)*(n - 1)/(N*(N - 1)*(N - 2)*(N - 3)))
    CTTerms <- c(0, 0, 0, m*n*(m - 1)/(N*(N - 1)*(N - 2)),
                 m*n*(n - 1)/(N*(N - 1)*(N - 2)),
                 m*n*(m - 1)*(n - 1)/(N*(N - 1)*(N - 2)*(N - 3)))
    expectedforvariance <- c(sum(casecounts*TTerms), sum(casecounts*TCTerms),
                             sum(casecounts*CTTerms), sum(casecounts*CTerms))
    expectedforvariance <- matrix(expectedforvariance, nrow = 2)
    variance <- expectedforvariance - expected %*% t(expected)
    rownames(variance) <- c("Test Win Sum", "Control Win Sum")
    colnames(variance) <- c("Test Win Sum", "Control Win Sum")

    out <- c(favorable = variance[1,1]/(m*n)^2,
             unfavorable = variance[2,2]/(m*n)^2,
             covariance = variance[1,2]/(m*n)^2,
             netBenefit = (variance[1,1] + variance[2,2] - 2*variance[1,2])/(m*n)^2,
             winRatio = as.numeric(NA))

    ## testrows <- which(treatment == level.treatment[2])
    ## comparisonmatrix <- winmatrix[testrows, -testrows]
    ## K <- pmax(comparisonmatrix, 0) # Test wins are the positive matrix entries; test losses and ties are 0
    ## L <- pmax(-comparisonmatrix, 0) # Control wins are the positive matrix entries; control losses and ties are 0
    ## winssum <- c("Test Win Sum" = sum(K), "Control Win Sum" = sum(L))
    ## stats <- c(winssum-expected,winssum[1]-winssum[2])/sqrt(c(variance[1,1],variance[2,2], variance[1,1]+variance[2,2]-2*variance[1,2]))
    ## 2*(1-pnorm(abs(stats)))

    return(out)
}
