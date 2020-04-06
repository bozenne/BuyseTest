## * inferenceResampling
## author Brice Ozenne
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
        if (!is.null(seed)) {set.seed(seed)} # set the seed

        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.resampling <- do.call(method.loop,
                                  args = list(X = 1:n.resampling,
                                              FUN = function(iB){
                                                  .BuyseTest(envir = envir,
                                                             iid = iid,
                                                             method.inference = method.inference,
                                                             pointEstimation = FALSE
                                                             )
                                              })
                                 )

        if(!is.null(seed)){rm(.Random.seed, envir=.GlobalEnv)} # restaure original seed
    }else { ## *** parallel resampling test

        ## define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.resampling, style = 3)          
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
        toExport <- c(".BuyseTest","calcPeron","calcSample")
        iB <- NULL ## [:forCRANcheck:] foreach        
        ls.resampling <- foreach::`%dopar%`(
                                      foreach::foreach(iB=1:n.resampling,
                                                       .export = toExport,
                                                       .packages = "data.table"),                                            
                                      {                                           
                                          if(trace>0){utils::setTxtProgressBar(pb, iB)}

                                           return(.BuyseTest(envir = envir,
                                                             iid = iid,
                                                             method.inference = method.inference,
                                                             pointEstimation = FALSE))
                      
                                       })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }

    ## ** post treatment
    test.resampling <- which(unlist(lapply(ls.resampling,is.null)) == FALSE)
    if(length(test.resampling) != n.resampling){
        n.failure <- n.resampling - length(test.resampling) 
        warning("The resampling procedure failed for ",n.failure," samples (",round(100*n.failure/n.resampling,2),"%)")
    }
    
    dim.delta <- c(n.resampling, D, 4, n.strata)
    dimnames.delta <- list(as.character(1:n.resampling), endpoint, c("favorable","unfavorable","netBenefit","winRatio"), level.strata)

    out <- list(deltaResampling = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                DeltaResampling = array(NA, dim = dim.delta[1:3], dimnames = dimnames.delta[1:3])
                )
    if(iid){
        out$covarianceResampling = array(NA, dim = c(n.resampling, D, 5))
    }else{
        out$covarianceResampling <- array(NA, dim = c(0,0,0))
    }
    
    for(iR in test.resampling){
        out$deltaResampling[iR,,,] <- ls.resampling[[iR]]$delta
        out$DeltaResampling[iR,,] <- ls.resampling[[iR]]$Delta
        
        if(iid){
            out$covarianceResampling[iR,,] <- ls.resampling[[iR]]$covariance
        }

    }

    ## ** export
    return(out)
}


## * inference U-statistic
## Implement the computation of the asymptotic variance via an Hajek projection
## used by BuysePower
inferenceUstatistic <- function(tablePairScore, order, weight, count.favorable, count.unfavorable,
                                n.pairs, n.C, n.T, level.strata, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    n.endpoint <- length(endpoint)
    
    ## ** merge tables
    ls.table <- wsumPairScore(tablePairScore, weight = weight, n.endpoint = n.endpoint)

    ## ** H-decomposition
    ## expectation
    U.favorable <- cumsum(count.favorable)/n.pairs
    U.unfavorable <- cumsum(count.unfavorable)/n.pairs

    ## apply(rbind(favorable = Upartial.favorable, unfavorable = Upartial.unfavorable),1,cumsum)
    ## lapply(ls.table, function(iDT){c(sum(iDT$favorable),sum(iDT$unfavorable))/NROW(ls.table[[1]])})
           
    ## storage
    out <- list(iidAverage_favorable = matrix(NA, nrow = n.T+n.C, ncol = n.endpoint, dimnames = list(NULL, endpoint)),
                iidAverage_unfavorable = matrix(NA, nrow = n.T+n.C, ncol = n.endpoint, dimnames = list(NULL, endpoint)),
                covariance = matrix(NA, nrow = n.endpoint, ncol = 5, dimnames = list(NULL,c("favorable","unfavorable","covariance","netBenefit","winRatio")))
                )
    if(order == 2){
        iidAverage_favorable2 <- matrix(NA, nrow = n.pairs, ncol = n.endpoint, dimnames = list(NULL, endpoint))
        iidAverage_unfavorable2 <- matrix(NA, nrow = n.pairs, ncol = n.endpoint, dimnames = list(NULL, endpoint))
    }

    ## loop
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        for(iE in 1:n.endpoint){ ## iE <- 1

            ## extract pairwise scores
            iTable <- ls.table[[iE]][ls.table[[iE]]$strata == level.strata[iStrata]]
            index2originalOrder.C <- iTable[!duplicated(iTable$index.C),
                                            stats::setNames(.SD$index.C,.SD$indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(iTable$index.T),
                                            stats::setNames(.SD$index.T,.SD$indexWithinStrata.T)]
            iN.strata <- NROW(iTable) ## number of pairs

            ## *** Hajek projection
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, list(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "index.T"]
            sumPair.T[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.T[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
            
            ## \E[X_i>=Y_j+\tau|Y_j] and \E[X_i+\tau<=Y_j|Y_j]
            sumPair.C <- iTable[, list(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "index.C"]
            sumPair.C[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.C[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]

            ## store
            out$iidAverage_favorable[index2originalOrder.C,iE] <- weight[iE] * (sumPair.C$E.favorable - U.favorable[iE]) / n.C
            out$iidAverage_favorable[index2originalOrder.T,iE] <- weight[iE] * (sumPair.T$E.favorable - U.favorable[iE]) / n.T
            out$iidAverage_unfavorable[index2originalOrder.C,iE] <- weight[iE] * (sumPair.C$E.unfavorable - U.unfavorable[iE]) / n.C
            out$iidAverage_unfavorable[index2originalOrder.T,iE] <- weight[iE] * (sumPair.T$E.unfavorable - U.unfavorable[iE]) / n.T
            
            ## *** second order
            if(order == 2){
                iidAverage_favorable2[,iE] <- weight[iE] * (iTable$favorable - out$iidAverage_favorable[iTable$index.C,iE] * n.C - out$iidAverage_favorable[iTable$index.T,iE] * n.T - U.favorable[iE])/n.pairs
                iidAverage_unfavorable2[,iE] <- weight[iE] * (iTable$unfavorable - out$iidAverage_unfavorable[iTable$index.C,iE] * n.C - out$iidAverage_unfavorable[iTable$index.T,iE] * n.T - U.unfavorable[iE])/n.pairs
            }
        }
    }

    ## ** compute Sigma
    if(n.endpoint==1){
        out$covariance[1,"favorable"] <- sum(out$iidAverage_favorable^2)
        out$covariance[1,"unfavorable"] <- sum(out$iidAverage_unfavorable^2)
        out$covariance[1,"covariance"] <- sum(out$iidAverage_favorable * out$iidAverage_unfavorable)

        if(order == 2){
            out$covariance[,"favorable"] <- out$covariance[,"favorable"] + sum(iidAverage_favorable2^2)
            out$covariance[,"unfavorable"] <- out$covariance[,"unfavorable"] + sum(iidAverage_unfavorable2^2)
            out$covariance[,"covariance"] <- out$covariance[,"covariance"] + sum(iidAverage_favorable2 * iidAverage_unfavorable2)
        }
    }else{
        ## cumsum because the iid decomposition is endpoint specific while the net benefit is the overall        
        out$covariance[,"favorable"] <- colSums(out$iidAverage_favorable^2)
        out$covariance[,"unfavorable"] <- colSums(out$iidAverage_unfavorable^2)
        out$covariance[,"covariance"] <- colSums(out$iidAverage_favorable * out$iidAverage_unfavorable)
        if(order == 2){
            out$covariance[,"favorable"] <- out$covariance[,"favorable"] + colSums(iidAverage_favorable2^2)
            out$covariance[,"unfavorable"] <- out$covariance[,"unfavorable"] + colSums(iidAverage_unfavorable2^2)
            out$covariance[,"covariance"] <- out$covariance[,"covariance"] + colSums(iidAverage_favorable2 * iidAverage_unfavorable2)
    
            ## // update variance:  Var(H2) = Var( 1/mn \sum_i,j H2_ij) = 1/mn Var(H2_ij)
            ## Mvar.col(0) += trans((H2_moments.row(2) - pow(H2_moments.row(0),2))/(ntot_pair));
            ## Mvar.col(1) += trans((H2_moments.row(3) - pow(H2_moments.row(1),2))/(ntot_pair));
            ## Mvar.col(2) += trans((H2_moments.row(4) - H2_moments.row(0) % H2_moments.row(1))/(ntot_pair));
        }
    }

    out$covariance[,"netBenefit"] <- out$covariance[,"favorable"] + out$covariance[,"unfavorable"] - 2 * out$covariance[,"covariance"]
    out$covariance[,"winRatio"] <- out$covariance[,"favorable"]/U.unfavorable^2 + out$covariance[,"unfavorable"] * U.favorable^2/U.unfavorable^4 - 2 * out$covariance[,"covariance"] * U.favorable/U.unfavorable^3
    
    ## ** export
    return(out)
}

## * inference U-statistic (Bebu et al 2015)
## Implement the computation of the asymptotic variance as described in
## Large sample inference for a win ratio analysis of a composite outcome based on prioritized components
## Biostatistics (2015), pp. 1–10 doi:10.1093/biostatistics/kxv032
## Give results equivalent to inferenceUstatistic

inferenceUstatisticBebu <- function(tablePairScore, order, weight, count.favorable, count.unfavorable,
                                    n.pairs, n.C, n.T, level.strata, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    n.endpoint <- length(endpoint)
    ntot.pairs <- sum(n.pairs)
        
    p1.favorable <- cumsum(count.favorable)/ntot.pairs
    p1.unfavorable <- cumsum(count.unfavorable)/ntot.pairs

    ## ** merge tables
    ls.table <- wsumPairScore(tablePairScore, weight = weight, n.endpoint = n.endpoint)

    ## ** compute variance component over strata
    M.cov <- matrix(0, nrow = n.endpoint, ncol = 3,
                    dimnames = list(endpoint, c("favorable","unfavorable","covariance")))
    
    strataSum <- matrix(NA, nrow = 6, ncol = n.endpoint,
                        dimnames = list(c("favorableT","favorableC","unfavorableT","unfavorableC","mixedC","mixedT"),
                                        endpoint))

    for(iStrata in 1:n.strata){ ## iStrata <- 1

        iN.strata <- n.pairs[iStrata]
        iDT.nCT <- ls.table[[1]][ls.table[[1]]$strata == level.strata[iStrata],
                                 .(n.C = length(unique(.SD$indexWithinStrata.C)),n.T = length(unique(.SD$index.T)))]
        iN.C <- iDT.nCT$n.C
        iN.T <- iDT.nCT$n.T
        
        for(iE in 1:n.endpoint){ ## iE <- 1

            iTable <- ls.table[[iE]][ls.table[[iE]]$strata == level.strata[iStrata]]
            
            index2originalOrder.C <- iTable[!duplicated(iTable$index.C),
                                            stats::setNames(.SD$index.C,.SD$indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(iTable$index.T),
                                            stats::setNames(.SD$index.T,.SD$indexWithinStrata.T)]
            ## *** Hajek projection
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "index.T"]
            sumPair.T[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.T[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
                
            ## \E[X_i>=Y_j+\tau|Y_j] and \E[X_i+\tau<=Y_j|Y_j]
            sumPair.C <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.C"]
            sumPair.C[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.C[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]

            ## *** variance 
            ## P[1(X_i,Y_j)1(X_i,Y_k)] = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) sum_k neq j 1(X_i,Y_j)
            ##                         = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) ( sum_k 1(X_i,Y_k) - 1(X_i,Y_j) )
            ## here we compute sum_k 1(X_i,Y_k) and m-1

            iTable[, c("sumFavorable.T") := sumPair.T$favorable[.SD$indexWithinStrata.T]]
            iTable[, c("sumUnfavorable.T") := sumPair.T$unfavorable[.SD$indexWithinStrata.T]]
            
            iTable[, c("sumFavorable.C") := sumPair.C$favorable[.SD$indexWithinStrata.C]]
            iTable[, c("sumUnfavorable.C") := sumPair.C$unfavorable[.SD$indexWithinStrata.C]]

            iN.setT <- NROW(sumPair.T)
            iN.setC <- NROW(sumPair.C)
            
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
        xi_10_11 <- strataSum["favorableT",] - p1.favorable^2
        ## P[X1>Y1 & X1'>Y1] - P[X1>Y1]^2
        xi_01_11 <- strataSum["favorableC",] - p1.favorable^2

        ## P[X1<Y1 & X1<Y1'] - P[X1<Y1]^2
        xi_10_22 <- strataSum["unfavorableT",] - p1.unfavorable^2
        ## P[X1<Y1 & X1'<Y1] - P[X1<Y1]^2
        xi_01_22 <- strataSum["unfavorableC",] - p1.unfavorable^2
    
        ## P[X1>Y1 & X1<Y1'] - P[X1>Y1]*P[X1<Y1]
        xi_10_12 <- strataSum["mixedC",] - p1.favorable * p1.unfavorable
        ## P[X1>Y1 & X1'<Y1] - P[X1>Y1]*P[X1<Y1]
        xi_01_12 <- strataSum["mixedT",] - p1.favorable * p1.unfavorable

        ## *** second order terms
        if(order == 2){
            Mfav <- do.call(cbind,lapply(ls.table,"[[","favorable"))
            Munfav <- do.call(cbind,lapply(ls.table,"[[","unfavorable"))
            
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
    iwFavorable <- cumsum(count.favorable * weight) / ntot.pairs
    iwUnfavorable <- cumsum(count.unfavorable * weight) / ntot.pairs
    return(list(Sigma = cbind(M.cov,
                              "netBenefit" = M.cov[,"favorable"] + M.cov[,"unfavorable"] - 2 * M.cov[,"covariance"],
                              "winRatio" =  M.cov[,"favorable"]/iwUnfavorable^2 + M.cov[,"unfavorable"]*iwFavorable^2/iwUnfavorable^4 - 2 * M.cov[,"covariance"]*iwFavorable/iwUnfavorable^3
                              ),
                iid1 = NULL,
                iid2 = NULL))
}

## * .iid2cov
.iid2cov <- function(A.iid, A2.iid, weight,
                     order, endpoint, n.endpoint){
    
    if(n.endpoint==1){
        M.cov <- cbind(favorable = sum(A.iid[,,"favorable"]^2),
                       unfavorable = sum(A.iid[,,"unfavorable"]^2),
                       covariance = sum(A.iid[,,"favorable"] * A.iid[,,"unfavorable"]))

        if(order == 2){
            M.cov[1,"favorable"] <- M.cov[1,"favorable"] + sum(A2.iid[,,"favorable"]^2)
            M.cov[1,"unfavorable"] <- M.cov[1,"unfavorable"] + sum(A2.iid[,,"unfavorable"]^2)
            M.cov[1,"covariance"] <- M.cov[1,"covariance"] + sum(A2.iid[,,"favorable"] * A2.iid[,,"unfavorable"])
        }
    }else{
        ## cumsum because the iid decomposition is endpoint specific while the net benefit is the overall        
        favorable.cumiid <- t(apply(A.iid[,,"favorable"],1,cumsum))
        unfavorable.cumiid <- t(apply(A.iid[,,"unfavorable"],1,cumsum))

        M.cov <- cbind(favorable = colSums(favorable.cumiid^2),
                       unfavorable = colSums(unfavorable.cumiid^2),
                       covariance = colSums(favorable.cumiid * unfavorable.cumiid))

        if(order == 2){
            favorable.cumiid2 <- t(apply(A2.iid[,,"favorable"],1,cumsum))
            unfavorable.cumiid2 <- t(apply(A2.iid[,,"unfavorable"],1,cumsum))

            M.cov[,"favorable"] <- M.cov[,"favorable"] + colSums(favorable.cumiid2^2)
            M.cov[,"unfavorable"] <- M.cov[,"unfavorable"] + colSums(unfavorable.cumiid2^2)
            M.cov[,"covariance"] <- M.cov[,"covariance"] + colSums(favorable.cumiid2 * unfavorable.cumiid2)
        }
    }
    rownames(M.cov) <- endpoint
    
    return(M.cov)
}



## * wsumPairScore
## cumulate over endpoint the scores
wsumPairScore <- function(pairScore, weight, n.endpoint){

    keep.col <- c("strata","index.C","index.T","index.pair","indexWithinStrata.C", "indexWithinStrata.T","favorableC","unfavorableC")
    old.col <- c("favorableC","unfavorableC")
    new.col <- c("favorable","unfavorable")

    out <- vector(mode = "list", length = n.endpoint)
    ## indexPair <- stats::setNames(1:NROW(pairScore[[1]]),pairScore[[1]]$index.pair)
    for(iE in 1:n.endpoint){ ## iE <- 2

        iTable <- data.table::copy(pairScore[[iE]][,.SD,.SDcols = keep.col])
        data.table::setnames(iTable, old = old.col, new = new.col)
        iTable[,c("favorable") := .SD$favorable * weight[iE]]
        iTable[,c("unfavorable") := .SD$unfavorable * weight[iE]]
        
        if(iE==1){
            out[[iE]] <- iTable
        }else{
            out[[iE]] <- data.table::copy(out[[iE-1]])
            indexMatch <- match(iTable$index.pair,out[[iE]]$index.pair) ## indexMatch - iTable$index.pair
            out[[iE]][indexMatch, c("favorable") := .SD$favorable + iTable$favorable]
            out[[iE]][indexMatch, c("unfavorable") := .SD$unfavorable + iTable$unfavorable]
        }
    }
    return(out)
}
