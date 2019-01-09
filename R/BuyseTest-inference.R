## * inferenceResampling
inferenceResampling <- function(envir){

    cpus <- envir$outArgs$cpus
    D <- envir$outArgs$D
    endpoint <- envir$outArgs$endpoint
    level.strata <- envir$outArgs$level.strata
    n.resampling <- envir$outArgs$n.resampling
    n.strata <- envir$outArgs$n.strata
    seed <- envir$outArgs$seed
    trace <- envir$outArgs$trace

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
                                                             keep.pairScore = FALSE,
                                                             method.inference = envir$outArgs$method.inference
                                                             )
                                              })
                                  )
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
        toExport <- c(".BuyseTest","initializeSurvival_Peron")

        iB <- NULL ## [:forCRANcheck:] foreach        
        ls.resampling <- foreach::`%dopar%`(
                                       foreach::foreach(iB=1:n.resampling,
                                                        .export = toExport),                                            
                                       {                                           
                                           if(trace>0){utils::setTxtProgressBar(pb, iB)}

                                           return(.BuyseTest(envir = envir,
                                                      keep.pairScore = FALSE,
                                                      method.inference = envir$outArgs$method.inference))
                      
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
    
    dim.delta <- c(n.strata, D, n.resampling)
    dimnames.delta <- list(level.strata, endpoint, as.character(1:n.resampling))

    out <- list(deltaResampling.netBenefit = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                deltaResampling.winRatio = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                DeltaResampling.netBenefit = matrix(NA, nrow = D, ncol = n.resampling,
                                                   dimnames = list(endpoint, as.character(1:n.resampling))),
                DeltaResampling.winRatio = matrix(NA, nrow = D, ncol = n.resampling,
                                                  dimnames = list(endpoint, as.character(1:n.resampling)))
                )

    for(iR in test.resampling){
        out$deltaResampling.netBenefit[,,iR] <- ls.resampling[[iR]]$delta_netBenefit
        out$deltaResampling.winRatio[,,iR] <- ls.resampling[[iR]]$delta_winRatio

        out$DeltaResampling.netBenefit[,iR] <- ls.resampling[[iR]]$Delta_netBenefit
        out$DeltaResampling.winRatio[,iR] <- ls.resampling[[iR]]$Delta_winRatio

    }

    ## ** export
    return(out)
}


## * inference U-statistic
## Implement the computation of the asymptotic variance via an Hajek projection
inferenceUstatistic <- function(tablePairScore, order, count.favorable, count.unfavorable,
                                n.pairs, n.C, n.T, level.strata, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    n.endpoint <- length(endpoint)
    
    ## ** merge tables
    keep.col <- c("strata","index.C","index.T","indexWithinStrata.C", "indexWithinStrata.T","favorableC","unfavorableC")
    old.col <- c("favorableC","unfavorableC")
    new.col <- c("favorable","unfavorable")

    ## first endpoint
    ls.table <- vector(mode = "list", length = n.endpoint)
    ls.table[[1]] <- tablePairScore[[1]][,.SD,.SDcols = keep.col]
    setnames(ls.table[[1]], old = old.col, new = new.col)
    
    if(n.endpoint>1){
        ##
        n.TCstrata <- tablePairScore[[1]][,length(unique(.SD$indexWithinStrata.C)), by = "strata"][[2]]

        for(iE in 2:n.endpoint){ ## iE <- 2
            iTable <- tablePairScore[[iE]][,.SD,.SDcols = keep.col]
            setnames(iTable, old = old.col, new = new.col)
            iTable[, c("indexTable") := (.SD$indexWithinStrata.T-1) * n.TCstrata[.GRP] + .SD$indexWithinStrata.C, by="strata"]

            ls.table[[iE]] <- data.table::copy(ls.table[[iE-1]])
            ls.table[[iE]][,c("favorable","unfavorable") := 0]
            ls.table[[iE]][iTable$indexTable, c("favorable") := iTable$favorable]
            ls.table[[iE]][iTable$indexTable, c("unfavorable") := iTable$unfavorable]
        }
    }
    
    ## ls.table[[1]][, mean(favorable)-mean(unfavorable)]
    ## ls.table[[2]][, mean(favorable)-mean(unfavorable)]
    ## ** H-decomposition
    ## expectation
    Upartial.favorable <- count.favorable/n.pairs
    Upartial.unfavorable <- count.unfavorable/n.pairs

    ## storage
    A.iid <- array(NA, dim = c(n.T+n.C, n.endpoint, 2), dimnames = list(NULL, endpoint, c("favorable","unfavorable")))
    if(order == 2){
        A2.iid <- array(NA, dim = c(n.pairs, n.endpoint, 2), dimnames = list(NULL, endpoint, c("favorable","unfavorable")))
    }else{
        A2.iid <- NULL
    }

    ## loop
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        for(iE in 1:n.endpoint){ ## iE <- 1

            ## extract pairwise scores
            iTable <- ls.table[[iE]][ls.table[[iE]]$strata == level.strata[iStrata]]
            index2originalOrder.C <- iTable[!duplicated(iTable$index.C),setNames(.SD$index.C,.SD$indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(iTable$index.T),setNames(.SD$index.T,.SD$indexWithinStrata.T)]
            iN.strata <- NROW(iTable) ## number of pairs

            ## *** Hajek projection
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.T"]
            sumPair.T[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.T[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
                
            ## \E[X_i>=Y_j+\tau|Y_j] and \E[X_i+\tau<=Y_j|Y_j]
            sumPair.C <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.C"]
            sumPair.C[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.C[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]

            ## store
            A.iid[index2originalOrder.C,iE,"favorable"] <- (sumPair.C$E.favorable - Upartial.favorable[iE]) / n.C
            A.iid[index2originalOrder.T,iE,"favorable"] <- (sumPair.T$E.favorable - Upartial.favorable[iE]) / n.T
            A.iid[index2originalOrder.C,iE,"unfavorable"] <- (sumPair.C$E.unfavorable - Upartial.unfavorable[iE]) / n.C
            A.iid[index2originalOrder.T,iE,"unfavorable"] <- (sumPair.T$E.unfavorable - Upartial.unfavorable[iE]) / n.T
            
            ## *** second order
            if(order == 2){
                A2.iid[,iE,"favorable"] <- (iTable$favorable - A.iid[iTable$index.C,iE,"favorable"] * n.C - A.iid[iTable$index.T,iE,"favorable"] * n.T - Upartial.favorable[iE])/n.pairs
                A2.iid[,iE,"unfavorable"] <- (iTable$unfavorable - A.iid[iTable$index.C,iE,"unfavorable"] * n.C - A.iid[iTable$index.T,iE,"unfavorable"] * n.T - Upartial.unfavorable[iE])/n.pairs                
            }
        }
    }

    ## ** compute Sigma
    M.cov <- .iid2cov(A.iid = A.iid, A2.iid = A2.iid,
                      order = order, endpoint = endpoint, n.endpoint = n.endpoint)
    
    ## ** export
    return(list(Sigma = M.cov,
                iid1 = A.iid,
                iid2 = A2.iid))
}

## * inference U-statistic (Bebu et al 2015)
## Implement the computation of the asymptotic variance as described in
## Large sample inference for a win ratio analysis of a composite outcome based on prioritized components
## Biostatistics (2015), pp. 1–10 doi:10.1093/biostatistics/kxv032
## Give results equivalent to inferenceUstatistic
inferenceUstatisticBebu <- function(tablePairScore, order, count.favorable, count.unfavorable,
                                    n.pairs, n.C, n.T, level.strata, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    n.endpoint <- length(endpoint)
    ntot.pairs <- sum(n.pairs)
        
    p1.favorable <- cumsum(count.favorable)/ntot.pairs
    p1.unfavorable <- cumsum(count.unfavorable)/ntot.pairs

    ## ** merge tables
    keep.col <- c("strata","index.C","index.T","indexWithinStrata.C", "indexWithinStrata.T","favorableC","unfavorableC")
    old.col <- c("favorableC","unfavorableC")
    new.col <- c("favorable","unfavorable")

    ## first endpoint
    ls.table <- vector(mode = "list", length = n.endpoint)
    ls.table[[1]] <- tablePairScore[[1]][,.SD,.SDcols = keep.col]
    setnames(ls.table[[1]], old = old.col, new = new.col)
    
    if(n.endpoint>1){
        ##
        n.TCstrata <- tablePairScore[[1]][,length(unique(.SD$indexWithinStrata.C)), by = "strata"][[2]]

        for(iE in 2:n.endpoint){ ## iE <- 2
            iTable <- tablePairScore[[iE]][,.SD,.SDcols = keep.col]
            setnames(iTable, old = old.col, new = new.col)
            iTable[, c("indexTable") := (.SD$indexWithinStrata.T-1) * n.TCstrata[.GRP] + .SD$indexWithinStrata.C, by="strata"]

            ls.table[[iE]] <- data.table::copy(ls.table[[iE-1]])
            ls.table[[iE]][iTable$indexTable, c("favorable") := .SD$favorable + iTable$favorable]
            ls.table[[iE]][iTable$indexTable, c("unfavorable") := .SD$unfavorable + iTable$unfavorable]
        }
    }
    ## ls.table[[1]][, mean(favorable)-mean(unfavorable)]
    ## ls.table[[2]][, mean(favorable)-mean(unfavorable)]


    ## ** compute variance component over strata
    M.cov <- matrix(0, nrow = n.endpoint, ncol = 3,
                    dimnames = list(endpoint, c("favorable","unfavorable","covariance")))
    
    strataSum <- matrix(NA, nrow = 6, ncol = n.endpoint,
                        dimnames = list(c("favorableT","favorableC","unfavorableT","unfavorableC","mixedC","mixedT"),
                                        endpoint))
    
    for(iStrata in 1:n.strata){ ## iStrata <- 1

        iN.strata <- n.pairs[iStrata]
        iDT.nCT <- ls.table[[1]][ls.table[[1]]$strata == level.strata[iStrata], .(n.C = length(unique(.SD$index.C)),n.T = length(unique(.SD$index.T)))]
        iN.C <- iDT.nCT$n.C
        iN.T <- iDT.nCT$n.T
        
        for(iE in 1:n.endpoint){ ## iE <- 1

            iTable <- ls.table[[iE]][ls.table[[iE]]$strata == level.strata[iStrata]]
            index2originalOrder.C <- iTable[!duplicated(iTable$index.C),setNames(.SD$index.C,.SD$indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(iTable$index.T),setNames(.SD$index.T,.SD$indexWithinStrata.T)]

            ## *** Hajek projection
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.T"]
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
            H2.favorable <- (p1.favorable*(1-p1.favorable) - xi_10_11 - xi_01_11)/(iN.strata)
            H2.unfavorable <- (p1.unfavorable*(1-p1.unfavorable) - xi_10_22 - xi_01_22)/(iN.strata)
            H2.covariance <- 0##(-xi_10_12)/(iN.strata)
            warning("No formula for the contribution of the second order term to the covariance. It is set to 0.")
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
    return(list(Sigma = M.cov,
                iid1 = NULL,
                iid2 = NULL))
}

## * .iid2cov
.iid2cov <- function(A.iid, A2.iid,
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
