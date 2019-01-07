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
inferenceUstatistic <- function(tablePairScore, count.favorable, count.unfavorable,
                                n.pairs, n.C, n.T, n.strata, n.endpoint, endpoint){
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
            ls.table[[iE]][iTable$indexTable, c("favorable") := .SD$favorable + iTable$favorable]
            ls.table[[iE]][iTable$indexTable, c("unfavorable") := .SD$unfavorable + iTable$unfavorable]
        }
    }
    ## ls.table[[1]][, mean(favorable)-mean(unfavorable)]
    ## ls.table[[2]][, mean(favorable)-mean(unfavorable)]
    
    ## ** compute sufficient statistics
    ## expectation
    p1.favorable <- cumsum(count.favorable)/n.pairs ## same as percentage favorable
    p1.unfavorable <- cumsum(count.unfavorable)/n.pairs ## same as percentage unfavorable

    ## variance T 
    strataSum.favorableT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                   dimnames = list(NULL, endpoint))
    strataSum.unfavorableT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                     dimnames = list(NULL, endpoint))
    strataSum.mixedT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                               dimnames = list(NULL, endpoint))

    ## variance C
    strataSum.favorableC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                   dimnames = list(NULL, endpoint))
    strataSum.unfavorableC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                     dimnames = list(NULL, endpoint))
    strataSum.mixedC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                               dimnames = list(NULL, endpoint))

    ## iid decomposition
    A.iid <- array(NA, dim = c(n.T+n.C, n.endpoint, 2), dimnames = list(NULL, endpoint, c("favorable","unfavorable")))

    ## loop
    for(iStrata in 1:n.strata){ ## iStrata <- 1
        for(iE in 1:n.endpoint){ ## iE <- 1

            iTable <- ls.table[[iE]][ls.table[[iE]]$strata == iStrata]
            index2originalOrder.C <- iTable[!duplicated(index.C),setNames(index.C,indexWithinStrata.C)]
            index2originalOrder.T <- iTable[!duplicated(index.T),setNames(index.T,indexWithinStrata.T)]
            iN.strata <- NROW(iTable) ## number of pairs

            ## ***  Hajek projection (iid decomposition)
            ## \E[X_i>=Y_j+\tau|X_i] and \E[X_i+\tau<=Y_j|X_i]
            sumPair.T <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.T"]
            sumPair.T[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.T[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]
                
            ## \E[X_i>=Y_j+\tau|Y_j] and \E[X_i+\tau<=Y_j|Y_j]
            sumPair.C <- iTable[, .(pairs  = .N, favorable = sum(.SD$favorable), unfavorable = sum(.SD$unfavorable)), by = "indexWithinStrata.C"]
            sumPair.C[, c("E.favorable") := .SD$favorable/.SD$pairs]
            sumPair.C[, c("E.unfavorable") := .SD$unfavorable/.SD$pairs]

            ## store
            A.iid[index2originalOrder.C,iE,"favorable"] <- (sumPair.C$E.favorable - (count.favorable[iE]/n.pairs)) / n.T
            A.iid[index2originalOrder.T,iE,"favorable"] <- (sumPair.T$E.favorable - (count.favorable[iE]/n.pairs)) / n.C
            A.iid[index2originalOrder.C,iE,"unfavorable"] <- (sumPair.C$E.unfavorable - (count.unfavorable[iE]/n.pairs)) / n.T
            A.iid[index2originalOrder.T,iE,"unfavorable"] <- (sumPair.T$E.unfavorable - (count.unfavorable[iE]/n.pairs)) / n.C

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
            
            if(NROW(sumPair.T) > 0){
                ## E[ 1(X_i>Y_j) 1(X_i>Y_k) ]
                strataSum.favorableT[iStrata,iE] <- iTable[,sum(.SD$sumFavorable.T * .SD$favorable)] / (iN.strata*iN.setT)
                ## E[ 1(X_i<Y_j) 1(X_i<Y_k) ]
                strataSum.unfavorableT[iStrata,iE] <- iTable[,sum(.SD$sumUnfavorable.T * .SD$unfavorable)] / (iN.strata*iN.setT)
                ## E[ 1(X_i>Y_j) 1(X_i<Y_k) ]
                strataSum.mixedT[iStrata,iE] <- iTable[,sum(.SD$sumUnfavorable.T * .SD$favorable)] / (iN.strata*iN.setT)
            }
            if(NCOL(sumPair.C) > 0){
                ## E[ 1(X_i>Y_j) 1(X_k>Y_j) ]
                strataSum.favorableC[iStrata,iE] <- iTable[,sum(.SD$sumFavorable.C * .SD$favorable)] / (iN.strata*iN.setC)
                ## E[ 1(X_i<Y_j) 1(X_k<Y_j) ]
                strataSum.unfavorableC[iStrata,iE] <- iTable[,sum(.SD$sumUnfavorable.C * .SD$unfavorable)] / (iN.strata*iN.setC)
                ## E[ 1(X_i>Y_j) 1(X_k<Y_j) ]
                strataSum.mixedC[iStrata,iE] <- iTable[,sum(.SD$sumUnfavorable.C * .SD$favorable)] / (iN.strata*iN.setC)
            }

        }
    }

    ## ** average over strata
    sum.favorableT <- colSums(strataSum.favorableT)
    sum.unfavorableT <- colSums(strataSum.unfavorableT)
    sum.mixedT <- colSums(strataSum.mixedT)

    sum.favorableC <- colSums(strataSum.favorableC)
    sum.unfavorableC <- colSums(strataSum.unfavorableC)
    sum.mixedC <- colSums(strataSum.mixedC)

    ## ** compute xi
    ## P[X1>Y1 & X1>Y1'] - P[X1>Y1]^2
    xi_10_11 <- cumsum(sum.favorableT) - p1.favorable^2
    ## P[X1>Y1 & X1'>Y1] - P[X1>Y1]^2
    xi_01_11 <- cumsum(sum.favorableC) - p1.favorable^2
    
    ## P[X1<Y1 & X1<Y1'] - P[X1<Y1]^2
    xi_10_22 <- cumsum(sum.unfavorableT) - p1.unfavorable^2
    ## P[X1<Y1 & X1'<Y1] - P[X1<Y1]^2
    xi_01_22 <- cumsum(sum.unfavorableC) - p1.unfavorable^2
    
    ## P[X1>Y1 & X1<Y1'] - P[X1>Y1]*P[X1<Y1]
    xi_10_12 <- cumsum(sum.mixedT) - p1.favorable * p1.unfavorable
    ## P[X1>Y1 & X1'<Y1] - P[X1>Y1]*P[X1<Y1]
    xi_01_12 <- cumsum(sum.mixedC) - p1.favorable * p1.unfavorable

    ## ** compute sigma
    ## N.TC = N.T+N.C
    ## asymptotic variance i.e. sqrt(N.TC)(Uhat - U) \sim N(0,Sigma)
    ## scaled asymptotic variance i.e. (Uhat - U) \sim N(0,Sigma/N.TC)
    M.cov <- cbind(favorable = 1/n.C * xi_10_11 + 1/n.T * xi_01_11,
                   unfavorable = 1/n.C * xi_10_22 + 1/n.T * xi_01_22,
                   covariance = 1/n.C * xi_10_12 + 1/n.T * xi_01_12)

    if(any((M.cov[,"favorable"] + M.cov[,"unfavorable"] - 2 * M.cov[,"covariance"]) <= 0)){
        warning("Non positive definite covariance matrix")
    }

    ## ** export
    return(list(Sigma = M.cov,
                iid = A.iid))
}
