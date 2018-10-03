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
    if (cpus == 1) { ## *** sequential permutation test
           
        if (!is.null(seed)) {set.seed(seed)} # set the seed

        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.permutation <- do.call(method.loop,
                                  args = list(X = 1:n.resampling,
                                              FUN = function(iB){
                                                  .BuyseTest(envir = envir,
                                                             keep.pairScore = FALSE,
                                                             method.inference = envir$outArgs$method.inference
                                                             )
                                              })
                                  )
    }else { ## *** parallel permutation test

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
        ls.permutation <- foreach::`%dopar%`(
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
    test.resampling <- which(unlist(lapply(ls.permutation,is.null)) == FALSE)
    if(length(test.resampling) != n.resampling){
        n.failure <- n.resampling - length(test.resampling) 
        warning("The resampling procedure failed for ",n.failure," samples (",round(100*n.failure/n.resampling,2),"%)")
    }
    
    dim.delta <- c(n.strata, D, n.resampling)
    dimnames.delta <- list(level.strata, endpoint, as.character(1:n.resampling))

    out <- list(deltaResampling.netChance = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                deltaResampling.winRatio = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                DeltaResampling.netChance = matrix(NA, nrow = D, ncol = n.resampling,
                                                   dimnames = list(endpoint, as.character(1:n.resampling))),
                DeltaResampling.winRatio = matrix(NA, nrow = D, ncol = n.resampling,
                                                  dimnames = list(endpoint, as.character(1:n.resampling)))
                )

    for(iR in test.resampling){
        out$deltaResampling.netChance[,,iR] <- ls.permutation[[iR]][paste0("delta.",1:n.strata),paste0("netChance.",1:D)]
        out$deltaResampling.winRatio[,,iR] <- ls.permutation[[iR]][paste0("delta.",1:n.strata),paste0("winRatio.",1:D)]

        out$DeltaResampling.netChance[,iR] <- ls.permutation[[iR]][paste0("Delta"),paste0("netChance.",1:D)]
        out$DeltaResampling.winRatio[,iR] <- ls.permutation[[iR]][paste0("Delta"),paste0("winRatio.",1:D)]

    }

    ## ** export
    return(out)
}


## * Inference U-statistic
inferenceUstatistic <- function(tablePairScore, count.favorable, count.unfavorable,
                                n.pairs, n.C, n.T, n.strata, n.endpoint, endpoint){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    n.endpoint <- length(endpoint)

    ## ** merge tables
    keep.col <- c("strata","index.C","index.T","indexWithinStrata.C", "indexWithinStrata.T","favorable.corrected","unfavorable.corrected")
    old.col <- c("favorable.corrected","unfavorable.corrected")
    new.col <- c("favorable","unfavorable")

    ## first endpoint
    ls.table <- vector(mode = "list", length = n.endpoint)
    ls.table[[1]] <- tablePairScore[[1]][,.SD,.SDcols = keep.col]
    setnames(ls.table[[1]], old = old.col, new = new.col)
    
    if(n.endpoint>1){

        ##
        n.TCstrata <- tablePairScore[[1]][,length(unique(.SD$indexWithinStrata.C)),by = "strata"][[2]]

        for(iE in 2:n.endpoint){ ## iE <- 2
            iTable <- tablePairScore[[iE]][,.SD,.SDcols = keep.col]
            setnames(iTable, old = old.col, new = new.col)
            iTable[, indexTable := (indexWithinStrata.T-1) * n.TCstrata[.GRP] + indexWithinStrata.C, by="strata"]

            ls.table[[iE]] <- data.table::copy(ls.table[[iE-1]])
            ls.table[[iE]][iTable$indexTable, favorable := favorable + iTable$favorable]
            ls.table[[iE]][iTable$indexTable, unfavorable := unfavorable + iTable$unfavorable]
        }
    }
    ## ls.table[[1]][, mean(favorable)-mean(unfavorable)]
    ## ls.table[[2]][, mean(favorable)-mean(unfavorable)]


    
    ## ** compute sufficient statistic for each endpoint
    
    ## *** compute sufficient statistics
    ## X1 Y1 X1 Y1' ie 10
    strataSum.favorableT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                  dimnames = list(NULL, endpoint))
    strataSum.unfavorableT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                    dimnames = list(NULL, endpoint))
    strataSum.mixedT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                               dimnames = list(NULL, endpoint))
    strataSum.setT <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                            dimnames = list(NULL, endpoint))

    ## X1 Y1 X1' Y1 ie 01
    strataSum.favorableC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                   dimnames = list(NULL, endpoint))
    strataSum.unfavorableC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                                     dimnames = list(NULL, endpoint))
    strataSum.mixedC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                               dimnames = list(NULL, endpoint))
    strataSum.setC <- matrix(NA, nrow = n.strata, ncol = n.endpoint,
                            dimnames = list(NULL, endpoint))

    colname.tempo <- c("set","favorable","unfavorable","mixed")
    M.iid.favorable <- matrix(NA, nrow = n.T+n.C, ncol = n.endpoint)
    M.iid.unfavorable <- matrix(NA, nrow = n.T+n.C, ncol = n.endpoint)
    
    for(iE in 1:n.endpoint){ ## iStrata <- 1
        for(iStrata in 1:n.strata){ ## iStrata <- 1

            ## P[1(X_i,Y_j)1(X_i,Y_k)] = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) sum_k neq j 1(X_i,Y_j)
            ##                         = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) ( sum_k 1(X_i,Y_k) - 1(X_i,Y_j) )
            ## here we compute sum_k 1(X_i,Y_k) and m-1
            iTable <- ls.table[[iE]][strata == iStrata]
            
            sumPair.T <- iTable[, .(favorable = sum(favorable), unfavorable = sum(unfavorable)), by = "indexWithinStrata.T"]
            sumPair.C <- iTable[strata == iStrata, .(favorable = sum(favorable), unfavorable = sum(unfavorable)), by = "indexWithinStrata.C"]

            iTable[, c("sumFavorable.T") := sumPair.T$favorable[.SD$indexWithinStrata.T]]
            iTable[, c("sumUnfavorable.T") := sumPair.T$unfavorable[.SD$indexWithinStrata.T]]
            
            iTable[, c("sumFavorable.C") := sumPair.C$favorable[.SD$indexWithinStrata.C]]
            iTable[, c("sumUnfavorable.C") := sumPair.C$unfavorable[.SD$indexWithinStrata.C]]

            iN.strata <- NROW(iTable)
            iN.setT <- NROW(sumPair.T) - 1
            iN.setC <- NROW(sumPair.C) - 1

            ## 1/iN.strata =  nm
            if(iN.setT > 0){
                ## E[ 1(X_i>Y_j) 1(X_i>Y_k) ]
                strataSum.favorableT[iStrata,iE] <- iTable[,sum((sumFavorable.T  - favorable) * favorable) / (iN.strata*iN.setT)]
                ## E[ 1(X_i<Y_j) 1(X_i<Y_k) ]
                strataSum.unfavorableT[iStrata,iE] <- iTable[,sum((sumUnfavorable.T  - unfavorable) * unfavorable) / (iN.strata*iN.setT)]
                ## E[ 1(X_i>Y_j) 1(X_i<Y_k) ]
                strataSum.mixedT[iStrata,iE] <- iTable[,sum((sumUnfavorable.T  - unfavorable) * favorable) / (iN.strata*iN.setT)]
                ## strataSum.mixedT[iStrata,iE] <- iTable[,sum((sumFavorable.T  - favorable) * unfavorable) / (iN.strata*iN.setT)]
            }
            if(iN.setC > 0){
                ## E[ 1(X_i>Y_j) 1(X_k>Y_j) ]
                strataSum.favorableC[iStrata,iE] <- iTable[,sum((sumFavorable.C  - favorable) * favorable) / (iN.strata*iN.setC)]
                ## E[ 1(X_i<Y_j) 1(X_k<Y_j) ]
                strataSum.unfavorableC[iStrata,iE] <- iTable[,sum((sumUnfavorable.C  - unfavorable) * unfavorable) / (iN.strata*iN.setC)]
                ## E[ 1(X_i>Y_j) 1(X_k<Y_j) ]
                strataSum.mixedC[iStrata,iE] <- iTable[,sum((sumUnfavorable.C  - unfavorable) * favorable) / (iN.strata*iN.setC)]
                ## strataSum.mixedC[iStrata,iE] <- iTable[,sum((sumFavorable.C  - favorable) * unfavorable) / (iN.strata*iN.setC)]
            }

            ## iid decomposition
            dt.iid.T <- iTable[,.(favorable = sum((sumFavorable.T  - favorable) * favorable) / (iN.strata*iN.setT),
                                  unfavorable = sum((sumUnfavorable.T  - unfavorable) * unfavorable) / (iN.strata*iN.setT)),by = "index.T"]
            dt.iid.C <- iTable[,.(favorable = sum((sumFavorable.C  - favorable) * favorable) / (iN.strata*iN.setC),
                                  unfavorable = sum((sumUnfavorable.C  - unfavorable) * unfavorable) / (iN.strata*iN.setC)),by = "index.C"]

            centering.tempo <- mean(c(dt.iid.T$favorable, dt.iid.C$favorable))
            M.iid.favorable[dt.iid.T$index.T,iE] <- dt.iid.T$favorable - centering.tempo
            M.iid.favorable[dt.iid.C$index.C,iE] <- dt.iid.C$favorable - centering.tempo

            centering.tempo <- mean(c(dt.iid.T$unfavorable, dt.iid.C$unfavorable))
            M.iid.unfavorable[dt.iid.T$index.T,iE] <- dt.iid.T$unfavorable - centering.tempo
            M.iid.unfavorable[dt.iid.C$index.C,iE] <- dt.iid.C$unfavorable - centering.tempo

            ## sum(M.iid.favorable^2)
            ## sum(M.iid.unfavorable^2)
        }
    }    

    sum.favorableT <- colSums(strataSum.favorableT)
    sum.unfavorableT <- colSums(strataSum.unfavorableT)
    sum.mixedT <- colSums(strataSum.mixedT)
    n.setT <- strataSum.setT[1,]

    sum.favorableC <- colSums(strataSum.favorableC)
    sum.unfavorableC <- colSums(strataSum.unfavorableC)
    sum.mixedC <- colSums(strataSum.mixedC)
    n.setC <- strataSum.setC[1,]

    ## ** compute xi
    ## pre-compute
    p1.favorable <- cumsum(count.favorable)/n.pairs ## same as percentage favorable
    p1.unfavorable <- cumsum(count.unfavorable)/n.pairs ## same as percentage unfavorable

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
    n <- n.T
    m <- n.C
    N <- n+m

    ## asymptotic variance i.e. sqrt(n+m)(Uhat - U) \sim N(0,Sigma)
    ## scaled asymptotic variance i.e. (Uhat - U) \sim N(0,Sigma/N)
    M.cov <- cbind(favorable = 1/m * xi_10_11 + 1/n * xi_01_11,
                   unfavorable = 1/m * xi_10_22 + 1/n * xi_01_22,
                   covariance = 1/m * xi_10_12 + 1/n * xi_01_12)
    ## crossprod(cbind(M.iid.favorable,M.iid.unfavorable))

    if((M.cov[1,"favorable"] + M.cov[1,"unfavorable"] - 2 * M.cov[1,"covariance"]) <= 0){
        warning("Non positive definite covariance matrix")
    }

    return(list(Sigma = M.cov,
                iid.favorable = M.iid.favorable,
                iid.unfavorable = M.iid.unfavorable))
}
