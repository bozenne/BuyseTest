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
    old.col <- c("strata","indexWithinStrata.C", "indexWithinStrata.T","favorable.corrected","unfavorable.corrected")
    new.col <- c("strata","index.C", "index.T","favorable","unfavorable")

    ## first endpoint
    ls.table <- vector(mode = "list", length = n.endpoint)
    ls.table[[1]] <- tablePairScore[[1]][,.SD,.SDcols = old.col]
    setnames(ls.table[[1]], old = old.col, new = new.col)
    
    if(n.endpoint>1){

        ##
        n.TCstrata <- tablePairScore[[1]][,length(unique(.SD$index.C)),by = "strata"][[2]]

        for(iE in 2:n.endpoint){ ## iE <- 2
            iTable <- tablePairScore[[iE]][,.SD,.SDcols = old.col]
            setnames(iTable, old = old.col, new = new.col)
            iTable[, indexTable := (index.T-1) * n.TCstrata[.GRP] + index.C, by="strata"]

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
    
    for(iE in 1:n.endpoint){ ## iStrata <- 1
        for(iStrata in 1:n.strata){ ## iStrata <- 1

            ## P[1(X_i,Y_j)1(X_i,Y_k)] = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) sum_k neq j 1(X_i,Y_j)
            ##                         = 1/nm(m-1) sum_i sum_j 1(X_i,Y_j) ( sum_k 1(X_i,Y_k) - 1(X_i,Y_j) )
            ## here we compute sum_k 1(X_i,Y_k) and m-1
            iTable <- ls.table[[iE]][strata == iStrata]
            
            sumPair.T <- iTable[, .(favorable = sum(favorable), unfavorable = sum(unfavorable)), by = "index.T"]
            sumPair.C <- iTable[strata == iStrata, .(favorable = sum(favorable), unfavorable = sum(unfavorable)), by = "index.C"]

            iTable[, c("sumFavorable.T") := sumPair.T$favorable[.SD$index.T]]
            iTable[, c("sumUnfavorable.T") := sumPair.T$unfavorable[.SD$index.T]]
            
            iTable[, c("sumFavorable.C") := sumPair.C$favorable[.SD$index.C]]
            iTable[, c("sumUnfavorable.C") := sumPair.C$unfavorable[.SD$index.C]]

            iN.strata <- NROW(iTable)
            iN.setT <- NROW(sumPair.T) - 1
            iN.setC <- NROW(sumPair.C) - 1

            if(iN.setT > 0){
                ## 1(X_i>Y_j)1(X_i>Y_k)
                strataSum.favorableT[iStrata,iE] <- iTable[,sum((sumFavorable.T  - favorable) * favorable) / (iN.strata*iN.setT)]
                ## 1(X_i<Y_j)1(X_i<Y_k)
                strataSum.unfavorableT[iStrata,iE] <- iTable[,sum((sumUnfavorable.T  - unfavorable) * unfavorable) / (iN.strata*iN.setT)]
                ## 1(X_i>Y_j)1(X_i<Y_k)
                strataSum.mixedT[iStrata,iE] <- iTable[,sum((sumUnfavorable.T  - unfavorable) * favorable) / (iN.strata*iN.setT)]
                ## strataSum.mixedT[iStrata,iE] <- iTable[,sum((sumFavorable.T  - favorable) * unfavorable) / (iN.strata*iN.setT)]
            }
            if(iN.setC > 0){
                ## 1(X_i>Y_j)1(X_k>Y_j)
                strataSum.favorableC[iStrata,iE] <- iTable[,sum((sumFavorable.C  - favorable) * favorable) / (iN.strata*iN.setC)]
                ## 1(X_i<Y_j)1(X_k<Y_j)
                strataSum.unfavorableC[iStrata,iE] <- iTable[,sum((sumUnfavorable.C  - unfavorable) * unfavorable) / (iN.strata*iN.setC)]
                ## 1(X_i>Y_j)1(X_k<Y_j)
                strataSum.mixedC[iStrata,iE] <- iTable[,sum((sumUnfavorable.C  - unfavorable) * favorable) / (iN.strata*iN.setC)]
                ## strataSum.mixedC[iStrata,iE] <- iTable[,sum((sumFavorable.C  - favorable) * unfavorable) / (iN.strata*iN.setC)]
            }
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

    ## ** compute sigma for each endpoint

    ## *** P[X1>Y1 & X1>Y1']
    p1.favorable <- cumsum(count.favorable)/n.pairs
    
    ## *** P[X1<Y1 & X1<Y1']
    p1.unfavorable <- cumsum(count.unfavorable)/n.pairs
    
    ## *** P[X1>Y1 & X1>Y1']
    p2.favorableC <- cumsum(sum.favorableC)/n.setC
    p2.favorableT <- cumsum(sum.favorableT)/n.setT

    ## *** P[X1<Y1 & X1<Y1']
    p2.unfavorableC <- cumsum(sum.unfavorableC)/n.setC
    p2.unfavorableT <- cumsum(sum.unfavorableT)/n.setT

    ## *** P[X1>Y1 & X1<Y1']
    p2.mixedC <- cumsum(sum.mixedC)/n.setC
    p2.mixedT <- cumsum(sum.mixedT)/n.setT

    ## cat("\n")
    ## cat("Treatment: fav.",p2.favorableT," unfav.",p2.unfavorableT," mixed",p2.mixedT,"\n")
    ## cat("Treatment: fav.",p2.favorableC," unfav.",p2.unfavorableC," mixed",p2.mixedC,"\n")
    
    ## *** compute xi
    xi_10_11 <- p2.favorableT - p1.favorable^2
    xi_01_11 <- p2.favorableC - p1.favorable^2
    
    xi_10_22 <- p2.unfavorableT - p1.unfavorable^2
    xi_01_22 <- p2.unfavorableC - p1.unfavorable^2
    
    xi_10_12 <- p2.mixedT - p1.favorable*p1.unfavorable
    xi_01_12 <- p2.mixedC - p1.favorable*p1.unfavorable

    ## ** compute sigma
    n <- n.T
    m <- n.C
    N <- n+m

    ## asymptotic variance i.e. sqrt(n+m)(Uhat - U) \sim N(0,Sigma)
    M.cov <- cbind(favorable = N/m * xi_10_11 + N/n *xi_01_11,
                   unfavorable = N/m * xi_10_22 + N/n *xi_01_22,
                   covariance = N/m * xi_10_12 + N/n * xi_01_12)

    if(M.cov[1,"favorable"]+M.cov[1,"unfavorable"]-2*M.cov[1,"covariance"]<=0){
        warning("Non positive definite covariance matrix")
    }
    ## scaled asymptotic variance i.e. (Uhat - U) \sim N(0,Sigma/(n+m))
    return(M.cov/(n+m))
}
