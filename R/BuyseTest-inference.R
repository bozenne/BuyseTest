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
            iTable <- tablePairScore[[iE]][strata == iStrata, .SD, .SDcols = c("strata","indexWithinStrata.C", "indexWithinStrata.T","favorable.corrected","unfavorable.corrected")]
            data.table::setnames(iTable,
                                 old = c("strata","indexWithinStrata.C", "indexWithinStrata.T","favorable.corrected","unfavorable.corrected"),
                                 new = c("strata","index.C", "index.T","favorable","unfavorable"))

            iN.pairs <- NROW(iTable)
            ls.indexC <- lapply(1:n.C,function(iC){which(iTable$index.C==iC)})
            ls.indexT <- lapply(1:n.T,function(iT){which(iTable$index.T==iT)})
            length.indexC <- sapply(ls.indexC, length) - 1
            length.indexT <- sapply(ls.indexT, length) - 1
            
            iResC <- data.table(matrix(0, nrow = 1, ncol = 4,
                                       dimnames = list(NULL, colname.tempo)))
            iResT <- data.table(matrix(0, nrow = 1, ncol = 4,
                                       dimnames = list(NULL, colname.tempo)))

            for(iPair in 1:iN.pairs){ ## iPair <- 1
                iIndex.C <- iTable[iPair, .SD$index.C]
                iIndex.T <- iTable[iPair, .SD$index.T]
                iFavorable <- iTable[iPair,.SD$favorable]  
                iUnfavorable <- iTable[iPair,.SD$unfavorable]  
              
                strataSum.setC[iStrata,iE] <- length.indexC[iIndex.C]
                strataSum.setT[iStrata,iE] <- length.indexT[iIndex.T]
                                
                if(strataSum.setC[iStrata,iE]>0 && iFavorable + iUnfavorable > 0){
                    iResC <- iResC + iTable[ls.indexC[[iIndex.C]], .(set = .N - 1,
                                                                     favorable = (sum(.SD$favorable)  - iFavorable) * iFavorable,
                                                                     unfavorable = (sum(.SD$unfavorable)  - iUnfavorable) * iUnfavorable,
                                                                     mixed = (sum(.SD$favorable)  - iFavorable) * iUnfavorable)] 
                }else{
                    iResC[, set := set + 1]
                }

                if(strataSum.setT[iStrata,iE]>0 && iFavorable + iUnfavorable > 0){
                    iResT <- iResT + iTable[ls.indexT[[iIndex.T]], .(set = .N - 1,
                                                                     favorable = (sum(.SD$favorable)  - iFavorable) * iFavorable,
                                                                     unfavorable = (sum(.SD$unfavorable)  - iUnfavorable) * iUnfavorable,
                                                                     mixed = (sum(.SD$unfavorable)  - iUnfavorable) * iFavorable)] 

                }else{
                    iResT[, set := set + 1]
                }
                
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
browser()   
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
