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
                                n.pairs, n.C, n.T,                                
                                correction.uninf){
    . <- NULL ## for CRAN test
    
    ## ** extract informations
    endpoint <- names(tablePairScore)
    D <- length(endpoint)
    col.id <- names(tablePairScore[[1]])[1:3]
    keep.col <- c(col.id,"favorable","unfavorable","neutral","uninformative")
    
    ## ** fct
    myFct_T <- function(favorable,unfavorable){
        iN.pairs <- length(favorable)
        if(iN.pairs>1){
            iCombination <- utils::combn(1:iN.pairs, 2)
            return(c(NCOL(iCombination),
                     sum( favorable[iCombination[1,]] * favorable[iCombination[2,]] ),  ## > >
                     sum( unfavorable[iCombination[1,]] * unfavorable[iCombination[2,]] ), ## < <
                     sum( favorable[iCombination[1,]] * unfavorable[iCombination[2,]] ) ## > <
                     ))
        }else{
            return(c(0,
                     0,  ## > >
                     0, ## < <
                     0 ## > <
                     ))
        }
    }

    myFct_C <- function(favorable,unfavorable){
        iN.pairs <- length(favorable)
        if(iN.pairs>1){
            iCombination <- utils::combn(1:iN.pairs, 2)
            return(c(NCOL(iCombination),
                     sum( favorable[iCombination[1,]] * favorable[iCombination[2,]] ),  ## > >
                     sum( unfavorable[iCombination[1,]] * unfavorable[iCombination[2,]] ), ## < <
                     sum( unfavorable[iCombination[1,]] * favorable[iCombination[2,]] ) ## > <
                     ))
        }else{
            return(c(0,
                     0,  ## > >
                     0, ## < <
                     0 ## > <
                     ))
        }
    }
    
    ## ** compute sufficient statistic for each endpoint
    M.sufficient <- matrix(NA, nrow = D, ncol = 8,
                           dimnames = list(endpoint, c("n.pairT","sum.favorableT","sum.unfavorableT","sum.mixedT",
                                                       "n.pairC","sum.favorableC","sum.unfavorableC","sum.mixedC")
                                           ))
    
    for(iE in 1:D){ ## iE <- 1

        iTable <- data.table::copy(tablePairScore[[iE]][,.SD,.SDcols = keep.col])
        data.table::setnames(iTable, old = col.id, new = c("strata","indexT","indexC"))

        ## *** perform correction
        if(correction.uninf){
            stop("Inference U statistic in presence of correction.uninf needs to be updated!!!")
            ## vec.tempo <- unlist(iTable[,.(favorable = sum(favorable),
            ##                               unfavorable = sum(unfavorable),
            ##                               neutral = sum(neutral),
            ##                               uninformative = sum(uninformative))])
            ## factor <- sum(vec.tempo)/sum(vec.tempo[1:3])
            
            ## iTable[, favorable := favorable * factor]
            ## iTable[, unfavorable := unfavorable * factor]
            ## iTable[, neutral := neutral * factor]
            ## iTable[, uninformative := 0]
        }

        ## *** compute sufficient statistics
        ls.count_T <- iTable[,.(list = .(myFct_T(.SD$favorable, .SD$unfavorable))), by = c("strata","indexT") ][["list"]]
        M.sufficient[iE,1:4] <- colSums(do.call(rbind,ls.count_T))
    
        ls.count_C <- iTable[,.(list = .(myFct_C(.SD$favorable, .SD$unfavorable))), by = c("strata","indexC") ][["list"]]
        M.sufficient[iE,5:8] <- colSums(do.call(rbind,ls.count_C))
        
    }

    ## ** compute sigma for each endpoint

    ## *** P[X1>Y1 & X1>Y1']
    p1.favorable <- cumsum(count.favorable)/n.pairs
    
    ## *** P[X1<Y1 & X1<Y1']
    p1.unfavorable <- cumsum(count.unfavorable)/n.pairs
    
    ## *** P[X1>Y1 & X1>Y1']
    p2.favorableT <- cumsum(M.sufficient[,"sum.favorableT"])/cumsum(M.sufficient[,"n.pairT"])
    p2.favorableC <- cumsum(M.sufficient[,"sum.favorableC"])/cumsum(M.sufficient[,"n.pairC"])

    ## *** P[X1<Y1 & X1<Y1']
    p2.unfavorableT <- cumsum(M.sufficient[,"sum.unfavorableT"])/cumsum(M.sufficient[,"n.pairT"])
    p2.unfavorableC <- cumsum(M.sufficient[,"sum.unfavorableC"])/cumsum(M.sufficient[,"n.pairC"])

    ## *** P[X1>Y1 & X1<Y1']
    p2.mixedT <- cumsum(M.sufficient[,"sum.mixedT"])/cumsum(M.sufficient[,"n.pairT"])
    p2.mixedC <- cumsum(M.sufficient[,"sum.mixedC"])/cumsum(M.sufficient[,"n.pairC"])

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
