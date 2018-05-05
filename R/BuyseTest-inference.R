## * inferenceResampling
inferenceResampling <- function(envir){

    cpus <- envir$cpus
    D <- envir$D
    endpoint <- envir$endpoint
    level.strata <- envir$level.strata
    n.resampling <- envir$n.resampling
    n.strata <- envir$n.strata
    seed <- envir$seed
    trace <- envir$trace
    
    ## ** computation
    if (cpus == 1) { ## *** sequential permutation test
           
            if (!is.null(seed)) {set.seed(seed)} # set the seed

        if(trace>0){
            if (trace > 1) {cat("Sequential permutation test \n")}
            requireNamespace("pbapply")
            ls.permutation <- pbapply::pblapply(1:n.resampling, function(iB){
                return(warperResampling(iB, envir = envir))
            })
        }else{
            ls.permutation <- lapply(1:n.resampling, function(iB){
                return(warperResampling(iB, envir = envir))
            })      
        }
            
    }else { ## *** parallel permutation test
            
        if (trace > 1) {cat("Parallel permutation test \n")}


        cl <- snow::makeSOCKcluster(cpus)
        doSNOW::registerDoSNOW(cl)
        
        if(trace){
            pb <- utils::txtProgressBar(min=1, max=n.resampling, style=3)
            ls.options <- list(progress = function(n){ utils::setTxtProgressBar(pb, n) })
        }else{
            ls.options <- NULL
        }

        toExport <- c("warperResampling","initializeSurvival_Peto","initializeSurvival_Peron")

        iB <- NULL ## [:forCRANcheck:] foreach        
        ls.permutation <- foreach::`%dopar%`(
                                       foreach::foreach(iB=1:n.resampling,
                                                        .packages = "BuyseTest",
                                                        .options.snow = ls.options,
                                                        .export = toExport),                                            
                                       {
                                           return(warperResampling(iB, envir = envir))
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

    out <- list(deltaResampling_netChance = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                deltaResampling_winRatio = array(NA, dim = dim.delta, dimnames = dimnames.delta),
                DeltaResampling_netChance = matrix(NA, nrow = D, ncol = n.resampling,
                                         dimnames = list(endpoint, as.character(1:n.resampling))),
                DeltaResampling_winRatio = matrix(NA, nrow = D, ncol = n.resampling,
                                        dimnames = list(endpoint, as.character(1:n.resampling)))
                )

        for(iR in test.resampling){
            
            out$deltaResampling_netChance[,,iR] <- ls.permutation[[iR]][paste0("delta.",1:n.strata),paste0("netChance.",1:D)]
            out$deltaResampling_winRatio[,,iR] <- ls.permutation[[iR]][paste0("delta.",1:n.strata),paste0("winRatio.",1:D)]

            out$DeltaResampling_netChance[,iR] <- ls.permutation[[iR]][paste0("Delta"),paste0("netChance.",1:D)]
            out$DeltaResampling_winRatio[,iR] <- ls.permutation[[iR]][paste0("Delta"),paste0("winRatio.",1:D)]

        }

    ## ** export
    return(out)
}
## * warperResampling
warperResampling <- function(x, envir){

    ## ** prepare
    Mnew.Treatment <- NULL
    Mnew.Control <- NULL
    Mnew.delta.Treatment <- NULL
    Mnew.delta.Control <- NULL    
    new.strataT <- list()
    new.strataC <- list()

  
    ## ** Data management
    new.nT <- 0
    new.nC <- 0
    if(envir$method.inference %in% c("bootstrap","permutation")){
        if(envir$method.inference == "bootstrap"){
            indexT <- sample.int(envir$n.T+envir$n.C,
                                 size = envir$n.T, 
                                 replace = TRUE)
        
            indexC <- sample.int(envir$n.T+envir$n.C,
                                 size = envir$n.C, 
                                 replace = TRUE)
        
        }else if(envir$method.inference == "permutation"){
            indexT <- sample.int(envir$n.T+envir$n.C,
                                 size = envir$n.T, 
                                 replace = FALSE)

            indexC <- setdiff(1:(envir$n.T+envir$n.C), indexT)
        }

        indexT.T <- intersect(indexT,1:envir$n.T)
        indexT.C <- setdiff(indexT,indexT.T) - envir$n.T
        indexC.T <- intersect(indexC,1:envir$n.T)
        indexC.C <- setdiff(indexC,indexC.T) - envir$n.T
     
    }

    for (iterS in 1:envir$n.strata) {  ## randomisation : new allocation of the treatment and control arm

        if(envir$method.inference %in% c("bootstrap","permutation")){
            ## +1 because index.strataX starts at 0 since it is sent to C++ in the punctual estimation
            iIndexT.T <- intersect(indexT.T, envir$index.strataT[[iterS]] + 1)
            iIndexT.C <- intersect(indexT.C, envir$index.strataC[[iterS]] + 1) 
            iIndexC.T <- intersect(indexC.T, envir$index.strataT[[iterS]] + 1)
            iIndexC.C <- intersect(indexC.C, envir$index.strataC[[iterS]] + 1)
            ## sort(c(iIndexT, iIndexC))
            
        }else{
            
            if(envir$method.inference == "stratified bootstrap"){
                iIndexT <- sample.int(envir$n.eachStrataT[iterS]+envir$n.eachStrataC[iterS],
                                      size = envir$n.eachStrataT[iterS], 
                                      replace = TRUE)
                iIndexC <- sample.int(envir$n.eachStrataT[iterS]+envir$n.eachStrataC[iterS],
                                      size = envir$n.eachStrataC[iterS], 
                                      replace = TRUE)
            }else if(envir$method.inference == "stratified permutation"){
                iIndexT <- sample.int(envir$n.eachStrataT[iterS]+envir$n.eachStrataC[iterS],
                                      size = envir$n.eachStrataT[iterS], 
                                      replace = FALSE)
                iIndexC <- setdiff(1:(envir$n.eachStrataT[iterS]+envir$n.eachStrataC[iterS]), iIndexT)
            }

            iIndexT.T <- iIndexT[iIndexT<=envir$n.eachStrataT[iterS]]
            iIndexT.C <- setdiff(iIndexT, iIndexT.T) - envir$n.eachStrataT[iterS]
            
            iIndexC.T <- iIndexC[iIndexC<=envir$n.eachStrataT[iterS]] 
            iIndexC.C <- setdiff(iIndexC, iIndexC.T) - envir$n.eachStrataT[iterS]

            ## convert back to global index (instead of within strata index)
            ## +1 because index.strataX starts at 0 since it is sent to C++ in the punctual estimation
            iIndexT.T <- envir$index.strataT[[iterS]][iIndexT.T] + 1
            iIndexT.C <- envir$index.strataC[[iterS]][iIndexT.C] + 1
            iIndexC.T <- envir$index.strataT[[iterS]][iIndexC.T] + 1
            iIndexC.C <- envir$index.strataC[[iterS]][iIndexC.C] + 1
        }

        ## sort(c(iIndexT.T,iIndexC.T))
        ## sort(c(iIndexT.C,iIndexC.C))
        nResampling.T <- length(iIndexT.T) + length(iIndexT.C) 
        nResampling.C <- length(iIndexC.T) + length(iIndexC.C)
        
        ## *** test whether there are observations in each group
        if (nResampling.T == 0 || nResampling.C == 0) {
            permutation.failure <- TRUE ; break ;
        }else{
            permutation.failure <- FALSE
        }
    
        ## *** update dataset
        Mnew.Treatment <- rbind(Mnew.Treatment, 
                                envir$M.Treatment[iIndexT.T,,drop = FALSE], 
                                envir$M.Control[iIndexT.C,,drop = FALSE])
        Mnew.Control <- rbind(Mnew.Control, 
                              envir$M.Treatment[iIndexC.T,,drop = FALSE],
                              envir$M.Control[iIndexC.C,,drop = FALSE]) 
    
        ## *** update strata - the new strata iIndex minus 1 (for C++ compatibility, vector begin at 0)
        new.strataT[[iterS]] <- seq(from = new.nT, by = 1, length = nResampling.T) 
        new.strataC[[iterS]] <- seq(from = new.nC, by = 1, length = nResampling.C)
    
        ## *** update censoring variables
        if (envir$D.TTE > 0) { 
            Mnew.delta.Treatment <- rbind(Mnew.delta.Treatment, 
                                          envir$M.delta.Treatment[iIndexT.T,, drop = FALSE], 
                                          envir$M.delta.Control[iIndexT.C,, drop = FALSE])
            Mnew.delta.Control <- rbind(Mnew.delta.Control, 
                                        envir$M.delta.Treatment[iIndexC.T,, drop = FALSE],
                                        envir$M.delta.Control[iIndexC.C,, drop = FALSE])
      
            if (envir$method.tte == 2) { # "Efron": set last event to non-censored
                Mnewstrata.Treatment <- Mnew.Treatment[new.strataT[[iterS]] + 1, which(envir$type == 3), drop = FALSE]
                Mnewstrata.Control <- Mnew.Control[new.strataC[[iterS]] + 1, which(envir$type == 3), drop = FALSE]
        
                for (iEndpoint.TTE in 1:envir$D.TTE) {
                    iIndexT_maxCensored <- which.max(Mnewstrata.Treatment[,iEndpoint.TTE])
                    Mnew.delta.Treatment[new.strataT[[iterS]][iIndexT_maxCensored] + 1,iEndpoint.TTE] <- 1
                    iIndexC_maxCensored <- which.max(Mnewstrata.Control[,iEndpoint.TTE])
                    Mnew.delta.Control[new.strataC[[iterS]][iIndexC_maxCensored] + 1,iEndpoint.TTE] <- 1
                }
            }
        }else {
            Mnew.delta.Treatment <- matrix(nrow = 0, ncol = 0)
            Mnew.delta.Control <- matrix(nrow = 0, ncol = 0)
        }
    
        new.nT <- new.nT + nResampling.T
        new.nC <- new.nC + nResampling.C
    }
  
    if (permutation.failure) {
        ##        return(matrix(NA, nrow = envir$n.strata + 1, ncol = 2*envir$D))
        return(NULL)
    }
    
    ## ** Update survival
    if(envir$method.tte == 0){ ## Gehan
        new.survivalT <- lapply(1:envir$D.TTE, matrix)
        new.survivalC <- lapply(1:envir$D.TTE, matrix)
    }else if(envir$method.tte == 1){ ## Peto
        outSurv <- initializeSurvival_Peto(M.Treatment=Mnew.Treatment,
                                           M.Control=Mnew.Control,
                                           M.delta.Treatment=Mnew.delta.Treatment,
                                           M.delta.Control=Mnew.delta.Control,
                                           endpoint=envir$endpoint,
                                           D.TTE=envir$D.TTE,
                                           type=envir$type,
                                           threshold=envir$threshold,
                                           index.strataT=new.strataT,
                                           index.strataC=new.strataC,
                                           n.strata=envir$n.strata)
        new.survivalT <- outSurv$list.survivalT
        new.survivalC <- outSurv$list.survivalC
    }else if(envir$method.tte %in% 2:3){
        outSurv <- initializeSurvival_Peron(M.Treatment=Mnew.Treatment,
                                            M.Control=Mnew.Control,
                                            M.delta.Treatment=Mnew.delta.Treatment,
                                            M.delta.Control=Mnew.delta.Control,
                                            endpoint=envir$endpoint,
                                            D.TTE=envir$D.TTE,
                                            type=envir$type,
                                            threshold=envir$threshold,
                                            index.strataT=new.strataT,
                                            index.strataC=new.strataC,
                                            n.strata=envir$n.strata)
        new.survivalT <- outSurv$list.survivalT
        new.survivalC <- outSurv$list.survivalC
    }
    
    ## ** Computation
    resBT <-   GPC_cpp(Treatment = Mnew.Treatment,
                       Control = Mnew.Control,
                       threshold = envir$threshold,
                       survEndpoint = (envir$type == 3),
                       delta_Treatment = Mnew.delta.Treatment,
                       delta_Control = Mnew.delta.Control,
                       D = envir$D,
                       returnIndex = FALSE,
                       strataT = new.strataT,
                       strataC = new.strataC,
                       n_strata = envir$n.strata,
                       n_TTE = envir$D.TTE,
                       Wscheme = envir$Wscheme,
                       index_survivalM1 = envir$index.survivalM1,
                       threshold_TTEM1 = envir$threshold.TTEM1,
                       list_survivalT = new.survivalT,
                       list_survivalC = new.survivalC,
                       methodTTE = envir$method.tte,
                       correctionTTE = envir$correctionTTE,
                       neutralAsUninf = envir$neutral.as.uninf,
                       keepComparison = FALSE
                       )

    ## ** export
    Mout <- cbind(rbind(resBT$delta_netChance, resBT$Delta_netChance),
                  rbind(resBT$delta_winRatio, resBT$Delta_winRatio))
    dimnames(Mout) <- list(c(paste0("delta.",1:envir$n.strata),"Delta"),
                           c(paste0("netChance.",1:envir$D),paste0("winRatio.",1:envir$D))
                           )
    return(Mout)
}







