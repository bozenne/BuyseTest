## * inferenceResampling
inferenceResampling <- function(envir){

    cpus <- envir$outArgs$cpus
    D <- envir$outArgs$D
    endpoint <- envir$outArgs$endpoint
    level.strata <- envir$outArgs$level.strata
    n.resampling <- envir$outArgs$n.resampling
    n.strata <- envir$outArgs$n.strata
    seed <- envir$outArgs$seed
    time <- envir$time
    trace <- envir$outArgs$trace

    ## ** display
    if (trace > 1) {
        if (time[3] == 0) {
            time.punctual <- "<0.001 s"
            time.permutation <- paste0("<",signif(0.001*n.resampling/cpus,4)," s")
        }else{
            time.punctual <- paste(time[3],"s")
            time.permutation <- paste(signif(time[3]*n.resampling/cpus,4),"s")
        }
        txt.type <- switch(envir$outArgs$method.inference,
                           "bootstrap" = "bootstrap resampling",
                           "stratified bootstrap" = "stratified bootstrap resampling",
                           "permutation" = "permutation test",
                           "stratified permutation" = "stratified permutation test")
            
            cat("Settings (",txt.type,"): \n",
                "   > requested time for one sample: ", time.punctual, "\n",
                "   > estimated time for ", n.resampling, " samples with ", cpus, " core", if (cpus > 1) {"s"}, ": ", time.permutation, "\n", sep = "")
            if (!is.null(seed)) {
                cat("   > seed", if (cpus > 1) {"s"}, ": ",paste(seq(seed,seed + cpus - 1), collapse = " "), sep = "")       
            }
            cat("\n")         
        }
    
    ## ** computation
    if (cpus == 1) { ## *** sequential permutation test
           
        if (!is.null(seed)) {set.seed(seed)} # set the seed

        if (trace > 0) {
            if(trace>1){
                cat("Sequential permutation test \n")
            }
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.permutation <- do.call(method.loop,
                                  args = list(X = 1:n.resampling,
                                              FUN = function(iB){
                                                  .BuyseTest(data = envir$outArgs$data,
                                                             censoring = envir$outArgs$censoring,
                                                             correction.tte = envir$outArgs$correction.tte,
                                                             D = envir$outArgs$D,
                                                             D.TTE = envir$outArgs$D.TTE,
                                                             endpoint = envir$outArgs$endpoint,
                                                             index.survivalM1 = envir$outArgs$index.survivalM1,                       
                                                             keep.comparison = FALSE,
                                                             level.treatment = envir$outArgs$level.treatment,
                                                             method.inference = envir$outArgs$method.inference,
                                                             method.tte = envir$outArgs$method.tte,
                                                             neutral.as.uninf = envir$outArgs$neutral.as.uninf,
                                                             n.strata = envir$outArgs$n.strata,
                                                             returnIndex = FALSE,
                                                             strata = envir$outArgs$allstrata,
                                                             threshold = envir$outArgs$threshold,
                                                             treatment = envir$outArgs$treatment,
                                                             threshold.TTEM1 = envir$outArgs$threshold.TTEM1,
                                                             type = envir$outArgs$type,
                                                             Wscheme = envir$outArgs$Wscheme)
                                              })
                                  )
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

        toExport <- c(".BuyseTest","initializeSurvival_Peto","initializeSurvival_Peron")

        iB <- NULL ## [:forCRANcheck:] foreach        
        ls.permutation <- foreach::`%dopar%`(
                                       foreach::foreach(iB=1:n.resampling,
                                                        .packages = "BuyseTest",
                                                        .options.snow = ls.options,
                                                        .export = toExport),                                            
                                       {
                                           .BuyseTest(data = envir$outArgs$data,
                                                      censoring = envir$outArgs$censoring,
                                                      correction.tte = envir$outArgs$correction.tte,
                                                      D = envir$outArgs$D,
                                                      D.TTE = envir$outArgs$D.TTE,
                                                      endpoint = envir$outArgs$endpoint,
                                                      index.survivalM1 = envir$outArgs$index.survivalM1,                       
                                                      keep.comparison = FALSE,
                                                      level.treatment = envir$outArgs$level.treatment,
                                                      method.inference = envir$outArgs$method.inference,
                                                      method.tte = envir$outArgs$method.tte,
                                                      neutral.as.uninf = envir$outArgs$neutral.as.uninf,
                                                      n.strata = envir$outArgs$n.strata,
                                                      returnIndex = FALSE,
                                                      strata = envir$outArgs$allstrata,
                                                      threshold = envir$outArgs$threshold,
                                                      treatment = envir$outArgs$treatment,
                                                      threshold.TTEM1 = envir$outArgs$threshold.TTEM1,
                                                      type = envir$outArgs$type,
                                                      Wscheme = envir$outArgs$Wscheme)
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
