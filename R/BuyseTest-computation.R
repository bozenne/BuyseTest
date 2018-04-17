## * Documentation
#' @name internal-computation
#' @rdname internal-computation
#' @title internal functions for BuyseTest - computations
#' @aliases calcPermutation calcCI
#' 
#' @description Functions called by \code{\link{BuyseTest}} to perform calculations.
#' 
#' @details
#' \code{calcPermutation} performs the permutation test in sequential or parallel mode. \cr
#' \code{calcCI} post process the results to obtain the p.value and the confidence interval. \cr
#' \code{warper_BTpermutation} resample the data for the cpp function. \cr
#' 
#' @keywords function internal BuyseTest

## * calcPermutation
#' @rdname internal-computation
calcPermutation <- function(envir){
  
  resPermutation <- array(NA, dim = c(envir$n.strata + 1, 2*envir$D,  envir$n.permutation))
  
    if (envir$cpus == 1) { ## ** sequential permutation test
        if (envir$trace > 1) {cat("Sequential permutation test \n")}
        if (envir$trace > 0) {
            envir$index.trace <- unique(round(seq(1, envir$n.permutation, length.out = 101)[-1])) # number of iterations corresponding to each percentage of progress in the permutation test computation
            envir$cpu_name <- "cpu 1 : " # only 1 cpu is used
            envir$title_pb <- paste(envir$cpu_name,"Permutation test",sep = "") # title of the bar displaying the progress
            label_pb <- envir$label_pb <- paste(envir$cpu_name,"0% done ",sep = "") # inital legend of the bar displaying the progress
            pb <- envir$pb <- tcltk::tkProgressBar(envir$title_pb,envir$label_pb, 0, 1, 0) # create and display the progress bar         
            on.exit(close(envir$pb)) # close the bar
        }else{
            envir$index.trace <- -1
        }
        if (!is.null(envir$seed)) {set.seed(envir$seed)} # set the seed 
        resPermutation[] <- vapply(1:envir$n.permutation, 
                                   function(x){warper_BTpermutation(x, envir = envir)},
                                   matrix(1.5, envir$n.strata + 1, 2*envir$D)) # perform the permutation test and return for each sample a (n.strata + 1)*(2*D) matrix
        if (envir$trace > 1) {cat("   > done \n")}
    }else { ## ** parallel permutation test
        if (envir$trace > 1) {cat("Parallel permutation test \n")}
    
        ## *** initsfLibrary
        if (envir$trace == 3) { # affect the cpus
            snowfall::sfInit(parallel = TRUE, cpus = envir$cpus)  
        }else if (envir$trace %in% c(1,2)) { 
            capture.output(suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = envir$cpus)))
        } else {
            capture.output(suppressWarnings(suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = envir$cpus)))) # warning si proxy dans commandArgs() (provient de la fonction searchCommandline)
        }
    
        pid <- unlist(snowfall::sfClusterEval(Sys.getpid())) # identifier of each R session
        snowfall::sfExport("pid")
        capture.output(suppressMessages(snowfall::sfLibrary("BuyseTest", character.only = TRUE)))
    
        ## *** wrapper function to be called by each cpu
        wrapper <- function(x){
            return(vapply(1:envir$nParallel.permutation,
                          function(x){envir$warper_BTpermutation(x, envir = envir)},
                          matrix(1.5, envir$n.strata + 1, 2 * envir$D))) # perform the permutation test and return for each sample a (n.strata+1)*(2*D) matrix
        }
    
        ## *** export the needed variables
        envir$nParallel.permutation <- round(envir$n.permutation/envir$cpus) # number of sample by cpu
        envir$warper_BTpermutation <- warper_BTpermutation
        envir$pid <- pid
        if (envir$trace > 0) {
            envir$index.trace <- unique(round(seq(1, envir$nParallel.permutation, length.out = 101)[-1])) # number of iterations corresponding to each percentage of progress of the percentage test
        }else{
            envir$index.trace <- -1
        }
        snowfall::sfExport("envir")
    
        ## *** seed
        if (!is.null(envir$seed)) { 
            invisible(snowfall::sfClusterEval(set.seed(envir$seed + which(pid %in% Sys.getpid()) - 1)))
        }
    
        ## *** trace
        if (envir$trace > 0) {
            capture.output(suppressMessages(snowfall::sfLibrary("tcltk", character.only = TRUE))) # export tcltk library (required for the progress bar) to each cpu
            invisible(snowfall::sfClusterEval(title_pb <- paste0("cpu ", which(pid %in% Sys.getpid()), ": Permutation test"))) # affect the title the progress bar to each cpu 
            invisible(snowfall::sfClusterEval(label_pb <- paste0("cpu ", which(pid %in% Sys.getpid()), ": 0% done"))) # affect the initial legend of the progress bar to each cpu 
            envir$pb <- snowfall::sfClusterEval(tcltk::tkProgressBar(title_pb, label_pb, 0, 1, 0)) # affect and display the progress bar to each cpu 
            snowfall::sfExport("envir")
        }
    
        ## *** permutation test
        resIter <- snowfall::sfClusterApply(rep(NA, envir$cpus), wrapper) # call the cpu to perform the permutation test
    
        if (envir$trace == 3) { # free the cpus
            snowfall::sfStop()
        }else{
            suppressMessages(snowfall::sfStop())
        }
    
        ## *** merge results
        for(iterCPU in 1:envir$cpus){
            resPermutation[,,seq((iterCPU - 1) * envir$nParallel.permutation + 1, iterCPU * envir$nParallel.permutation)] <- resIter[[iterCPU]]
        }
        if (envir$trace > 1) {cat("   > done \n")}
    }
  
  #### export
  envir$delta_permutation$netChance[] <- resPermutation[1:envir$n.strata,1:envir$D,]
  envir$Delta_permutation$netChance[] <- resPermutation[envir$n.strata + 1,1:envir$D,]
  
  envir$delta_permutation$winRatio[] <- resPermutation[1:envir$n.strata,(envir$D + 1):(2 * envir$D),]
  envir$Delta_permutation$winRatio[] <- resPermutation[envir$n.strata + 1,(envir$D + 1):(2 * envir$D),]
  
}

## * calcCI
#' @rdname internal-computation
calcCI <- function(Delta,Delta_permutation,
                   endpoint,D,alternative,alpha,
                   n.permutation,cpus,trace){

    option <- BuyseTest.options()
        
    ## Confidence interval
    ## Delta_permutation <- apply(delta_permutation,c(2,3),sum) # sum of the proportion in favor of treatment over the strata for the permutation test.
    ## if (D > 1) {Delta_permutation <- apply(Delta_permutation,2,cumsum)} # cumulative sum over the endpoints to obtaine the cumulative proportion in favor of treatment for the permutation test
    n.permutation_real <- lapply(Delta_permutation, function(x){n.permutation - apply(is.na(x), 1, sum)})
  
    if (trace && any(unlist(n.permutation_real) < n.permutation)) {
        message <- "BuyseTest : for some samples, the permutation procedure failed and returned NA \n the sample will be ignore in the finale analysis \n"
        indexD <- lapply(n.permutation_real, function(x){which(x<n.permutation)})
    
        for (iter_endpoint in sort(unique(unlist(indexD)))) {
            message <- paste(message,
                             paste("endpoint ",iter_endpoint," (\"",endpoint[iter_endpoint],"\") : ",
                                   sum(is.na(Delta_permutation[[option$statistic]][iter_endpoint,])),
                                   " (",100*sum(is.na(Delta_permutation[[option$statistic]][iter_endpoint,]))/n.permutation,"%) failures  \n",sep = ""),
                             sep = " ")
        }
    
        if (cpus > 1 && ((n.permutation/cpus) %% 1 != 0)) {
            message <- paste(message,"probable cause : some iterations were not launched \n",
                             "because n.permutation (=",n.permutation,") is not a multiple of cpus ",cpus," \n",sep = "") 
        }else{
            message <- paste(message,"possible cause : resampling leaded to no control or no case \n",sep = "")
        }
    
        warning(message)
    
    }
    
    ## ** compute the quantiles for each endpoint of the cumulative proportions in favor of treatment (in the sample)
    Delta_quantile <- lapply(1:2, function(type){
        apply(Delta_permutation[[type]],1,function(x){stats::quantile(x,probs = c(alpha/2,1 - alpha/2),na.rm = TRUE)})
    })
    names(Delta_quantile) <- names(Delta)
  
    ## ** p. value (we don't really need to differientiate netChance and winRatio)
    testp.value <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                          "two.sided" = lapply(1:2, function(type){
                              apply(Delta_permutation[[type]], 2,function(x){abs(Delta[[type]] - c(0,0.5)[type]) > abs(x - c(0,0.5)[type])})
                          }),
                          "less" = lapply(1:2, function(type){
                              apply(Delta_permutation[[type]], 2,function(x){(Delta[[type]] - c(0,0.5)[type]) < (x - c(0,0.5)[type])})
                          }),
                          "more" =  lapply(1:2, function(type){
                              apply(Delta_permutation[[type]], 2,function(x){(Delta[[type]] - c(0,0.5)[type]) > (x - c(0,0.5)[type])})
                          })
                          )
  
    ## ** mean over the samples by endpoint
    ## (i.e. proportion of samples in which the statistic obtained after permutation is more extreme than the observed one)
    if (D == 1) {
        p.value <- lapply(testp.value, function(x){1 - mean(x, na.rm = TRUE)})
    }else{
        p.value <- lapply(testp.value, function(x){1 - rowMeans(x, na.rm = TRUE)})
    }
    names(p.value) <- names(Delta)
  
    ## ** export  
    res <- list()
    res$p.value <- p.value
    res$Delta_quantile <- Delta_quantile
    res$n_permutation <- n.permutation_real
    return(res)
}

## * warper_BTpermutation
#' @rdname internal-computation
warper_BTpermutation <- function(x,envir){
  
    ## ** prepare
    Mnew.Treatment <- NULL
    Mnew.Control <- NULL
    Mnew.delta_Treatment <- NULL
    Mnew.delta_Control <- NULL
    new.survivalT <- if(envir$method == 0){lapply(1:envir$D.TTE, matrix)}else{vector(length = envir$D.TTE, mode = "list")}
    new.survivalC <- if(envir$method == 0){lapply(1:envir$D.TTE, matrix)}else{vector(length = envir$D.TTE, mode = "list")}
    new.strataT <- list()
    new.strataC <- list()
  
    ## ** for data management
    new.nT <- 0
    new.nC <- 0
  
    if (envir$stratified == FALSE) {
        groupT.Tall <- which(stats::rbinom(envir$n.Treatment, size = 1, prob = envir$prob.alloc) == 1)
        groupT.Call <- which(stats::rbinom(envir$n.Control, size = 1, prob = envir$prob.alloc) == 1)
    }
  
    for (iterS in 1:envir$n.strata) {  ## randomisation : new allocation of the treatment and control arm
    
        indexStrataT <- seq(envir$nCumSum.strataTreatment[iterS],envir$nCumSum.strataTreatment[iterS + 1] - 1)
        indexStrataC <- seq(envir$nCumSum.strataControl[iterS],envir$nCumSum.strataControl[iterS + 1] - 1)
    
        if (envir$stratified) {
            indexT.T <- sample.int(indexStrataT, 
                                   size = round(envir$n.eachStrataT[iterS]*envir$prob.alloc), 
                                   replace = FALSE)
            indexT.C <- sample.int(indexStrataC, 
                                   size = round(envir$n.eachStrataC[iterS]*envir$prob.alloc), 
                                   replace = FALSE)
        }else{
            indexT.T <- intersect(groupT.Tall, indexStrataT)
            indexT.C <- intersect(groupT.Call, indexStrataC)
        }
    
        indexC.T <- setdiff(indexStrataT,indexT.T)
        indexC.C <- setdiff(indexStrataC,indexT.C)
    
        nPermutation.T <- length(indexT.T) + length(indexT.C)
        nPermutation.C <- length(indexC.T) + length(indexC.C)
    
        ## *** test whether there are observations in each group
        if (nPermutation.T == 0 || nPermutation.C == 0) {
            permutation.failure <- TRUE ; break ;
        }else{
            permutation.failure <- FALSE
        }
    
        ## *** update dataset
        Mnew.Treatment <- rbind(Mnew.Treatment, 
                                envir$M.Treatment[indexT.T,,drop = FALSE], 
                                envir$M.Control[indexT.C,,drop = FALSE])
        Mnew.Control <- rbind(Mnew.Control, 
                              envir$M.Treatment[indexC.T,,drop = FALSE],
                              envir$M.Control[indexC.C,,drop = FALSE]) 
    
        ## *** update strata - the new strata index minus 1 (for C++ compatibility, vector begin at 0)
        new.strataT[[iterS]] <- seq(from = new.nT, by = 1, length = nPermutation.T) 
        new.strataC[[iterS]] <- seq(from = new.nC, by = 1, length = nPermutation.C)
    
        ## *** update censoring variables
        if (envir$D.TTE > 0) { 
            Mnew.delta_Treatment <- rbind(Mnew.delta_Treatment, 
                                          envir$M.delta_Treatment[indexT.T,, drop = FALSE], 
                                          envir$M.delta_Control[indexT.C,, drop = FALSE])
            Mnew.delta_Control <- rbind(Mnew.delta_Control, 
                                        envir$M.delta_Treatment[indexC.T,, drop = FALSE],
                                        envir$M.delta_Control[indexC.C,, drop = FALSE])
      
            if (envir$method == 1) { # "Peto"
                for (iter_endpointTTE in 1:envir$D.TTE) {
                    new.survivalT[[iter_endpointTTE]] <- rbind(new.survivalT[[iter_endpointTTE]], 
                                                               envir$list_survivalT[[iter_endpointTTE]][indexT.T,], 
                                                               envir$list_survivalC[[iter_endpointTTE]][indexT.C,])
                    new.survivalC[[iter_endpointTTE]] <- rbind(new.survivalC[[iter_endpointTTE]], 
                                                               envir$list_survivalT[[iter_endpointTTE]][indexC.T,], 
                                                               envir$list_survivalC[[iter_endpointTTE]][indexC.C,])
                }
            }else if (envir$method == 2) { # "Efron" # set last event to non-censored
                Mnewstrata.Treatment <- Mnew.Treatment[new.strataT[[iterS]] + 1,which(envir$type == 3), drop = FALSE]
                Mnewstrata.Control <- Mnew.Control[new.strataC[[iterS]] + 1,which(envir$type == 3), drop = FALSE]
        
                for (iter_endpointTTE in 1:envir$D.TTE) {
                    indexT_maxCensored <- which.max(Mnewstrata.Treatment[,iter_endpointTTE])
                    Mnew.delta_Treatment[new.strataT[[iterS]][indexT_maxCensored] + 1,iter_endpointTTE] <- 1
                    indexC_maxCensored <- which.max(Mnewstrata.Control[,iter_endpointTTE])
                    Mnew.delta_Control[new.strataC[[iterS]][indexC_maxCensored] + 1,iter_endpointTTE] <- 1
                }
            }
        }else {
            Mnew.delta_Treatment <- matrix(-1,1,1)
            Mnew.delta_Control <- matrix(-1,1,1)
        }
    
    new.nT <- new.nT + nPermutation.T
    new.nC <- new.nC + nPermutation.C
  }
  
    if (permutation.failure == FALSE) {
        if (envir$method %in% 2:3) { # c("Efron","Peron") # update survival
            res_init <- initSurvival(M.Treatment = Mnew.Treatment,M.Control = Mnew.Control,
                                     M.delta_Treatment = Mnew.delta_Treatment,M.delta_Control = Mnew.delta_Control,
                                     endpoint = envir$endpoint,D.TTE = envir$D.TTE,type = envir$type,threshold = envir$threshold,
                                     index.strataT = new.strataT,index.strataC = new.strataC,n.strata = envir$n.strata,
                                     method = envir$method)
        }
    
        ## ** Computation
        resBT <-   GPC_cpp(Treatment = Mnew.Treatment,
                           Control = Mnew.Control,
                           threshold = envir$threshold,
                           survEndpoint = (envir$type == 3),
                           delta_Treatment = Mnew.delta_Treatment,
                           delta_Control = Mnew.delta_Control,
                           D = envir$D,
                           returnIndex = FALSE,
                           strataT = new.strataT,
                           strataC = new.strataC,
                           n_strata = envir$n.strata,
                           n_TTE = envir$D.TTE,
                           Wscheme = envir$Wscheme,
                           index_survivalM1 = envir$index_survivalM1,
                           threshold_TTEM1 = envir$threshold_TTEM1,
                           list_survivalT = if (envir$method %in% 2:3) {res_init$list_survivalT} else {new.survivalT},
                           list_survivalC = if (envir$method %in% 2:3) {res_init$list_survivalC} else {new.survivalC},
                           methodTTE = envir$method,
                           correctionTTE = envir$correctionTTE,
                           neutral.as.uninf = envir$neutral.as.uninf,
                           keepComparison = FALSE
                           )
    
    
      
      resWarper <- cbind(rbind(resBT$delta_netChance, resBT$Delta_netChance),
                         rbind(resBT$delta_winRatio, resBT$Delta_winRatio))
    
    }else{
        resWarper <- matrix(NA, nrow = envir$n.strata + 1, ncol = 2*envir$D)
    }
  
  
    if (x == envir$index.trace[1]) { ## ** display
        ## if the following percentage of iterations is reached
        ## the label on the bar is updated and the progress bar is displayed
        ## then the first element of index.trace to test for the next percentage of iterations
    
        if (envir$cpus == 1) {
            pc_done <- envir$index.trace[1]/envir$n.permutation
            label_pb <- paste(envir$cpu_name,round(100 * pc_done), "% done", sep = "")
            tcltk::setTkProgressBar(pb = envir$pb, value = pc_done, title = paste(envir$title_pb, "(", round(100 * pc_done), "%)", sep = ""), label = label_pb)
            envir$index.trace <- envir$index.trace[-1]                 
        }else{
            pc_done <- envir$index.trace[1]/envir$nParallel.permutation
      
            cpu_index <- which(envir$pid %in% Sys.getpid())
            title_pb <- paste0("cpu ", cpu_index, ": Permutation test")
            label_pb <- paste0("cpu ", cpu_index, ": ",round(100 * pc_done), "% done")
      
            tcltk::setTkProgressBar(pb = envir$pb[[cpu_index]], value = pc_done,title = paste(title_pb,"(",round(100 * pc_done),"%)",sep = ""), label = label_pb)
            envir$index.trace <- envir$index.trace[-1]                 
        } 
    }
  
    ## ** export
    return(resWarper)
}






