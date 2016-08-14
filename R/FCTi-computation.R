#' @name internal-computation
#' @rdname internal-computation
#' @title internal functions for BuyseTest - computations
#' @aliases calcBootstrap calcCI
#' 
#' @description Functions called by \code{\link{BuyseTest}} to perform calculations.
#' 
#' @details
#' \code{calcBootstrap} performs the bootstrap iterations in sequential or parallel mode. \cr
#' \code{calcCI} post process the results to obtain the p.value and the confidence interval. \cr
#' \code{warper_BTboot} resample the data for the cpp function. \cr
#' 
#' @keywords function internal BuyseTest

#' @rdname internal-computation
calcBootstrap <- function(envir){
  
  if (envir$cpus == 1) { ## sequential boostrap
    if (envir$trace > 1) {cat("Sequential boostrap \n")}
    if (envir$trace > 0) {
      envir$index.trace <- unique(round(seq(1, envir$n.bootstrap, length.out = 101)[-1])) # number of iterations corresponding to each percentage of progress in the bootstrap computation
      envir$cpu_name <- "cpu 1 : " # only 1 cpu is used
      envir$title_pb <- paste(envir$cpu_name,"bootstrap iterations",sep = "") # title of the bar displaying the progress
      label_pb <- envir$label_pb <- paste(envir$cpu_name,"0% done ",sep = "") # inital legend of the bar displaying the progress
      pb <- envir$pb <- tcltk::tkProgressBar(envir$title_pb,envir$label_pb, 0, 1, 0) # create and display the progress bar         
      on.exit(close(envir$pb)) # close the bar
    }else{
      envir$index.trace <- -1
    }
    if (!is.null(envir$seed)) {set.seed(envir$seed)} # set the seed 
    envir$delta_boot[] <- vapply(1:envir$n.bootstrap, 
                                 function(x){warper_BTboot(x, envir = envir)},
                                 matrix(1.5, envir$n.strata, envir$D)) # perform the bootstrap and return for each bootstrap sample a n.strata*D matrix
    
  }else {## parallel boostrap
     #envir$cpus <- 1 #browser()
    if (envir$trace > 1) {cat("Parallel boostrap \n")}
    
    ## init
    if (envir$trace == 3) { # affect the cpus
      snowfall::sfInit(parallel = TRUE, cpus = envir$cpus)  
    }else if (envir$trace %in% c(1,2)) { 
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = envir$cpus))
    } else {
      suppressWarnings(suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = envir$cpus))) # warning si proxy dans commandArgs() (provient de la fonction searchCommandline)
    }
    pid <- unlist(snowfall::sfClusterEval(Sys.getpid())) # identifier of each R session
    snowfall::sfExport("pid")
    
    ## wrapper function to be called by each cpu
    wrapper <- function(x){
      return(vapply(1:envir$nParallel.bootstrap,
                    function(x){envir$warper_BTboot(x, envir = envir)},
                    matrix(1.5, envir$n.strata, envir$D))) # perform the bootstrap and return for each bootstrap sample a n.strata*D matrix
    }
    
    ## export the needed variables
    envir$nParallel.bootstrap <- round(envir$n.bootstrap/envir$cpus) # number of bootstrap sample by cpu
    envir$warper_BTboot <- warper_BTboot
    envir$pid <- pid
    if (envir$trace > 0) {
      envir$index.trace <- unique(round(seq(1, envir$nParallel.bootstrap, length.out = 101)[-1])) # number of iterations corresponding to each percentage of progress in the bootstrap computation
    }else{
      envir$index.trace <- -1
    }
    snowfall::sfExport("envir")
    
    ## seed
    if (!is.null(envir$seed)) { 
      invisible(snowfall::sfClusterEval(set.seed(envir$seed + which(pid %in% Sys.getpid()) - 1)))
    }
    
    ## trace
    if (envir$trace > 0) {
      snowfall::sfLibrary("tcltk", character.only = TRUE) # export tcltk library (required for the progress bar) to each cpu
      invisible(snowfall::sfClusterEval(title_pb <- paste0("cpu ", which(pid %in% Sys.getpid()), ": bootstrap iterations"))) # affect the title the progress bar to each cpu 
      invisible(snowfall::sfClusterEval(label_pb <- paste0("cpu ", which(pid %in% Sys.getpid()), ": 0% done"))) # affect the initial legend of the progress bar to each cpu 
      envir$pb <- snowfall::sfClusterEval(tcltk::tkProgressBar(title_pb, label_pb, 0, 1, 0)) # affect and display the progress bar to each cpu 
      snowfall::sfExport("envir")
    }
  
    ## bootstrap    
    resIter <- snowfall::sfClusterApply(rep(NA, envir$cpus), wrapper) # call the cpu to perform bootstrap computations
    
    if (envir$trace == 3) { # free the cpus
      snowfall::sfStop()
    }else{
      suppressMessages(snowfall::sfStop())
    }
    
    # merge results
    sapply(1:envir$cpus,function(iter_CPU){
      envir$delta_boot[,,seq((iter_CPU - 1) * envir$nParallel.bootstrap + 1, iter_CPU * envir$nParallel.bootstrap)] <- resIter[[iter_CPU]]
      return(1)
    })
    
    
  }
  
}

#' @rdname internal-computation
calcCI <- function(delta,delta_boot,endpoint,D,alternative,alpha,
                   n.bootstrap,cpus,trace){
  
  ## Confidence interval
  Delta_boot <- apply(delta_boot,c(2,3),sum) # sum of the proportion in favor of treatment over the strata for the bootstrap samples
  if (D > 1) {Delta_boot <- apply(Delta_boot,2,cumsum)} # cumulative sum over the endpoints to obtaine the cumulative proportion in favor of treatment for the bootstrap samples
  
  n.bootstrap_real <- n.bootstrap - sum(is.na(Delta_boot[nrow(Delta_boot),]))
  
  if (trace && sum(is.na(Delta_boot)) > 0) {
    
    message <- "BuyseTest : bootsrap failed for some iterations \n"
    for (iter_endpoint in 1:D) {
      message <- paste(message,
                       paste("endpoint ",iter_endpoint," (\"",endpoint[iter_endpoint],"\") : ",
                             sum(is.na(Delta_boot[iter_endpoint,]))," (",100*sum(is.na(Delta_boot[iter_endpoint,]))/n.bootstrap,"%) failures  \n",sep = ""),
                       sep = " ")
    }
    
    if (cpus > 1 && ((n.bootstrap/cpus) %% 1 != 0)) {
      message <- paste(message,"probable cause : some iterations were not launched \n",
                       "because n.boostrap (=",n.bootstrap,") is not a multiple of cpus ",cpus," \n",sep = "") 
    }else{
      message <- paste(message,"probable cause : resampling leaded to no control or no case \n",sep = "")
    }
    
    warning(message)
    
  }
  
  # compute the quantiles for each endpoint of the cumulative proportions in favor of treatment (in the bootstrap sample)
  Delta_quantile <- apply(Delta_boot,1,function(x){stats::quantile(x,probs = c(alpha/2,1 - alpha/2),na.rm = TRUE)})
  
  ## p. value
  cum_Delta <- apply(delta,2,sum) # sum of the proportion in favor of treatment over the strata for the punctual estimation
  if (D > 1) {cum_Delta <- cumsum(cum_Delta)} # cumulative proportions in favor of treatment for the punctual estimation (sum over the strata followed by cumulative sum over the endpoints)
  
  testp.value <- switch(alternative, # test whether each bootstrap samples is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                        "two.sided" = apply(Delta_boot, 2,function(x){abs(cum_Delta) > abs(x)}),
                        "less" = apply(Delta_boot, 2,function(x){cum_Delta < x}),
                        "more" = apply(Delta_boot, 2,function(x){cum_Delta > x})
  )
  
  # mean over the bootstrap samples by endpoint (i.e. proportion of bootstrap sample that is more extreme)
  if (D == 1) {
    p.value <- 1 - mean(testp.value, na.rm = TRUE)
  }else{
    p.value <- 1 - rowMeans(testp.value, na.rm = TRUE)
  }
  
  #### export ####
  res <- list()
  res$p.value <- p.value
  res$Delta_quantile <- Delta_quantile
  res$n.bootstrap_real <- n.bootstrap_real
  return(res)
}


#' @rdname internal-computation
warper_BTboot <- function(x,envir){
  
  #### prepare ####
  Mnew.Treatment <- NULL
  Mnew.Control <- NULL
  if (envir$D.TTE > 0) { 
    Mnew.delta_Treatment <- NULL
    Mnew.delta_Control <- NULL
    if (envir$method == "Peto") {
      new.survivalT <- vector(length = envir$D.TTE, mode = "list")
      new.survivalC <- vector(length = envir$D.TTE, mode = "list")
    }
  }
  new.strataT <- list()
  new.strataC <- list()
  
  ## for data management
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
    
    nboot.T <- length(indexT.T) + length(indexT.C)
    nboot.C <- length(indexC.T) + length(indexC.C)
    
    ## test whether there are observations in each group
    if (nboot.T == 0 || nboot.C == 0) {
      boot.failure <- TRUE ; break ;
    }else{
      boot.failure <- FALSE
    }
    
    ## update dataset
    Mnew.Treatment <- rbind(Mnew.Treatment, 
                            envir$M.Treatment[indexT.T,,drop = FALSE], 
                            envir$M.Control[indexT.C,,drop = FALSE])
    Mnew.Control <- rbind(Mnew.Control, 
                          envir$M.Treatment[indexC.T,,drop = FALSE],
                          envir$M.Control[indexC.C,,drop = FALSE]) 
    
    ## update strata - the new strata index minus 1 (for C++ compatibility, vector begin at 0)
    new.strataT[[iterS]] <- seq(from = new.nT, by = 1, length = nboot.T) 
    new.strataC[[iterS]] <- seq(from = new.nC, by = 1, length = nboot.C)
    
    ## update censoring variables
    if (envir$D.TTE > 0) { 
      Mnew.delta_Treatment <- rbind(Mnew.delta_Treatment, 
                                    envir$M.delta_Treatment[indexT.T,, drop = FALSE], 
                                    envir$M.delta_Control[indexT.C,, drop = FALSE])
      Mnew.delta_Control <- rbind(Mnew.delta_Control, 
                                  envir$M.delta_Treatment[indexC.T,, drop = FALSE],
                                  envir$M.delta_Control[indexC.C,, drop = FALSE])
      
      if (envir$method == "Peto") {
        for (iter_endpointTTE in 1:envir$D.TTE) {
          new.survivalT[[iter_endpointTTE]] <- rbind(new.survivalT[[iter_endpointTTE]], 
                                                     envir$list_survivalT[[iter_endpointTTE]][indexT.T,], 
                                                     envir$list_survivalC[[iter_endpointTTE]][indexT.C,])
          new.survivalC[[iter_endpointTTE]] <- rbind(new.survivalC[[iter_endpointTTE]], 
                                                     envir$list_survivalT[[iter_endpointTTE]][indexC.T,], 
                                                     envir$list_survivalC[[iter_endpointTTE]][indexC.C,])
        }
      }else if (envir$method == "Efron") { # set last event to non-censored
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
    
    new.nT <- new.nT + nboot.T
    new.nC <- new.nC + nboot.C
  }
  
  if (boot.failure == FALSE) {
    if (envir$method %in% c("Efron","Peron")) { # update survival
      res_init <- BuyseTest::initSurvival(M.Treatment = Mnew.Treatment,M.Control = Mnew.Control,
                                          M.delta_Treatment = Mnew.delta_Treatment,M.delta_Control = Mnew.delta_Control,
                                          endpoint = envir$endpoint,D.TTE = envir$D.TTE,type = envir$type,threshold = envir$threshold,
                                          index.strataT = new.strataT,index.strataC = new.strataC,n.strata = envir$n.strata,
                                          method = envir$method)
    }
    
    #### Computation
    if (envir$method %in% c("Peto","Efron","Peron")) {
      delta_boot <-   BuyseTest::BuyseTest_PetoEfronPeron_cpp(Treatment = Mnew.Treatment,Control = Mnew.Control,threshold = envir$threshold,type = envir$type,
                                                              delta_Treatment = Mnew.delta_Treatment,delta_Control = Mnew.delta_Control,
                                                              D = envir$D,returnIndex = FALSE,
                                                              strataT = new.strataT,strataC = new.strataC,n_strata = envir$n.strata,n_TTE = envir$D.TTE,
                                                              Wscheme = envir$Wscheme,index_survivalM1 = envir$index_survivalM1,threshold_TTEM1 = envir$threshold_TTEM1,
                                                              list_survivalT = if (envir$method %in% c("Efron","Peron")) {res_init$list_survivalT} else {new.survivalT},
                                                              list_survivalC = if (envir$method %in% c("Efron","Peron")) {res_init$list_survivalC} else {new.survivalC},
                                                              PEP = which(c("Peto","Efron","Peron") == envir$method))$delta
      
    }else if (envir$method == "Gehan") {
      
      delta_boot <- BuyseTest::BuyseTest_Gehan_cpp(Treatment = Mnew.Treatment,Control = Mnew.Control,threshold = envir$threshold,type = envir$type,
                                                   delta_Treatment = Mnew.delta_Treatment,delta_Control = Mnew.delta_Control,
                                                   D = envir$D,returnIndex = FALSE,
                                                   strataT = new.strataT,strataC = new.strataC,n_strata = envir$n.strata,n_TTE = envir$D.TTE)$delta
     
    }
  }else{
    delta_boot <- matrix(NA, nrow = envir$n.strata, ncol = envir$D)
  }
  
  
  if (x == envir$index.trace[1]) { ## display
    # if the following percentage of iterations is reached
    # the label on the bar is updated and the progress bar is displayed
    # then the first element of index.trace to test for the next percentage of iterations
    
    if (envir$cpus == 1) {
      pc_done <- envir$index.trace[1]/envir$n.bootstrap
      label_pb <- paste(envir$cpu_name,round(100 * pc_done), "% done", sep = "")
      tcltk::setTkProgressBar(pb = envir$pb, value = pc_done, title = paste(envir$title_pb, "(", round(100 * pc_done), "%)", sep = ""), label = label_pb)
      envir$index.trace <- envir$index.trace[-1]                 
    }else{
      pc_done <- envir$index.trace[1]/envir$nParallel.bootstrap
      
      cpu_index <- which(envir$pid %in% Sys.getpid())
      title_pb <- paste0("cpu ", cpu_index, ": bootstrap iterations")
      label_pb <- paste0("cpu ", cpu_index, ": ",round(100 * pc_done), "% done")
      
      tcltk::setTkProgressBar(pb = envir$pb[[cpu_index]], value = pc_done,title = paste(title_pb,"(",round(100 * pc_done),"%)",sep = ""), label = label_pb)
      envir$index.trace <- envir$index.trace[-1]                 
    } 
  }
  
  ## export
  return(delta_boot)
}






