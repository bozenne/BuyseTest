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
#' 
#' @keywords function internal BuyseTest

#' @rdname internal-computation
calcBootstrap <- function(envir){
  
  if(envir$cpus==1){ ## sequential boostrap
    if(envir$trace>1){cat("Sequential boostrap \n")}
    if(envir$trace>0){
      envir$index.trace <- unique(round(seq(1,envir$n.bootstrap,length.out=101)[-1])) # number of iterations corresponding to each percentage of progress in the bootstrap computation
      envir$cpu_name <- "cpu 1 : " # 1 seul cpu est utilise
      envir$title_pb <- paste(envir$cpu_name,"bootstrap iterations",sep="") # title of the bar displaying the progress
      label_pb <- envir$label_pb <- paste(envir$cpu_name,"0% done ",sep="") # inital legend of the bar displaying the progress
      pb <- envir$pb <- tcltk::tkProgressBar(envir$title_pb,envir$label_pb,0, 1, 0) # create and display the progress bar         
    }
    if(!is.null(envir$seed)){set.seed(envir$seed)} # set the seed 
    envir$delta_boot[] <- vapply(1:envir$n.bootstrap,function(x){envir$boot_ci(x,envir=envir)},matrix(1.5,envir$n.strata,envir$D)) # perform the bootstrap and return for each bootstrap sample a n.strata*D matrix
    if(envir$trace>0){close(envir$pb)}  # close the bar
    
    
  }else{ ## parallel boostrap
    #envir$cpus <- 1 #browser()
    if(envir$trace>1){cat("Parallel boostrap \n")}
    envir$nParallel.bootstrap <- round(envir$n.bootstrap/envir$cpus) # number of bootstrap sample by cpu
    
    # init
    if(envir$trace==3){ # affect the cpus
      snowfall::sfInit(parallel=TRUE, cpus=envir$cpus)  
    }else if(envir$trace %in% c(1,2)){ 
      suppressMessages(snowfall::sfInit(parallel=TRUE, cpus=envir$cpus))
    } else {
      suppressWarnings(suppressMessages(snowfall::sfInit(parallel=TRUE, cpus=envir$cpus))) # warning si proxy dans commandArgs() (provient de la fonction searchCommandline)
    }
    
    ## wrapper function to be called by each cpu
    wrapper <- function(x){
      return(vapply(1:envir$nParallel.bootstrap,function(x){envir$boot_ci(x,envir=envir)},matrix(1.5,envir$n.strata,envir$D))) # perform the bootstrap and return for each bootstrap sample a n.strata*D matrix
    }
    
    pid <- snowfall::sfClusterEval(Sys.getpid()) # identifier of each R session
    snowfall::sfExport("pid")
    invisible(snowfall::sfClusterEval(cpu_index <- which(pid %in% Sys.getpid()))) # identify to index of the session of each cpu 
    
    if(envir$trace>0){
      envir$index.trace <- unique(round(seq(1,envir$nParallel.bootstrap,length.out=101)[-1])) # number of iterations corresponding to each percentage of progress in the bootstrap computation
      
      snowfall::sfLibrary("tcltk", character.only=TRUE) # export tcltk library (required for the progress bar) to each cpu
      invisible(snowfall::sfClusterEval(cpu_name <- paste("cpu ",cpu_index," : ",sep=""))) # identify the name of each cpu session
      invisible(snowfall::sfClusterEval(title_pb <- paste(cpu_name,"bootstrap iterations",sep=""))) # affect the title the progress bar to each cpu 
      invisible(snowfall::sfClusterEval(label_pb <- paste(cpu_name," 0% done",sep=""))) # affect the initial legend of the progress bar to each cpu 
      invisible(snowfall::sfClusterEval(pb <- tcltk::tkProgressBar(title_pb,label_pb,0, 1, 0))) # affect and display the progress bar to each cpu 
    }
    
    ## export the needed variables
    snowfall::sfExport("envir")
    
    ## bootstrap    
    if(!is.null(envir$seed)){      
      snowfall::sfClusterEval(set.seed(envir$seed+cpu_index-1))
    }
    
    resIter <- snowfall::sfClusterApply(rep(NA,envir$cpus),wrapper) # call the cpu to perform bootstrap computations
    
    if(envir$trace==3){ # free the cpus
      snowfall::sfStop()
    }else{
      suppressMessages(snowfall::sfStop())
    }
    
    # merge results
    sapply(1:envir$cpus,function(iter_CPU){
      envir$delta_boot[,,seq((iter_CPU-1)*envir$nParallel.bootstrap+1,iter_CPU*envir$nParallel.bootstrap)] <- resIter[[iter_CPU]]
      return(1)
    })
    
    
  }
  
}

#' @rdname internal-computation
calcCI <- function(delta,delta_boot,endpoint,D,alternative,alpha,
                   n.bootstrap,cpus,trace){
  
  ## Confidence interval
  Delta_boot <- apply(delta_boot,c(2,3),sum) # sum of the proportion in favor of treatment over the strata for the bootstrap samples
  if(D>1){Delta_boot <- apply(Delta_boot,2,cumsum)} # cumulative sum over the endpoints to obtaine the cumulative proportion in favor of treatment for the bootstrap samples
  
  n.bootstrap_real <- n.bootstrap - sum(is.na(Delta_boot[nrow(Delta_boot),]))
  
  if(trace && sum(is.na(Delta_boot))>0){
    
    message <- "BuyseTest : bootsrap failed for some iterations \n"
    for(iter_endpoint in 1:D){
      message <- paste(message,
                       paste("endpoint ",iter_endpoint," (\"",endpoint[iter_endpoint],"\") : ",
                             sum(is.na(Delta_boot[iter_endpoint,]))," (",100*sum(is.na(Delta_boot[iter_endpoint,]))/n.bootstrap,"%) failures  \n",sep=""),
                       sep=" ")
    }
    
    if(cpus>1 && ((n.bootstrap/cpus) %% 1 != 0)){
      message <- paste(message,"probable cause : some iterations were not launched \n",
                       "because n.boostrap (=",n.bootstrap,") is not a multiple of cpus ",cpus," \n",sep="") 
    }else{
      message <- paste(message,"probable cause : resampling leaded to no control or no case \n",sep="")
    }
    
    warning(message)
    
  }
  
  # compute the quantiles for each endpoint of the cumulative proportions in favor of treatment (in the bootstrap sample)
  Delta_quantile <- apply(Delta_boot,1,function(x){stats::quantile(x,probs=c(alpha/2,1-alpha/2),na.rm=TRUE)})
  
  ## p. value
  cum_Delta <- apply(delta,2,sum) # sum of the proportion in favor of treatment over the strata for the punctual estimation
  if(D>1){cum_Delta <- cumsum(cum_Delta)} # cumulative proportions in favor of treatment for the punctual estimation (sum over the strata followed by cumulative sum over the endpoints)
  
  testp.value <- switch(alternative, # test whether each bootstrap samples is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                        "two.sided"=apply(Delta_boot,2,function(x){abs(cum_Delta)>abs(x)}),
                        "less"=apply(Delta_boot,2,function(x){cum_Delta<x}),
                        "more"=apply(Delta_boot,2,function(x){cum_Delta>x})
  )
  
  # mean over the bootstrap samples by endpoint (i.e. proportion of bootstrap sample that is more extreme)
  if(D==1){
    p.value <- 1-mean(testp.value,na.rm=TRUE)
  }else{
    p.value <- 1-rowMeans(testp.value,na.rm=TRUE)
  }
  
  #### export ####
  res <- list()
  res$p.value <- p.value
  res$Delta_quantile <- Delta_quantile
  res$n.bootstrap_real <- n.bootstrap_real
  return(res)
}



