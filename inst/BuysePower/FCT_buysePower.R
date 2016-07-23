#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% sample size computation 
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BuysePower <- function(n,hypothesis,type, threshold=NULL,endpoint=NULL,alpha=0.05, callMethod="Gehan",
                       n.bootstrap=1000,prob.alloc=0.5,n.simul=1000,
                       alternative = "two.sided",seed=10,cpus=1,trace=3){
  
  #### preliminary tests ####
  
  ## n
  n <- initN(n,callMethod="BuysePower")
  n.n <- length(n)
  
  ## type
  type <- initType(type=type,callMethod="BuyseTest")
  D <- length(type) # number of endpoints
  D.TTE <- sum(type==3) # number of time to event endpoints
  
  if(method %in% c("Peto","Efron") && D.TTE==0){ 
    if(trace>0){cat("NOTE : Kaplan Meier imputation is requested but there is no survival endpoint \n")}
    method <- "Gehan"
  }
  
  ## parameters
  hypothesis <- initHypothesis(hypothesis=hypothesis,type=type,D=D,callMethod="BuysePower")
  
  ## endpoint
  if(is.null(endpoint)){
    endpoint <- paste(c("bin","cont","TTE")[type],1:length(type),sep="_")
  }else{
    if(length(endpoint)!=D){
      stop("BuysePower : incorrect specification of \'endpoint\' \n",
           "\'endpoint\' must have length : ",D," \n",
           "length(endpoint) : ",length(endpoint),"\n")
    }
  }
  
  ## threshold
  threshold <- initThreshold(threshold=threshold,type=type,D=D,
                             method=method,endpoint=endpoint,
                             callMethod="BuysePower")
  
  ## miscellaneous
  if(alternative %in% c("two.sided", "less", "greater")==FALSE){
    stop("BuysePower : incorrect specification of \'alternative\' \n",
         "valid \'alternative\' : \"two.sided\" \"less\" \"greater\" \n",
         "proposed \'alternative\' : ",alternative,"\n")
  } 
  
  if(!is.numeric(n.bootstrap) || n.bootstrap < 0 || (n.bootstrap %% 1 !=0)){
    stop("BuysePower : wrong specification of \'n.bootstrap\' \n",
         "\'n.bootstrap\' must be a positive integer \n",
         "proposed \'n.bootstrap\' : ",n.bootstrap,"\n")
  }
  
  if(!is.numeric(n.simul) || n.simul < 0 || (n.simul %% 1 !=0)){
    stop("BuysePower : wrong specification of \'n.simul\' \n",
         "\'n.simul\' must be a positive integer \n",
         "proposed \'n.simul\' : ",n.simul,"\n")
  }
  
  ## proba
  prob.alloc <- initAlloc(prob.alloc,n.Treatment=n.Treatment,n.Control=n.Control,callMethod="BuysePower")
  n.Treatment <- round(n*prob.alloc)
  n.Control <- round(n*prob.alloc)
  n.Integer <- n.Treatment+n.Control
  
  ## seed
  if(!is.null(seed) && (!is.numeric(seed) || seed <= 0 || (seed %% 1 !=0))){
    stop("BuysePower : wrong specification of \'seed\' \n",
         "\'seed\' must be a strictly positive integer \n",
         "proposed \'seed\' : ",seed,"\n")
  }
  
  ## cpu
  cpus <- BuyseTest::initCPU(cpus=cpus,callMethod="BuysePower")
  export_names <- c("nParallel.bootstrap","n.Treatment","prob.alloc",
                    "M.Treatment","M.Control","threshold","type","D","D.TTE",
                    "n.n","n.Treatment","n.Control")
  
  if(D.TTE>0){export_names <- c(export_names,"M.delta_Treatment","M.delta_Control")}
  if(method=="Peto"){export_names <- c(export_names,"Wscheme","list_survivalT","list_survivalC")}
  
  
  ## method
  if(method %in% c("Peto","Efron")){
    res_init <- initWscheme(endpoint.TTE=endpoint[type==3],D.TTE=sum(type==3),threshold=threshold,type=type)
    Wscheme <- res_init$Wscheme  
    index_TTEM1 <- res_init$index_TTEM1 
    threshold_TTEM1 <- res_init$threshold_TTEM1 
  }else{
    Wscheme <- NULL  
    index_TTEM1 <- NULL
  }
  
  #### display
  if(trace>1){
    n_nchar <- nchar(n)
    
    n.Integer_space <- initSpace(n_nchar-nchar(n.Integer))
    n.Treatment_space <- initSpace(n_nchar-nchar(n.Treatment))
    n.Control_space <- initSpace(n_nchar-nchar(n.Control))
    
    cat("Settings \n")
    
    cat("   # first species risk              : ",alpha," (",alternative,") \n",sep="")    
    cat("   # proposed allocation probability : ",prob.alloc," \n")    
    cat("   # proposed sample size            : ",paste(n,collapse=" ")," \n")    
    cat("   # real sample size                : ",paste(n.Integer,n.Integer_space,collapse=""),"\n")
    cat("   # real size of the treatment arm  : ",paste(n.Treatment,n.Treatment_space,collapse=""),"\n")
    cat("   # real size of the control arm    : ",paste(n.Control,n.Control_space,collapse=""),"\n")
    
    #     iter_mTTE <- 0
    for(iter_m in 1:D){
      cat("      >  endpoint ",iter_m," = \"",endpoint[iter_m],"\" - type = \"",c("binary","continuous","timeToEvent")[type[iter_m]],"\"",sep="")
      if(type[iter_m] %in% c(2,3)){cat(" , threshold =",threshold[iter_m])}
      if(type[iter_m]==3){
        p.censoringT <- exp(-hypothesis[[iter_m]]["lambda_t"]*hypothesis[[iter_m]]["T_followUp"])*(1-exp(-hypothesis[[iter_m]]["lambda_t"]*hypothesis[[iter_m]]["T_inclusion"]))/(hypothesis[[iter_m]]["lambda_t"]*hypothesis[[iter_m]]["T_inclusion"])
        p.censoringC <- exp(-hypothesis[[iter_m]]["lambda_c"]*hypothesis[[iter_m]]["T_followUp"])*(1-exp(-hypothesis[[iter_m]]["lambda_c"]*hypothesis[[iter_m]]["T_inclusion"]))/(hypothesis[[iter_m]]["lambda_c"]*hypothesis[[iter_m]]["T_inclusion"])
        
        cat(" | expected censoring = ",round(100*p.censoringT,2),"% (treatment) and ",round(100*p.censoringC,2),"% (control)",sep="");
      }
      cat("\n")
    }
    cat("   # management of censored survival pairs : ")
    switch(method,
           "Gehan"=cat("uninformative pairs \n"),
           "Peto"=cat("imputation using one survival curve estimated on all patients \n"),
           "Efron"=cat("imputation using different survival curve for control and treatment patients \n"))
    if(method %in% c("Peto","Efron")){
      cat("   # weights of the pairs relatively to the enpoints : \n")
      print(Wscheme)
      index_TTEM1_display <-  rbind(paste("[",round(c(+Inf,threshold[type==3])[index_TTEM1+1],4)," ; ",round(threshold[type==3][index_TTEM1+1],4),"] ",sep=""))
      colnames(index_TTEM1_display) <- endpoint[type==3]      
      rownames(index_TTEM1_display) <- "threshold interval"
      cat("   # intervals thresholds for survival endpoints : \n")
      print(index_TTEM1_display)
    }
    cat("   # n.simul = ",n.simul," | n.bootstrap = ",n.bootstrap," \n",sep="")    
    if(n.bootstrap>0){
      cat("   # cpus = ",cpus,sep="")
      if(!is.null(seed)){
        cat(" (seeds : ",paste(seq(seed,seed+cpus-1),collapse=" "),")",sep="")       
      }
      cat("\n")
    }
    
    
  }
  
  #   trace.simul <- trace
  #   trace <- FALSE
  
  #### boot_ci function #################################################################################### 
  
  fctText <- "function(x,envir){ \n"
  
  fctText <- paste(fctText,
                   "groupT <- rbinom(envir$n.Treatment[envir$n.n],size=1,prob=envir$prob.alloc) 
                   groupC <- rbinom(envir$n.Control[envir$n.n],size=1,prob=1-envir$prob.alloc)
                   indexT.T <- which(groupT==1)
                   indexT.C <- which(groupC==1)
                   indexC.T <- which(groupT==0)
                   indexC.C <- which(groupC==0) \n \n",sep="")
  
  fctText <- paste(fctText,
                   "delta_boot <- ",if(method %in% c("Peto","Efron")){"BuyseTest::pairSampleMainStrataKM_cpp"}else{"BuyseTest::pairSampleMainStrata_cpp"},"(
                   Treatment=rbind(envir$M.Treatment[indexT.T,,drop=FALSE],envir$M.Control[indexT.C,,drop=FALSE]),
                   Control=rbind(envir$M.Treatment[indexC.T,,drop=FALSE],envir$M.Control[indexC.C,,drop=FALSE]),
                   threshold=envir$threshold,type=envir$type,
                   delta_Treatment=",if(D.TTE>0){"rbind(envir$M.delta_Treatment[indexT.T,,drop=FALSE],envir$M.delta_Control[indexT.C,,drop=FALSE])"}else{"matrix(-1,1,1)"},",
                   delta_Control=",if(D.TTE>0){"rbind(envir$M.delta_Treatment[indexC.T,,drop=FALSE],envir$M.delta_Control[indexC.C,,drop=FALSE])"}else{"matrix(-1,1,1)"},",
                   D=envir$D,n_TTE=envir$D.TTE",if(method %in% c("Peto","Efron")){",Wscheme=envir$Wscheme,
                                                           list_survivalT=lapply(1:envir$D.TTE,function(x){rbind(envir$list_survivalT[[x]][indexT.T,],envir$list_survivalC[[x]][indexT.C,])}),
                                                           list_survivalC=lapply(1:envir$D.TTE,function(x){rbind(envir$list_survivalT[[x]][indexC.T,],envir$list_survivalC[[x]][indexC.C,])})"},",
                   numT=c(indexT.T-1,indexT.C-1),numC=c(indexC.T-1,indexC.C-1),sizeT=envir$n.Treatment,sizeC=envir$n.Control,n_size=envir$n.n)$delta  \n \n",sep="")
  
  fctText <- paste(fctText,"return(delta_boot)  \n }")  
  
  boot_ci <- eval(parse(text=paste(fctText))) # to evaluate the text which lead to assign the boot_ci function
  
  ########################################################################################################## 
  
  
  #### preparation ####
  index_endpoint <- which(duplicated(endpoint[type==3])==FALSE)
  if(D.TTE>0){
    censoring.rate <- array(NA,dim=c(n.simul,length(index_endpoint),n.n,2),
                            dimnames=list(1:n.simul,endpoint[type==3][index_endpoint],n.Integer,c("treatment","control")))
  }
  
  p.value <- array(NA,dim=c(n.simul,D,n.n),
                   dimnames=list(1:n.simul,endpoint,n.Integer))
  
  delta <- array(NA,dim=c(n.simul,D,n.n),
                 dimnames=list(1:n.simul,endpoint,n.Integer))
  
  DeltaCI <- array(NA,dim=c(n.simul,4,n.n),
                   dimnames=list(1:n.simul,c("Delta","Delta_inf","Delta_sup","n.bootstrap"),n.Integer))
  
  #### iterations ####
  
  wrapper_power <- function(iter_simul,envir_base,title_pb,pb){
    
    ## loading
    n.Treatment <- envir_base$n.Treatment
    n.Control <- envir_base$n.Control
    prob.alloc <- envir_base$prob.alloc
    n.n <- envir_base$n.n
    type <- envir_base$type
    D <- envir_base$D
    D.TTE <- envir_base$D.TTE
    threshold <- envir_base$threshold
    Wscheme <- envir_base$Wscheme
    
    ## storage
    DeltaCI <- matrix(NA,nrow=4,ncol=n.n)
    rownames(DeltaCI) <- c("Delta","Delta_inf","Delta_sup","n.bootstrap")
    p.value <- matrix(NA,nrow=D,ncol=n.n)
    
    #### simulation ####
    res <- BuyseTest::calcSample(n.Treatment=n.Treatment[n.n],n.Control=n.Control[n.n],
                                 hypothesis=envir_base$hypothesis,endpoint=envir_base$endpoint,
                                 type=type,D=D,D.TTE=D.TTE,prob.alloc=prob.alloc)
    
    res <- BuyseTest::initData(data=res$data,Ind.Treatment=res$Ind.Treatment,type=type,endpoint=envir_base$endpoint,censoring=res$censoring,
                               method=envir_base$method,D.TTE=D.TTE,threshold=threshold,Wscheme=Wscheme,
                               trace=FALSE,test=TRUE,callMethod="BuysePower")
    
    M.Treatment <- res$M.Treatment
    M.Control <- res$M.Control
    M.delta_Treatment <- res$M.delta_Treatment
    M.delta_Control <- res$M.delta_Control
    list_survivalT <- res$list_survivalT
    list_survivalC <- res$list_survivalC
    
    ## censoring rate
    if(D.TTE>0){
      censoring.rateT <- sapply(n.Treatment,function(x){1-apply(M.delta_Treatment[1:x,index_endpoint,drop=FALSE],2,mean)})
      censoring.rateC <- sapply(n.Control,function(x){1-apply(M.delta_Control[1:x,index_endpoint,drop=FALSE],2,mean)})
    }else{
      censoring.rateT <- NULL
      censoring.rateC <- NULL
    }
    
    #### poncual estimation #####
    if(envir_base$method %in% c("Peto","Efron")){       
      delta <-   BuyseTest::pairSampleMainStrataKM_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                                       delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                                       D=D,n_TTE=D.TTE,
                                                       Wscheme=Wscheme,list_survivalT=list_survivalT,list_survivalC=list_survivalC,
                                                       numT=0:(n.Treatment[n.n]-1),numC=0:(n.Control[n.n]-1),sizeT=n.Treatment,sizeC=n.Control,n_size=n.n)$delta
    }else{
      delta <-   BuyseTest::pairSampleMainStrata_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                                     delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                                     D=D,n_TTE=D.TTE,
                                                     numT=0:(n.Treatment[n.n]-1),numC=0:(n.Control[n.n]-1),sizeT=n.Treatment,sizeC=n.Control,n_size=n.n)$delta    
    }
    
    
    if(D==1){
      DeltaCI["Delta",] <- delta
    }else{
      DeltaCI["Delta",] <- colSums(delta)
    }
    
    #### boostrap estimation #####
    envir <- environment()
    
    delta_boot <- array(NA,dim=c(envir_base$D,envir_base$n.n,envir_base$n.bootstrap))
    delta_boot[] <- vapply(1:envir_base$n.bootstrap,function(x){envir_base$boot_ci(x,envir=envir)},matrix(1.5,envir_base$D,envir_base$n.n))
    delta_boot <- aperm(delta_boot,c(2,1,3))
    
    for(iter_n in 1:n.n){
      res <- BuyseTest::calcCI(delta=t(delta[,iter_n]),delta_boot=delta_boot[iter_n,,,drop=FALSE],
                               endpoint=endpoint,D=D,alternative=alternative,alpha=alpha,
                               n.bootstrap=n.bootstrap,cpus=cpus,trace=FALSE)
      
      p.value[,iter_n] <-  res$p.value
      
      DeltaCI[c("Delta_inf","Delta_sup","n.bootstrap"),iter_n] <- c(DeltaCI[c("Delta"),iter_n]+res$Delta_quantile[,D],res$n.bootstrap_real)
    }
    
    
    #### display ####
    if(envir_base$trace>0){
      
      if(iter_simul == envir_base$index.trace[1]){ 
        
        if(envir_base$cpus==1){n.tot <- envir_base$n.simul}else{n.tot <- envir_base$nParallel.simul}  
        
        envir_base$label_pb <- paste(round(100*envir_base$index.trace[1]/n.tot),"% done ",sep="")
        tcltk::setTkProgressBar(pb=pb, value=envir_base$index.trace[1]/n.tot, title=paste(title_pb,"(",round(100*envir_base$index.trace[1]/n.tot),"%)",sep=""), label=envir_base$label_pb)
        envir_base$index.trace <- envir_base$index.trace[-1]      
        
      }
      
    }
    
    #### export ####
    return(list(delta=delta,
                DeltaCI=DeltaCI,
                p.value=p.value,
                censoringT=censoring.rateT,
                censoringC=censoring.rateC #,seedBoot=seedBoot                
    ))
  }
  
  calcSimulation(environment())
  
  #### transfomration into BuyseRes object ####
  power <- apply(p.value,c(2,3),function(x){mean(x<alpha,na.rm=TRUE)})
  
  BuyseRes.object <- new("BuyseSample",
                         alpha=alpha,
                         alternative=alternative,
                         method=method,
                         n.bootstrap=n.bootstrap,
                         type=type,
                         threshold=threshold,
                         power=power,
                         delta = delta, 
                         censoring.rate = if(D.TTE>0){censoring.rate}else{array(NA)}, 
                         p.value = p.value,    
                         DeltaCI =DeltaCI
  )
  
  return(BuyseRes.object)
  
}

#### 2- calc functions ####

# calcCI : cf FCT_buyseTest.R

calcSample <- function(n.Treatment,n.Control,hypothesis,endpoint,
                       type,D,D.TTE,prob.alloc){
  
  ## Ind.Treatment
  Ind.Treatment <- c(rep(1,n.Treatment),
                     rep(0,n.Control)
  )
  
  ## censoring
  
  
  ## data
  endpoint.unique <- unique(endpoint)
  D.unique <- length(endpoint.unique)  
  D.TTE.unique <- length(unique(endpoint[type==3]))
  if(D.TTE>0){
    censoring <- paste("censoring",endpoint[type==3],sep="_")
    iter_D.TTE <- 0
  }else{
    censoring <- NULL  
  }
  
  data <- matrix(NA,nrow=n.Treatment+n.Control,ncol=D.unique+1+D.TTE.unique)
  
  
  for(iter_endpoint.unique in 1:D.unique){
    
    iter_endpoint <- which(endpoint == endpoint.unique[iter_endpoint.unique])[1]
    
    if(type[iter_endpoint]==1){
      data[1:n.Treatment,iter_endpoint.unique+1] <- rbinom(n.Treatment,size=1,prob=hypothesis[[iter_endpoint]]["proba_t"])
      data[(n.Treatment+1):(n.Treatment+n.Control),iter_endpoint.unique+1] <- rbinom(n.Treatment,size=1,prob=hypothesis[[iter_endpoint]]["proba_c"])
    }
    if(type[iter_endpoint]==2){
      data[1:n.Treatment,iter_endpoint.unique+1] <- rnorm(n.Treatment,mean=hypothesis[[iter_endpoint]]["mu_t"],sd=hypothesis[[iter_endpoint]]["sigma"])
      data[(n.Treatment+1):(n.Treatment+n.Control),iter_endpoint.unique+1] <- rnorm(n.Treatment,mean=hypothesis[[iter_endpoint]]["mu_c"],sd=hypothesis[[iter_endpoint]]["sigma"])
    }
    if(type[iter_endpoint]==3){
      iter_D.TTE <- iter_D.TTE + 1 
      
      # Treatment arm
      Time.ObservationT <- hypothesis[[iter_endpoint]][["T_followUp"]]+runif(n.Treatment,min=0,max=hypothesis[[iter_endpoint]][["T_inclusion"]])
      Time.EventT <- rexp(n.Treatment,rate=hypothesis[[iter_endpoint]][["lambda_t"]])
      
      index_censored <- which(Time.EventT>Time.ObservationT)
      index_event <- which(Time.EventT<=Time.ObservationT)
      if(length(index_censored)>0){
        data[1:n.Treatment,D.unique+iter_D.TTE+1][index_censored] <- 0
        data[1:n.Treatment,iter_endpoint.unique+1][index_censored] <- Time.ObservationT[index_censored]
      }
      if(length(index_event)>0){
        data[1:n.Treatment,D.unique+iter_D.TTE+1][index_event] <- 1
        data[1:n.Treatment,iter_endpoint.unique+1][index_event] <- Time.EventT[index_event]
      }
      
      # Control arm
      Time.ObservationC <- hypothesis[[iter_endpoint]][["T_followUp"]]+runif(n.Treatment,min=0,max=hypothesis[[iter_endpoint]][["T_inclusion"]])
      Time.EventC <- rexp(n.Treatment,rate=hypothesis[[iter_endpoint]][["lambda_c"]])
      
      index_censored <- which(Time.EventC>Time.ObservationC)
      index_event <- which(Time.EventC<=Time.ObservationC)
      if(length(index_censored)>0){
        data[(n.Treatment+1):(n.Treatment+n.Control),D.unique+iter_D.TTE+1][index_censored] <- 0
        data[(n.Treatment+1):(n.Treatment+n.Control),iter_endpoint.unique+1][index_censored] <- Time.ObservationC[index_censored]
      }
      if(length(index_event)>0){
        data[(n.Treatment+1):(n.Treatment+n.Control),D.unique+iter_D.TTE+1][index_event] <- 1
        data[(n.Treatment+1):(n.Treatment+n.Control),iter_endpoint.unique+1][index_event] <- Time.EventC[index_event]
      }
    }
    
  }
  
  data <- as.data.frame(data)
  names(data) <- c("treatment",unique(endpoint),unique(censoring))
  data$treatment <- Ind.Treatment
  
  #### export ####
  res <- list()
  res$data <- data
  res$censoring <- censoring
  res$Ind.Treatment <- Ind.Treatment
  
  return(res)
}

calcSimulation <- function(envir){
  
  if(envir$cpus==1){ ## sequential boostrap
    if(envir$trace>1){cat("Sequential simulation \n")}
    if(envir$trace>0){
      envir$index.trace <- unique(round(seq(1,envir$n.simul,length.out=101)[-1])) # number of iterations corresponding to each percentage of progress in the simulation computation
      envir$cpu_name <- "cpu 1 : " # 1 seul cpu est utilise
      envir$title_pb <- paste(envir$cpu_name,"simulation iterations",sep="") # title of the bar displaying the progress
      envir$label_pb <- paste(envir$cpu_name,"0% done ",sep="") # inital legend of the bar displaying the progress
      envir$pb <- tcltk::tkProgressBar(envir$title_pb,envir$label_pb,0, 1, 0) # create and display the progress bar   
    }
    
    # envir$seedBoot <- rep(NA,envir$n.simul)
    
    # perform the bootstrap and return for each bootstrap sample a n.strata*D matrix
    if(!is.null(envir$seed)){set.seed(envir$seed)}
    
    res <- sapply(1:envir$n.simul,function(x){
      res <- envir$wrapper_power(x,envir=envir,title_pb=envir$title_pb,pb=envir$pb)
      
      envir$delta[x,,] <- res$delta
      envir$DeltaCI[x,,] <- res$DeltaCI
      envir$p.value[x,,] <-  res$p.value
      if(envir$D.TTE>0){
        envir$censoring.rate[x,,,"treatment"] <- res$censoringT
        envir$censoring.rate[x,,,"control"] <- res$censoringC
      }
    })
    
    #print(envir$seedBoot)
    if(envir$trace>0){close(envir$pb)}  # close the bar
    
  }else{ ## parallel boostrap
    if(envir$trace>1){cat("Parallel boostrap \n")}
    envir$nParallel.simul <- round(envir$n.simul/envir$cpus) # number of simualtion sample by cpu
    
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
      
      # store results
      envirCPU <- environment()
      if(envir$D.TTE>0){envirCPU$censoring.rate <- array(NA,dim=c(envir$nParallel.simul,length(envir$index_endpoint),envir$n.n,2))}else{envirCPU$censoring.rate <- NULL}
      envirCPU$p.value <- array(NA,dim=c(envir$nParallel.simul,envir$D,envir$n.n))
      envirCPU$delta <- array(NA,dim=c(envir$nParallel.simul,envir$D,envir$n.n))      
      envirCPU$DeltaCI <- array(NA,dim=c(envir$nParallel.simul,4,envir$n.n))
      
      # apply simul
      res <- sapply(1:envir$nParallel.simul,function(iter){
        res <- envir$wrapper_power(iter,envir=envir,title_pb=title_pb,pb=pb)
        
        envirCPU$delta[iter,,] <- res$delta
        envirCPU$DeltaCI[iter,,] <- res$DeltaCI
        envirCPU$p.value[iter,,] <-  res$p.value
        if(envir$D.TTE>0){
          envirCPU$censoring.rate[iter,,,1] <- res$censoringT
          envirCPU$censoring.rate[iter,,,2] <- res$censoringC
        }
      })
      
      # export
      return(list(delta=envirCPU$delta,
                  DeltaCI=envirCPU$DeltaCI,
                  p.value=envirCPU$p.value,
                  censoring.rate=envirCPU$censoring.rate #,seedBoot=envirCPU$seedBoot
      ))
      
    }
    
    ## init display    
    pid <- snowfall::sfClusterEval(Sys.getpid()) # identifier of each R session
    snowfall::sfExport("pid")
    invisible(snowfall::sfClusterEval(cpu_index <- which(pid %in% Sys.getpid()))) # identify to index of the session of each cpu 
    if(envir$trace>0){
      envir$index.trace <- unique(round(seq(1,envir$nParallel.simul,length.out=101)[-1])) # number of iterations corresponding to each percentage of progress in the bootstrap computation
      snowfall::sfLibrary("tcltk", character.only=TRUE) # export tcltk library (required for the progress bar) to each cpu
      invisible(snowfall::sfClusterEval(cpu_name <- paste("cpu ",cpu_index," : ",sep=""))) # name of each cpu session
      invisible(snowfall::sfClusterEval(title_pb <- paste(cpu_name,"simulation iterations",sep=""))) # affect the title the progress bar to each cpu 
      invisible(snowfall::sfClusterEval(label_pb <- paste(cpu_name," 0% done",sep=""))) # affect the initial legend of the progress bar to each cpu 
      invisible(snowfall::sfClusterEval(pb <- tcltk::tkProgressBar(title_pb,label_pb,0, 1, 0))) # affect and display the progress bar to each cpu 
    }
    
    ## export variables 
    snowfall::sfExport("envir")
    
    ## bootstrap
    if(!is.null(envir$seed)){
      invisible(snowfall::sfClusterEval( set.seed(envir$seed+cpu_index-1) ) )  
    }
    resIter <- snowfall::sfClusterApply(rep(NA,envir$cpus),wrapper) # call the cpu to perform bootstrap computations
    
    if(envir$trace==3){ # free the cpus
      snowfall::sfStop()
    }else{
      suppressMessages(snowfall::sfStop())
    }
    
    # merge results
    sapply(1:envir$cpus,function(iter_CPU){
      index_iterations <- seq((iter_CPU-1)*envir$nParallel.simul+1,iter_CPU*envir$nParallel.simul)
      if(envir$D.TTE>0){envir$censoring.rate[index_iterations,,,] <- resIter[[iter_CPU]]$censoring.rate}
      envir$delta[index_iterations,,] <- resIter[[iter_CPU]]$delta
      envir$DeltaCI[index_iterations,,] <- resIter[[iter_CPU]]$DeltaCI
      envir$p.value[index_iterations,,] <- resIter[[iter_CPU]]$p.value
      return(1)
    })
  }
  
}