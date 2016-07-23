#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% fonctions implementant le test de (Buyse, 2010) sous R
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BuyseTest <- function(data,treatment,endpoint,threshold=NULL,strata=NULL,censoring=NULL,type, method="Peron",
                      n.bootstrap=0,prob.alloc=NULL, alternative = "two.sided",seed=10,cpus=1,trace=3){
  
  
  #### preliminary tests ####
  Buysecall <- match.call()
  
  ## type
  type <- initType(type=type,callMethod="BuyseTest")
  D <- length(type) # number of endpoints
  D.TTE <- sum(type==3) # number of time to event endpoints
  
  if(D.TTE==0){ 
    if("method" %in% names(Buysecall) && trace>0){
      cat("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
    }
    method <- "Gehan"
  }
  
  ## endpoint
  initEndpoint(endpoint=endpoint,data=data,D=D,
               treatment=treatment,censoring=censoring,strata=strata,callMethod="BuyseTest")
  
  ## censoring
  censoring <- initCensoring(censoring=censoring,data=data,endpoint=endpoint,type=type,D=D,D.TTE=D.TTE,
                             treatment=treatment,strata=strata,callMethod="BuyseTest")
  
  ## treatment
  res <- initTreatment(data=data,treatment=treatment,
                       endpoint=endpoint,censoring=censoring,strata=strata,
                       callMethod="BuyseTest")
  
  Ind.Treatment <- res$Ind.Treatment
  reference <- res$reference
  n.Treatment <- res$n.Treatment
  n.Control <- res$n.Control
  
  ## strata
  res <- initStrata(strata=strata,data=data,Ind.Treatment=Ind.Treatment,n.Treatment=n.Treatment,n.Control=n.Control,
                    endpoint=endpoint,treatment=treatment,censoring=censoring,callMethod="BuyseTest")
  
  strataT <- res$strataT 
  strataC <- res$strataC 
  index.strataT <- res$index.strataT 
  index.strataC <- res$index.strataC 
  n.strata <- res$n.strata 
  levels.strata <- res$levels.strata 
  
  ## threshold
  threshold <- initThreshold(threshold=threshold,type=type,D=D,
                             method=method,endpoint=endpoint,
                             callMethod="BuyseTest")
  
  ## miscellaneous
  if(alternative %in% c("two.sided", "less", "greater")==FALSE){
    stop("BuyseTest : incorrect specification of \'alternative\' \n",
         "valid \'alternative\' : \"two.sided\" \"less\" \"greater\" \n",
         "proposed \'alternative\' : ",alternative,"\n")
  } 
  
  if(!is.numeric(n.bootstrap) || n.bootstrap < 0 || (n.bootstrap %% 1 !=0)){
    stop("BuyseTest : wrong specification of \'n.bootstrap\' \n",
         "\'n.bootstrap\' must be a positive integer \n",
         "proposed \'n.bootstrap\' : ",n.bootstrap,"\n")
  }
  
  ## proba
  prob.alloc <- initAlloc(prob.alloc,n.Treatment=n.Treatment,n.Control=n.Control,callMethod="BuyseTest")
  
  ## seed
  if(!is.null(seed) && (!is.numeric(seed) || seed <= 0 || (seed %% 1 !=0))){
    stop("BuyseTest : wrong specification of \'seed\' \n",
         "\'seed\' must be a stricly positive integer \n",
         "proposed \'seed\' : ",seed,"\n")
  }
  
  # cpu
  if(cpus!=1){
    #### test package

    test.package <- requireNamespace("snow",quietly=TRUE)
    if("package:snow" %in% search()==FALSE){
      test.package <- try(attachNamespace("snow"), silent = TRUE)
      if("package:snow" %in% search()==FALSE){
        stop("BuyseTest : this function with argument cpus>1 requires to have installed the snow package to work \n")
      }
    }
    # Error in attachNamespace(\"snow\") : namespace is already attached
    # Error in loadNamespace(name) : there is no package called "snow"\n"
    
    if("package:parallel" %in% search()==FALSE){
      test.package <- try(attachNamespace("parallel"), silent = TRUE)
      if("package:parallel" %in% search()==FALSE){
        stop("BuyseTest : this function with argument cpus>1 requires to have installed the parallel package to work \n")
      }
    }
    
    test.package <- requireNamespace("snowfall",quietly=TRUE)
    if(test.package==FALSE){
      stop("BuyseTest : this function with argument cpus>1 requires to have installed the snowfall package to work \n")
    }
    
    test.package <- requireNamespace("parallel",quietly=TRUE)
    if(test.package==FALSE){
      stop("BuyseTest : this function with argument cpus>1 requires to have installed the parallel package to work \n")
    }
    
    test.package <- requireNamespace("tcltk",quietly=TRUE)
    if(test.package==FALSE){
      stop("BuyseTest : this function with argument cpus>1 and trace>0 requires to have installed the tcltk package to work \n")
    }    
    
  }
  
  cpus <- initCPU(cpus=cpus,callMethod="BuyseTest")
  
  #### dataset ####

  res <- initData(data=data,Ind.Treatment=Ind.Treatment,type=type,endpoint=endpoint,D=D,censoring=censoring,
                  index.strataT=index.strataT,index.strataC=index.strataC,n.strata=n.strata,                  
                  method=method,D.TTE=D.TTE,threshold=threshold,Wscheme=NULL,
                  test=TRUE,trace=trace,callMethod="BuyseTest")
  
  M.Treatment <- res$M.Treatment
  M.Control <- res$M.Control 
  M.delta_Treatment <- res$M.delta_Treatment
  M.delta_Control <- res$M.delta_Control
  Wscheme <- res$Wscheme
  threshold_TTEM1 <- res$threshold_TTEM1
  index_survivalM1 <- res$index_survivalM1
  list_survivalT <- res$list_survivalT
  list_survivalC <- res$list_survivalC
  
  ## display
  if(trace>1){
    cat("Settings \n")
    cat("   # chosen reference : Control = ",reference[1]," and Treatment = ",reference[2],"\n")
    cat("   # number of endpoints : ",D," \n")
    iter_mTTE <- 0
    for(iter_m in 1:D){
      cat("      >  endpoint ",iter_m," = \"",endpoint[iter_m],"\" - type = \"",c("binary","continuous","timeToEvent")[type[iter_m]],"\"",sep="")
      if(type[iter_m] %in% c(2,3)){cat(" | threshold =",threshold[iter_m])}
      if(type[iter_m]==3){iter_mTTE <- iter_mTTE + 1 ; cat(" | censoring = \"",censoring[iter_mTTE],"\"",sep="");}
      cat("\n")
    }
    cat("   # n.strata = ",n.strata," : ",paste(levels.strata,collapse=" "),"\n")
    cat("   # n.bootstrap = ",n.bootstrap," | prob.alloc = ",prob.alloc," \n",sep="")    
    cat("   # management of censored survival pairs : ")
    switch(method,
           "Gehan"=cat("uninformative pairs \n"),
           "Peto"=cat("imputation using one survival curve estimated on all patients \n"),
           "Efron"=cat("imputation using different survival curve for control and treatment patients \n"),
           "Peron"=cat("imputation using different survival curve for control and treatment patients \n")
    ) 
    if(method %in% c("Peto","Efron","Peron")){
  
      cat("   # weights of the pairs relatively to the enpoints : \n")
      print(Wscheme)
 
      cat("   # intervals thresholds for survival endpoints : \n")
      
      if(length(threshold_TTEM1)>0){
      threshold_TTEM1.display <- threshold_TTEM1
      threshold_TTEM1.display[threshold_TTEM1.display<0] <- +Inf
      }else{
      threshold_TTEM1.display <- +Inf
      }
      
      threshold.display <- rbind(sapply(1:D.TTE,function(x){paste(c("[",round(threshold_TTEM1.display[x],4)," ; ",round(threshold[type==3][x],4),"] "),collapse="")}))
      colnames(threshold.display) <- endpoint[type==3]      
      rownames(threshold.display) <- "threshold interval"
      print(threshold.display)
    }
    if(n.bootstrap>0){
      cat("   # cpus = ",cpus,sep="")
      if(!is.null(seed)){
        cat(" (seeds : ",paste(seq(seed,seed+cpus-1),collapse=" "),")",sep="")       
      }
      cat("\n")
    }
  }
  
  #### Punctual estimation ####
  if(trace>1){cat("Punctual estimation \n")}
 
  if(method %in% c("Peto","Efron","Peron")){  
    resPonctual <-   BuyseTest::BuyseTest_PetoEfronPeron_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                                             delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                                             D=D,returnIndex=TRUE,
                                                             strataT=index.strataT,strataC=index.strataC,n_strata=n.strata,n_TTE=D.TTE,
                                                             Wscheme=Wscheme,index_survivalM1=index_survivalM1,threshold_TTEM1=threshold_TTEM1,list_survivalT=list_survivalT,list_survivalC=list_survivalC,
                                                             PEP=which(c("Peto","Efron","Peron") == method)
    )
  }else if(method=="Gehan"){
    resPonctual <-   BuyseTest::BuyseTest_Gehan_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                                    delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                                    D=D,returnIndex=TRUE,
                                                    strataT=index.strataT,strataC=index.strataC,n_strata=n.strata,n_TTE=D.TTE)    
  }
  
  if(n.bootstrap==0){
    if(trace>0){cat("*** only ponctual estimation requested ***\n",
                    "set \'n.bootstrap\' argument to a strictly postive value for confidence interval and p.value \n"
    )}
    #### transfomration into BuyseRes object ####
    BuyseRes.object <- new("BuyseRes",
                           delta = resPonctual$delta, 
                           count_favorable = resPonctual$count_favorable,      
                           count_unfavorable = resPonctual$count_unfavorable,
                           count_neutral = resPonctual$count_neutral,    
                           count_uninf = resPonctual$count_uninf,
                           index_neutralT= resPonctual$index_neutralT,
                           index_neutralC = resPonctual$index_neutralC,
                           index_uninfT= resPonctual$index_uninfT,
                           index_uninfC= resPonctual$index_uninfC,
                           n_pairs = resPonctual$n_pairs,
                           delta_boot = array(NA,dim=c(n.strata,D,1)), 
                           p.value = rep(NA,D),    
                           Delta_quantile = matrix(NA,nrow=2,ncol=D,dimnames=list(c("2.5%","97.5%"))),
                           endpoints=paste(endpoint,threshold,sep="_")
    )
    return(BuyseRes.object)
  }
  
  
  #### boot_ci function #################################################################################### 
  # the function boot_ci will create a bootstrap dataset and apply the Buyse Test and export the obtained delta (proportion in favor of treatment)
  # instead of being directly written it is created as text to allow adapt the function to the situation (presence of survival oucome, trace, method)
  # it avoid the use of if statement inside the function which may slow down computations.
  fctDisplay <- paste("if(x == envir$index.trace[1]){
                     pc_done <- envir$index.trace[1]/",if(cpus==1){"envir$n.bootstrap"}else{"envir$nParallel.bootstrap"},"
                     label_pb <- paste(",if(cpus==1){"envir$"},"cpu_name,round(100*pc_done),\"% done\",sep=\"\")
                     tcltk::setTkProgressBar(pb=",if(cpus==1){"envir$"},"pb, value=pc_done,title=paste(",if(cpus==1){"envir$"},"title_pb,\"(\",round(100*pc_done),\"%)\",sep=\"\"), label=label_pb)
                     envir$index.trace <- envir$index.trace[-1]                 
                      } \n \n",sep="")
  
  fctText <- "function(x,envir){ \n"  # browser cat(x,\" \")
  ## randomisation : new allocation of the treatment and control arm
  # groupT and group C contain the sampled indicator variables of the new allocation : treatment (1) and control (0) respectively for the previously treatment and control patients.
  # indexT.T contains the index of the new treatment patients in the (old) treatment arm
  # indexT.C contains the index of the new treatment patients in the (old) control arm
  # indexC.T contains the index of the new control patients in the (old) treatment arm
  # indexC.C contains the index of the new control patients in the (old) control arm
  fctText <- paste(fctText,
                   "groupT <- rbinom(envir$n.Treatment,size=1,prob=envir$prob.alloc) \n",
                   "groupC <- rbinom(envir$n.Control,size=1,prob=envir$prob.alloc) \n",
                   "indexT.T <- which(groupT==1) \n",
                   "indexT.C <- which(groupC==1) \n",
                   "indexC.T <- which(groupT==0) \n",
                   "indexC.C <- which(groupC==0) \n",
                   "strataT_boot <- c(envir$strataT[indexT.T],envir$strataC[indexT.C]) \n",
                   "strataC_boot <- c(envir$strataT[indexC.T],envir$strataC[indexC.C])  \n",
                   "new.strataT <- lapply(1:envir$n.strata,function(x){which(strataT_boot==x)-1}) \n",
                   "new.strataC <- lapply(1:envir$n.strata,function(x){which(strataC_boot==x)-1}) \n \n",
                   sep="")
  
  ## sample with no treatment or control in one strata are ignored
  fctText <- paste(fctText,
                   "if(any(unlist(lapply(new.strataT,length))==0) || any(unlist(lapply(new.strataC,length))==0) ){ \n",
                   if(trace>0){fctDisplay},
                   "return(matrix(NA,nrow=envir$n.strata,ncol=envir$D)) \n",
                   "} \n \n")
  
  ## call the C++ pairCompMainStrata(KM)_cpp for performing the test and select delta from the output
  # the new Treatment matrix is composed of the endpoints of the new treatements in the (old) treatment arm and of the endpoints of the new treatements in the (old) control arm 
  # the new Control matrix is composed of the endpoints of the new controls in the (old) treatment arm and of the endpoints of the new controls in the (old) control arm 
  # threshold and type have not changed
  # the new delta_Treatment matrix is composed of the event indicator of the new treatements in the (old) treatment arm and of the event indicator of the new treatements in the (old) control arm 
  # the new delta_Control matrix is composed of the event indicator of the new controls in the (old) treatment arm and of the event indicator of the new controls in the (old) control arm 
  # D has not changed. 
  # there is no need to return the index of the neutral and uninformative pairs
  # the new strata index for the treatments contains for each strata the index corresponding to the strata variable minus 1 (for C++ compatibility, vector begin at 0)
  # the new strata index for the control contains for each strata the index corresponding to the strata variable minus 1 (for C++ compatibility, vector begin at 0)
  # n_strata and D.TTE have not changed
  # if KM imputation is requested Wscheme must be specified as well as list_survivalT and list_survivalC
  # list_survivalT contains for each survival endpoint a matrix of the survival times of the new treatment patients 
  # list_survivalC contains for each survival endpoint a matrix of the survival times of the new control patients 
  
  fctText <- paste(fctText,
                   if(method %in% c("Efron","Peron")){"Mnew.Treatment <- rbind(envir$M.Treatment[indexT.T,,drop=FALSE],envir$M.Control[indexT.C,,drop=FALSE]) \n"},
                   if(method %in% c("Efron","Peron")){"Mnew.Control <- rbind(envir$M.Treatment[indexC.T,,drop=FALSE],envir$M.Control[indexC.C,,drop=FALSE]) \n"},
                   if(method %in% c("Efron","Peron")){"Mnew.delta_Treatment <- rbind(envir$M.delta_Treatment[indexT.T,,drop=FALSE],envir$M.delta_Control[indexT.C,,drop=FALSE]) \n"},
                   if(method %in% c("Efron","Peron")){"Mnew.delta_Control <- rbind(envir$M.delta_Treatment[indexC.T,,drop=FALSE],envir$M.delta_Control[indexC.C,,drop=FALSE]) \n"},  
                   if(method == "Efron"){" 
                            for(iter_strata in 1:n.strata){
                                    Mnewstrata.Treatment <- Mnew.Treatment[new.strataT[[iter_strata]]+1,,drop=FALSE]
                                    Mnewstrata.Control <- Mnew.Control[new.strataC[[iter_strata]]+1,,drop=FALSE]
                                    for(iter_endpointTTE in 1:envir$D.TTE){
                                        indexT_maxCensored <- which(Mnewstrata.Treatment[,which(type==3)[iter_endpointTTE]]==max(Mnewstrata.Treatment[,which(type==3)[iter_endpointTTE]]))
                                        Mnew.delta_Treatment[new.strataT[[iter_strata]][indexT_maxCensored]+1,iter_endpointTTE] <- 1
                                        indexC_maxCensored <- which(Mnewstrata.Control[,which(type==3)[iter_endpointTTE]]==max(Mnewstrata.Control[,which(type==3)[iter_endpointTTE]]))
                                        Mnew.delta_Control[new.strataC[[iter_strata]][indexC_maxCensored]+1,iter_endpointTTE] <- 1
                                      }
                             } \n"},
                   if(method %in% c("Efron","Peron")){"res_init <- BuyseTest::initSurvival(M.Treatment=Mnew.Treatment,M.Control=Mnew.Control,
                                                    M.delta_Treatment=Mnew.delta_Treatment,M.delta_Control=Mnew.delta_Control,
                                                    endpoint=envir$endpoint,D.TTE=envir$D.TTE,type=envir$type,threshold=envir$threshold,
                                                    index.strataT=new.strataT,index.strataC=new.strataC,n.strata=envir$n.strata,
                                                    method=envir$method,callMethod=\"delta_boot\") \n \n"},
                   sep="")
  
  fctText <- paste(fctText,"delta_boot <- ",
                   if(method=="Gehan"){"BuyseTest::BuyseTest_Gehan_cpp("}else{"BuyseTest::BuyseTest_PetoEfronPeron_cpp("},
                   if(method %in% c("Efron","Peron")){"Treatment=Mnew.Treatment,"}else{"Treatment=rbind(envir$M.Treatment[indexT.T,,drop=FALSE],envir$M.Control[indexT.C,,drop=FALSE]), \n"},
                   if(method %in% c("Efron","Peron")){"Control=Mnew.Control,"}else{"Control=rbind(envir$M.Treatment[indexC.T,,drop=FALSE],envir$M.Control[indexC.C,,drop=FALSE]), \n"},
                   "threshold=envir$threshold,type=envir$type, \n",
                   "delta_Treatment=",if(D.TTE==0){"matrix(-1,1,1), \n"}else if(method %in% c("Efron","Peron")){"Mnew.delta_Treatment, \n"}else{"rbind(envir$M.delta_Treatment[indexT.T,,drop=FALSE],envir$M.delta_Control[indexT.C,,drop=FALSE]), \n"},
                   "delta_Control=",if(D.TTE==0){"matrix(-1,1,1), \n"}else if(method %in% c("Efron","Peron")){"Mnew.delta_Control, \n"}else{"rbind(envir$M.delta_Treatment[indexC.T,,drop=FALSE],envir$M.delta_Control[indexC.C,,drop=FALSE]), \n"},
                   "D=envir$D,returnIndex=FALSE, \n",
                   "strataT=new.strataT, \n",
                   "strataC=new.strataC, \n",
                   "n_strata=envir$n.strata,n_TTE=envir$D.TTE",
                   if(method %in% c("Peto","Efron","Peron")){"\n ,Wscheme=envir$Wscheme,index_survivalM1=index_survivalM1,threshold_TTEM1=envir$threshold_TTEM1, \n"},
                   if(method == "Peto"){"list_survivalT=lapply(1:envir$D.TTE,function(x){rbind(envir$list_survivalT[[x]][indexT.T,],envir$list_survivalC[[x]][indexT.C,])}), \n"}else if(method %in% c("Efron","Peron")){"list_survivalT=res_init$list_survivalT, \n"},
                   if(method == "Peto"){"list_survivalC=lapply(1:envir$D.TTE,function(x){rbind(envir$list_survivalT[[x]][indexC.T,],envir$list_survivalC[[x]][indexC.C,])}), \n"}else if(method %in% c("Efron","Peron")){"list_survivalC=res_init$list_survivalC, \n"},
                   if(method != "Gehan"){"PEP=which(c(\"Peto\",\"Efron\",\"Peron\") == method)"},
                   ")$delta  \n \n",sep="")
  
  ## display 
  # if the following percentage of iterations is reached
  # the label on the bar is updated and the progress bar is displayed
  # then the first element of index.trace to test for the next percentage of iterations
  
  fctText <- paste(fctText,if(trace>0){fctDisplay},sep="")
  
  fctText <- paste(fctText,"return(delta_boot)  \n }")
  
  boot_ci <- eval(parse(text=paste(fctText))) # to evaluate the text which lead to assign the boot_ci function

  ########################################################################################################## 
  
  export_names <- c("nParallel.bootstrap","n.Treatment","prob.alloc","strataT","strataC","M.Treatment","M.Control","threshold","type","D","n.strata","D.TTE")
  if(D.TTE>0){export_names <- c(export_names,"M.delta_Treatment","M.delta_Control")}
  if(method %in% c("Peto","Efron","Peron")){export_names <- c(export_names,"Wscheme","list_survivalT","list_survivalC")}
 
  delta_boot <- array(NA,dim=c(n.strata,D,n.bootstrap))
  calcBootstrap(environment())
  
  if(trace>1){cat("Post-Treatment \n")}
  res <- calcCI(delta=resPonctual$delta,delta_boot=delta_boot,endpoint=endpoint,D=D,alternative=alternative,alpha=0.05,
                n.bootstrap=n.bootstrap,cpus=cpus,trace=TRUE)
  p.value <- res$p.value  
  Delta_quantile <-  res$Delta_quantile
  
  #### transfomration into BuyseRes object ####
  BuyseRes.object <- new("BuyseRes",
                         delta = resPonctual$delta, 
                         count_favorable = resPonctual$count_favorable,      
                         count_unfavorable = resPonctual$count_unfavorable,
                         count_neutral = resPonctual$count_neutral,    
                         count_uninf = resPonctual$count_uninf,
                         index_neutralT= resPonctual$index_neutralT,
                         index_neutralC = resPonctual$index_neutralC,
                         index_uninfT= resPonctual$index_uninfT,
                         index_uninfC= resPonctual$index_uninfC,
                         n_pairs = resPonctual$n_pairs,
                         delta_boot = delta_boot, 
                         p.value = p.value,    
                         Delta_quantile =Delta_quantile,
                         endpoints=paste(endpoint,threshold,sep="_")
  )
  return(BuyseRes.object)
}


#### associated functions ####

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
