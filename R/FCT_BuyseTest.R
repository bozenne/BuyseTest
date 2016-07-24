#' @name BuyseTest
#' @title Generalized Pairwise Comparisons
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event outcomes.
#' @param data A \code{data.frame} containing the variables.
#' @param treatment the name of the treatment variable identifying the control and the experimental group. \emph{character}.
#' @param endpoint the name of the endpoint variable(s). \emph{character vector}.
#' @param threshold the thresholds, one for each endpoint variable. \emph{numeric vector}. Default is \code{NULL} indicating no threshold.
#' @param the name of the strata variable(s). \emph{numeric vector}. Default is \code{NULL} indicating only one strata.
#' @param censoring the name of the censoring variable(s), one for each endpoint. \emph{character vector}. Default is \code{NULL}.
#' @param type the type of each endpoint. \emph{character vector}. Can be \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.
#' @param method Is defined when at least one time-to-event outcome is analyzed. Defines the method used to handle pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations.  Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}. See details. 
#' @param n.bootstrap the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}. Default is \code{0} meaning no bootstrap (and thus only ponctual estimation).
#' @param prob.alloc the resampling probability for assignement to the experimental group in the bootstrap samples. \emph{double}. Default is \code{NULL} indicating to use the proportion of patients in the experimental group.
#' @param alternative a \emph{character} specifying the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. Default is \code{"two.sided"}.
#' @param seed the seed to consider for the bootstrap. \emph{integer}. Default is \code{10}.
#' @param cpus the number of CPU to use. \emph{integer}. Default is \code{1}.
#' @param trace Should the execution of the function be traced ? \emph{integer}. Default is \code{3}.
#' 
#' @details 
#' \bold{Treatment:} The variable corresponding to \code{treatment} in data must have only two levels (e.g. \code{0} and \code{1}). \cr
#' \bold{Endpoint, threshold, censoring, and type:} Arguments \code{endpoint}, \code{threshold}, \code{censoring}  and \code{type} must have the same length. \cr
#' \code{threshold} must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. \cr
#' \code{censoring} must be \code{NA} for binary or continuous endpoints and indicate a variable in data for time to event endpoints. 
#' Short forms for endpoint \code{type} are \code{"bin"} (binary endpoint), \code{"cont"} (continuous endpoint), \code{"TTE"} (time-to-event endpoint). 
#' 
#' \bold{Bootstrap:} The number of bootstrap replications (argument \code{n.bootstrap}) must be specified to enable the computation of the confidence intervals and the p.value. 
#' A large number of bootstrap samples (e.g. \code{n.bootstrap=10000})  are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 
#' 
#' \bold{Trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the bootstrap,  and \code{0} remains silent.
#' 
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus. The parallelization relies on the \emph{snowfall} package (function \emph{sfClusterApplyLB}). The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevend from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.
#' 
#' \bold{Method:} Pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}). Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions. \code{method="Peto"} estimate these probabilities using the common Kaplan-Meier estimator of the survival function for treated and control patients. \code{method="Efron"}, and \code{method="Peron"} estimate these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. When the largest observation is censored, it is not possible to estimate the survival probability by the Kaplan-Meier estimator beyond this time point.  \code{method="Efron"} treats the largest observations in each patient group as if it were uncensored. \code{method="Peron"} treats the probability of survival beyond the last observation as NA, resulting in a non null probability that the pair is uninformative    
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references 
#' Marc Buyse (2010) Generalized pairwise comparisons of prioritized endpoints in the two-sample problem. \emph{Statistics in Medicine} \bold{vol. 29} 3245-3257 \cr
#' Efron B (1967) The two sample problem with censored data \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} \bold{vol. 4} 831-583 \cr
#' Peto R, Peto J (1972) Asymptotically efficient rank invariant test procedures \emph{J R Stat Soc A} \bold{vol. 135(2)} 185-198 \cr
#' Gehan EA (1965) A generalized two-sample Wilcoxon test for doubly censored data \emph{Biometrika} \bold{vol. 52(3)} 650-653 \cr
#'
#' @seealso 
#' \code{\link{summary,BuyseRes-method}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' 
#' @examples 
#'    
#'   #### real example : Veteran dataset of the survival package ####
#'   #### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#'   #### method = "Gehan"
#'   \dontrun{
#'     data(veteran,package="survival")
#'     require(BuyseTest)
#'     BuyseTest_veteran_Gehan <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
#'                                          strata="celltype",type="timeToEvent",censoring="status",threshold=0,
#'                                          n.bootstrap=1000,method="Gehan",cpus="all")
#'     
#'     summary_veteran_Gehan <- summary(BuyseTest_veteran_Gehan)
#'     
#'     #### method = "Peron"
#'     
#'     BuyseTest_veteran_Peron <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
#'                                          strata="celltype",type="timeToEvent",censoring="status",threshold=0,
#'                                          n.bootstrap=1000,method="Peron",cpus="all")
#'     
#'     summary_veteran_Peron <- summary(BuyseTest_veteran_Peron)
#'   } 
#'     
#'     
#'     #### Several endpoints :
#'     #######Survival, a time-to-event endpoint
#'     #######Toxicity, a continuous/ordinal endpoint : 6 grades of maximal adverse event 
#'     
#'     n.Treatment <- 100
#'     n.Control<- 100
#'     prob.Treatment_TOX <- c(0.5,0.25,0.10,0.075,0.05,0.025)
#'     prob.Control_TOX <- c(0.7,0.15,0.05,0.05,0.025,0.025)
#'     
#'     lambda.Treatment_TTE <- 0.6
#'     lambda.Control_TTE <- 1
#'     
#'     
#'     set.seed(10)
#'     data_test <- data.frame(treatment=c(rep(1,n.Treatment),
#'                                         rep(0,n.Control)  ))
#'     data_test$toxicity <- c(apply(rmultinom(n.Treatment,size=1,
#'                                             prob=prob.Treatment_TOX)==1,2,which),
#'                             apply(rmultinom(n.Control,size=1,
#'                                             prob=prob.Control_TOX)==1,2,which))
#'     
#'     data_test$toxicityInv <-6-data_test$toxicity
#'     
#'     data_test$EventTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
#'                              rexp(n.Control,rate=lambda.Control_TTE))
#'     data_test$CensoringTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
#'                                  rexp(n.Control,rate=lambda.Control_TTE))
#'     data_test$CensoringTime[data_test$CensoringTime>4] <- 4
#'     
#'     data_test$Survival <- apply(data_test[,c("EventTime","CensoringTime")],1,min)
#'     data_test$event <- as.numeric(apply(data_test[,c("EventTime","CensoringTime")],
#'                                         1,which.min)==1)
#'     
#'     resKM_tempo <- survfit(Surv(data_test[,"Survival"],data_test[,"event"])~data_test$treatment)
#'     plot(resKM_tempo)
#'     
#'     #### method = "Gehan". 
#'     
#'     BuyseTest_severalendpoint_Gehan <- BuyseTest(data=data_test,method="Gehan",
#'                                                  endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
#'                                                  treatment="treatment",
#'                                                  censoring=c("event",NA,"event1",NA,"event",NA),
#'                                                  type=c("TTE","cont","TTE","cont","TTE","cont"),
#'                                                  threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
#'     summary(BuyseTest_severalendpoint_Gehan)
#'     
#'     #### method = "Peron". 
#'     
#'     BuyseTest_severalendpoint_Peron <- BuyseTest(data=data_test,method="Peron",
#'                                                  endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
#'                                                  treatment="treatment",
#'                                                  censoring=c("event",NA,"event",NA,"event",NA),
#'                                                  type=c("TTE","cont","TTE","cont","TTE","cont"),
#'                                                  threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
#'     summary(BuyseTest_severalendpoint_Peron)
#'     
#' @keywords function BuyseTest
#' @export
BuyseTest <- function(data,treatment,endpoint,type,threshold=NULL,strata=NULL,censoring=NULL,method="Peron",
                      n.bootstrap=0,prob.alloc=NULL, alternative = "two.sided",seed=10,cpus=1,trace=3){
  
  
  #### 1- data management + tests ####
  Buysecall <- match.call()
  
  ## Treatment: extract the 2 levels
  validCharacter(treatment, validLength = 1, method = "BuyseTest")
  validNames(data, requiredValues = treatment, validLength = NULL, method = "BuyseTest")
  levels.treatment <- levels(as.factor(data[[treatment]])) # extraction of the levels of the treatment variable
  if(length(levels.treatment) != 2){
    stop("BuyseTest : wrong specification of \'treatment\' \n",
         "the corresponding column in \'data\' must have exactly 2 levels \n",
         "proposed levels : ",paste(levels.treatment,collapse=" "),"\n")
  }
  
  ## endpoint
  validNames(data, requiredValues = endpoint, validLength = NULL, method = "BuyseTest")
  D <- length(endpoint) # number of endpoints
  
  ## type: convert type to numeric and count the number of endpoints
  validCharacter(type, validValues = c("bin","binary","cont","continuous","TTE","timeToEvent"), validLength = D, method = "BuyseTest")
  type[type %in% c("binary","bin")] <- "1" 
  type[type %in% c("continuous","cont")] <- "2"
  type[type %in% c("timeToEvent","TTE")] <- "3"
  type <- as.numeric(type) # type is an integer equal to 1 (binary endpoint), 2 (continuous endpoint) or 3 (time to event endpoint)
  
  D.TTE <- sum(type==3) # number of time to event endpoints
  
  ## censoring: (see .Rd internal-intilisation, section details for details)
  censoring <- initCensoring(censoring=censoring,endpoint=endpoint,type=type,D=D,D.TTE=D.TTE,
                             treatment=treatment,strata=strata)
  validNames(data, requiredValues = censoring, validLength = NULL, refuse.NULL = FALSE, method = "BuyseTest")
  
  ## data: split the data according to the two levels
  if(!is.null(strata)){
  validNames(data, requiredValues = strata, validLength = NULL, method = "BuyseTest")
  }
  
  if(is.data.table(data)){
    dataT <- data[data[[treatment]]==levels.treatment[1], c(endpoint,strata,censoring),with=FALSE]
    dataC <- data[data[[treatment]]==levels.treatment[2], c(endpoint,strata,censoring),with=FALSE]
  }else if(is.data.frame(data)){
    dataT <- as.data.table(data[data[[treatment]]==levels.treatment[1],c(endpoint,strata,censoring),drop=FALSE])
    dataC <- as.data.table(data[data[[treatment]]==levels.treatment[2],c(endpoint,strata,censoring),drop=FALSE])
  }else{
    stop("BuyseTest : data must be a data.frame or a data.table \n")
  }
  n.Treatment <- NROW(dataT) # number of patient in the treatment arm
  n.Control <- NROW(dataC) # number of patient in the control arm
  
  ## threshold: (see .Rd internal-intilisation, section details for details)
  threshold <- initThreshold(threshold=threshold,type=type,D=D,
                             method=method,endpoint=endpoint)
  
  ## strata: (see .Rd internal-intilisation, section details for details)
  res <- initStrata(strata=strata,
                    dataT=dataT,dataC=dataC,n.Treatment=n.Treatment,n.Control=n.Control,
                    endpoint=endpoint,censoring=censoring)
  
  index.strataT <- res$index.strataT 
  index.strataC <- res$index.strataC 
  n.strata <- res$n.strata 
  levels.strata <- res$levels.strata 
 
  ## method
  validCharacter(method, validLength = 1, validValues = c("Gehan","Peto","Efron","Peron"), method = "BuyseTest")
  if(D.TTE==0){
    method <- "Gehan"
    if("method" %in% names(Buysecall) && trace>0){
      message("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
    }
  }
  
  ## alternative
  validCharacter(alternative, validLength = 1, validValues = c("two.sided", "less", "greater"), method = "BuyseTest")
  
  ## n.bootstrap
  validInteger(n.bootstrap, validLength = 1, min = 0, method = "BuyseTest")
  
  ## proba
  if(is.null(prob.alloc)){ # if prob.alloc is not set by the user
    prob.alloc <- n.Treatment/(n.Treatment+n.Control)  # it is set to the proportion of patients in the treatment arm.
  }else{
    validNumeric(prob.alloc, validLength = 1, min = 0, max = 1, method = "BuyseTest")
  }
  
  ## seed
  validInteger(seed, validLength = 1, refuse.NULL = FALSE, min = 1, method = "BuyseTest")
  
  ## cpu
  if(cpus!=1){
    #### test package
    test.package <- requireNamespace("snow",quietly=TRUE)
    if("package:snow" %in% search()==FALSE){
      test.package <- try(attachNamespace("snow"), silent = TRUE)
      if("package:snow" %in% search()==FALSE){
        stop("BuyseTest : this function with argument cpus>1 requires to have installed the snow package to work \n")
      }
    }

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
  
  if(cpus=="all"){ 
    cpus <- parallel::detectCores() # this function detect the number of CPU cores 
  }else{  # if several cpus are intended to be used, check this correspond to a valid number of CPU cores
    validInteger(cpus, validLength = 1, validValues = 1:parallel::detectCores(), method = "BuyseTest")
  }
 
  ## data
  res <- initData(dataT=dataT, dataC=dataC,type=type,endpoint=endpoint,D=D,censoring=censoring,
                  index.strataT=index.strataT,index.strataC=index.strataC,n.strata=n.strata,                  
                  method=method,D.TTE=D.TTE,threshold=threshold,Wscheme=NULL,
                  test=TRUE,trace=trace)
  
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
    cat("   # chosen reference : Control = ",levels.treatment[1]," and Treatment = ",levels.treatment[2],"\n")
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
  
  #### 2- Punctual estimation ####
  if(trace>1){cat("Punctual estimation \n")}
  
  if(method %in% c("Peto","Efron","Peron")){  
    resPonctual <-   BuyseTest_PetoEfronPeron_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                                  delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                                  D=D,returnIndex=TRUE,
                                                  strataT=index.strataT,strataC=index.strataC,n_strata=n.strata,n_TTE=D.TTE,
                                                  Wscheme=Wscheme,index_survivalM1=index_survivalM1,threshold_TTEM1=threshold_TTEM1,list_survivalT=list_survivalT,list_survivalC=list_survivalC,
                                                  PEP=which(c("Peto","Efron","Peron") == method)
    )
  }else if(method=="Gehan"){
    resPonctual <-   BuyseTest_Gehan_cpp(Treatment=M.Treatment,Control=M.Control,threshold=threshold,type=type,
                                         delta_Treatment=M.delta_Treatment,delta_Control=M.delta_Control,
                                         D=D,returnIndex=TRUE,
                                         strataT=index.strataT,strataC=index.strataC,n_strata=n.strata,n_TTE=D.TTE)    
  }
  
  #### transfomration into BuyseRes object ####
  BuyseRes.object <- BuyseRes(
    delta = resPonctual$delta, 
    count_favorable = resPonctual$count_favorable,      
    count_unfavorable = resPonctual$count_unfavorable,
    count_neutral = resPonctual$count_neutral,    
    count_uninf = resPonctual$count_uninf,
    index_neutralT = resPonctual$index_neutralT,
    index_neutralC = resPonctual$index_neutralC,
    index_uninfT = resPonctual$index_uninfT,
    index_uninfC = resPonctual$index_uninfC,
    n_pairs = resPonctual$n_pairs,
    delta_boot = array(NA,dim = c(n.strata,D,1)), 
    p.value = rep(NA,D),    
    Delta_quantile = matrix(NA,nrow = 2, ncol = D, dimnames = list(c("2.5%","97.5%"))),
    endpoint = endpoint,
    threshold = threshold,
    strata = levels.strata
  )
  
  if(n.bootstrap==0){
    if(trace>0){cat("*** only ponctual estimation requested ***\n",
                    "set \'n.bootstrap\' argument to a strictly postive value for confidence interval and p.value \n"
    )}
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
                   if(method %in% c("Efron","Peron")){"res_init <- initSurvival(M.Treatment=Mnew.Treatment,M.Control=Mnew.Control,
                                                    M.delta_Treatment=Mnew.delta_Treatment,M.delta_Control=Mnew.delta_Control,
                                                    endpoint=envir$endpoint,D.TTE=envir$D.TTE,type=envir$type,threshold=envir$threshold,
                                                    index.strataT=new.strataT,index.strataC=new.strataC,n.strata=envir$n.strata,
                                                    method=envir$method,callMethod=\"delta_boot\") \n \n"},
                   sep="")
  
  fctText <- paste(fctText,"delta_boot <- ",
                   if(method=="Gehan"){"BuyseTest_Gehan_cpp("}else{"BuyseTest_PetoEfronPeron_cpp("},
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
  
  #### update BuyseRes object ####
  BuyseRes.object@delta_boot <- delta_boot
  BuyseRes.object@p.value <- p.value
  BuyseRes.object@Delta_quantile <- Delta_quantile
   
  #### export ####
  return(BuyseRes.object)
}