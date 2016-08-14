#' @name internal-intilisation
#' @rdname internal-intilisation
#' @title internal functions for BuyseTest - intilisation
#' @aliases initCensoring initData initN initHypothesis initSpace initStrata initSurvival initThreshold initWscheme
#' 
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' 
#' @details
#' All functions performs test to check the validity of the argument. In addition: \cr
#' 
#' \code{initCensoring}: If no TTE and censoring is full of NA, set it to NULL. Else when censoring contains NA for the non-TTE endpoints, remove them. \cr
#' \code{initStrata}: merge the strata into one with the interaction function and compute the index of the observations relative to each strata. \cr
#' \code{initThreshold}: set default threshold to 1e-12 expect for binary variable where it is set to 1. 
#' Also set null threshold to a small positive value (0->1e-12) 
#' - the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12. When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' 
#' 
#' 
#' @keywords function internal BuyseTest

#' @rdname internal-intilisation
initCensoring <- function(censoring,endpoint,type,D,D.TTE,
                          treatment,strata){
  
  if(D.TTE==0){
    
    if(!is.null(censoring) && all(is.na(censoring))){
      censoring <- NULL
    }
    
    if(!is.null(censoring)){
      stop("BuyseTest : \'censoring\' must be NULL when there are no TTE endpoints \n",
           "propose value : ",paste(censoring,collapse=" "),"\n")
    }
    
  }else{
    
    if(length(censoring) == D){
      
      if(any(!is.na(censoring[type!=3])) ){
        stop("BuyseTest : wrong specification of \'censoring\' \n",
             "\'censoring\' must be NA for binary or continuous endpoints \n",
             "binary or continuous endoints : ",paste(endpoint[type!=3],collapse=" "),"\n",
             "proposed \'censoring\' for these endoints : ",paste(censoring[type!=3],collapse=" "),"\n",
        )
      }
      
      if( any(is.na(censoring[type==3])) ){
        stop("BuyseTest : wrong specification of \'censoring\' \n",
             "\'censoring\' must indicate a variable in data for TTE endpoints \n",
             "TTE endoints : ",paste(endpoint[type==3],collapse=" "),"\n",
             "proposed \'censoring\' for these endoints : ",paste(censoring[type==3],collapse=" "),"\n",
        )
      }
      
      censoring <- censoring[type==3] # from now, censoring contains for each time to event endpoint the name of variable indicating censoring (0) or event (1)
    }else if(length(censoring) != D.TTE){
      stop("BuyseTest : \'censoring\' does not match \'endpoint\' size \n",
           "length(censoring) : ",length(censoring),"\n",
           "length(endpoint) : ",D,"\n")
      
    }
  }
  
  #### export #### 
  return(censoring)
}

#' @rdname internal-intilisation
initData <- function(dataT,dataC,type,endpoint,D,censoring,
                     index.strataT,index.strataC,n.strata,          
                     method,D.TTE,threshold,Wscheme=NULL,
                     trace,test=TRUE){
  
  ## check NA
  if(test && any(type==3)){
    test.naT <- colSums(is.na(dataT[,c(censoring,endpoint[type==3]),with=FALSE]))
    test.naC <- colSums(is.na(dataC[,c(censoring,endpoint[type==3]),with=FALSE]))
    test.na <- test.naT+test.naC
    if(any(test.na>0)){ 
      stop("BuyseTest : TTE endpoint / censoring variables must not contain NA \n",
           "number of NA in \'data\' : ",sum(test.na>0),"\n",
           "columns with \'data\' : ",paste(c(endpoint,censoring[type==3])[test.na>0],collapse=" "),"\n")
    }
  }
  
  if(test && trace>0 && any(type!=3)){
    test.naT <- colSums(is.na(dataT[,endpoint[type!=3],with=FALSE]))
    test.naC <- colSums(is.na(dataC[,endpoint[type!=3],with=FALSE]))
    test.na <- test.naT+test.naC
    if(any(test.na>0)){ 
      warning("BuyseTest : some binary or continuous endpoints contains NA \n",
              "number of NA in \'data\' : ",sum(test.na>0),"\n",
              "columns with \'data\' : ",paste(c(endpoint,censoring[type==3])[test.na>0],collapse=" "),"\n")
    }
  }
  
  ## endpoint checking : binary type
  indexY <- which(type==1)
  if(test && length(indexY)>0){
    for(iterY in indexY){
      validNumeric(dataT[[endpoint[iterY]]], name1 = endpoint[iterY], validValues = c(NA,0,1), validLength = NULL, method = "BuyseTest")
      validNumeric(dataC[[endpoint[iterY]]], name1 = endpoint[iterY], validValues = c(NA,0,1), validLength = NULL, method = "BuyseTest")
    }
  }
  
  ## endpoint checking : continuous type
  indexY <- which(type==2)
  if(test && length(indexY)>0){
    for(iterY in indexY){
      validNumeric(dataT[[endpoint[iterY]]], name1 = endpoint[iterY], validLength = NULL, method = "BuyseTest")
      validNumeric(dataC[[endpoint[iterY]]], name1 = endpoint[iterY], validLength = NULL, method = "BuyseTest")
    }
  }
  
  ## endpoint checking : time to event
  indexY <- which(type==3)
  if(test && length(indexY)>0){
    for(iterY in indexY){
      validNumeric(dataT[[endpoint[iterY]]], name1 = endpoint[iterY], validLength = NULL, method = "BuyseTest")
      validNumeric(dataC[[endpoint[iterY]]], name1 = endpoint[iterY], validLength = NULL, method = "BuyseTest")
      validNumeric(dataT[[censoring[which(indexY == iterY)]]], name1 = censoring[which(indexY == iterY)], validValues = c(0,1), validLength = NULL, method = "BuyseTest")
      validNumeric(dataC[[censoring[which(indexY == iterY)]]], name1 = censoring[which(indexY == iterY)], validValues = c(0,1), validLength = NULL, method = "BuyseTest")
    }
  }
  
  ## transformation
  M.Treatment <- as.matrix(dataT[,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
  M.Control <- as.matrix(dataC[,endpoint,with=FALSE]) # matrix of endpoints for the control arm
  
  if(!is.null(censoring)){
    M.delta_Treatment <- as.matrix(dataT[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    M.delta_Control <- as.matrix(dataC[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    
    if(method=="Efron"){ 
      for(iter_strata in 1:n.strata){
        
        Mstrata.Treatment <- M.Treatment[index.strataT[[iter_strata]]+1,,drop=FALSE]
        Mstrata.Control <- M.Control[index.strataC[[iter_strata]]+1,,drop=FALSE]
        
        # set last observation for each TTE endpoint to non-censored
        for(iter_endpointTTE in 1:D.TTE){
          indexT_maxCensored <- which(Mstrata.Treatment[,which(type==3)[iter_endpointTTE]]==max(Mstrata.Treatment[,which(type==3)[iter_endpointTTE]]))
          M.delta_Treatment[index.strataT[[iter_strata]][indexT_maxCensored]+1,iter_endpointTTE] <- 1
          indexC_maxCensored <- which(Mstrata.Control[,which(type==3)[iter_endpointTTE]]==max(Mstrata.Control[,which(type==3)[iter_endpointTTE]]))
          M.delta_Control[index.strataC[[iter_strata]][indexC_maxCensored]+1,iter_endpointTTE] <- 1
        }
      }
    }
    
  }else{ # if the is no time to event variables
    M.delta_Treatment <- matrix(-1,1,1) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    M.delta_Control <- matrix(-1,1,1) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
  }
  
  ####  KM imputation  #### 
  if(method %in% c("Peto","Efron","Peron")){
    
    endpoint.TTE <- endpoint[type==3] # vector of variable names of the TTE endpoints
    
    ## design matrix for the weights
    if(is.null(Wscheme)){
      res_init <- initWscheme(D=D,endpoint=endpoint,endpoint.TTE=endpoint.TTE,D.TTE=D.TTE,threshold=threshold,type=type)
      Wscheme <- res_init$Wscheme  
      index_survivalM1 <- res_init$index_survivalM1
      threshold_TTEM1 <- res_init$threshold_TTEM1       
    }
    
    ## Survival estimate using Kaplan Meier    
    res_init <- initSurvival(M.Treatment=M.Treatment,M.Control=M.Control,M.delta_Treatment=M.delta_Treatment,M.delta_Control=M.delta_Control,
                             endpoint=endpoint,D.TTE=D.TTE,type=type,threshold=threshold,
                             index.strataT=index.strataT,index.strataC=index.strataC,n.strata=n.strata,   
                             method=method)
    
    list_survivalT <- res_init$list_survivalT
    list_survivalC <- res_init$list_survivalC
    
  }else{
    Wscheme <- matrix()  # factice design matrix for the weights. Will be sent to the C++ arguments to fill the argument but not used by the function.
    list_survivalT <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
    list_survivalC <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
    index_survivalM1 <- numeric(0) # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    threshold_TTEM1 <- numeric(0)  # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }
  
  #### export ####
  res <- list()
  res$M.Treatment <-  M.Treatment
  res$M.Control <-  M.Control
  res$M.delta_Treatment <-  M.delta_Treatment
  res$M.delta_Control <-  M.delta_Control
  res$Wscheme <- Wscheme
  res$threshold_TTEM1 <- threshold_TTEM1
  res$index_survivalM1 <- index_survivalM1
  res$list_survivalT <- list_survivalT
  res$list_survivalC <- list_survivalC
  
  return(res)
}

#' @rdname internal-intilisation
initN <- function(n){
  
  #### tests ####
  
  if(length(n)!=3){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must have length 3 \n",
         "length(n) : ",length(n),"\n")
  }
  
  if(any(n %% 1>0) || any(n<=0)){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must contains strictly positive integers \n",
         "proposed n : ",paste(n,collpase=" "),"\n")
  }
  
  if(is.null(names(n))){names(n) <- c("from","to","by")}
  
  if(any(names(n) %in% c("from","to","by","length.out")==FALSE)){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "elements of \'n\' must be named with \"from\", \"to\", \"by\" or \"length.out\" \n",
         "invalid names : ",paste(names(n)[names(n) %in% c("from","to","by","length.out")==FALSE] ,collpase=" "),"\n")    
  }
  
  if(length(names(n)) != length(unique(names(n)))){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "elements in \'n\' have the same names \n")    
  }
  
  if("from" %in% names(n) ==FALSE){
    stop("BuyseTest : wrong specification of \'n\' \n",
         "\'n\' must contain an element named \"from\" \n",
         "names(n) : ",paste(names(n) ,collpase=" "),"\n")    
  }
  
  #### computation ####
  
  test_to <- "to" %in% names(n)
  test_by <- "by" %in% names(n)
  test_length.out <- "length.out" %in% names(n)
  if( test_to && test_by){
    n <- seq(from=n["from"],to=n["to"],by=n["by"])
  }
  if( test_to && test_length.out){
    n <- seq(from=n["from"],to=n["to"],length.out=n["length.out"])
  }
  if( test_by && test_length.out){
    n <- seq(from=n["from"],by=n["by"],length.out=n["length.out"])
  }
  
  #### export ####
  return(n)
}

#' @rdname internal-intilisation
initHypothesis <- function(hypothesis,type,D){
  
  if(!is.list(hypothesis) && D==1){
    hypothesis <- list(hypothesis)
  }
  
  if(!is.list(hypothesis) || length(hypothesis)!=D){
    stop("BuyseTest : proposed \'hypothesis\' does not match \'type\' \n",
         "\'hypothesis\' must be a list with ",D," elements \n",
         "is(hypothesis) : ",paste(is(hypothesis),collapse=" "),"\n",
         "length(hypothesis) : ",length(hypothesis),"\n"
    )    
  }
  
  valid_parameters <- list(binary=c("proba_t","proba_c"),
                           continuous=c("mu_t","mu_c","sigma"),
                           TTE=c("lambda_t","lambda_c","T_inclusion","T_followUp")
  )
  
  valid_range <- list(proba_t=c(0,1),
                      proba_c=c(0,1),
                      sigma=c(0,Inf),
                      lambda_t=c(0,Inf),
                      lambda_c=c(0,Inf),
                      T_inclusion=c(0,Inf),
                      T_followUp=c(0,Inf))
  
  for(iter_outcome in 1:D){
    
    parameters_tempo <- valid_parameters[[type[iter_outcome]]]
    
    if(length(hypothesis[[iter_outcome]]) != length(parameters_tempo)){
      stop("BuyseTest : proposed \'hypothesis\' does not match \'type\' \n",
           "type of outcome ",iter_outcome," : ",names(valid_parameters)[type[iter_outcome]]," \n",
           "number of parameters requested for this type of outcome : ",length(parameters_tempo)," \n",
           "length(hypothesis[[",iter_outcome,"]]) : ",length(hypothesis[[iter_outcome]])," \n")   
    }
    
    if(is.null(names(hypothesis[[iter_outcome]]))){
      names(hypothesis[[iter_outcome]]) <- parameters_tempo
    }
    
    if(any(parameters_tempo %in% names(hypothesis[[iter_outcome]]) == FALSE)){
      stop("BuyseTest : wrong specification of \'hypothesis\' \n",
           "type of outcome ",iter_outcome," : ",names(valid_parameters)[type[iter_outcome]]," \n",
           "parameters requested for this type of outcome : \"",paste(parameters_tempo,collapse="\" \""),"\" \n",
           "names(hypothesis[[",iter_outcome,"]]) : \"",paste(names(hypothesis[[iter_outcome]])[names(hypothesis[[iter_outcome]]) %in% parameters_tempo == FALSE],collapse="\" \""),"\" \n")   
    }
    
    for(iter_param in 1:length(parameters_tempo)){
      if(parameters_tempo[iter_param] %in% names(valid_range)){
        test.inf <- valid_range[[parameters_tempo[iter_param]]][1]>=hypothesis[[iter_outcome]][parameters_tempo[iter_param]]
        test.sup <- valid_range[[parameters_tempo[iter_param]]][2]<=hypothesis[[iter_outcome]][parameters_tempo[iter_param]]
        if(test.inf || test.sup){
          stop("BuyseTest : unvalid parameter value in \'hypothesis\' \n",
               "possible range of values : ",paste(valid_range[[parameters_tempo[iter_param]]],collapse=" ")," \n",
               "proposed value of hypothesis[[",iter_outcome,"]][",parameters_tempo[iter_param],"] :  ",hypothesis[[iter_outcome]][parameters_tempo[iter_param]]," \n"
          )  
        }
      }
    }
    
    
  }
  
  #### export ####
  return(hypothesis)
  
}

#' @rdname internal-intilisation
initSpace  <- function(nchar){  
  return(sapply(nchar,function(x){if(x>0){do.call(paste,c(as.list(rep(" ",x)),sep=""))}else{""}}))  
}

#' @rdname internal-intilisation
initStrata <- function(strata,
                       dataT,dataC,n.Treatment,n.Control,
                       endpoint,censoring){
  
  if(!is.null(strata)){
    strataT <- interaction(dataT[,strata,with=FALSE],drop = TRUE,lex.order=FALSE,sep=".") # strata variable for the treatment arm transformed into factor
    strataC <- interaction(dataC[,strata,with=FALSE],drop = TRUE,lex.order=FALSE,sep=".") # strata variable for the control arm transformed into factor
    levels.strata <- levels(strataT) # extraction of the levels of the strata variables
    
    if(length(levels.strata) != length(levels(strataC)) || any(levels.strata != levels(strataC))){
      stop("BuyseTest : wrong specification of \'strata\' \n",
           "different levels between Control and Treatment \n",
           "levels(strataT) : ",paste(levels(strataT),collapse=" "),"\n",
           "levels(strataC) : ",paste(levels(strataC),collapse=" "),"\n")
    }
    
    strataT <- as.numeric(strataT) # strata variable for the treatment arm converted into numeric 
    strataC <- as.numeric(strataC) # strata variable for the control arm converted into numeric
    
  }else{ # if there is no strata variable the same strata is used for all patient
    strataT <- rep(1,n.Treatment) # all patient of the treatment arm are in the same strata : strata 1
    strataC <- rep(1,n.Control)  # all patient of the control arm are in the same strata : strata 1
    levels.strata <- 1 # the only strata is strata 1
  }
  
  n.strata <- length(levels.strata) # number of strata
  index.strataT <- lapply(1:n.strata,function(x){which(x==strataT)-1}) # for each strata, the index of the patients belonging to this strata in the treatment arm is stored in an element of the list. 
  index.strataC <- lapply(1:n.strata,function(x){which(x==strataC)-1}) # for each strata, the index of the patients belonging to this strata in the control arm is stored in an element of the list. 
  # Because of the minus one after the which, both indexes begin at 0. This is compulsory for C++ compatibility
  
  # check there is at least one element per arm in each strata
  if(any(unlist(lapply(index.strataT,length))==0)){
    stop("BuyseTest : some strata contain no patient in the treatment arm \n",
         "strata : ",paste(levels.strata[which(unlist(lapply(index.strataT,length))==0)],collapse=" ")," \n")
  }

  if(any(unlist(lapply(index.strataC,length))==0)){
    stop("BuyseTest : some strata contain no patient in the control arm \n",
         "strata : ",paste(levels.strata[which(unlist(lapply(index.strataC,length))==0)],collapse=" ")," \n")
  }
  
  #### export ####
  res <- list()
  res$index.strataT <- index.strataT
  res$index.strataC <- index.strataC
  res$n.strata <- n.strata
  res$levels.strata <- levels.strata
  
  return(res)
}

#' @rdname internal-intilisation
#' @export
initSurvival <- function(M.Treatment,M.Control,M.delta_Treatment,M.delta_Control,
                         endpoint,D.TTE,type,threshold,
                         index.strataT,index.strataC,n.strata,
                         method){
  
  #### conversion to R index
  index.strataT <- lapply(index.strataT,function(x){x+1})
  index.strataC <- lapply(index.strataC,function(x){x+1})
  
  #### preparing 
  if(method %in% c("Efron","Peron")){
    colnamesT <- c(paste("SurvivalT_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalC_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalT_TimeT_0(ordered)","SurvivalC_TimeT-threshold(ordered)","SurvivalC_TimeT+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
    
    colnamesC <- c(paste("SurvivalC_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalT_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalC_TimeC_0(ordered)","SurvivalT_TimeC-threshold(ordered)","SurvivalT_TimeC+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
    
    list_survivalT <- rep(list(matrix(NA,nrow=nrow(M.Treatment),ncol=11,dimnames=list(1:nrow(M.Treatment),colnamesT))),times=D.TTE) # list of matrix of survival data for each endpoint in the treatment arm
    list_survivalC <- rep(list(matrix(NA,nrow=nrow(M.Control),ncol=11,dimnames=list(1:nrow(M.Control),colnamesC))),times=D.TTE) # list of matrix of survival data for each endpoint in the control arm
  }else{ # method == "Peto"
    colnamesT <- paste("Survival_TimeT",c("-threshold","_0","+threshold"),sep="")
    colnamesC <- paste("Survival_TimeC",c("-threshold","_0","+threshold"),sep="")
    
    list_survivalT <- rep(list(matrix(NA,nrow=nrow(M.Treatment),ncol=3,dimnames=list(1:nrow(M.Treatment),colnamesT))),times=D.TTE) # list of matrix of survival data for each endpoint in the treatment arm
    list_survivalC <- rep(list(matrix(NA,nrow=nrow(M.Control),ncol=3,dimnames=list(1:nrow(M.Control),colnamesC))),times=D.TTE) # list of matrix of survival data for each endpoint in the control arm
  }
  
  #### loop 
  for(iter_strata in 1:n.strata){

    Mstrata.Treatment <- M.Treatment[index.strataT[[iter_strata]],,drop=FALSE]
    Mstrata.Control <- M.Control[index.strataC[[iter_strata]],,drop=FALSE]
    Mstrata.delta_Treatment <- M.delta_Treatment[index.strataT[[iter_strata]],,drop=FALSE]
    Mstrata.delta_Control <- M.delta_Control[index.strataC[[iter_strata]],,drop=FALSE]
    
    for(iter_endpointTTE in 1:D.TTE){
     
      if(method %in% c("Efron","Peron")){
   
        ## treatment
        time_treatment <- Mstrata.Treatment[,endpoint[type==3][iter_endpointTTE]] # store the event times of the treatment arm for a given TTE endpoint
        resKMT_tempo <- survival::survfit(survival::Surv(time_treatment,Mstrata.delta_Treatment[,iter_endpointTTE])~1) # compute the survival over the controls with common Kaplan Meier estimator.  
        # prefer sort(time_treatment) to resKMT_tempo$time to avoid rounding error
        
        if(all(Mstrata.delta_Treatment[which(time_treatment==max(time_treatment)),iter_endpointTTE]==1)){ # step interpolation of the survival function (last = event)
          
          survestKMT_tempo <- stats::approxfun(x=unique(sort(time_treatment)),y=resKMT_tempo$surv,
                                               yleft=1,yright=0,f=0, method = "constant") # 0 at infinity if last event is a death      
          
        }else{ # step interpolation of the survival function (last = censoring)
          
          survestKMT_tempo <- stats::approxfun(x=unique(sort(time_treatment)),y=resKMT_tempo$surv,
                                               yleft=1,yright=NA,f=0, method = "constant") # NA at infinity if last event is censored
          
        }
        
        ## control
        time_control <- Mstrata.Control[,endpoint[type==3][iter_endpointTTE]] # store the event times of the control arm for a given TTE endpoint
        resKMC_tempo <- survival::survfit(survival::Surv(time_control,Mstrata.delta_Control[,iter_endpointTTE])~1) # compute the survival over the treatments with common Kaplan Meier estimator.  
        
        if(all(Mstrata.delta_Control[which(time_control==max(time_control)),iter_endpointTTE]==1)){  # step interpolation of the survival function (last = event)
          
          survestKMC_tempo <- stats::approxfun(x=unique(sort(time_control)), y=resKMC_tempo$surv,
                                               yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
          
        }else{ # step interpolation of the survival function (last = censoring)
          
          survestKMC_tempo <- stats::approxfun(x=unique(sort(time_control)),y=resKMC_tempo$surv,
                                               yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is censored      
          
        }
        
        #### survival
        ## treatment
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],1] <- survestKMT_tempo(time_treatment-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],2] <- survestKMT_tempo(time_treatment) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],3] <- survestKMT_tempo(time_treatment+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],4] <- survestKMC_tempo(time_treatment-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],5] <- survestKMC_tempo(time_treatment) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],6] <- survestKMC_tempo(time_treatment+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        
        order_time_treatment <- order(time_treatment)
        
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],7:9] <- list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],c(2,4,6),drop=FALSE][order_time_treatment,,drop=FALSE] # St[t_treatment], Sc[t_treatment-threshold,t_treatment+threshold]
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],10] <-  time_treatment[order_time_treatment]
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],11] <-   Mstrata.delta_Treatment[order_time_treatment,iter_endpointTTE] # to compute the integral (Efron)
        
        rownames(list_survivalT[[iter_endpointTTE]])[index.strataT[[iter_strata]]] <- time_treatment
        
        ## control
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],1] <- survestKMC_tempo(time_control-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],2] <- survestKMC_tempo(time_control) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],3] <- survestKMC_tempo(time_control+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],4] <- survestKMT_tempo(time_control-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],5] <- survestKMT_tempo(time_control) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],6] <- survestKMT_tempo(time_control+threshold[type==3][iter_endpointTTE]) # survival at t + tau                                              
        
        order_time_control <- order(time_control)
        
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],7:9] <- list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],c(2,4,6),drop=FALSE][order_time_control,,drop=FALSE] # Sc[t_control], St[t_control-threshold,t_control+threshold]
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],10] <- time_control[order_time_control]
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],11] <- Mstrata.delta_Control[order_time_control,iter_endpointTTE] # to compute the integral (Efron)
        
        rownames(list_survivalC[[iter_endpointTTE]])[index.strataC[[iter_strata]]] <- time_control    
    
      }else if(method=="Peto"){ # Peto
        time_treatment <- Mstrata.Treatment[,endpoint[type==3][iter_endpointTTE]] # store the event times of the treatment arm for a given TTE endpoint
        time_control <- Mstrata.Control[,endpoint[type==3][iter_endpointTTE]] # store the event times of the control arm for a given TTE endpoint
        
        # prefer c(time_treatment,time_control) to time_all to avoid rounding error
        censoring_all <- rbind(Mstrata.delta_Treatment,Mstrata.delta_Control)[,iter_endpointTTE]
        resKM_tempo <- survival::survfit(survival::Surv(c(time_treatment,time_control),censoring_all)~1) # compute the survival over the entire cohort (Control and Treatment) with common Kaplan Meier estimator.  
        max_time <- max(resKM_tempo$time)+10^{-12}
        
        if(all(censoring_all[which(c(time_treatment,time_control)==max(c(time_treatment,time_control)))]==1)){ # step interpolation of the survival function
          survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                              yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
          
        }else{
          survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                              yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is a death
        }
        
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],1] <- survestKM_tempo(time_treatment-threshold[type==3][iter_endpointTTE]) # survival at t - tau, tau being the threshold
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],2] <- survestKM_tempo(time_treatment) # survival at t
        list_survivalT[[iter_endpointTTE]][index.strataT[[iter_strata]],3] <- survestKM_tempo(time_treatment+threshold[type==3][iter_endpointTTE]) # survival at t + tau
       
        rownames(list_survivalT[[iter_endpointTTE]])[index.strataT[[iter_strata]]] <- time_treatment
        
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],1] <- survestKM_tempo(time_control-threshold[type==3][iter_endpointTTE]) # survival at t - tau
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],2] <- survestKM_tempo(time_control) # survival at t
        list_survivalC[[iter_endpointTTE]][index.strataC[[iter_strata]],3] <- survestKM_tempo(time_control+threshold[type==3][iter_endpointTTE]) # survival at t + tau
        
        rownames(list_survivalC[[iter_endpointTTE]])[index.strataC[[iter_strata]]] <- time_control    
        
      }
      
    }
  }
  
  return(list(list_survivalT=list_survivalT,
              list_survivalC=list_survivalC))
  
}

#' @rdname internal-intilisation
initThreshold <- function(threshold,type,D,method,endpoint){
  
  ## threshold
  if(is.null(threshold)){
    threshold <- rep(10^{-12},D)  # if no treshold is proposed all threshold are by default set to 10^{-12}
    if(any(type==1)){threshold[type==1] <- NA} # except for threshold corresponding to binary endpoints that are set to NA.
  }else  if(any(stats::na.omit(threshold)==0)){
    threshold[stats::na.omit(which(threshold==0))] <- 10^{-12}
  }
 
  validNumeric(threshold, validLength = D, refuse.NA = FALSE , method = "BuyseTest")
  
  if(any(!is.na(threshold[type==1])) ){
    stop("BuyseTest : wrong specification of \'threshold\' \n",
         "\'threshold\' must be NA for binary endpoints \n",
         "proposed \'threshold\' : ",paste(threshold[type==1],collapse=" "),"\n",
         "binary endpoints : ",paste(endpoint[type==1],collapse=" "),"\n")
  }
  
  ## duplicates 
  test.duplicated <- duplicated(paste(endpoint,threshold,sep=""))
  if(any(test.duplicated)){
    display_duplicated <- sapply(unique(endpoint[test.duplicated]),function(x){paste("\"",x,"\" (thresholds = ",paste(threshold[x==endpoint],collapse=" "),")",sep="")})
    stop("BuyseTest : wrong specification of \'endpoint\' or \'threshold\' \n",
         "there are duplicated endpoints : ",paste(display_duplicated,
                                                   collapse="\n                                 "),"\n")    
  }
  
  ## test that the thresholds corresponding to TTE endpoints are decreasing 
  if(method %in% c("Peto","Efron","Peron")){
    endpoint.TTE <- endpoint[type==3]
    
    test.decreasing <- sapply(unique(endpoint.TTE),function(x){
      thresholdTTE <- threshold[endpoint == x]
      sum(abs(thresholdTTE-sort(thresholdTTE,decreasing=TRUE)<0))+sum(diff(thresholdTTE)==0) # compare the thresolds to the thresholds in decreasing order
    })
    
    if(any(test.decreasing>0)){
      stop("BuyseTest : wrong specification of \'endpoint\' or \'threshold\' \n",
           "survival endpoints must be used with decreasing threshold when re-used with lower priority \n",
           "survival endpoints             : ",paste(endpoint[type==3][test.decreasing>0],
                                                     initSpace(nchar(threshold[type==3][test.decreasing>0])-nchar(endpoint[type==3][test.decreasing>0])),
                                                     sep="",collapse=" "),"\n",
           "have non-decreasing thresholds : ",paste(threshold[type==3][test.decreasing>0],
                                                     initSpace(nchar(endpoint[type==3][test.decreasing>0])-nchar(threshold[type==3][test.decreasing>0])),
                                                     sep="",collapse=" "),"\n")      
    }
  }
  
  #### export ####
  return(threshold)
}


#' @rdname internal-intilisation
initWscheme <- function(endpoint,D,endpoint.TTE,D.TTE,threshold,type){

  if(D>1){
        
    Wscheme <- matrix(NA,nrow=D.TTE,ncol=D-1) # design matrix indicating to combine the weights obtained at differents TTE endpoints
    rownames(Wscheme) <- paste("weigth of ",endpoint.TTE,"(",threshold[type==3],")",sep="")
    colnames(Wscheme) <- paste("for ",endpoint[-1],"(",threshold[-1],")",sep="")
    
    index_survivalM1 <- rep(-1,D.TTE) # index of previous TTE endpoint (-1 if no previous TTE endpoint i.e. endpoint has not already been used)
    threshold_TTEM1 <- rep(-1,D.TTE) # previous threshold (-1 if no previous threshold i.e. endpoint has not already been used)
    
    iter_endpoint.TTE <- if(type[1]==3){1}else{0} #  index_survivalM1 and  threshold_TTEM1 are -1 even if the first endoint is a survival endpoint
    
    for(iter_endpoint in 2:D){   
      
      if(type[iter_endpoint]==3){ 
        iter_endpoint.TTE <- iter_endpoint.TTE + 1
      }
      
      # select valid rows
      index_rowTTE <- which(paste(endpoint.TTE,threshold[type==3],sep="_") %in% paste(endpoint,threshold,sep="_")[1:(iter_endpoint-1)])
      if(length(index_rowTTE)==0){next} # not yet TTE endpoints (no valid rows)
      Wscheme[index_rowTTE,iter_endpoint-1] <- 0 # potential weights
      
      # keep only the last repeated endpoints
      index_rowTTE <- sapply(unique(endpoint.TTE[index_rowTTE]),
                                   function(x){utils::tail(index_rowTTE[endpoint.TTE[index_rowTTE]==x],1)})
      
      # if survival endpoint remove similar endpoints
      if(endpoint[iter_endpoint] %in% endpoint.TTE[index_rowTTE]){
        index_survivalM1[iter_endpoint.TTE] <- index_rowTTE[endpoint.TTE[index_rowTTE]==endpoint[iter_endpoint]]-1
        threshold_TTEM1[iter_endpoint.TTE] <- threshold[type==3][index_survivalM1[iter_endpoint.TTE]+1]
        index_rowTTE <- setdiff(index_rowTTE,index_survivalM1[iter_endpoint.TTE]+1)
      }
      
      # update Wscheme
      if(length(index_rowTTE)>0){
        Wscheme[index_rowTTE,iter_endpoint-1] <- 1
      }
      
    }

  }else{
    Wscheme <- matrix(nrow=0,ncol=0)
    index_survivalM1 <- numeric(0)
    threshold_TTEM1 <- numeric(0)
  }
  
  #### export ####
  res <- list()
  res$Wscheme <- Wscheme
  res$index_survivalM1 <- index_survivalM1
  res$threshold_TTEM1 <- threshold_TTEM1
  
  return(res)
  
}