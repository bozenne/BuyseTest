#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test the computation of the survival
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# See initSurvival function file BuyseInit.R (approx. line 433)
#

## no teste : saut au moment du dernier evenement.

#### local 
# path <- "E:/Creation_package/Package_BuyseTest/BuyseTest" # path to the uncompressed tar.gz file
# Rcpp::sourceCpp(file.path(path,"src/FCT_BuyseTest.cpp"),rebuild=TRUE)
# source(file.path(path,"R/FCT_buyseTest.R"))
# source(file.path(path,"R/OBJET_buyseTest.R"))
# source(file.path(path,"R/FCT_buyseInit.R"))

options(error=function() traceback(2)) 
options(max.print=10000)

#### package
require(BuyseTest)

#### save
save <- FALSE #   TRUE #  
precision <- 10^{-7}

#### 1- No strata ####
M.Treatment <- cbind(time=1:5)
M.Control <-  cbind(time=c(1:5-0.1,5,5))
threshold <- 0.001 

if(save){
  results_initSurvival <- list()
  results_initSurvival$NoStrata <- list()
  results_initSurvival$Strata <- list()
}else{
  load("BuyseTest/tests/test_initSurvival-results_initSurvival.RData")
}

for(iter_dataset in 1:3){
  
  if(save){
    results_initSurvival$NoStrata[[iter_dataset]] <- list()
  }
  
  if(iter_dataset==1){
    M.delta_Treatment <- cbind(time=c(1,0,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,0,1,0,0,0))
  }else if(iter_dataset==2){
    M.delta_Treatment <- cbind(time=c(1,0,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,0,1,0,1,1))
  }else if(iter_dataset==3){
    M.delta_Treatment <- cbind(time=c(1,1,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,1,1,1,1,1))
  }
  
#### Peto ####

SurvivalPeto_noStrata <- initSurvival(M.Treatment=M.Treatment, 
                                      M.Control=M.Control, 
                                      M.delta_Treatment=M.delta_Treatment, 
                                      M.delta_Control=M.delta_Control, 
                                      endpoint="time",
                                      index.strataT=list(seq(0,nrow(M.Treatment)-1)),index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                                      D.TTE=1, type=3, threshold=threshold, method="Peto", 
                                      callMethod = "initSurvival")
if(save){
  results_initSurvival$NoStrata[[iter_dataset]]$Peto <- SurvivalPeto_noStrata
}else{
  
  # mark last events
  last.time <- which.max(c(M.Treatment,M.Control))
  indexC.last <- which(M.Control==last.time)
  indexT.last <- which(M.Treatment==last.time)
  
  #### jump at event (expept last event)
  index.Tevent <- setdiff(which(M.delta_Treatment[,1]==1),indexT.last)
  test.Tevent_before <- SurvivalPeto_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT-threshold"] -  SurvivalPeto_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT_0"] 
  test.Tevent_after <- SurvivalPeto_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT_0"] -  SurvivalPeto_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT+threshold"] 

  if(any(test.Tevent_before<=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Tevent_after>0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  index.Cevent <- setdiff(which(M.delta_Control[,1]==1),indexC.last)
  test.Cevent_before <- SurvivalPeto_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC-threshold"] -  SurvivalPeto_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC_0"] 
  test.Cevent_after <- SurvivalPeto_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC_0"] -  SurvivalPeto_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC+threshold"] 
  
  if(any(test.Cevent_before<=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Cevent_after>0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  #### jump at censoring 
  index.Tcensoring <- setdiff(which(M.delta_Treatment[,1]==0),indexT.last)
  test.Tcensoring_before <- SurvivalPeto_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT-threshold"] -  SurvivalPeto_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT_0"] 
  test.Tcensoring_after <- SurvivalPeto_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT_0"] -  SurvivalPeto_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT+threshold"] 
  
  if(any(test.Tcensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Tcensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  index.Ccensoring <- setdiff(which(M.delta_Control[,1]==0),indexC.last)
  test.Ccensoring_before <- SurvivalPeto_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC-threshold"] -  SurvivalPeto_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC_0"] 
  test.Ccensoring_after <- SurvivalPeto_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC_0"] -  SurvivalPeto_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC+threshold"] 
  
  if(any(test.Ccensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Ccensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  #### last events
  test.censoringLast <- any(c(M.delta_Treatment[indexT.last],M.delta_Control[indexC.last])==0)
  
  test.Tafter <-  SurvivalPeto_noStrata$list_survivalT[[1]][indexT.last,"Survival_TimeT+threshold"]
  test.Cafter <-  SurvivalPeto_noStrata$list_survivalC[[1]][indexC.last,"Survival_TimeC+threshold"]
  
  if(test.censoringLast && (!is.na(test.Tafter) || !is.na(test.Cafter)) ){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect \n",
         "may be due to bad management of survival after the last event when the last observation is censored \n",)
  }
  
  if(test.censoringLast==FALSE && (any(test.Tafter!=0) || any(test.Cafter!=0)) ){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect \n",
         "may be due to bad management of survival after the last event when the last observation is censored \n",)
  }
  
  #### overall 
  test.T <- results_initSurvival$NoStrata[[iter_dataset]]$Peto$list_survivalT[[1]] - SurvivalPeto_noStrata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Peto$list_survivalT[[1]]) - is.na(SurvivalPeto_noStrata$list_survivalT[[1]])
  test.C <- results_initSurvival$NoStrata[[iter_dataset]]$Peto$list_survivalT[[1]] - SurvivalPeto_noStrata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Peto$list_survivalC[[1]]) - is.na(SurvivalPeto_noStrata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(last.time,indexC.last,indexT.last,
     index.Tevent,test.Tevent_before,test.Tevent_after,index.Cevent,test.Cevent_before,test.Cevent_after,
     index.Tcensoring,test.Tcensoring_before,test.Tcensoring_after,index.Ccensoring,test.Ccensoring_before,test.Ccensoring_after,
     test.censoringLast,test.Tafter,test.Cafter,
     test.T,test.TNA,test.C,test.CNA)
}

#### display

timesT_th <- c(M.Treatment[,1]-threshold,M.Treatment[,1],M.Treatment[,1]+threshold)
plot(timesT_th[order(timesT_th)],as.vector(SurvivalPeto_noStrata$list_survivalT[[1]])[order(timesT_th)],
     ylab="Survival Peto (Treatement)",xlab="time",ylim=c(0,1),
     type="o")
#     type="p")

SurvivalPeto_noStrata$list_survivalT[[1]][,"Survival_TimeT-threshold"] # survival before jump points
SurvivalPeto_noStrata$list_survivalT[[1]][,"Survival_TimeT_0"] # survival at jump points
SurvivalPeto_noStrata$list_survivalT[[1]][,"Survival_TimeT+threshold"] # survival after jump points

timesC_th <- c(M.Control[,1]-threshold,M.Control[,1],M.Control[,1]+threshold)
plot(timesC_th[order(timesC_th)],as.vector(SurvivalPeto_noStrata$list_survivalC[[1]])[order(timesC_th)],
     ylab="Survival Peto (Control)",xlab="time",ylim=c(0,1),
     type="o")
# type="p")

SurvivalPeto_noStrata$list_survivalC[[1]][,"Survival_TimeC-threshold"] # survival before jump points
SurvivalPeto_noStrata$list_survivalC[[1]][,"Survival_TimeC_0"] # survival at jump points
SurvivalPeto_noStrata$list_survivalC[[1]][,"Survival_TimeC+threshold"] # survival after jump points

#### Efron ####	
SurvivalEfron_noStrata <- BuyseTest:::initData(data=data.frame(time=c(M.Treatment,M.Control),censoring=c(M.delta_Treatment,M.delta_Control)), D=1,
                     Ind.Treatment=c(rep(1,nrow(M.Treatment)),rep(0,nrow(M.Control))),
                     type=3, endpoint="time", censoring="censoring", method="Efron", 
                     index.strataT=list(seq(0,nrow(M.Treatment)-1)),index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                     D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE, callMethod = "initData") 

if(save){
  results_initSurvival$NoStrata[[iter_dataset]]$Efron <- SurvivalEfron_noStrata
}else{
  # mark last events
  last.time <- which.max(c(M.Treatment,M.Control))
  indexC.last <- which(M.Control==last.time)
  indexT.last <- which(M.Treatment==last.time)
  
  #### jump at event (expept last event)
  index.Tevent <- setdiff(which(M.delta_Treatment[,1]==1),indexT.last)
  test.Tevent_before <- SurvivalEfron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT-threshold"] -  SurvivalEfron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] 
  test.Tevent_after <- SurvivalEfron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] -  SurvivalEfron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT+threshold"] 
  
  if(any(test.Tevent_before<=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Tevent_after>0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  index.Cevent <- setdiff(which(M.delta_Control[,1]==1),indexC.last)
  test.Cevent_before <- SurvivalEfron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC-threshold"] -  SurvivalEfron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] 
  test.Cevent_after <- SurvivalEfron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] -  SurvivalEfron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC+threshold"] 
  
  if(any(test.Cevent_before<=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Cevent_after>0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  #### jump at censoring 
  index.Tcensoring <- setdiff(which(M.delta_Treatment[,1]==0),indexT.last)
  test.Tcensoring_before <- SurvivalEfron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT-threshold"] -  SurvivalEfron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] 
  test.Tcensoring_after <- SurvivalEfron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] -  SurvivalEfron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT+threshold"] 
  
  if(any(test.Tcensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Tcensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  index.Ccensoring <- setdiff(which(M.delta_Control[,1]==0),indexC.last)
  test.Ccensoring_before <- SurvivalEfron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC-threshold"] -  SurvivalEfron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] 
  test.Ccensoring_after <- SurvivalEfron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] -  SurvivalEfron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC+threshold"] 
  
  if(any(test.Ccensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Ccensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  #### last events   
  test.Tafter <-  SurvivalEfron_noStrata$list_survivalT[[1]][indexT.last,"SurvivalT_TimeT+threshold"]
  test.Cafter <-  SurvivalEfron_noStrata$list_survivalC[[1]][indexC.last,"SurvivalC_TimeC+threshold"]
  
  if(any(test.Tafter!=0) || any(test.Cafter!=0) ){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect \n",
         "may be due to bad management of survival after the last event\n",)
  }
  
  #### overall 
  test.T <- results_initSurvival$NoStrata[[iter_dataset]]$Efron$list_survivalT[[1]] - SurvivalEfron_noStrata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Efron$list_survivalT[[1]]) - is.na(SurvivalEfron_noStrata$list_survivalT[[1]])
  test.C <- results_initSurvival$NoStrata[[iter_dataset]]$Efron$list_survivalT[[1]] - SurvivalEfron_noStrata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Efron$list_survivalC[[1]]) - is.na(SurvivalEfron_noStrata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>0 || sum(abs(test.TNA))>0 || sum(abs(test.C),na.rm=T)>0 || sum(abs(test.CNA))>0){
    stop("test_initSurvival[NoStrata - Efron] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(last.time,indexC.last,indexT.last,
     index.Tevent,test.Tevent_before,test.Tevent_after,index.Cevent,test.Cevent_before,test.Cevent_after,
     index.Tcensoring,test.Tcensoring_before,test.Tcensoring_after,index.Ccensoring,test.Ccensoring_before,test.Ccensoring_after,
     test.Tafter,test.Cafter,
     test.T,test.TNA,test.C,test.CNA)
}

## do not work because initData appply a correction if last observation is censored (it transform it into an event)
# SurvivalEfron_noStrata <- BuyseTest:::initSurvival(M.Treatment=M.Treatment, 
#                                               M.Control=M.Control, 
#                                               M.delta_Treatment=M.delta_Treatment, 
#                                               M.delta_Control=M.delta_Control, 
#                                               endpoint="time",
#                                               D.TTE=1, type=3, threshold=threshold, method="Efron", 
#                                               callMethod = "initSurvival")

#### display
SurvivalEfron_noStrata$M.delta_Treatment
SurvivalEfron_noStrata$M.delta_Control

timesT_th <- c(M.Treatment[,1]-threshold,M.Treatment[,1],M.Treatment[,1]+threshold)
plot(timesT_th[order(timesT_th)],as.vector(SurvivalEfron_noStrata$list_survivalT[[1]])[order(timesT_th)],
     ylab="Survival Efron (Treatement)",xlab="time",ylim=c(0,1),
     type="o")
#     type="p")

SurvivalEfron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT-threshold"] # survival before jump points
SurvivalEfron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT_0"] # survival at jump points
SurvivalEfron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT+threshold"] # survival after jump points

timesC_th <- c(M.Control[,1]-threshold,M.Control[,1],M.Control[,1]+threshold)
plot(timesC_th[order(timesC_th)],as.vector(SurvivalEfron_noStrata$list_survivalC[[1]])[order(timesC_th)],
     ylab="Survival Efron (Control)",xlab="time",ylim=c(0,1),
     type="o")
# type="p")

SurvivalEfron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC-threshold"] # survival before jump points
SurvivalEfron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC_0"] # survival at jump points
SurvivalEfron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC+threshold"] # survival after jump points

#### Peron ####
					
SurvivalPeron_noStrata <- BuyseTest:::initSurvival(M.Treatment=M.Treatment, 
                                              M.Control=M.Control, 
                                              M.delta_Treatment=M.delta_Treatment, 
                                              M.delta_Control=M.delta_Control, 
                                              endpoint="time",
                                              index.strataT=list(seq(0,nrow(M.Treatment)-1)),index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                                              D.TTE=1, type=3, threshold=threshold, method="Peron", 
											                        callMethod = "initSurvival")

if(save){
  results_initSurvival$NoStrata[[iter_dataset]]$Peron <- SurvivalPeron_noStrata
}else{
  # mark last events
  last.time <- which.max(c(M.Treatment,M.Control))
  indexC.last <- which(M.Control==last.time)
  indexT.last <- which(M.Treatment==last.time)
  
  #### jump at event (expept last event)
  index.Tevent <- setdiff(which(M.delta_Treatment[,1]==1),indexT.last)
  test.Tevent_before <- SurvivalPeron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT-threshold"] -  SurvivalPeron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] 
  test.Tevent_after <- SurvivalPeron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] -  SurvivalPeron_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT+threshold"] 
  
  if(any(test.Tevent_before<=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Tevent_after>0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  index.Cevent <- setdiff(which(M.delta_Control[,1]==1),indexC.last)
  test.Cevent_before <- SurvivalPeron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC-threshold"] -  SurvivalPeron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] 
  test.Cevent_after <- SurvivalPeron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] -  SurvivalPeron_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC+threshold"] 
  
  if(any(test.Cevent_before<=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an event) \n")
  }
  
  if(any(test.Cevent_after>0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an event) \n")
  }
  
  #### jump at censoring 
  index.Tcensoring <- setdiff(which(M.delta_Treatment[,1]==0),indexT.last)
  test.Tcensoring_before <- SurvivalPeron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT-threshold"] -  SurvivalPeron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] 
  test.Tcensoring_after <- SurvivalPeron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] -  SurvivalPeron_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT+threshold"] 
  
  if(any(test.Tcensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Tcensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Treated) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  
  index.Ccensoring <- setdiff(which(M.delta_Control[,1]==0),indexC.last)
  test.Ccensoring_before <- SurvivalPeron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC-threshold"] -  SurvivalPeron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] 
  test.Ccensoring_after <- SurvivalPeron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] -  SurvivalPeron_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC+threshold"] 
  
  if(any(test.Ccensoring_before!=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t- (t is an censoring) \n")
  }
  
  if(any(test.Ccensoring_after!=0)){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect (Control) \n",
         "may be due to bad management of survival at t+ (t is an censoring) \n")
  }
  
  #### last events
  test.Tafter <-  SurvivalPeron_noStrata$list_survivalT[[1]][indexT.last,"SurvivalT_TimeT+threshold"]
  test.Cafter <-  SurvivalPeron_noStrata$list_survivalC[[1]][indexC.last,"SurvivalC_TimeC+threshold"]
  
  test.TcensoringLast <- any(M.delta_Treatment[indexT.last]==0)
  test.CcensoringLast <- any(M.delta_Control[indexC.last]==0)
 
  if( (test.TcensoringLast && !is.na(test.Tafter)) || (test.CcensoringLast && !is.na(test.Cafter)) ){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect \n",
         "may be due to bad management of survival after the last event when the last observation is censored \n",)
  }
  
  if( (test.TcensoringLast==FALSE && any(test.Tafter!=0)) || (test.CcensoringLast==FALSE && any(test.Cafter!=0)) ){
    stop("test_initSurvival[NoStrata - Peron] Survival computation is incorrect \n",
         "may be due to bad management of survival after the last event when the last observation is censored \n",)
  }
  
  #### overall
  test.T <- results_initSurvival$NoStrata[[iter_dataset]]$Peron$list_survivalT[[1]] - SurvivalPeron_noStrata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Peron$list_survivalT[[1]]) - is.na(SurvivalPeron_noStrata$list_survivalT[[1]])
  test.C <- results_initSurvival$NoStrata[[iter_dataset]]$Peron$list_survivalT[[1]] - SurvivalPeron_noStrata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$NoStrata[[iter_dataset]]$Peron$list_survivalC[[1]]) - is.na(SurvivalPeron_noStrata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[NoStrata - Peto] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(last.time,indexC.last,indexT.last,
     index.Tevent,test.Tevent_before,test.Tevent_after,index.Cevent,test.Cevent_before,test.Cevent_after,
     index.Tcensoring,test.Tcensoring_before,test.Tcensoring_after,index.Ccensoring,test.Ccensoring_before,test.Ccensoring_after,
     test.TcensoringLast,test.CcensoringLast,test.Tafter,test.Cafter,
     test.T,test.TNA,test.C,test.CNA)
  
}

#### Display
timesT_th <- c(M.Treatment[,1]-threshold,M.Treatment[,1],M.Treatment[,1]+threshold)
plot(timesT_th[order(timesT_th)],as.vector(SurvivalPeron_noStrata$list_survivalT[[1]])[order(timesT_th)],
     ylab="Survival Peron (Treatement)",xlab="time",ylim=c(0,1),
     type="o")
#     type="p")

SurvivalPeron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT-threshold"] # survival before jump points
SurvivalPeron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT_0"] # survival at jump points
SurvivalPeron_noStrata$list_survivalT[[1]][,"SurvivalT_TimeT+threshold"] # survival after jump points

timesC_th <- c(M.Control[,1]-threshold,M.Control[,1],M.Control[,1]+threshold)
plot(timesC_th[order(timesC_th)],as.vector(SurvivalPeron_noStrata$list_survivalC[[1]])[order(timesC_th)],
     ylab="Survival Peron (Control)",xlab="time",ylim=c(0,1),
     type="o")
# type="p")

SurvivalPeron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC-threshold"] # survival before jump points
SurvivalPeron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC_0"] # survival at jump points
SurvivalPeron_noStrata$list_survivalC[[1]][,"SurvivalC_TimeC+threshold"] # survival after jump points

} # end iter_dataset

#### 2- Avec strates ####
n.strata <- 2
M.Treatment <- cbind(time=rep(1:5,n.strata))
M.Control <-  cbind(time=rep(c(1:5-0.1,5,5),n.strata))
M.delta_Treatment <- cbind(time=rep(c(1,0,1,1,1),n.strata))
M.delta_Control <- cbind(time=rep(c(1,1,0,1,0,0,0),n.strata))
threshold <- 0.001 
index.StrataT <- list(0:4,5:9)
index.StrataC <- list(0:6,7:13)

#### Peto ####

SurvivalPeto_Strata <- BuyseTest:::initSurvival(M.Treatment=M.Treatment, 
                                             M.Control=M.Control, 
                                             M.delta_Treatment=M.delta_Treatment, 
                                             M.delta_Control=M.delta_Control, 
                                             endpoint="time",
                                             index.strataT=index.StrataT,index.strataC=index.StrataC,n.strata=2,
                                             D.TTE=1, type=3, threshold=threshold, method="Peto", 
                                             callMethod = "initSurvival")

if(save){
  results_initSurvival$Strata$Peto <- SurvivalPeto_Strata
}else{
  
  #### identical between strata
  test.T <- SurvivalPeto_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,] - SurvivalPeto_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,]
  test.TNA <- is.na(SurvivalPeto_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,]) - is.na(SurvivalPeto_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
  test.C <- SurvivalPeto_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,] - SurvivalPeto_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,]
  test.CNA <- is.na(SurvivalPeto_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,]) - is.na(SurvivalPeto_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[Strata - Peto] Survival computation is incorrect \n",
         "incoherency between strata \n")
  }
  
  #### overalll
  test.T <- results_initSurvival$Strata$Peto$list_survivalT[[1]] - SurvivalPeto_Strata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$Strata$Peto$list_survivalT[[1]]) - is.na(SurvivalPeto_Strata$list_survivalT[[1]])
  test.C <- results_initSurvival$Strata$Peto$list_survivalT[[1]] - SurvivalPeto_Strata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$Strata$Peto$list_survivalC[[1]]) - is.na(SurvivalPeto_Strata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[Strata - Peto] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(test.T,test.TNA,test.C,test.CNA)
}

#### Efron #### 

SurvivalEfron_Strata <- BuyseTest:::initData(data=data.frame(time=c(M.Treatment,M.Control),censoring=c(M.delta_Treatment,M.delta_Control)), 
                                          Ind.Treatment=c(rep(1,nrow(M.Treatment)),rep(0,nrow(M.Control))),D=1,
                                          type=3, endpoint="time", censoring="censoring", method="Efron", 
                                          index.strataT=index.StrataT,index.strataC=index.StrataC,n.strata=2,
                                          D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE, callMethod = "initData") 

if(save){
  results_initSurvival$Strata$Efron <- SurvivalEfron_Strata
}else{
  
  #### identical between strata
  test.T <- SurvivalEfron_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,] - SurvivalEfron_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,]
  test.TNA <- is.na(SurvivalEfron_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,]) - is.na(SurvivalEfron_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
  test.C <- SurvivalEfron_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,] - SurvivalEfron_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,]
  test.CNA <- is.na(SurvivalEfron_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,]) - is.na(SurvivalEfron_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[Strata - Efron] Survival computation is incorrect \n",
         "incoherency between strata \n")
  }
  
  #### overalll
  
  test.T <- results_initSurvival$Strata$Efron$list_survivalT[[1]] - SurvivalEfron_Strata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$Strata$Efron$list_survivalT[[1]]) - is.na(SurvivalEfron_Strata$list_survivalT[[1]])
  test.C <- results_initSurvival$Strata$Efron$list_survivalT[[1]] - SurvivalEfron_Strata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$Strata$Efron$list_survivalC[[1]]) - is.na(SurvivalEfron_Strata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>0 || sum(abs(test.TNA))>0 || sum(abs(test.C),na.rm=T)>0 || sum(abs(test.CNA))>0){
    stop("test_initSurvival[Strata - Efron] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(test.T,test.TNA,test.C,test.CNA)
}

#### Peron #### 

SurvivalPeron_Strata <- BuyseTest:::initData(data=data.frame(time=c(M.Treatment,M.Control),censoring=c(M.delta_Treatment,M.delta_Control)), 
                                                Ind.Treatment=c(rep(1,nrow(M.Treatment)),rep(0,nrow(M.Control))),
                                                type=3, endpoint="time", censoring="censoring", method="Peron",D=1, 
                                                index.strataT=index.StrataT,index.strataC=index.StrataC,n.strata=2,
                                                D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE, callMethod = "initData") 
 
if(save){
  results_initSurvival$Strata$Peron <- SurvivalPeron_Strata
}else{
  
  #### identical between strata
  test.T <- SurvivalPeron_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,] - SurvivalPeron_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,]
  test.TNA <- is.na(SurvivalPeron_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,]) - is.na(SurvivalPeron_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
  test.C <- SurvivalPeron_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,] - SurvivalPeron_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,]
  test.CNA <- is.na(SurvivalPeron_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,]) - is.na(SurvivalPeron_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
  
  if(sum(abs(test.T),na.rm=T)>precision || sum(abs(test.TNA))>precision || sum(abs(test.C),na.rm=T)>precision || sum(abs(test.CNA))>precision){
    stop("test_initSurvival[Strata - Peron] Survival computation is incorrect \n",
         "incoherency between strata \n")
  }
  
  #### overalll
  test.T <- results_initSurvival$Strata$Peron$list_survivalT[[1]] - SurvivalPeron_Strata$list_survivalT[[1]]
  test.TNA <- is.na(results_initSurvival$Strata$Peron$list_survivalT[[1]]) - is.na(SurvivalPeron_Strata$list_survivalT[[1]])
  test.C <- results_initSurvival$Strata$Peron$list_survivalT[[1]] - SurvivalPeron_Strata$list_survivalT[[1]]
  test.CNA <- is.na(results_initSurvival$Strata$Peron$list_survivalC[[1]]) - is.na(SurvivalPeron_Strata$list_survivalC[[1]])
  
  if(sum(abs(test.T),na.rm=T)>0 || sum(abs(test.TNA))>0 || sum(abs(test.C),na.rm=T)>0 || sum(abs(test.CNA))>0){
    stop("test_initSurvival[Strata - Peron] Survival computation is incorrect \n")
  }
  
  #### cleaning
  rm(test.T,test.TNA,test.C,test.CNA)
}


#### Z- export ####

if(save){
  save(results_initSurvival,file="test_initSurvival-results_initSurvival.RData")
}
