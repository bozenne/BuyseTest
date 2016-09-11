verboseContext("Check KM computation")

#### additional spec
n.patients <- 200
precision <- 10^{-7}
results_initSurvival <- list()
#### 1- No strata ####
M.Treatment <- cbind(time=1:5)
M.Control <-  cbind(time=c(1:5-0.1,5,5))
threshold <- 0.001 # threshold smaller than the interval between two events
results_initSurvival$NoStrata <- list()

for(iter_dataset in 1:3){
  results_initSurvival$NoStrata[[iter_dataset]] <- list()
  
  if(iter_dataset==1){ # death and censoring at the last event
    M.delta_Treatment <- cbind(time=c(1,0,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,0,1,0,0,0))
  }else if(iter_dataset==2){ # only death at the last event
    M.delta_Treatment <- cbind(time=c(1,0,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,0,1,0,1,1))
  }else if(iter_dataset==3){ # only death at the last event
    M.delta_Treatment <- cbind(time=c(1,1,1,1,1))
    M.delta_Control <- cbind(time=c(1,1,1,1,1,1,1))
  }
  
  for(method in c("Peto","Efron","Peron")){
    
    #### Compute survival
    if(method == "Efron"){ # necessary because efron forces the last observation to be an event
      Survival_noStrata <- BuyseTest:::initData(dataT=setNames(data.table(M.Treatment, M.delta_Treatment), c("time","censoring")),
                                  dataC=setNames(data.table(M.Control, M.delta_Control), c("time","censoring")), 
                                  D=1,
                                  type=3, endpoint="time", censoring="censoring", method="Efron", 
                                  index.strataT=list(seq(0,nrow(M.Treatment)-1)), 
                                  index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                                  D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE)
     }else{
      Survival_noStrata <- BuyseTest:::initSurvival(M.Treatment=M.Treatment, 
                                        M.Control=M.Control, 
                                        M.delta_Treatment=M.delta_Treatment, 
                                        M.delta_Control=M.delta_Control, 
                                        endpoint="time",
                                        index.strataT=list(seq(0,nrow(M.Treatment)-1)),index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                                        D.TTE=1, type=3, threshold=threshold, method=method)  
    }
    
    
    
    results_initSurvival$NoStrata[[iter_dataset]][[method]] <- Survival_noStrata
    
    # mark last events
    last.time <- which.max(c(M.Treatment,M.Control))
    indexC.last <- which(M.Control==last.time)
    indexT.last <- which(M.Treatment==last.time)
    
    #### jump at event (except last event)
    index.Tevent <- setdiff(which(M.delta_Treatment[,1]==1),indexT.last)
    index.Cevent <- setdiff(which(M.delta_Control[,1]==1),indexC.last)
    if(method == "Peto"){
      test.Tevent_before <- Survival_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT-threshold"] -  Survival_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT_0"] 
      test.Tevent_after <- Survival_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT_0"] -  Survival_noStrata$list_survivalT[[1]][index.Tevent,"Survival_TimeT+threshold"] 
      
      test.Cevent_before <- Survival_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC-threshold"] -  Survival_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC_0"] 
      test.Cevent_after <- Survival_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC_0"] -  Survival_noStrata$list_survivalC[[1]][index.Cevent,"Survival_TimeC+threshold"] 
    }else{
      test.Tevent_before <- Survival_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT-threshold"] -  Survival_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] 
      test.Tevent_after <- Survival_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT_0"] -  Survival_noStrata$list_survivalT[[1]][index.Tevent,"SurvivalT_TimeT+threshold"] 
      
      test.Cevent_before <- Survival_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC-threshold"] -  Survival_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] 
      test.Cevent_after <- Survival_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC_0"] -  Survival_noStrata$list_survivalC[[1]][index.Cevent,"SurvivalC_TimeC+threshold"] 
    }
    
    
    test_that(paste0("survival at t- noStrata, event, ",method), { ## jump at event time
      Vexpect_more_than(test.Tevent_before, 0)
      Vexpect_more_than(test.Cevent_before, 0)
    })
    
    test_that(paste0("survival at t+ noStrata, event, ",method), { ## no jump just after event time
      Vexpect_equal(test.Tevent_after, 0)
      Vexpect_equal(test.Cevent_after, 0)
    })
    
    #### jump at censoring 
    index.Tcensoring <- setdiff(which(M.delta_Treatment[,1]==0),indexT.last)
    index.Ccensoring <- setdiff(which(M.delta_Control[,1]==0),indexC.last)
    if(method == "Peto"){
      test.Tcensoring_before <- Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT-threshold"] -  Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT_0"] 
      test.Tcensoring_after <- Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT_0"] -  Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"Survival_TimeT+threshold"] 
      
      test.Ccensoring_before <- Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC-threshold"] -  Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC_0"] 
      test.Ccensoring_after <- Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC_0"] -  Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"Survival_TimeC+threshold"] 
    }else{
      test.Tcensoring_before <- Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT-threshold"] -  Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] 
      test.Tcensoring_after <- Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT_0"] -  Survival_noStrata$list_survivalT[[1]][index.Tcensoring,"SurvivalT_TimeT+threshold"] 
      
      test.Ccensoring_before <- Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC-threshold"] -  Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] 
      test.Ccensoring_after <- Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC_0"] -  Survival_noStrata$list_survivalC[[1]][index.Ccensoring,"SurvivalC_TimeC+threshold"] 
    }
    
     
    test_that(paste0("survival at t- noStrata, censoring, ",method), { ## no jump at event time
      Vexpect_equal(test.Tcensoring_before, 0)
      Vexpect_equal(test.Ccensoring_before, 0)
    })
    
    test_that(paste0("survival at t+ noStrata, censoring,  ",method), { ## no jump just after event time
      Vexpect_equal(test.Tcensoring_after, 0)
      Vexpect_equal(test.Ccensoring_after, 0)
    })
    
    #### last events
    if(method == "Peto"){
      test.censoringLastT <- test.censoringLastC <- any(c(M.delta_Treatment[indexT.last],M.delta_Control[indexC.last])==0)
      test.Tafter <-  Survival_noStrata$list_survivalT[[1]][indexT.last,"Survival_TimeT+threshold"]
      test.Cafter <-  Survival_noStrata$list_survivalC[[1]][indexC.last,"Survival_TimeC+threshold"]
    }else{
      test.censoringLastT <- any(M.delta_Treatment[indexT.last]==0)
      test.censoringLastC <- any(M.delta_Control[indexC.last]==0)
      
      test.Tafter <-  Survival_noStrata$list_survivalT[[1]][indexT.last,"SurvivalT_TimeT+threshold"]
      test.Cafter <-  Survival_noStrata$list_survivalC[[1]][indexC.last,"SurvivalC_TimeC+threshold"]
    }
    
    
      test_that(paste0("NA after last event if death - noStrata, ",method), {
        if(test.censoringLastT && method!="Efron"){Vexpect_NA(test.Tafter)}
        if(test.censoringLastC && method!="Efron"){Vexpect_NA(test.Cafter)}
      })
      
      test_that(paste0("0 after last event if censoing - noStrata, ",method), {# Efron always put 0 after last event
        if(test.censoringLastT == FALSE || method=="Efron"){Vexpect_equal(test.Tafter,0)}
        if(test.censoringLastC == FALSE || method=="Efron"){Vexpect_equal(test.Cafter,0)}
      })
      
  }
}
  
#### display
test <- FALSE
if(test){
  timesT_th <- c(M.Treatment[,1]-threshold,M.Treatment[,1],M.Treatment[,1]+threshold)
  plot(timesT_th[order(timesT_th)],as.vector(Survival_noStrata$list_survivalT[[1]])[order(timesT_th)],
       ylab="Survival Peto (Treatement)",xlab="time",ylim=c(0,1),
       type="o")
  #     type="p")
  
  Survival_noStrata$list_survivalT[[1]][,"Survival_TimeT-threshold"] # survival before jump points
  Survival_noStrata$list_survivalT[[1]][,"Survival_TimeT_0"] # survival at jump points
  Survival_noStrata$list_survivalT[[1]][,"Survival_TimeT+threshold"] # survival after jump points
  
  timesC_th <- c(M.Control[,1]-threshold,M.Control[,1],M.Control[,1]+threshold)
  plot(timesC_th[order(timesC_th)],as.vector(Survival_noStrata$list_survivalC[[1]])[order(timesC_th)],
       ylab="Survival Peto (Control)",xlab="time",ylim=c(0,1),
       type="o")
  # type="p")
  
  Survival_noStrata$list_survivalC[[1]][,"Survival_TimeC-threshold"] # survival before jump points
  Survival_noStrata$list_survivalC[[1]][,"Survival_TimeC_0"] # survival at jump points
  Survival_noStrata$list_survivalC[[1]][,"Survival_TimeC+threshold"] # survival after jump points
}

#### 2- With strata ####
n.strata <- 2
M.Treatment <- cbind(time=rep(1:5,n.strata))
M.Control <-  cbind(time=rep(c(1:5-0.1,5,5),n.strata))
M.delta_Treatment <- cbind(time=rep(c(1,0,1,1,1),n.strata))
M.delta_Control <- cbind(time=rep(c(1,1,0,1,0,0,0),n.strata))
threshold <- 0.001 
index.StrataT <- list(0:4,5:9)
index.StrataC <- list(0:6,7:13)

results_initSurvival$Strata <- list()

for(method in c("Peto","Efron","Peron")){

  if(method == "Efron"){ # necessary because efron forces the last observation to be an event
    Survival_Strata <- BuyseTest:::initData(dataT=setNames(data.table(M.Treatment, M.delta_Treatment), c("time","censoring")),
                                  dataC=setNames(data.table(M.Control, M.delta_Control), c("time","censoring")), 
                                  D=1,
                                  type=3, endpoint="time", censoring="censoring", method="Efron", 
                                  index.strataT=index.StrataT, 
                                  index.strataC=index.StrataC,n.strata=2,
                                  D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE)
  }else{
    Survival_Strata <- BuyseTest:::initSurvival(M.Treatment=M.Treatment, 
                                      M.Control=M.Control, 
                                      M.delta_Treatment=M.delta_Treatment, 
                                      M.delta_Control=M.delta_Control, 
                                      endpoint="time",
                                      index.strataT=index.StrataT,index.strataC=index.StrataC,n.strata=2,
                                      D.TTE=1, type=3, threshold=threshold, method=method)  
  }

    results_initSurvival$Strata[[method]] <- Survival_Strata
    
    #### identical between strata
    test_that(paste0("identical between strata, ",method),{
      expect_equal(Survival_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,],
                   Survival_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
      expect_equal(Survival_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,],
                   Survival_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
    })
   
}

#### Z- export ####
if(identical(save, TRUE)){
  saveRDS(results_initSurvival,file=file.path(dirSave,"test-initSurvival.rds"))
}else {
  GS <- readRDS(file=file.path(dirSave,"test-initSurvival.rds"))
  # expect_equal(results_initSurvival, GS)
}

# expect_equal(results_initSurvival$NoStrata[[2]]$Efron$M.Treatment, 
#              GS$NoStrata[[2]]$Efron$M.Treatment)
