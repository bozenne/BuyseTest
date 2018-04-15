if(FALSE){
    library(testthat)
    library(BuyseTest)
}

context("Check KM computation")

## * functions
initSurvival <- BuyseTest:::initSurvival
initData <- BuyseTest:::initData

## * No strata

## ** settings
M.Treatment <- cbind(time=1:5)
M.Control <-  cbind(time=c(1:5-0.1,5,5))
threshold <- 0.001 # threshold smaller than the interval between two events

## ** tests
for(iter_dataset in 1:3){

  ## *** define dataset
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
    method.num <- which(c("Gehan","Peto","Efron","Peron") == method) - 1
    
    ## *** Compute survival
    ## Cannot use initSurvival instead of initData
    ## because for method="Efron", initData needs to force the last observation to be an event
    iSurv <- initData(dataT=setNames(data.table(M.Treatment, M.delta_Treatment), c("time","censoring")),
                      dataC=setNames(data.table(M.Control, M.delta_Control), c("time","censoring")), 
                      D=1,
                      type=3, endpoint="time", censoring="censoring", method=method.num, 
                      index.strataT=list(seq(0,nrow(M.Treatment)-1)), 
                      index.strataC=list(seq(0,nrow(M.Control)-1)),n.strata=1,
                      D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE)
    
    ## *** index relative to each group and the last events
    if(method == "Peto"){
        col.Surv <- c(T = "Survival_TimeT", C = "Survival_TimeC")
    }else{
        col.Surv <- c(T = "SurvivalT_TimeT", C = "SurvivalC_TimeC")
    }

    last.time <- which.max(c(M.Treatment,M.Control))
    indexC.last <- which(M.Control==last.time)
    indexT.last <- which(M.Treatment==last.time)

    index.Tevent <- setdiff(which(M.delta_Treatment[,1]==1),indexT.last)
    index.Cevent <- setdiff(which(M.delta_Control[,1]==1),indexC.last)

    index.Tcensoring <- setdiff(which(M.delta_Treatment[,1]==0),indexT.last)
    index.Ccensoring <- setdiff(which(M.delta_Control[,1]==0),indexC.last)

    ## *** [Test] jump at event (except last event)
    Surv_before <- list(T = iSurv$list_survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"-threshold")],
                        C = iSurv$list_survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"-threshold")])
    Surv_at <- list(T = iSurv$list_survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"_0")],
                    C = iSurv$list_survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"_0")])
    Surv_after <- list(T = iSurv$list_survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"+threshold")],
                       C = iSurv$list_survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"+threshold")])
    
    test_that(paste0("jump at event time (t-, noStrata), ",method), { ## jump at event time
        expect_true(all(Surv_before[["T"]] > Surv_at[["T"]]))
        expect_true(all(Surv_before[["C"]] > Surv_at[["C"]]))
    })
    
    test_that(paste0("no jump after event time (t+, noStrata),  ",method), { ## no jump just after event time
        expect_true(all(Surv_at[["T"]] == Surv_after[["T"]]))
        expect_true(all(Surv_at[["C"]] == Surv_after[["C"]]))
    })

    test_that(paste0("check survival at event time (noStrata),  ",method), { 
        if(method == "Peto"){
            GS <- switch(as.character(iter_dataset),
                         "1" = list(T = c(0.8333333, 0.6428571, 0.4285714), C = c(0.9166667, 0.7500000, 0.5357143)),
                         "2" = list(T = c(0.8333333, 0.6428571, 0.4285714), C = c(0.9166667, 0.7500000, 0.5357143)),
                         "3" = list(T = c(0.8333333, 0.6666667, 0.5000000, 0.3333333), C = c(0.9166667, 0.7500000, 0.5833333, 0.4166667, 0.2500000)))
            expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
            expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
        }else if(method %in% c("Efron","Peron")){
            GS <- switch(as.character(iter_dataset),
                         "1" = list(T = c(0.8000000, 0.5333333, 0.2666667), C = c(0.8571429, 0.7142857, 0.5357143)),
                         "2" = list(T = c(0.8000000, 0.5333333, 0.2666667), C = c(0.8571429, 0.7142857, 0.5357143)),
                         "3" = list(T = c(0.8, 0.6, 0.4, 0.2), C = c(0.8571429, 0.7142857, 0.5714286, 0.4285714, 0.2857143)))
            expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
            expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
        }
    })

    ## *** [Test] jump at censoring 
    Surv_before <- list(T = iSurv$list_survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"-threshold")],
                        C = iSurv$list_survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"-threshold")])
    Surv_at <- list(T = iSurv$list_survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"_0")],
                    C = iSurv$list_survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"_0")])
    Surv_after <- list(T = iSurv$list_survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"+threshold")],
                       C = iSurv$list_survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"+threshold")])
    
    test_that(paste0("no jump at censoring time (t-, noStrata), ",method), { ## no jump at censoring time
      expect_true(all(Surv_before[["T"]] == Surv_at[["T"]]))
      expect_true(all(Surv_before[["C"]] == Surv_at[["C"]]))
    })
    
    test_that(paste0("no jump after censoring time (t+, noStrata),  ",method), { ## no jump just after censoring time
      expect_true(all(Surv_at[["T"]] == Surv_after[["T"]]))
      expect_true(all(Surv_at[["C"]] == Surv_after[["C"]]))
    })
    
    test_that(paste0("check survival at censoring time (noStrata),  ",method), { 
        if(method == "Peto"){
            GS <- switch(as.character(iter_dataset),
                         "1" = list(T = 0.75, C = c(0.75, 0.4285714)),
                         "2" = list(T = 0.75, C = c(0.75, 0.4285714)),
                         "3" = list(T = numeric(0), C = numeric(0)))
            expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
            expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
        }else if(method %in% c("Efron","Peron")){
            GS <- switch(as.character(iter_dataset),
                         "1" = list(T = 0.8, C = c(0.7142857, 0.5357143)),
                         "2" = list(T = 0.8, C = c(0.7142857, 0.5357143)),
                         "3" = list(T = numeric(0), C = numeric(0)))
            expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
            expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
        }
    })
    
    ## *** last events
    Surv_after <- list(T = iSurv$list_survivalT[[1]][indexT.last,paste0(col.Surv["T"],"+threshold")],
                       C = iSurv$list_survivalC[[1]][indexC.last,paste0(col.Surv["C"],"+threshold")])
    
    if(method == "Peto"){ ## same survival curves
        test.censoringLastT <- test.censoringLastC <- any(c(M.delta_Treatment[indexT.last],M.delta_Control[indexC.last])==0)
    }else{ ##  stratified survival curves
        test.censoringLastT <- any(M.delta_Treatment[indexT.last]==0)
        test.censoringLastC <- any(M.delta_Control[indexC.last]==0)      
    }
        
    test_that(paste0("NA after last event if death - noStrata, ",method), {        
        if(test.censoringLastT && method!="Efron"){expect_true(all(is.na(Surv_after[["T"]])))}
        if(test.censoringLastC && method!="Efron"){expect_true(all(is.na(Surv_after[["C"]])))}
    })
      
    test_that(paste0("0 after last event if censoing - noStrata, ",method), {# Efron always put 0 after last event
        if(test.censoringLastT == FALSE || method=="Efron"){expect_true(all(Surv_after[["T"]]==0))}
        if(test.censoringLastC == FALSE || method=="Efron"){expect_true(all(Surv_after[["C"]]==0))}
    })
      
  }
}
  
## * With strata

## ** settings
n.strata <- 2
M.Treatment <- cbind(time=rep(1:5,n.strata))
M.Control <-  cbind(time=rep(c(1:5-0.1,5,5),n.strata))
M.delta_Treatment <- cbind(time=rep(c(1,0,1,1,1),n.strata))
M.delta_Control <- cbind(time=rep(c(1,1,0,1,0,0,0),n.strata))
threshold <- 0.001 
index.StrataT <- list(0:4,5:9)
index.StrataC <- list(0:6,7:13)

## ** tests
for(method in c("Peto","Efron","Peron")){
    method.num <- which(c("Gehan","Peto","Efron","Peron") == method) - 1

  ## *** Compute survival
  ## Cannot use initSurvival instead of initData
  ## because for method="Efron", initData needs to force the last observation to be an event
  Survival_Strata <- initData(dataT=setNames(data.table(M.Treatment, M.delta_Treatment), c("time","censoring")),
                              dataC=setNames(data.table(M.Control, M.delta_Control), c("time","censoring")), 
                              D=1,
                              type=3, endpoint="time", censoring="censoring", method=method.num, 
                              index.strataT=index.StrataT, 
                              index.strataC=index.StrataC,n.strata=2,
                              D.TTE=1, threshold=threshold, Wscheme = NULL, trace=TRUE, test = TRUE)

  ## *** [test] identical survival between strata
  test_that(paste0("identical between strata, ",method),{
      expect_equal(Survival_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,],
                   Survival_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
      expect_equal(Survival_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,],
                   Survival_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
  })

  ## *** [test] check survival
    test_that(paste0("check survival (strata),  ",method), {
        ## to import the expected values:
        ## rownames(Survival_Strata$list_survivalT[[1]]) <- NULL
        ## butils::object2script(Survival_Strata$list_survivalT[[1]], digit = 6)
        ## rownames(Survival_Strata$list_survivalC[[1]]) <- NULL
        ## butils::object2script(Survival_Strata$list_survivalC[[1]], digit = 6)
      if(method == "Peto"){
          GS.list_survivalT <- cbind(c("Survival_TimeT-threshold" = 0.916667, 0.75, 0.75, 0.535714, 0.428571, 0.916667, 0.75, 0.75, 0.535714, 0.428571), 
                                     c("Survival_TimeT_0" = 0.833333, 0.75, 0.642857, 0.428571, 0.285714, 0.833333, 0.75, 0.642857, 0.428571, 0.285714), 
                                     c("Survival_TimeT+threshold" = 0.833333, 0.75, 0.642857, 0.428571, NA, 0.833333, 0.75, 0.642857, 0.428571, NA)
                                     )
          GS.list_survivalC <- cbind(c("Survival_TimeC-threshold" = 1, 0.833333, 0.75, 0.642857, 0.428571, 0.428571, 0.428571, 1, 0.833333, 0.75, 0.642857, 0.428571, 0.428571, 0.428571), 
                                     c("Survival_TimeC_0" = 0.916667, 0.75, 0.75, 0.535714, 0.428571, 0.285714, 0.285714, 0.916667, 0.75, 0.75, 0.535714, 0.428571, 0.285714, 0.285714), 
                                     c("Survival_TimeC+threshold" = 0.916667, 0.75, 0.75, 0.535714, 0.428571, NA, NA, 0.916667, 0.75, 0.75, 0.535714, 0.428571, NA, NA)
                                     )
         
      }else if(method == "Efron"){
          GS.list_survivalT <- cbind("SurvivalT_TimeT-threshold" = c(1, 0.8, 0.8, 0.53333, 0.26667, 1, 0.8, 0.8, 0.53333, 0.26667), 
                                     "SurvivalT_TimeT_0" = c(0.8, 0.8, 0.53333, 0.26667, 0, 0.8, 0.8, 0.53333, 0.26667, 0), 
                                     "SurvivalT_TimeT+threshold" = c(0.8, 0.8, 0.53333, 0.26667, 0, 0.8, 0.8, 0.53333, 0.26667, 0), 
                                     "SurvivalC_TimeT-threshold" = c(0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571), 
                                     "SurvivalC_TimeT_0" = c(0.85714, 0.71429, 0.71429, 0.53571, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0), 
                                     "SurvivalC_TimeT+threshold" = c(0.85714, 0.71429, 0.71429, 0.53571, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0), 
                                     "SurvivalT_TimeT_0(ordered)" = c(0.8, 0.8, 0.53333, 0.26667, 0, 0.8, 0.8, 0.53333, 0.26667, 0), 
                                     "SurvivalC_TimeT-threshold(ordered)" = c(0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571), 
                                     "SurvivalC_TimeT+threshold(ordered)" = c(0.85714, 0.71429, 0.71429, 0.53571, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0), 
                                     "time_control(ordered)" = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
                                     "events(ordered)" = c(1, 0, 1, 1, 1, 1, 0, 1, 1, 1)
                                     )
          GS.list_survivalC <- cbind("SurvivalC_TimeC-threshold" = c(1, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0.53571, 1, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0.53571), 
                                     "SurvivalC_TimeC_0" = c(0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0), 
                                     "SurvivalC_TimeC+threshold" = c(0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0), 
                                     "SurvivalT_TimeC-threshold" = c(1, 0.8, 0.8, 0.53333, 0.26667, 0.26667, 0.26667, 1, 0.8, 0.8, 0.53333, 0.26667, 0.26667, 0.26667), 
                                     "SurvivalT_TimeC_0" = c(1, 0.8, 0.8, 0.53333, 0.26667, 0, 0, 1, 0.8, 0.8, 0.53333, 0.26667, 0, 0), 
                                     "SurvivalT_TimeC+threshold" = c(1, 0.8, 0.8, 0.53333, 0.26667, 0, 0, 1, 0.8, 0.8, 0.53333, 0.26667, 0, 0), 
                                     "SurvivalC_TimeC_0(ordered)" = c(0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0, 0.85714, 0.71429, 0.71429, 0.53571, 0.53571, 0, 0), 
                                     "SurvivalT_TimeC-threshold(ordered)" = c(1, 0.8, 0.8, 0.53333, 0.26667, 0.26667, 0.26667, 1, 0.8, 0.8, 0.53333, 0.26667, 0.26667, 0.26667), 
                                     "SurvivalT_TimeC+threshold(ordered)" = c(1, 0.8, 0.8, 0.53333, 0.26667, 0, 0, 1, 0.8, 0.8, 0.53333, 0.26667, 0, 0), 
                                     "time_control(ordered)" = c(0.9, 1.9, 2.9, 3.9, 4.9, 5, 5, 0.9, 1.9, 2.9, 3.9, 4.9, 5, 5), 
                                     "events(ordered)" = c(1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1)
                                     ) 
      }else if(method == "Peron"){
          GS.list_survivalT <-  cbind(c("SurvivalT_TimeT-threshold" = 1, 0.8, 0.8, 0.533333, 0.266667, 1, 0.8, 0.8, 0.533333, 0.266667), 
                                      c("SurvivalT_TimeT_0" = 0.8, 0.8, 0.533333, 0.266667, 0, 0.8, 0.8, 0.533333, 0.266667, 0), 
                                      c("SurvivalT_TimeT+threshold" = 0.8, 0.8, 0.533333, 0.266667, 0, 0.8, 0.8, 0.533333, 0.266667, 0), 
                                      c("SurvivalC_TimeT-threshold" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714), 
                                      c("SurvivalC_TimeT_0" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714), 
                                      c("SurvivalC_TimeT+threshold" = 0.857143, 0.714286, 0.714286, 0.535714, NA, 0.857143, 0.714286, 0.714286, 0.535714, NA), 
                                      c("SurvivalT_TimeT_0(ordered)" = 0.8, 0.8, 0.533333, 0.266667, 0, 0.8, 0.8, 0.533333, 0.266667, 0), 
                                      c("SurvivalC_TimeT-threshold(ordered)" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714), 
                                      c("SurvivalC_TimeT+threshold(ordered)" = 0.857143, 0.714286, 0.714286, 0.535714, NA, 0.857143, 0.714286, 0.714286, 0.535714, NA), 
                                      c("time_control(ordered)" = 1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
                                      c("events(ordered)" = 1, 0, 1, 1, 1, 1, 0, 1, 1, 1)
                                      )

          GS.list_survivalC <- cbind(c("SurvivalC_TimeC-threshold" = 1, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714, 1, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714), 
                                     c("SurvivalC_TimeC_0" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714, 0.535714, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714, 0.535714), 
                                     c("SurvivalC_TimeC+threshold" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, NA, NA, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, NA, NA), 
                                     c("SurvivalT_TimeC-threshold" = 1, 0.8, 0.8, 0.533333, 0.266667, 0.266667, 0.266667, 1, 0.8, 0.8, 0.533333, 0.266667, 0.266667, 0.266667), 
                                     c("SurvivalT_TimeC_0" = 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0, 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0), 
                                     c("SurvivalT_TimeC+threshold" = 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0, 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0), 
                                     c("SurvivalC_TimeC_0(ordered)" = 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714, 0.535714, 0.857143, 0.714286, 0.714286, 0.535714, 0.535714, 0.535714, 0.535714), 
                                     c("SurvivalT_TimeC-threshold(ordered)" = 1, 0.8, 0.8, 0.533333, 0.266667, 0.266667, 0.266667, 1, 0.8, 0.8, 0.533333, 0.266667, 0.266667, 0.266667), 
                                     c("SurvivalT_TimeC+threshold(ordered)" = 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0, 1, 0.8, 0.8, 0.533333, 0.266667, 0, 0), 
                                     c("time_control(ordered)" = 0.9, 1.9, 2.9, 3.9, 4.9, 5, 5, 0.9, 1.9, 2.9, 3.9, 4.9, 5, 5), 
                                     c("events(ordered)" = 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0)
                                     )
      }      
      expect_equal(unname(Survival_Strata$list_survivalT[[1]]), unname(GS.list_survivalT), tol = 1e-5)
      expect_equal(unname(Survival_Strata$list_survivalC[[1]]), unname(GS.list_survivalC), tol = 1e-5)
  })

}

