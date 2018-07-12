if(FALSE){
    library(testthat)
    library(BuyseTest)
    library(data.table)
}

context("Check KM computation")

## * functions
initializeData <- BuyseTest:::initializeData
initializeSurvival_Peron <- BuyseTest:::initializeSurvival_Peron

calcSurvTest <- function(data, treatment, endpoint, strata, censoring,
                         type, D.TTE, threshold,
                         method.tte){
    level.treatment <- levels(as.factor(data[[treatment]])) ## added
    if(!is.null(strata)){
        n.strata <- length(unique(data[[strata]])) ## added
    } else {
        n.strata <- 1
    }
################ START: CODE COPIED FROM .BuyseTest.R #######################
    ## *** data: split the data according to the two levels
    indexT <- which(data[[treatment]] == level.treatment[2])
    indexC <- which(data[[treatment]] == level.treatment[1])

    ## *** data: extract endpoint 
    M.Treatment <- as.matrix(data[indexT,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
    M.Control <- as.matrix(data[indexC,endpoint,with=FALSE]) # matrix of endpoints for the control arm

    ## *** strata
    if(!is.null(strata)){
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
         index.strataT <- lapply(1:n.strata,function(iS){
            which(data[indexT,strata[1],with=FALSE] == iS) - 1
        })
        index.strataC <- lapply(1:n.strata,function(iS){
            which(data[indexC,strata[1],with=FALSE] == iS) - 1
        })
                                        
    }else{ # if there is no strata variable the same strata is used for all patient
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
        index.strataT <- list(0:(length(indexT)-1))
        index.strataC <- list(0:(length(indexC)-1))
    }

    ## *** data: censoring
    if(!is.null(censoring)){
        M.delta.Treatment <- as.matrix(data[indexT,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
        M.delta.Control <- as.matrix(data[indexC,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    }else{ # if the is no time to event variables
        M.delta.Treatment <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
        M.delta.Control <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }

    ## *** Update survival
    if(method.tte == 0){ ## Gehan
        outSurv <- list(list.survivalT = lapply(1:D.TTE, matrix),
                        list.survivalC = lapply(1:D.TTE, matrix))        
    }else {
        outSurv <- initializeSurvival_Peron(M.Treatment = M.Treatment,
                                            M.Control = M.Control,
                                            M.delta.Treatment = M.delta.Treatment,
                                            M.delta.Control = M.delta.Control,
                                            endpoint = endpoint,
                                            D.TTE = D.TTE,
                                            type = type,
                                            threshold = threshold,
                                            index.strataT = index.strataT,
                                            index.strataC = index.strataC,
                                            n.strata = n.strata)
    }
    ################ END: CODE COPIED FROM .BuyseTest.R #######################

    outSurv$index.strataT <- index.strataT
    outSurv$index.strataC <- index.strataC
    return(outSurv)
}


## * No strata

## ** settings
dataT <- data.table(time = 1:5,
                    treatment = "T",
                    status1 = c(1,0,1,1,1),
                    status2 = c(1,0,1,1,1),
                    status3 = c(1,1,1,1,1))
dataC <- data.table(time = c(1:5-0.1,5,5),
                    treatment = "C",
                    status1 = c(1,1,0,1,0,0,0),
                    status2 = c(1,1,0,1,0,1,1),
                    status3 = c(1,1,1,1,1,1,1))
data <- rbind(dataC, dataT)
threshold <- 0.001 # threshold smaller than the interval between two events

## ** tests
for(iData in 1:3){ ## iData <- 1
    
      ## *** Compute survival
      ## Cannot use initSurvival instead of initData because for method="Efron",
      ## initData needs to force the last observation to be an event
      outSurv <- calcSurvTest(data,
                              treatment = "treatment",
                              endpoint = "time",
                              strata = NULL,
                              censoring = paste0("status",iData),
                              type = 3,
                              D.TTE = 1,
                              threshold = 1e-12,
                              method.tte = 1
                              )
    
    ## *** index relative to each group and the last events
    col.Surv <- c(T = "SurvivalT_TimeT", C = "SurvivalC_TimeC")

    last.time <- max(data$time)
    indexT.last <- data[treatment == "T",which(time==last.time)]
    indexC.last <- data[treatment == "C",which(time==last.time)]

    index.Tevent <- setdiff(which(data[treatment == "T"][[paste0("status",iData)]]==1),indexT.last)
    index.Cevent <- setdiff(which(data[treatment == "C"][[paste0("status",iData)]]==1),indexC.last)

    index.Tcensoring <- setdiff(which(data[treatment == "T"][[paste0("status",iData)]]==0),indexT.last)
    index.Ccensoring <- setdiff(which(data[treatment == "C"][[paste0("status",iData)]]==0),indexC.last)

    ## *** [Test] jump at event (except last event)
    Surv_before <- list(T = outSurv$list.survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"-threshold")],
                        C = outSurv$list.survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"-threshold")])
    Surv_at <- list(T = outSurv$list.survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"_0")],
                    C = outSurv$list.survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"_0")])
    Surv_after <- list(T = outSurv$list.survivalT[[1]][index.Tevent,paste0(col.Surv["T"],"+threshold")],
                       C = outSurv$list.survivalC[[1]][index.Cevent,paste0(col.Surv["C"],"+threshold")])
    
    test_that("jump at event time (t-, noStrata), Peron", { ## jump at event time
        expect_true(all(Surv_before[["T"]] > Surv_at[["T"]]))
        expect_true(all(Surv_before[["C"]] > Surv_at[["C"]]))
    })
    
    test_that("no jump after event time (t+, noStrata), Peron", { ## no jump just after event time
        expect_true(all(Surv_at[["T"]] == Surv_after[["T"]]))
        expect_true(all(Surv_at[["C"]] == Surv_after[["C"]]))
    })

    test_that("check survival at event time (noStrata), Peron", { 
        GS <- switch(as.character(iData),
                     "1" = list(T = c(0.8000000, 0.5333333, 0.2666667), C = c(0.8571429, 0.7142857, 0.5357143)),
                     "2" = list(T = c(0.8000000, 0.5333333, 0.2666667), C = c(0.8571429, 0.7142857, 0.5357143)),
                     "3" = list(T = c(0.8, 0.6, 0.4, 0.2), C = c(0.8571429, 0.7142857, 0.5714286, 0.4285714, 0.2857143)))
        expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
        expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
    })

    ## *** [Test] jump at censoring 
    Surv_before <- list(T = outSurv$list.survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"-threshold")],
                        C = outSurv$list.survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"-threshold")])
    Surv_at <- list(T = outSurv$list.survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"_0")],
                    C = outSurv$list.survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"_0")])
    Surv_after <- list(T = outSurv$list.survivalT[[1]][index.Tcensoring,paste0(col.Surv["T"],"+threshold")],
                       C = outSurv$list.survivalC[[1]][index.Ccensoring,paste0(col.Surv["C"],"+threshold")])
    
    test_that("no jump at censoring time (t-, noStrata), Peron", { ## no jump at censoring time
      expect_true(all(Surv_before[["T"]] == Surv_at[["T"]]))
      expect_true(all(Surv_before[["C"]] == Surv_at[["C"]]))
    })
    
    test_that("no jump after censoring time (t+, noStrata), Peron", { ## no jump just after censoring time
      expect_true(all(Surv_at[["T"]] == Surv_after[["T"]]))
      expect_true(all(Surv_at[["C"]] == Surv_after[["C"]]))
    })
    
    test_that("check survival at censoring time (noStrata), Peron", { 
            GS <- switch(as.character(iData),
                         "1" = list(T = 0.8, C = c(0.7142857, 0.5357143)),
                         "2" = list(T = 0.8, C = c(0.7142857, 0.5357143)),
                         "3" = list(T = numeric(0), C = numeric(0)))
            expect_equal(unname(Surv_at[["T"]]), GS[["T"]], tol = 1e-6)
            expect_equal(unname(Surv_at[["C"]]), GS[["C"]], tol = 1e-6)
    })
    
    ## *** last events
    lastData <- data[time == max(time)]
    lastDataT <- lastData[treatment == "T"]
    lastDataC <- lastData[treatment == "C"]
    Surv_after <- list(T = outSurv$list.survivalT[[1]][indexT.last,paste0(col.Surv["T"],"+threshold")],
                       C = outSurv$list.survivalC[[1]][indexC.last,paste0(col.Surv["C"],"+threshold")])
    
    test.censoringLastT <- any(lastDataT[[paste0("status",iData)]]==0)
    test.censoringLastC <- any(lastDataC[[paste0("status",iData)]]==0)
        
    test_that("NA after last event if censoring - noStrata, Peron", {        
        if(test.censoringLastT){expect_true(all(is.na(Surv_after[["T"]])))}
        if(test.censoringLastC){expect_true(all(is.na(Surv_after[["C"]])))}
    })
      
    test_that("0 after last event if death - noStrata, Peron", {
        if(test.censoringLastT == FALSE){expect_true(all(Surv_after[["T"]]==0))}
        if(test.censoringLastC == FALSE){expect_true(all(Surv_after[["C"]]==0))}
    })
      
}
    
## * With strata

## ** settings
threshold <- 0.001

dataT <- data.table(time = 1:5,
                    treatment = "T",
                    status = c(1,0,1,1,1)
                    )
dataC <- data.table(time = c(1:5-0.1,5,5),
                    treatment = "C",
                    status = c(1,1,0,1,0,0,0)
                    )
data <- rbind(dataC, dataT)
dataStrata <- rbind(cbind(data, strata = 1),
                    cbind(data, strata = 2)
                    )

## ** tests
## *** Compute survival
## Cannot use initSurvival instead of initData because for method="Efron",
## initData needs to force the last observation to be an event
outSurv <- calcSurvTest(dataStrata,
                        treatment = "treatment",
                        endpoint = "time",
                        strata = "strata",
                        censoring = "status",
                        type = 3,
                        D.TTE = 1,
                        threshold = 1e-12,
                        method.tte = 1
                        )

## *** [test] identical survival between strata
test_that("identical between strata, Peron",{
    expect_equal(outSurv$list.survivalT[[1]][outSurv$index.strataT[[1]]+1,],
                 outSurv$list.survivalT[[1]][outSurv$index.strataT[[2]]+1,])
    expect_equal(outSurv$list.survivalC[[1]][outSurv$index.strataC[[1]]+1,],
                 outSurv$list.survivalC[[1]][outSurv$index.strataC[[2]]+1,])
})

## *** [test] check survival
test_that("check survival (strata),  Peron", {
    ## to import the expected values:
    ## rownames(Survival_Strata$list_survivalT[[1]]) <- NULL
    ## butils::object2script(Survival_Strata$list_survivalT[[1]], digit = 6)
    ## rownames(Survival_Strata$list_survivalC[[1]]) <- NULL
    ## butils::object2script(Survival_Strata$list_survivalC[[1]], digit = 6)
        
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
           
    expect_equal(unname(outSurv$list.survivalT[[1]]), unname(GS.list_survivalT), tol = 1e-5)
    expect_equal(unname(outSurv$list.survivalC[[1]]), unname(GS.list_survivalC), tol = 1e-5)
})




