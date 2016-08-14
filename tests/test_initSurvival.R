#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test the computation of the survival
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# no strata: Check jump at t- and no jump at t+ in the estimated survival by KM
# no strata: Check survival=NA after last event if censored else death (rational being that we do not know what it would be for censored but we do know that the survival if 0 if everybody is dead)
# no strata and strata: check identical with previous version
#

## no test : saut au moment du dernier evenement.

#### spec
library(BuyseTest)
library(testthat)
library(data.table)
library(lava)

precision <- 10^{-7}
n.patients <- 200
save <- NULL # TRUE to save results, FALSE to test, NULL to ignore

#### data ####
version <- packageVersion("BuyseTest")
dir <- paste0("tests/Results-version",version)

if(identical(save, TRUE)){
  results_initSurvival <- list()
  results_initSurvival$NoStrata <- list()
  results_initSurvival$Strata <- list()
}else if(identical(save, FALSE)){
  load(file = file.path(dir,"test_initSurvival-results_initSurvival.RData"))
}

#### 0- function - file FCT/FCT_check.R ####
validPairs <- function(BuyseRes, type = c("strata","sum")){
  
  BuyseSummary <- summary(BuyseRes, show = NULL)
  enpoint_threshold <- paste(BuyseSummary$nb$endpoint,BuyseSummary$nb$threshold, sep = "_")
  endpoints <- unique(enpoint_threshold)
  D <- length(endpoints)
  index_strata <- which(BuyseSummary$nb$strata!="global")
  levels.strata <- unique(BuyseSummary$nb$strata[index_strata])
  n.strata <- length(levels.strata) 
  
  diff <- NULL
  
  if("strata" %in% type){
    diff.strata <- matrix(NA,nrow=(1+n.strata)*D,ncol=5)
    colnames(diff.strata) <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
    rownames(diff.strata) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    for(iter_endpoint in 1:D){
      index_endpoint <- which(enpoint_threshold==endpoints[iter_endpoint])  
      res_strata <- as.matrix(BuyseSummary$nb[intersect(index_strata,index_endpoint),c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")])
      diff.strata[intersect(index_strata,index_endpoint),] <- t(apply(res_strata,1,function(x){x-res_strata[1,,drop=TRUE]}))
    }
    diff <- cbind(diff, diff.strata)
  }
  if("sum" %in% type){
    diff.sum <- cbind(sum = BuyseSummary$nb$n.total - (BuyseSummary$nb$n.favorable+BuyseSummary$nb$n.unfavorable+BuyseSummary$nb$n.neutral+BuyseSummary$nb$n.uninf))
    rownames(diff.sum) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    diff <- cbind(diff, diff.sum)
  }
  return(diff)
}

expect_equalPairsBT <- function(BuyseRes1, BuyseRes2){
  count1 <- getCount(BuyseRes1) 
  count2 <- getCount(BuyseRes2)
  expect_equal(count1, count2)
}

expect_equalBT <- function(BuyseRes1, BuyseRes2, slots = NULL, trace = 1){
  if(is.null(slots)){slots <- setdiff(intersect(names(attributes(BuyseRes1)), names(attributes(BuyseRes2))), "class")}
  
  test.strata1 <- !("try-error" %in% class(try(BuyseRes1@strata, silent = TRUE)))
  test.strata2 <- !("try-error" %in% class(try(BuyseRes2@strata, silent = TRUE)))
  
  if (test.strata1 && length(BuyseRes1@strata) > 1) {
    if(test.strata2 == FALSE){
      
      if(!identical(BuyseRes1@delta, BuyseRes2@delta)){
        if(trace>0){cat("* reorder according to strata in BuyseRes 1 and 2 \n")}
        index.match <- match(BuyseRes1@delta[,1], BuyseRes2@delta[,1]) # could be wrong to consider only the first outcome
        BuyseRes2@delta <- BuyseRes2@delta[index.match,,drop=FALSE]
        BuyseRes2@count_favorable <- BuyseRes2@count_favorable[index.match,,drop=FALSE]
        BuyseRes2@count_unfavorable <- BuyseRes2@count_unfavorable[index.match,,drop=FALSE]
        BuyseRes2@count_neutral <- BuyseRes2@count_neutral[index.match,,drop=FALSE]
        BuyseRes2@count_uninf <- BuyseRes2@count_uninf[index.match,,drop=FALSE]
        BuyseRes2@delta_boot <- BuyseRes2@delta_boot[index.match,,,drop = FALSE]
        
        BuyseRes2@index_neutralT <- sort(BuyseRes2@index_neutralT)
        BuyseRes2@index_neutralC <- sort(BuyseRes2@index_neutralC)
        BuyseRes2@index_uninfT <- sort(BuyseRes2@index_uninfT)
        BuyseRes2@index_uninfC <- sort(BuyseRes2@index_uninfC)
        
        BuyseRes1@index_neutralT <- sort(BuyseRes1@index_neutralT)
        BuyseRes1@index_neutralC <- sort(BuyseRes1@index_neutralC)
        BuyseRes1@index_uninfT <- sort(BuyseRes1@index_uninfT)
        BuyseRes1@index_uninfC <- sort(BuyseRes1@index_uninfC)
      }
    }
    
  }
  
  test.error <- FALSE
  for(iterSlot in slots){
    res <- try(expect_equal(slot(BuyseRes1,  iterSlot), slot(BuyseRes2,  iterSlot)), silent = TRUE)
    if("try-error" %in% class(res) ){
      test.error <- TRUE
      cat("Differences in slot: ",iterSlot,"\n")
      if(trace == 2){cat(res[1])}
    }
  }
  if(test.error){stop("difference between the slots of the two BuyseRes objects \n")}
}


Vexpect_less_than <- function(x,y,...){
  sapply(x, function(X){expect_less_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_more_than <- function(x,y,...){
  sapply(x, function(X){expect_more_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_equal <- function(x,y,...){
  sapply(x, function(X){expect_equal(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_NA <- function(x,...){
  sapply(x, function(X){expect_true(is.na(X),...)})
  return(invisible(TRUE))
}

#### 1- No strata ####
M.Treatment <- cbind(time=1:5)
M.Control <-  cbind(time=c(1:5-0.1,5,5))
threshold <- 0.001 # threshold smaller than the interval between two events

for(iter_dataset in 1:3){
  if(identical(save,TRUE)){results_initSurvival$NoStrata[[iter_dataset]] <- list()}
  
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
    
    
    
    if(identical(save, TRUE)){results_initSurvival$NoStrata[[iter_dataset]][[method]] <- Survival_noStrata ; break}
    
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
    
    #### match previous implementation
      if(!is.null(save)){
        test_that("comparison with the previous version", {
          expect_equal(results_initSurvival$NoStrata[[iter_dataset]][[method]]$list_survivalT, Survival_noStrata$list_survivalT, tolerance = precision)
          expect_equal(results_initSurvival$NoStrata[[iter_dataset]][[method]]$list_survivalC, Survival_noStrata$list_survivalC, tolerance = precision)
        })
      }
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

    if(identical(save, TRUE)){results_initSurvival$Strata[[method]] <- Survival_Strata ; break;}
    
    #### identical between strata
    test_that(paste0("identical between strata, ",method),{
      expect_equal(Survival_Strata$list_survivalT[[1]][index.StrataT[[1]]+1,],
                   Survival_Strata$list_survivalT[[1]][index.StrataT[[2]]+1,])
      expect_equal(Survival_Strata$list_survivalC[[1]][index.StrataC[[1]]+1,],
                   Survival_Strata$list_survivalC[[1]][index.StrataC[[2]]+1,])
    })
    
    #### match previous implementation
    if(!is.null(save)){
      test_that(paste0("comparison with the previous version ",method), {
        expect_equal(results_initSurvival$Strata[[method]]$list_survivalT, 
                     Survival_Strata$list_survivalT, tolerance = precision)
        expect_equal(results_initSurvival$Strata[[method]]$list_survivalC, 
                     Survival_Strata$list_survivalC, tolerance = precision)
      })
    }
}
#### Z- export ####

if(identical(save, TRUE)){
  save(results_initSurvival,file="test_initSurvival-results_initSurvival.RData")
}
