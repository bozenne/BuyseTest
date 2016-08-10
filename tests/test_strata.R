#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with strata
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# check whether the number of pairs displayed in the summary is ok (e.g. total = favorable + unfavorable + neutral + uninf)
# check whether when using identical values across strata, the result is the same across strata
# check consistency with results from previous version (on common slots)

# PB with the bootstrap !!!!!!

#### spec
# library(BuyseTest) # butils:::package.source("BuyseTest", Rcode = TRUE, Ccode = TRUE)
library(testthat)
library(lava)
library(data.table)
precision <- 10^{-7}
n.patients <- 200
n.bootstrap <- 10
save <- FALSE


#### function
source("tests/FCT_check.R")

#### 1- Check number of pairs ####
set.seed(10)
n.patients <- 10
argsBin <- list(p.T = c(0.5,0.75))
argsCont <- list(mu.T = 1:3, sigma.T = rep(1,3))
argsTTE <- list(rates.T = 1:3, rates.Censor = rep(1,3))

#### binary #### 
data_Bin <- simulBT(n.patients, argsBin = argsBin, argsCont = NULL, argsTTE = NULL)
data_BinS <- rbind(cbind(data_Bin, strata = 1),
                       cbind(data_Bin, strata = 2),
                       cbind(data_Bin, strata = 3))


BT_Bin1 <- BuyseTest(data=data_BinS,endpoint=c("Y_bin1","Y_bin2"),
                                treatment="Treatment", type=c("bin","bin"),strata="strata",
                                n.bootstrap=n.bootstrap,trace=0)
summary(BT_Bin1)

summary(BT_Bin1)

test_that("count pairs summary - Binary",{
  valTest <- as.double(validPairs(BT_Bin1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Binary",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Bin1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### continuous ####
set.seed(10)
data_Cont <- simulBT(n.patients, argsBin = NULL, argsCont = argsCont, argsTTE = NULL)
data_ContS <- rbind(cbind(data_Cont, strata = 1),
                        cbind(data_Cont, strata = 2),
                        cbind(data_Cont, strata = 3))

BT_Cont1 <- BuyseTest(data=data_ContS,endpoint=c("Y_cont1","Y_cont2"),
                                 treatment="Treatment", type=c("cont","cont"),threshold=c(1,1),strata="strata",
                                 n.bootstrap=n.bootstrap,trace=0)

summary(BT_Cont1)

test_that("count pairs summary - Continuous",{
  valTest <- as.double(validPairs(BT_Cont1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Continuous",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Cont1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### TTE ####
set.seed(10)
data_TTE <- simulBT(n.patients, argsBin = NULL, argsCont = NULL, argsTTE = argsTTE)
data_TTES <- rbind(cbind(data_TTE, strata = 1),
                   cbind(data_TTE, strata = 2),
                   cbind(data_TTE, strata = 3))

BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE1[[method]] <- BuyseTest(data=data_TTES,endpoint=c("Y_TTE1","Y_TTE2","Y_TTE3"),method=method,
                                 treatment="Treatment",censoring=c("event1","event2","event1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),
                                 n.bootstrap = n.bootstrap, trace=0)
  (BT_TTE1[[method]])
  
  test_that("count pairs summary - TTE",{
    valTest <- as.double(validPairs(BT_TTE1[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - TTE",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_TTE1[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### mixed outcomes ####
set.seed(10)
data_Mix <- simulBT(n.patients, argsBin = argsBin, argsCont = argsCont, argsTTE = argsTTE)
data_MixS <- rbind(cbind(data_Mix, strata = 1),
                   cbind(data_Mix, strata = 2),
                   cbind(data_Mix, strata = 3))

#### without strata ####
BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_Mix[[method]] <- BuyseTest(data=data_MixS,
                                endpoint=c("Y_TTE1","Y_cont1","Y_bin1","Y_TTE1","Y_cont1"),method=method,
                                treatment="Treatment",censoring=c("event1",NA,NA,"event1",NA),strata="strata",
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                                threshold=c(0.5,1,NA,0.25,0.5),
                                n.bootstrap = n.bootstrap, trace=0)
  
  SBT_Mix <- summary(BT_Mix[[method]])
  
  test_that("count pairs summary - mixed",{
    valTest <- as.double(validPairs(BT_Mix[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - mixed",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_Mix[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}


##### export
version <- packageVersion("BuyseTest")
dir <- paste0("tests/Results-version",version)

if(!is.null(save) && save == TRUE){
  results_strata <- list(OutcomeBin = list(data = data_BinS, BT = BT_Bin1),
                         OutcomeCont = list(data = data_ContS, BT = BT_Cont1),
                         OutcomeTTE = list(data = data_TTES, BT = BT_TTE1), 
                         OutcomeMix = list(data = data_MixS, BT = BT_Mix))
  if(dir.exists(dir) == FALSE){dir.create(dir)}
  save(results_strata, file = file.path(dir,"test_strata.RData"))
}else if(!is.null(save) && save == FALSE){
  load(file = file.path(dir,"test_strata.RData"))
  
  test_that("comparison with the previous version", {
    expect_equalBT(BT_Bin1, results_strata$OutcomeBin$BT)
    expect_equalBT(BT_Cont1, results_strata$OutcomeCont$BT)
    expect_equalBT(BT_TTE1$Gehan, results_strata$OutcomeTTE$BT$Gehan)
    expect_equalBT(BT_TTE1$Peto, results_strata$OutcomeTTE$BT$Peto)
    expect_equalBT(BT_TTE1$Efron, results_strata$OutcomeTTE$BT$Efron)
    expect_equalBT(BT_TTE1$Peron, results_strata$OutcomeTTE$BT$Peron)
    
    expect_equalBT(BT_Mix$Gehan, results_strata$OutcomeMix$BT$Gehan)
    expect_equalBT(BT_Mix$Peto, results_strata$OutcomeMix$BT$Peto)
    expect_equalBT(BT_Mix$Efron, results_strata$OutcomeMix$BT$Efron)
    expect_equalBT(BT_Mix$Peron, results_strata$OutcomeMix$BT$Peron)
  })
}