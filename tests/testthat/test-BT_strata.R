verboseContext("Check BuyseTest with strata")

#### additional spec
n.patients <- c(90,100)
n.bootstrap <- 10

#### 1- Simulated data ####
set.seed(10)
argsBin <- list(p.T = c(0.5,0.75))
argsCont <- list(mu.T = 1:3, sigma.T = rep(1,3))
argsTTE <- list(rates.T = 1:3, rates.Censor = rep(1,3))

#### binary #### 
cat("* binary endpoint \n")
data_Bin <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = NULL, argsTTE = NULL)
data_BinS <- rbind(cbind(data_Bin, strata = 1),
                   cbind(data_Bin, strata = 2),
                   cbind(data_Bin, strata = 3))
if(conv2df){data_BinS <- as.data.frame(data_BinS)}


BT_Bin1 <- BuyseTest(data=data_BinS,endpoint=c("toxicity1","toxicity2"),
                     treatment="Treatment", type=c("bin","bin"),strata="strata",
                     n.bootstrap=n.bootstrap)

test_that("count pairs summary - Binary",{
  valTest <- as.double(validPairs(BT_Bin1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Binary",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Bin1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### continuous ####
cat("* continuous endpoint \n")
set.seed(10)
data_Cont <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = argsCont, argsTTE = NULL)
data_ContS <- rbind(cbind(data_Cont, strata = 1),
                        cbind(data_Cont, strata = 2),
                        cbind(data_Cont, strata = 3))
if(conv2df){data_ContS <- as.data.frame(data_ContS)}

BT_Cont1 <- BuyseTest(data=data_ContS,endpoint=c("score1","score2"),
                                 treatment="Treatment", type=c("cont","cont"),threshold=c(1,1),strata="strata",
                                 n.bootstrap=n.bootstrap)

test_that("count pairs summary - Continuous",{
  valTest <- as.double(validPairs(BT_Cont1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Continuous",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Cont1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### TTE ####
cat("* TTE endpoint \n")
set.seed(10)
data_TTE <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = NULL, argsTTE = argsTTE)
data_TTES <- rbind(cbind(data_TTE, strata = 1),
                   cbind(data_TTE, strata = 2),
                   cbind(data_TTE, strata = 3))
if(conv2df){data_TTES <- as.data.frame(data_TTES)}

## different endpoints
BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE1[[method]] <- BuyseTest(data=data_TTES,endpoint=c("eventtime1","eventtime2","eventtime3"),method=method,
                                 treatment="Treatment",censoring=c("status1","status2","status3"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),
                                 n.bootstrap = n.bootstrap)
  
  test_that("count pairs summary - TTE",{
    valTest <- as.double(validPairs(BT_TTE1[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - TTE",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_TTE1[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

## same endpoint
BT_TTE2 <- vector(length = 4, mode = "list")
names(BT_TTE2) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE2[[method]] <- BuyseTest(data=data_TTES,endpoint=c("eventtime1","eventtime1","eventtime1"),method=method,
                                 treatment="Treatment",censoring=c("status1","status1","status1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.5,0.25),
                                 n.bootstrap = n.bootstrap)
  
  test_that("count pairs summary - TTE",{
    valTest <- as.double(validPairs(BT_TTE2[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - TTE",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_TTE2[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### mixed endpoints ####
cat("* mixed endpoints \n")
set.seed(10)
data_Mix <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = argsCont, argsTTE = argsTTE)
data_MixS <- rbind(cbind(data_Mix, strata = 1),
                   cbind(data_Mix, strata = 2),
                   cbind(data_Mix, strata = 3))
if(conv2df){data_MixS <- as.data.frame(data_MixS)}

BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_Mix[[method]] <- BuyseTest(data=data_MixS,
                                endpoint=c("eventtime1","score1","toxicity1","eventtime1","score1"),method=method,
                                treatment="Treatment",censoring=c("status1",NA,NA,"status1",NA),strata="strata",
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                                threshold=c(0.5,1,NA,0.25,0.5),
                                n.bootstrap = n.bootstrap)
  
  test_that("count pairs summary - mixed",{
    valTest <- as.double(validPairs(BT_Mix[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - mixed",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_Mix[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### 2- Real data ####
cat("* Veteran \n")
data(veteran, package = "survival")

BT_veteran <- BuyseTest(data = veteran, endpoint = "time", treatment = "trt", strata = "celltype",
                        type = "timeToEvent", censoring = "status",threshold = 0, 
                        n.bootstrap = n.bootstrap)
BT_veteran

##### export
if(identical(save, TRUE)){
  results_strata <- list(OutcomeBin = list(data = data_BinS, BT = BT_Bin1),
                         OutcomeCont = list(data = data_ContS, BT = BT_Cont1),
                         OutcomeTTE = list(data = data_TTES, BT1 = BT_TTE1, BT2 = BT_TTE2), 
                         OutcomeMix = list(data = data_MixS, BT = BT_Mix),
                         veteran = list(data = veteran, BT = BT_veteran))
  saveRDS(results_strata, file = file.path(dirSave,"test-strata.rds"))
}else if(identical(save, FALSE)){
  cat("* Previous version \n")
  GS <- readRDS(file = file.path(dirSave,"test-strata.rds"))
  
  test_that("comparison with the previous version", {
    expect_equalBT(BT_Bin1, GS$OutcomeBin$BT)
    expect_equalBT(BT_Cont1, GS$OutcomeCont$BT)
    expect_equalBT(BT_TTE1$Gehan, GS$OutcomeTTE$BT1$Gehan)
    expect_equalBT(BT_TTE1$Peto, GS$OutcomeTTE$BT1$Peto)
    expect_equalBT(BT_TTE1$Efron, GS$OutcomeTTE$BT1$Efron)
    expect_equalBT(BT_TTE1$Peron, GS$OutcomeTTE$BT1$Peron)
    expect_equalBT(BT_TTE2$Gehan, GS$OutcomeTTE$BT2$Gehan)
    expect_equalBT(BT_TTE2$Peto, GS$OutcomeTTE$BT2$Peto)
    expect_equalBT(BT_TTE2$Efron, GS$OutcomeTTE$BT2$Efron)
    expect_equalBT(BT_TTE2$Peron, GS$OutcomeTTE$BT2$Peron)
    
    expect_equalBT(BT_Mix$Gehan, GS$OutcomeMix$BT$Gehan)
    expect_equalBT(BT_Mix$Peto, GS$OutcomeMix$BT$Peto)
    expect_equalBT(BT_Mix$Efron, GS$OutcomeMix$BT$Efron)
    expect_equalBT(BT_Mix$Peron, GS$OutcomeMix$BT$Peron)
  
    expect_equalBT(BT_veteran, GS$veteran$BT)
  })
}
