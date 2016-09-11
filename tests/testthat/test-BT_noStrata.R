verboseContext("Check BuyseTest without strata")

#### additional spec
n.patients <- c(90,100)
n.bootstrap <- 1000

#### 1- Check number of pairs ####
set.seed(10)
argsBin <- list(p.T = c(0.5,0.75))
argsCont <- list(mu.T = 1:3, sigma.T = rep(1,3))
argsTTE <- list(rates.T = 1:3, rates.Censor = rep(1,3))

#### binary #### 
data_Bin <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = NULL, argsTTE = NULL)
if(conv2df){data_Bin <- as.data.frame(data_Bin)}

BT_Bin1 <- BuyseTest(data=data_Bin,endpoint=c("Y_bin1","Y_bin2"),
                     treatment="Treatment", type=c("bin","bin"),
                     n.bootstrap=n.bootstrap)

test_that("bootstrap approximately matches fisher test - Binary",{
  fisherP <- fisher.test(table(data_Bin$Y_bin1,data_Bin$Treatment))$p.value
  expect_equal(fisherP,BT_Bin1@p.value[1],1/sqrt(n.bootstrap)) # 1/sqrt(n.bootstrap) is an approximation of the convergence rate of the bootstrap - no theorical justification just a guess
})

test_that("count pairs summary - Binary",{
  valTest <- as.double(validPairs(BT_Bin1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### continuous ####
set.seed(10)
data_Cont <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = argsCont, argsTTE = NULL)
if(conv2df){data_Cont <- as.data.frame(data_Cont)}

BT_Cont1 <- BuyseTest(data=data_Cont,endpoint=c("Y_cont1","Y_cont2"),
                      treatment="Treatment", type=c("cont","cont"),threshold=c(0,1),
                      n.bootstrap=n.bootstrap)

test_that("bootstrap approximately matches wilcoxon test - Continuous",{
  wilcoxRes <- wilcox.test(data_Cont[Treatment==0,Y_cont1],
                           data_Cont[Treatment==1,Y_cont1], 
                           conf.int = TRUE)
  expect_equal(wilcoxRes$p.value,BT_Cont1@p.value[1],1/sqrt(n.bootstrap)) # 1/sqrt(n.bootstrap) is an approximation of the convergence rate of the bootstrap - no theorical justification just a guess
  # wilcoxRes$estimate - BT_Cont1@delta
})
test_that("count pairs summary - Continuous",{
  valTest <- as.double(validPairs(BT_Cont1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### TTE ####
set.seed(10)
data_TTE <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = NULL, argsTTE = argsTTE)
if(conv2df){data_TTE <- as.data.frame(data_TTE)}

## different endpoints
BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  
  n.bootstrapTempo <- if(method=="Gehan"){n.bootstrap}else{10}
  
  BT_TTE1[[method]] <- BuyseTest(data=data_TTE,endpoint=c("Y_TTE1","Y_TTE2","Y_TTE3"),method=method,
                                 treatment="Treatment",censoring=c("event1","event2","event1"),
                                 type=c("TTE","TTE","TTE"),threshold=c(0,0.5,0.25),
                                 n.bootstrap = n.bootstrapTempo)
  
  test_that("count pairs summary - TTE different endpoint",{
    valTest <- as.double(validPairs(BT_TTE1[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  
}

## same endpoint
BT_TTE2 <- vector(length = 4, mode = "list")
names(BT_TTE2) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  
  BT_TTE2[[method]] <- BuyseTest(data=data_TTE,endpoint=c("Y_TTE1","Y_TTE1","Y_TTE1"),method=method,
                                 treatment="Treatment",censoring=c("event1","event1","event1"),
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.5,0.25),
                                 n.bootstrap = 10)
  
  test_that("count pairs summary - TTE same endpoint",{
    valTest <- as.double(validPairs(BT_TTE2[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  
}

#### mixed outcomes ####
set.seed(10)
data_Mix <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = argsCont, argsTTE = argsTTE)
if(conv2df){data_Mix <- as.data.frame(data_Mix)}

BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_Mix[[method]] <- BuyseTest(data=data_Mix,
                                endpoint=c("Y_TTE1","Y_cont1","Y_bin1","Y_TTE1","Y_cont1"),method=method,
                                treatment="Treatment",censoring=c("event1",NA,NA,"event1",NA),
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                                threshold=c(0.5,1,NA,0.25,0.5),
                                n.bootstrap = 10)
  
  test_that("count pairs summary - mixed",{
    valTest <- as.double(validPairs(BT_Mix[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### Real data ####
data(veteran, package = "survival")

BT_veteran <- BuyseTest(data = veteran, endpoint = "time", treatment = "trt", 
                        type = "timeToEvent", censoring = "status",threshold = 0, 
                        n.bootstrap = 10)
BT_veteran

##### export
if(identical(save, TRUE)){
  results_noStrata <- list(OutcomeBin = list(data = data_BinS, BT = BT_Bin1),
                           OutcomeCont = list(data = data_ContS, BT = BT_Cont1),
                           OutcomeTTE = list(data = data_TTES, BT1 = BT_TTE1, BT2 = BT_TTE2), 
                           OutcomeMix = list(data = data_MixS, BT = BT_Mix),
                           veteran = list(data = veteran, BT = BT_veteran))
  saveRDS(results_noStrata, file = file.path(dirSave,"test-noStrata.rds"))
}else if(identical(save, FALSE)){
  GS <- readRDS(file = file.path(dirSave,"test-noStrata.rds"))
  
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
