if(FALSE){
    library(testthat)
    library(BuyseTest)
}

context("Check BuyseTest without strata")

## * settings
n.patients <- c(90,100)
format <- "data.table"
BuyseTest.options(n.permutation = 0, trace = 0, keepComparison = TRUE)

## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))
## butils::object2script(dt.sim)

## * binary endpoint 

## ** run BuyseTest

test_that("BuyseTest - binary ", {
    BT.bin <- BuyseTest(Treatment ~ bin(toxicity1),
                        data = dt.sim)
    
    BT2 <- BuyseTest(data = dt.sim,
                     endpoint = "toxicity1",
                     treatment = "Treatment",
                     type = "bin")
    
    ## test against fixed value    
    expect_equal(as.double(BT.bin@count_favorable),1739)
    expect_equal(as.double(BT.bin@count_unfavorable),2809)
    expect_equal(as.double(BT.bin@count_neutral),4452)
    expect_equal(as.double(BT.bin@count_uninf),0)
    expect_equal(as.double(BT.bin@Delta$netChance),-0.1188889,tol=1e-6)
    expect_equal(as.double(BT.bin@Delta$winRatio),0.619081,tol=1e-6)

    expect_equal(BT.bin,BT2)

    ## fisherP <- fisher.test(table(dt.sim$toxicity1,dt.sim$Treatment))

    ## count pairs
    tableS <- summary(BT.bin, show = FALSE, percentage = FALSE)
    expect_equal(tableS[,c("n.total")],
                 unname(rowSums(tableS[,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])))
})

## ** continuous endpoint
## test_that("BuyseTest - binary ", {
##     BT.cont <- BuyseTest(Treatment ~ bin(score1),
##                          data = dt.sim)
    
##     BT2 <- BuyseTest(data = dt.sim,
##                      endpoint = "toxicity1",
##                      treatment = "Treatment",
##                      type = "bin")
    
##     ## test against fixed value    
##     expect_equal(as.double(BT.bin@count_favorable),1739)
##     expect_equal(as.double(BT.bin@count_unfavorable),2809)
##     expect_equal(as.double(BT.bin@count_neutral),4452)
##     expect_equal(as.double(BT.bin@count_uninf),0)
##     expect_equal(as.double(BT.bin@Delta$netChance),-0.1188889,tol=1e-6)
##     expect_equal(as.double(BT.bin@Delta$winRatio),0.619081,tol=1e-6)

##     expect_equal(BT.bin,BT2)

##     ## fisherP <- fisher.test(table(dt.sim$toxicity1,dt.sim$Treatment))

##     ## count pairs
##     tableS <- summary(BT.bin, show = FALSE, percentage = FALSE)
##     expect_equal(tableS[,c("n.total")],
##                  unname(rowSums(tableS[,c("n.favorable","n.unfavorable","n.neutral","n.uninf")])))
## })


BT_Cont1 <- BuyseTest(data=data_Cont,endpoint=c("score1","score2"),
                      treatment="Treatment", type=c("cont","cont"),threshold=c(0,1),
                      n.bootstrap=n.bootstrap)

test_that("bootstrap approximately matches wilcoxon test - Continuous",{
  wilcoxRes <- wilcox.test(data_Cont[Treatment==0,score1],
                           data_Cont[Treatment==1,score1], 
                           conf.int = TRUE)
  expect_equal(wilcoxRes$p.value,BT_Cont1@p.value$netChance[1],1/sqrt(n.bootstrap)) # 1/sqrt(n.bootstrap) is an approximation of the convergence rate of the bootstrap - no theorical justification just a guess
  # wilcoxRes$estimate - BT_Cont1@delta
})
test_that("count pairs summary - Continuous",{
  valTest <- as.double(validPairs(BT_Cont1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

## ** Time to event endpoint
cat("* TTE endpoint \n")
set.seed(10)
data_TTE <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = NULL, argsTTE = argsTTE)
if(conv2df){data_TTE <- as.data.frame(data_TTE)}

## *** different endpoints
BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){ # method <- "Gehan"
  
  n.bootstrapTempo <- if(method=="Gehan"){n.bootstrap}else{10}

  BT_TTE1[[method]] <- BuyseTest(data=data_TTE,endpoint=c("eventtime1","eventtime2","eventtime3"),method=method,
                                 treatment="Treatment",censoring=c("status1","status2","status3"),
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.5,0.25),
                                 n.bootstrap = n.bootstrapTempo)
  
  test_that("count pairs summary - TTE different endpoint",{
    valTest <- as.double(validPairs(BT_TTE1[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  
}

## *** same endpoint
BT_TTE2 <- vector(length = 4, mode = "list")
names(BT_TTE2) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  
  BT_TTE2[[method]] <- BuyseTest(data=data_TTE,endpoint=c("eventtime1","eventtime1","eventtime1"),method=method,
                                 treatment="Treatment",censoring=c("status1","status1","status1"),
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.5,0.25),
                                 n.bootstrap = 10)
  
  test_that("count pairs summary - TTE same endpoint",{
    valTest <- as.double(validPairs(BT_TTE2[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  
}

## ** mixed endpoints 
cat("* mixed endpoints \n")
set.seed(10)
data_Mix <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = argsCont, argsTTE = argsTTE)
if(conv2df){data_Mix <- as.data.frame(data_Mix)}

BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_Mix[[method]] <- BuyseTest(data=data_Mix,
                                endpoint=c("eventtime1","score1","toxicity1","eventtime1","score1"),method=method,
                                treatment="Treatment",censoring=c("status1",NA,NA,"status1",NA),
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                                threshold=c(0.5,1,NA,0.25,0.5),
                                n.bootstrap = 10)
  
  test_that("count pairs summary - mixed",{
    valTest <- as.double(validPairs(BT_Mix[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

## * Real data
cat("* Veteran \n")
data(veteran, package = "survival")

BT_veteran <- BuyseTest(data = veteran, endpoint = "time", treatment = "trt", 
                        type = "timeToEvent", censoring = "status",threshold = 0, 
                        n.bootstrap = 10)
BT_veteran




## ## ** with many pairs (lambda.C = 0.5)
## set.seed(10)
## TpsFin <- 1 # 0.75
## lambda.T <- 0.5
## n.Treatment <- 10
## n.Control <- 10
## n <- n.Treatment + n.Control
## group <- c(rep(1, n.Treatment),rep(0, n.Control))

## lambda.C <- 0.5
## TimeEvent <- c(rexp(n.Treatment,rate=lambda.T),
##              rexp(n.Control,rate=lambda.C))
## Time.Cens <- runif(n,0,TpsFin)
## Time <-pmin(Time.Cens,TimeEvent)
## Event <- Time == TimeEvent
## Event <- as.numeric(Event)
## tab <- data.frame(group,Time,Event)

## ## *** Gehan
## test_that("lambdaC = 0.5 - Gehan",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Gehan",seed=11)
  
##   expect_equal(as.double(BT@count_favorable),6)
##   expect_equal(as.double(BT@count_unfavorable),8)
##   expect_equal(as.double(BT@count_neutral),0)
##   expect_equal(as.double(BT@count_uninf),86)
  
##   expect_equal(as.double(BT@delta$netChance),-0.02)
##   expect_equal(as.double(BT@delta$winRatio),15/20)
## })

## ## *** Peto
## test_that("lambdaC = 0.5 - Peto",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peto",seed=11)
  
##   expect_equal(as.double(BT@count_favorable),40.95, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable),41.67, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral),0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf),17.38, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),-0.00712, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),0.9829167, tolerance = 1e-3)
## })

## ## *** Efron
## test_that("lambdaC = 0.5 - Efron",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Efron",seed=11)
  
##   expect_equal(as.double(BT@count_favorable),11.11, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable),11.11, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral),1, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf),76.78, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),0, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),1, tolerance = 1e-3)
  
## })

## ## *** Peron
## test_that("lambdaC = 0.5 - Peron",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peron",seed=11)
 
##   expect_equal(as.double(BT@count_favorable), 11.11, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable), 11.11, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral), 0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf), 77.78, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),0, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),1, tolerance = 1e-3)
## })

## ## ** with many pairs (lambda.C = 1)
## lambda.C <- 1
## TimeEvent<-c(rexp(n.Treatment,rate=lambda.T),
##              rexp(n.Control,rate=lambda.C))
## Time.Cens<-runif(n,0,TpsFin)
## Time<-pmin(Time.Cens,TimeEvent)
## Event<-Time==TimeEvent
## Event<-as.numeric(Event)
## tab<-data.frame(group,Time,Event)

## ## *** Gehan
## test_that("lambdaC = 1 - Gehan",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Gehan",seed=11)
  
##   expect_equal(as.double(BT@count_favorable),0)
##   expect_equal(as.double(BT@count_unfavorable),9)
##   expect_equal(as.double(BT@count_neutral),0)
##   expect_equal(as.double(BT@count_uninf),91)
  
##   expect_equal(as.double(BT@delta$netChance),-0.09)
##   expect_equal(as.double(BT@delta$winRatio),0)
## })

## ## *** Peto
## test_that("lambdaC = 1 - Peto",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peto",seed=11)
  
##   expect_equal(as.double(BT@count_favorable),36, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable),42, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral),0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf),22, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),-0.06, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),0.8571429, tolerance = 1e-3)
## })

## ## *** Efron
## test_that("lambdaC = 1 - Efron",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Efron",seed=11)
##   expect_equal(as.double(BT@count_favorable),0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable),55, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral),1, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf),44, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),-0.55, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),0, tolerance = 1e-3)
  
## })

## ## *** Peron
## test_that("lambdaC = 1 - Peron",{
##   BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
##                   type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peron",seed=11)
  
##   expect_equal(as.double(BT@count_favorable), 0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_unfavorable), 10, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_neutral), 0, tolerance = 1e-3)
##   expect_equal(as.double(BT@count_uninf), 90, tolerance = 1e-3)
  
##   expect_equal(as.double(BT@delta$netChance),-0.1, tolerance = 1e-3)
##   expect_equal(as.double(BT@delta$winRatio),0, tolerance = 1e-3)
## })
