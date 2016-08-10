#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test the consistency of the number of pairs over outcomes
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options(error=function() traceback(2)) 
options(max.print=10000)

#### package
library(BuyseTest)
library(testthat)
precision <- 10^{-7}

#### function
valid.summary <- function(BuyseSummary,precision=10^{-7}){
  index_global <- which(BuyseSummary$nb$strata=="global")
  diff <- BuyseSummary$nb[index_global[-1],"n.total"] - rowSums(BuyseSummary$nb[index_global[-length(index_global)],c("n.neutral","n.uninf")])
  return(diff)
}

#### 1- binary endpoint ####
n.Treatment_testBin <- 500
n.Control_testBin <- 500
prob.Treatment_testBin <- c(0.5,0.75)
prob.Control_testBin <- c(0.5,0.25)

set.seed(10)
data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)

data_testBin$endpoint3 <- 1

#### test count
BuyseTest_testBin1 <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("bin","bin"),
                                n.bootstrap=0,trace=0)

SBuyseTest_testBin1 <- summary(BuyseTest_testBin1)

if(valid.summary(SBuyseTest_testBin1)$test){
  stop("test_countPairs[Bin] incorrect summary \n",
       "pair management may be wrong \n")
}

#### test endpoint
BuyseTest_testBin2 <- BuyseTest(data=data_testBin,endpoint=c("endpoint3","endpoint1"),
                                treatment="treatment", type=c("bin","bin"),
                                n.bootstrap=0,trace=0)

SBuyseTest_testBin2 <- summary(BuyseTest_testBin2)

test.count <- SBuyseTest_testBin1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testBin2$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count>precision))){
  stop("test_countPairs[Bin] incorrect summary \n",
       "problem when using various outcomes \n")
}

#### 2- continous endpoint #####
n.Treatment_testCont <- 500
n.Control_testCont <- 500
mu.Treatment_testCont <- c(2,4)
mu.Control_testCont <- c(2,0)

set.seed(10)
data_testCont <- data.frame(treatment=c(rep(1,n.Treatment_testCont),rep(0,n.Treatment_testCont)))
data_testCont$endpoint1 <- c(rnorm(n.Treatment_testCont,mean=mu.Treatment_testCont[1]),
                             rnorm(n.Control_testCont,mean=mu.Control_testCont[1]))
data_testCont$endpoint2 <- c(rnorm(n.Control_testCont,mean=mu.Treatment_testCont[2]),
                             rnorm(n.Control_testCont,mean=mu.Control_testCont[2]))
data_testCont$strata <- rbinom(n.Treatment_testCont+n.Control_testCont,size=4,prob=0.5)

data_testCont$endpoint3 <- 1

#### test count
BuyseTest_testCont1 <- BuyseTest(data=data_testCont,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("cont","cont"),threshold=c(1,1),
                                n.bootstrap=0,trace=0)

SBuyseTest_testCont1 <- summary(BuyseTest_testCont1)

if(valid.summary(SBuyseTest_testCont1)$test){
  stop("test_countPairs[Cont] incorrect summary \n",
       "pair management may be wrong \n")
}

#### test endpoint
BuyseTest_testCont2 <- BuyseTest(data=data_testCont,endpoint=c("endpoint3","endpoint1"),
                                treatment="treatment", type=c("cont","cont"),threshold=c(1,1),
                                n.bootstrap=0,trace=0)

SBuyseTest_testCont2 <- summary(BuyseTest_testCont2)

test.count <- SBuyseTest_testCont1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testCont2$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count>precision))){
  stop("test_countPairs[Cont] incorrect summary \n",
       "problem when using various outcomes \n")
}


#### 3- survival outcome ####

#### data ####

n.Treatment_testTTE <- 500
n.Control_testTTE <- 500
lambda.Treatment_testTTE <- c(0.75,0.5)
lambda.Control_testTTE <- c(0.75,5)
lambda.Censoring_testTTE <- c(0.5,0.5)

set.seed(10)
data_testTTE <- data.frame(treatment=c(rep(1,n.Treatment_testTTE),rep(0,n.Treatment_testTTE)))
data_testTTE$EventTime1 <- c(rexp(n.Treatment_testTTE,rate=lambda.Treatment_testTTE[1]),
                             rexp(n.Control_testTTE,rate=lambda.Control_testTTE[1]))
data_testTTE$EventTime2 <- c(rexp(n.Control_testTTE,rate=lambda.Treatment_testTTE[2]),
                             rexp(n.Control_testTTE,rate=lambda.Control_testTTE[2]))
data_testTTE$CensoringTime1 <- c(rexp(n.Treatment_testTTE,rate=lambda.Censoring_testTTE[1]),
                                 rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[1]))
data_testTTE$CensoringTime2 <- c(rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[2]),
                                 rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[2]))
data_testTTE$strata <- rbinom(n.Treatment_testTTE+n.Control_testTTE,size=4,prob=0.5)

data_testTTE$endpoint1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,min)
data_testTTE$event1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,which.min)==1
data_testTTE$event1 <- as.numeric(data_testTTE$event1)
data_testTTE$endpoint2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,min)
data_testTTE$event2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,which.min)==1
data_testTTE$event2 <- as.numeric(data_testTTE$event2)

data_testTTE$endpoint3 <- 1
data_testTTE$event3 <- 1

#### Gehan ####

#### a) test count Gehan - same TTE endpoint

BuyseTest_testGehan1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Gehan",
                                   treatment="treatment",censoring=c("event1","event1","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testGehan1 <- summary(BuyseTest_testGehan1)


if(valid.summary(SBuyseTest_testGehan1)$test){
  stop("test_countPairs[TTE - Gehan] incorrect summary \n",
       "pair management may be wrong with repeted TTE endpoint \n")
}

#### b) test count Gehan - different TTE endpoint

BuyseTest_testGehan2 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1","endpoint2","endpoint1","endpoint2"),method="Gehan",
                                  treatment="treatment",censoring=c("event1","event2","event1","event2","event1","event2"),
                                  type=c("TTE","TTE","TTE","TTE","TTE","TTE"),threshold=c(0.75,0.75,0.5,0.5,0.25,0.25),trace=0)
SBuyseTest_testGehan2 <- summary(BuyseTest_testGehan2)


if(valid.summary(SBuyseTest_testGehan2)$test){
  stop("test_countPairs[TTE - Gehan] incorrect summary \n",
       "pair management may be wrong with different TTE endpoint \n")
}

#### c) test endpoint Gehan
BuyseTest_testGehan3 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint3","endpoint1","endpoint1"),method="Gehan",
                                  treatment="treatment",censoring=c("event3","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(1,0.75,0.5),trace=0)
SBuyseTest_testGehan3 <- summary(BuyseTest_testGehan3)


test.count1 <- SBuyseTest_testGehan1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testGehan3$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
test.count2 <- SBuyseTest_testGehan1$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testGehan3$nb[5,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count1)>precision) || any(abs(test.count2)>precision)){
  stop("test_countPairs[TTE - Gehan] incorrect summary \n",
       "problem when using various outcomes with first uninformative \n")
}

#### Peto ####
#### a) test count Peto - same TTE endpoint

BuyseTest_testPeto1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Peto",
                                  treatment="treatment",censoring=c("event1","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testPeto1 <- summary(BuyseTest_testPeto1)

if(valid.summary(SBuyseTest_testPeto1)$test){
  stop("test_countPairs[TTE - Peto] incorrect summary \n",
       "pair management may be wrong with repeted TTE endpoint \n")
}

#### b) test count Peto - different TTE endpoint

BuyseTest_testPeto2 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1","endpoint2","endpoint1","endpoint2"),method="Peto",
                                  treatment="treatment",censoring=c("event1","event2","event1","event2","event1","event2"),
                                  type=c("TTE","TTE","TTE","TTE","TTE","TTE"),threshold=c(0.75,0.75,0.5,0.5,0.25,0.25),trace=0)
SBuyseTest_testPeto2 <- summary(BuyseTest_testPeto2)

if(valid.summary(SBuyseTest_testPeto2)$test){
  stop("test_countPairs[TTE - Peto] incorrect summary \n",
       "pair management may be wrong with different TTE endpoint \n")
}

#### c) test endpoint Peto

BuyseTest_testPeto3 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint3","endpoint1","endpoint1"),method="Peto",
                                  treatment="treatment",censoring=c("event3","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(1,0.75,0.5),trace=0)
SBuyseTest_testPeto3 <- summary(BuyseTest_testPeto3)


test.count1 <- SBuyseTest_testPeto1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testPeto3$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
test.count2 <- SBuyseTest_testPeto1$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testPeto3$nb[5,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count1)>precision) || any(abs(test.count2)>precision)){
  stop("test_countPairs[TTE - Peto] incorrect summary \n",
       "problem when using various outcomes \n")
}

#### Efron ####
#### a) test count Efron - same TTE endpoint

BuyseTest_testEfron1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Efron",
                                 treatment="treatment",censoring=c("event1","event1","event1"),
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testEfron1 <- summary(BuyseTest_testEfron1)

if(valid.summary(SBuyseTest_testEfron1)$test){
  stop("test_countPairs[TTE - Efron] incorrect summary \n",
       "pair management may be wrong with repeted TTE endpoint \n")
}

#### b) test count Efron - different TTE endpoint

BuyseTest_testEfron2 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1","endpoint2","endpoint1","endpoint2"),method="Efron",
                                 treatment="treatment",censoring=c("event1","event2","event1","event2","event1","event2"),
                                 type=c("TTE","TTE","TTE","TTE","TTE","TTE"),threshold=c(0.75,0.75,0.5,0.5,0.25,0.25),trace=0)
SBuyseTest_testEfron2 <- summary(BuyseTest_testEfron2)

if(valid.summary(SBuyseTest_testEfron2)$test){
  stop("test_countPairs[TTE - Efron] incorrect summary \n",
       "pair management may be wrong with different TTE endpoint \n")
}

#### c) test endpoint Efron

BuyseTest_testEfron3 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint3","endpoint1","endpoint1"),method="Efron",
                                 treatment="treatment",censoring=c("event3","event1","event1"),
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.75,0.5),trace=0)
SBuyseTest_testEfron3 <- summary(BuyseTest_testEfron3)


test.count1 <- SBuyseTest_testEfron1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testEfron3$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
test.count2 <- SBuyseTest_testEfron1$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testEfron3$nb[5,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count1)>precision) || any(abs(test.count2)>precision)){
  stop("test_countPairs[TTE - Efron] incorrect summary \n",
       "problem when using various outcomes \n")
}

#### Peron ####
#### a) test count Peron - same TTE endpoint

BuyseTest_testPeron1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Peron",
                                  treatment="treatment",censoring=c("event1","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testPeron1 <- summary(BuyseTest_testPeron1)

if(valid.summary(SBuyseTest_testPeron1)$test){
  stop("test_countPairs[TTE - Peron] incorrect summary \n",
       "pair management may be wrong with repeted TTE endpoint \n")
}

#### b) test count Peron - different TTE endpoint

BuyseTest_testPeron2 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1","endpoint2","endpoint1","endpoint2"),method="Peron",
                                  treatment="treatment",censoring=c("event1","event2","event1","event2","event1","event2"),
                                  type=c("TTE","TTE","TTE","TTE","TTE","TTE"),threshold=c(0.75,0.75,0.5,0.5,0.25,0.25),trace=0)
SBuyseTest_testPeron2 <- summary(BuyseTest_testPeron2)

if(valid.summary(SBuyseTest_testPeron2)$test){
  stop("test_countPairs[TTE - Peron] incorrect summary \n",
       "pair management may be wrong with different TTE endpoint \n")
}

#### c) test endpoint Peron

BuyseTest_testPeron3 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint3","endpoint1","endpoint1"),method="Peron",
                                  treatment="treatment",censoring=c("event3","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(1,0.75,0.5),trace=0)
SBuyseTest_testPeron3 <- summary(BuyseTest_testPeron3)


test.count1 <- SBuyseTest_testPeron1$nb[1,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testPeron3$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
test.count2 <- SBuyseTest_testPeron1$nb[3,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]-SBuyseTest_testPeron3$nb[5,c("n.favorable","n.unfavorable","n.neutral","n.uninf")]
if(any(abs(test.count1)>precision) || any(abs(test.count2)>precision)){
  stop("test_countPairs[TTE - Peron] incorrect summary \n",
       "problem when using various outcomes \n")
}


#### 3- mixed outcome ####

data_testMixed <- data.frame(data_testTTE[,c("treatment","strata","event1")],
                             endpointBin1=data_testBin[,"endpoint1"],
                             endpointCont1=data_testCont[,"endpoint1"],
                             endpointTTE1=data_testTTE[,"endpoint1"])

#### without strata ####
for(method in c("Gehan","Peto","Efron","Peron")){
  BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                   endpoint=c("endpointTTE1","endpointCont1","endpointBin1","endpointTTE1","endpointCont1"),
                                   treatment="treatment",censoring=c("event1",NA,NA,"event1",NA),
                                   type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),method=method,
                                   threshold=c(0.5,1,NA,0.25,0.5),trace=0)
  
  SBuyseTest_testMixte <- summary(BuyseTest_testMixte)
  
  
  if(valid.summary(SBuyseTest_testMixte)$test){
    stop("test_countPairs[mixte1 - ",method,"] incorrect summary \n",
         "pair management may be wrong with mixed endpoints \n")
  }
  
  BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                   endpoint=c("endpointCont1","endpointTTE1","endpointBin1","endpointTTE1","endpointCont1"),
                                   treatment="treatment",censoring=c(NA,"event1",NA,"event1",NA),
                                   type=c("continuous","timeToEvent","binary","timeToEvent","continuous"),method=method,
                                   threshold=c(1,0.5,NA,0.25,0.5),trace=0)
  
  SBuyseTest_testMixte <- summary(BuyseTest_testMixte)
  
  
  if(valid.summary(SBuyseTest_testMixte)$test){
    stop("test_countPairs[mixte2 - ",method,"] incorrect summary \n",
         "pair management may be wrong with mixed endpoints \n")
  }
}

#### with strata ####
for(method in c("Gehan","Peto","Efron","Peron")){
  BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                   endpoint=c("endpointTTE1","endpointCont1","endpointBin1","endpointTTE1","endpointCont1"),
                                   treatment="treatment",censoring=c("event1",NA,NA,"event1",NA),strata="strata",
                                   type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),method=method,
                                   threshold=c(0.5,1,NA,0.25,0.5),trace=0)
  
  SBuyseTest_testMixte <- summary(BuyseTest_testMixte)
  
  
  if(valid.summary(SBuyseTest_testMixte)$test){
    stop("test_countPairs[mixte1Strata - ",method,"] incorrect summary \n",
         "pair management may be wrong with mixed endpoints \n")
  }

  BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                   endpoint=c("endpointCont1","endpointTTE1","endpointBin1","endpointTTE1","endpointCont1"),
                                   treatment="treatment",censoring=c(NA,"event1",NA,"event1",NA),strata="strata",
                                   type=c("continuous","timeToEvent","binary","timeToEvent","continuous"),method=method,
                                   threshold=c(1,0.5,NA,0.25,0.5),trace=0)
  
  SBuyseTest_testMixte <- summary(BuyseTest_testMixte)
  
  
  if(valid.summary(SBuyseTest_testMixte)$test){
    stop("test_countPairs[mixte2Strata - ",method,"] incorrect summary \n",
         "pair management may be wrong with mixed endpoints \n")
  }
}
