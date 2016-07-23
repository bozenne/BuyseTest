#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with strata
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
precision <- 10^{-7}

#### function
validIdentical.strata <- function(BuyseSummary,precision=10^{-7}){
 
  endpoints <- unique(BuyseSummary$nb$endpoint_threshold)
  D <- length(endpoints)
  index_strata <- which(BuyseSummary$nb$strata!="global")
  levels.strata <- unique(BuyseSummary$nb$strata[index_strata])
  n.strata <- length(levels.strata) 
  
  diff <- matrix(NA,nrow=n.strata*D,ncol=5)
  for(iter_endpoint in 1:D){
   index_endpoint <- which(BuyseSummary$nb$endpoint_threshold==endpoints[iter_endpoint])  
   res_strata <- as.matrix(BuyseSummary$nb[intersect(index_strata,index_endpoint),c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")])
   diff[(iter_endpoint-1)*n.strata+1:n.strata,] <- t(apply(res_strata,1,function(x){x-res_strata[1,,drop=TRUE]}))
  }
  
  return(list(diff=diff,
              test=any(abs(diff)>precision)
              )
  )
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
data_testBin <- rbind(cbind(data_testBin,strata=0),
                      cbind(data_testBin,strata=1),
                      cbind(data_testBin,strata=2))

#### test count
BuyseTest_testBin1 <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("bin","bin"),strata="strata",
                                n.bootstrap=0,trace=0)

SBuyseTest_testBin1 <- summary(BuyseTest_testBin1)

if(validIdentical.strata(SBuyseTest_testBin1)$test){
  stop("test_strata[Bin] incorrect summary \n",
       "strata management may be wrong \n")
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
data_testCont <- rbind(cbind(data_testCont,strata=0),
                       cbind(data_testCont,strata=1),
                       cbind(data_testCont,strata=2))

#### test count
BuyseTest_testCont1 <- BuyseTest(data=data_testCont,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("cont","cont"),threshold=c(1,1),strata="strata",
                                n.bootstrap=0,trace=0)

SBuyseTest_testCont1 <- summary(BuyseTest_testCont1)

if(validIdentical.strata(SBuyseTest_testCont1)$test){
  stop("test_strata[Cont] incorrect summary \n",
       "strata management may be wrong \n")
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
data_testTTE$endpoint1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,min)
data_testTTE$event1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,which.min)==1
data_testTTE$event1 <- as.numeric(data_testTTE$event1)
data_testTTE$endpoint2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,min)
data_testTTE$event2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,which.min)==1
data_testTTE$event2 <- as.numeric(data_testTTE$event2)


data_testTTE <- rbind(cbind(data_testTTE,strata=0),
                      cbind(data_testTTE,strata=1),
                      cbind(data_testTTE,strata=2))

#### Gehan ####

#### a) test count Gehan - different TTE endpoint 

BuyseTest_testGehan1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1"),method="Gehan",
                                   treatment="treatment",censoring=c("event1","event2","event1"),strata="strata",
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testGehan1 <- summary(BuyseTest_testGehan1)


if(validIdentical.strata(SBuyseTest_testGehan1)$test){
  stop("test_strata[TTE - Gehan] incorrect summary \n",
       "strata management may be wrong \n")
}


#### b) test count Peto - different TTE endpoint

BuyseTest_testPeto1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1"),method="Peto",
                                  treatment="treatment",censoring=c("event1","event2","event1"),strata="strata",
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testPeto1 <- summary(BuyseTest_testPeto1)


if(validIdentical.strata(SBuyseTest_testPeto1)$test){
  stop("test_strata[TTE - Peto] incorrect summary \n",
       "strata management may be wrong \n")
}

#### c) test endpoint Efron - different TTE endpoint
BuyseTest_testEfron1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1"),method="Efron",
                                 treatment="treatment",censoring=c("event1","event2","event1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testEfron1 <- summary(BuyseTest_testEfron1)


if(validIdentical.strata(SBuyseTest_testEfron1)$test){
  stop("test_strata[TTE - Efron] incorrect summary \n",
       "strata management may be wrong \n")
}

#### c) test endpoint Peron - different TTE endpoint
BuyseTest_testPeron1 <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1"),method="Peron",
                                  treatment="treatment",censoring=c("event1","event2","event1"),strata="strata",
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
SBuyseTest_testPeron1 <- summary(BuyseTest_testPeron1)


if(validIdentical.strata(SBuyseTest_testPeron1)$test){
  stop("test_strata[TTE - Peron] incorrect summary \n",
       "strata management may be wrong \n")
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
                                   treatment="treatment",censoring=c("event1",NA,NA,"event1",NA),strata="strata",
                                   type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),method=method,
                                   threshold=c(0.5,1,NA,0.25,0.5),trace=0)
  
  SBuyseTest_testMixte <- summary(BuyseTest_testMixte)
  
  
  if(validIdentical.strata(SBuyseTest_testMixte)$test){
    stop("test_strata[mixte - ",method,"] incorrect summary \n",
         "strata management may be wrong \n")
  }
  
}
