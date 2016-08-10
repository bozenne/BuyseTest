#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest 
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options(error=function() traceback(2)) 
options(max.print=10000)

require(BuyseTest)

#### real example ####
data(veteran,package="survival")

BuyseTest_veteran <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
                     type="timeToEvent",censoring="status",threshold=0,n.bootstrap=110)

summary_BuyseTest_veteran <- summary(BuyseTest_veteran)

BuyseTest_veteranStrata <- BuyseTest(data=veteran,endpoint="time",treatment="trt",strata="celltype",
                           type="timeToEvent",censoring="status",threshold=0,n.bootstrap=10,trace=0)

summary_BuyseTest_veteranStrata <- summary(BuyseTest_veteranStrata)


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

#### no strata, n.bootstrap=0
BuyseTest_testBin <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                               treatment="treatment", type=c("bin","bin"))
summary_BuyseTest_testBin <- summary(BuyseTest_testBin)

#### no strata, n.bootstrap=10000
BuyseTest_testBin <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                                  treatment="treatment", type=c("bin","bin"),
                                  n.bootstrap=10,trace=0)

summary_BuyseTest_testBin <- summary(BuyseTest_testBin)

fisher.test(data_testBin[data_testBin$treatment==1,"endpoint1"],
	        data_testBin[data_testBin$treatment==1,"endpoint2"]
)

#### strata, n.bootstrap=0
BuyseTest_testBin <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                               treatment="treatment",strata="strata",
                               type=c("bin","bin"))
summary_BuyseTest_testBin <- summary(BuyseTest_testBin)

#### 2- continuous endpoint ####

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

#### no strata, n.bootstrap=0
BuyseTest_testCont <- BuyseTest(data=data_testCont,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("cont","cont"),threshold=c(1,1))
summary_BuyseTest_testCont <- summary(BuyseTest_testCont)

#### no strata, n.bootstrap=10000
BuyseTest_testCont <- BuyseTest(data=data_testCont,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment", type=c("cont","cont"),threshold=c(1,1),
                                n.bootstrap=10,trace=0)

summary_BuyseTest_testCont <- summary(BuyseTest_testCont)

wilcox.test(data_testCont[data_testCont$treatment==1,"endpoint1"],
            data_testCont[data_testCont$treatment==0,"endpoint1"],conf.int=TRUE,correct=FALSE)
wilcox.test(data_testCont[data_testCont$treatment==1,"endpoint2"],
            data_testCont[data_testCont$treatment==0,"endpoint2"],conf.int=TRUE,correct=FALSE)

#### strata, n.bootstrap=0
BuyseTest_testCont <- BuyseTest(data=data_testCont,endpoint=c("endpoint1","endpoint2"),
                                treatment="treatment",strata="strata",
                                type=c("cont","cont"),threshold=c(1,1))
summary_BuyseTest_testCont <- summary(BuyseTest_testCont)

#### 3- survival endpoint ####

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

#### no strata, n.bootstrap=0
BuyseTest_testTTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2"),
                               treatment="treatment",censoring=c("event1","event2"),
                               type=c("TTE","TTE"),threshold=c(0.25,0.25))
summary_BuyseTest_testTTE <- summary(BuyseTest_testTTE)
survdiff(Surv(endpoint1,event1)~treatment,data=data_testTTE)

#### no strata, n.bootstrap=10000
BuyseTest_testTTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2"),
                                  treatment="treatment",censoring=c("event1","event2"),
                                  type=c("TTE","TTE"),threshold=c(0.25,0.25),n.bootstrap=10,trace=0)


summary_BuyseTest_testTTE <- summary(BuyseTest_testTTE)

#### strata, n.bootstrap=0
BuyseTest_testTTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2"),
                            treatment="treatment",censoring=c("event1","event2"),strata="strata",
                            type=c("TTE","TTE"),threshold=c(0.25,0.25))
summary_BuyseTest_testTTE <- summary(BuyseTest_testTTE)


#### 4- mixed endpoint ####
data_testMixed <- data.frame(data_testTTE[,c("treatment","strata","event1")],
                             endpointBin1=data_testBin[,"endpoint1"],
                             endpointCont1=data_testCont[,"endpoint1"],
                             endpointTTE1=data_testTTE[,"endpoint1"])

BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointBin1","endpointCont1","endpointTTE1"),
                                 treatment="treatment",censoring=c(NA,NA,"event1"),strata="strata",
                                 type=c("binary","continuous","timeToEvent"),
                                 threshold=c(NA,1,0.5))


summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)



#### 5- KM imputation ####
n.Treatment_testMixte4 <- 500
n.Control_testMixte4 <- 500
prob.Treatment_testBin4 <- c(0.5,0.25,0.10,0.075,0.05,0.025)
prob.Control_testBin4 <- c(0.7,0.15,0.05,0.05,0.025,0.025)

lambda.Treatment_testTTE4 <- 0.75
lambda.Control_testTTE4 <- 1


set.seed(10)
data_testMixte <- data.frame(treatment=c(rep(1,n.Treatment_testMixte4),
                                         rep(0,n.Control_testMixte4)  ))
data_testMixte$toxicity <- c(apply(rmultinom(n.Treatment_testMixte4,size=1,
                             prob=prob.Treatment_testBin4)==1,2,which),
                             apply(rmultinom(n.Control_testMixte4,size=1,
							 prob=prob.Control_testBin4)==1,2,which))

data_testMixte$toxicity5 <- as.numeric(data_testMixte$toxicity<5)
data_testMixte$toxicity4 <- as.numeric(data_testMixte$toxicity<4)
data_testMixte$toxicity3 <- as.numeric(data_testMixte$toxicity<3)

data_testMixte$EventTime1 <- c(rexp(n.Treatment_testMixte4,rate=lambda.Treatment_testTTE4),
                               rexp(n.Control_testMixte4,rate=lambda.Control_testTTE4))
data_testMixte$CensoringTime1 <- c(rexp(n.Treatment_testMixte4,rate=lambda.Treatment_testTTE4),
                                  rexp(n.Control_testMixte4,rate=lambda.Control_testTTE4))
data_testMixte$CensoringTime1[data_testMixte$CensoringTime1>4] <- 4

data_testMixte$Survival <- apply(data_testMixte[,c("EventTime1","CensoringTime1")],1,min)
data_testMixte$event1 <- as.numeric(apply(data_testMixte[,c("EventTime1","CensoringTime1")],
                                          1,which.min)==1)

# boxplot(data_testMixte$EventTime1 ~ data_testMixte$treatment)
# boxplot(data_testMixte$Survival ~ data_testMixte$treatment)
resKM_tempo <- survfit(Surv(data_testMixte[,"Survival"],data_testMixte[,"event1"])~1)
plot(resKM_tempo)

BuyseTest4.1_testMixte <- BuyseTest(data=data_testMixte,method="Gehan",
               endpoint=c("Survival","toxicity5","Survival","toxicity4","Survival","toxicity3"),
               treatment="treatment",
               censoring=c("event1",NA,"event1",NA,"event1",NA),
               type=c("TTE","bin","TTE","bin","TTE","bin"),
               threshold=c(1.5,NA,0.75,NA,0.25,NA),n.bootstrap=10,trace=0)

summary_BuyseTest1.2_testTTE1 <- summary(BuyseTest4.1_testMixte)

BuyseTest4.1_testMixteKM <- BuyseTest(data=data_testMixte,method="Peto",
               endpoint=c("Survival","toxicity5","Survival","toxicity4","Survival","toxicity3"),
               treatment="treatment",
               censoring=c("event1",NA,"event1",NA,"event1",NA),
               type=c("TTE","bin","TTE","bin","TTE","bin"),
               threshold=c(1.5,NA,0.75,NA,0.25,NA),n.bootstrap=10,trace=2)
			   
summary_BuyseTest1.2_testTTE1KM <- summary(BuyseTest4.1_testMixteKM)
summary_BuyseTest1.2_testTTE1KM <- summary(BuyseTest4.1_testMixteKM,"nb")


#### 6- parallel computing ####
BuyseTest_testBin <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                                  treatment="treatment", type=c("bin","bin"),
                                  n.bootstrap=10000,cpus="all")

summary_BuyseTest_testBin <- summary(BuyseTest_testBin)

