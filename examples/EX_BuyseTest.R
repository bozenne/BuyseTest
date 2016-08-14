#### real example : Veteran dataset of the survival package ####
#### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#### method = "Gehan"

\dontrun{
  data(veteran,package="survival")
  library(BuyseTest)
  BuyseTest_veteran_Gehan <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
                                       strata="celltype",type="timeToEvent",censoring="status",threshold=0,
                                       n.bootstrap=10000,method="Gehan",cpus="all")
  
  summary_veteran_Gehan <- summary(BuyseTest_veteran_Gehan)
  
  #### method = "Peron"
  
  BuyseTest_veteran_Peron <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
                                       strata="celltype",type="timeToEvent",censoring="status",threshold=0,
                                       n.bootstrap=1000,method="Peron",cpus="all")
  
  summary_veteran_Peron <- summary(BuyseTest_veteran_Peron)
}

#### Several endpoints :
#######Survival, a time-to-event endpoint
#######Toxicity, a continuous/ordinal endpoint : 6 grades of maximal adverse event 

set.seed(10)

n.Treatment <- 100
n.Control <- 100
prob.Treatment_TOX <- c(0.5,0.25,0.10,0.075,0.05,0.025)
prob.Control_TOX <- c(0.7,0.15,0.05,0.05,0.025,0.025)

lambda.Treatment_TTE <- 0.6
lambda.Control_TTE <- 1

data_test <- data.frame(treatment=c(rep(1,n.Treatment),
                                    rep(0,n.Control)  ))
data_test$toxicity <- c(apply(rmultinom(n.Treatment,size=1,
                                        prob=prob.Treatment_TOX)==1,2,which),
                        apply(rmultinom(n.Control,size=1,
                                        prob=prob.Control_TOX)==1,2,which))

data_test$toxicityInv <-6-data_test$toxicity

data_test$EventTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
                         rexp(n.Control,rate=lambda.Control_TTE))
data_test$CensoringTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
                             rexp(n.Control,rate=lambda.Control_TTE))
data_test$CensoringTime[data_test$CensoringTime>4] <- 4

data_test$Survival <- apply(data_test[,c("EventTime","CensoringTime")],1,min)
data_test$event <- as.numeric(apply(data_test[,c("EventTime","CensoringTime")],
                                    1,which.min)==1)

resKM_tempo <- survfit(Surv(data_test[,"Survival"],data_test[,"event"])~data_test$treatment)
plot(resKM_tempo)

#### method = "Gehan". 

BuyseTest_severalendpoint_Gehan <- BuyseTest(data=data_test,method="Gehan",
                                             endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
                                             treatment="treatment",
                                             censoring=c("event",NA,"event",NA,"event",NA),
                                             type=c("TTE","cont","TTE","cont","TTE","cont"),
                                             threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
summary(BuyseTest_severalendpoint_Gehan)

#### method = "Peron". 

BuyseTest_severalendpoint_Peron <- BuyseTest(data=data_test,method="Peron",
                                             endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
                                             treatment="treatment",
                                             censoring=c("event",NA,"event",NA,"event",NA),
                                             type=c("TTE","cont","TTE","cont","TTE","cont"),
                                             threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
summary(BuyseTest_severalendpoint_Peron)
