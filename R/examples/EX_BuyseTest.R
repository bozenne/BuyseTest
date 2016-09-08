#### real example : Veteran dataset of the survival package ####
#### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#### method = "Gehan"

\dontrun{
  data(veteran,package="survival")
  library(BuyseTest)
  BT_Gehan <- BuyseTest(data=veteran,endpoint="time",treatment="trt",strata="celltype",
                        type="timeToEvent",censoring="status",threshold=0,
                        n.bootstrap=10000,method="Gehan",cpus="all")
  
  summary_Gehan <- summary(BT_Gehan)
  
  #### method = "Peron"
  
  BT_Peron <- BuyseTest(data=veteran,endpoint="time",treatment="trt",strata="celltype",
                        type="timeToEvent",censoring="status",threshold=0,
                        n.bootstrap=1000,method="Peron",cpus="all")
  
  summary_Peron <- summary(BT_Peron)
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

all.endpoints <- c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv")

#### method = "Gehan". 
\dontrun{
BT_Gehan <- BuyseTest(data=data_test,method="Gehan",
                      endpoint=all.endpoints,
                      treatment="treatment",
                      censoring=c("event",NA,"event",NA,"event",NA),
                      type=c("TTE","cont","TTE","cont","TTE","cont"),
                      threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=100,trace=2,cpus="all")
}
\dontshow{
  BT_Gehan <- BuyseTest(data=data_test,method="Gehan",
                        endpoint=all.endpoints,
                        treatment="treatment",
                        censoring=c("event",NA,"event",NA,"event",NA),
                        type=c("TTE","cont","TTE","cont","TTE","cont"),
                        threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1,trace=0,cpus=1)
  
}
summary(BT_Gehan)

#### method = "Peron". 

\dontrun{
  BT_Peron <- BuyseTest(data=data_test,method="Peron",
                      endpoint=all.endpoints,
                      treatment="treatment",
                      censoring=c("event",NA,"event",NA,"event",NA),
                      type=c("TTE","cont","TTE","cont","TTE","cont"),
                      threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=100,trace=2,cpus="all")
}
\dontshow{
  BT_Peron <- BuyseTest(data=data_test,method="Peron",
                        endpoint=all.endpoints,
                        treatment="treatment",
                        censoring=c("event",NA,"event",NA,"event",NA),
                        type=c("TTE","cont","TTE","cont","TTE","cont"),
                        threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1,trace=0,cpus=1)
}
summary(BT_Peron)
