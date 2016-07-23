#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with parallel
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

require(BuyseTest)


#### data ####
res = NULL #start with an empty dataset

n.Treatment <- 5
n.Control <- 5
n<-n.Treatment+n.Control
group<-c(rep(1, n.Treatment),rep(0, n.Control))
lambda.T <- 0.5
lambda.C <- 1
TpsFin <-1.5936
TimeEvent<-c(rexp(n.Treatment,rate=lambda.T),
             rexp(n.Control,rate=lambda.C))
Time.Cens<-runif(n,0,TpsFin)
Time<-pmin(Time.Cens,TimeEvent)
Event<-Time==TimeEvent
Event<-as.numeric(Event)
tab<-data.frame(group,Time,Event)
#plot(survfit(Surv(time=Time, event=Event) ~ group, data=tab))
seedGPC<-15

#### model ####
n.bootstrap <- 1000

BuyseresG<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Gehan",seed=seedGPC)

BuyseresG_parallel<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Gehan",seed=seedGPC,cpus=2)

BuyseresP<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peto",seed=seedGPC)

BuyseresP_parallel<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peto",seed=seedGPC,cpus=2)

BuyseresE<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Efron",seed=seedGPC)

BuyseresE_parallel<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Efron",seed=seedGPC,cpus=2)

BuyseresPer<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peron",seed=seedGPC)

BuyseresPer_parallel<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peron",seed=seedGPC,cpus=2)

summary(BuyseresG)
summary(BuyseresG_parallel)

summary(BuyseresP)
summary(BuyseresP_parallel)

summary(BuyseresE)
summary(BuyseresE_parallel)

summary(BuyseresPer)
summary(BuyseresPer_parallel)

