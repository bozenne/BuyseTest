library(butils)
??butils

## fct
dirR  <- "/home/brice/Téléchargements/BuyseTest-master/R/"
dirC  <- "/home/brice/Téléchargements/BuyseTest-master/src/"
sapply(file.path(dirR,setdiff(list.files(dirR),"RcppExports.R")),source)
sapply(file.path(dirC,"FCT_buyseTest.cpp"),sourceCpp)

#### test_strata ####

## data
dir <- paste0("/home/brice/GitHub/BuyseTest/tests/Results-version1.0")
dir.exists(dir)
load(file = file.path(dir,"test_strata.RData"))

## bin
data_BinS <- as.data.frame(results_strata$OutcomeBin$data)

BT_Bin1 <- BuyseTest(data=data_BinS,endpoint=c("Y_bin1","Y_bin2"),
                     treatment="Treatment", type=c("bin","bin"),strata="strata",
                     n.bootstrap=0,trace=0)

## cont
data_ContS <- as.data.frame(results_strata$OutcomeCont$data)

BT_Cont1 <- BuyseTest(data=data_ContS,endpoint=c("Y_cont1","Y_cont2"),
                      treatment="Treatment", type=c("cont","cont"),threshold=c(1,1),strata="strata",
                      n.bootstrap=0,trace=0)


## TTE
data_TTES <- as.data.frame(results_strata$OutcomeTTE$data)

BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE1[[method]] <- BuyseTest(data=data_TTES,endpoint=c("Y_TTE1","Y_TTE2","Y_TTE3"),method=method,
                                 treatment="Treatment",censoring=c("event1","event2","event1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),trace=0)
}

## mix 

data_MixS <- as.data.frame(results_strata$OutcomeMix$data)

BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){ # method <- "Gehan"
  BT_Mix[[method]] <- BuyseTest(data=data_MixS,
                                endpoint=c("Y_TTE1","Y_cont1","Y_bin1","Y_TTE1","Y_cont1"),
                                treatment="Treatment",censoring=c("event1",NA,NA,"event1",NA),strata="strata",
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),method=method,
                                threshold=c(0.5,1,NA,0.25,0.5),trace=0)
}


results_strata <- list(OutcomeBin = list(data = data_BinS, BT = BT_Bin1),
                       OutcomeCont = list(data = data_ContS, BT = BT_Cont1),
                       OutcomeTTE = list(data = data_TTES, BT = BT_TTE1), 
                       OutcomeMix = list(data = data_MixS, BT = BT_Mix))
if(dir.exists(dir) == FALSE){dir.create(dir)}
save(results_strata, file = file.path(dir,"test_strata.RData"))


#### test_parallel ####

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
