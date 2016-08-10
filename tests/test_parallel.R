#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with parallel
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### spec
library(BuyseTest)
library(testthat)
library(lava)
precision <- 10^{-7}
n.patients <- 500
save <- FALSE
# butils:::package.source("BuyseTest", Rcode = TRUE, Ccode = TRUE)

#### data ####
set.seed(10)
dt.BT <- simulBT(n.patients)

#### model ####
n.bootstrap <- 10
seedGPC <- 15

#### pb avec Wscheme dans initData
ls.times <- list()
ls.sum <- list()
for(method in c("Gehan","Peto","Efron","Peron")){
  cat("method: ",method,"\n")
  
  ls.times[[method]] <- system.time(
    ls.sum[[method]] <- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                        type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=0,method=method,seed=seedGPC)
  )
  print(ls.times[[method]])
  summary(ls.sum[[method]])
}

BuyseresG
BuyseresP<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=100,trace=2,method="Efron",seed=seedGPC)
BuyseresP
BuyseresE<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=100,trace=2,method="Efron",seed=seedGPC)
BuyseresE
BuyseresE<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=100,trace=2,method="Efron",seed=seedGPC)
BuyseresE


BuyseresG<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Gehan",seed=seedGPC)

BuyseresG_parallel<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Gehan",seed=seedGPC,cpus=2)

BuyseresP<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peto",seed=seedGPC)

BuyseresP_parallel<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peto",seed=seedGPC,cpus=2)

BuyseresE<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Efron",seed=seedGPC)

BuyseresE_parallel<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Efron",seed=seedGPC,cpus=2)

BuyseresPer<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peron",seed=seedGPC)

BuyseresPer_parallel<- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                      type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=2,method="Peron",seed=seedGPC,cpus=2)

summary(BuyseresG)
summary(BuyseresG_parallel)

summary(BuyseresP)
summary(BuyseresP_parallel)

summary(BuyseresE)
summary(BuyseresE_parallel)

summary(BuyseresPer)
summary(BuyseresPer_parallel)

