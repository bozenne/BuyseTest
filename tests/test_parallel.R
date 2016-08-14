#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with parallel
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### spec
library(BuyseTest) # butils:::package.source("BuyseTest", Rcode = TRUE, Ccode = TRUE)
library(testthat)
library(data.table)
library(lava)
precision <- 10^{-7}
n.patients <- 100
n.bootstrap <- 100


#### data ####
set.seed(10)
dt.BT <- simulBT(n.patients)

#### sequential computation ####
ls.times <- list()
ls.sum <- list()
for(method in c("Gehan","Peto","Efron","Peron")){
  cat("method: ",method,"\n")
  
  ls.times[[method]] <- system.time(
    ls.sum[[method]] <- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                        type="TTE",censoring="event1",threshold=0,n.bootstrap=n.bootstrap,trace=0,method=method)
  )
}

#### parallel computation ####
ls.times <- list()
ls.sum <- list()
for(method in c("Gehan","Peto","Efron","Peron")){ # method <- "Gehan"
  cat("method: ",method,"\n")
  
  ls.times[[method]] <- system.time(
    ls.sum[[method]] <- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                                  type="TTE",censoring="event1",threshold=0, cpus = "all",
                                  n.bootstrap=n.bootstrap,trace=1,method=method)
  )
}
