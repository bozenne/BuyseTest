#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with parallel
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#context("Parallel")

library(testthat)
library(BuyseTest)

library(lava)
library(data.table)
library(survival)

#### additional spec
n.patients <- 100
n.bootstrap <- 100
precision <- 10^{-7}
save <- NULL # TRUE to save results, FALSE to test, NULL to ignore
conv2df <- FALSE

#### data ####
set.seed(10)
dt.BT <- simulBT(n.patients)

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


##### export
version <- packageVersion("BuyseTest")
dir <- paste0("tests/Results-version",version)

if(identical(save, TRUE)){
  results <- list(data = dt.BT,
                  Gehan = ls.sum$Gehan,
                  Peto = ls.sum$Peto,
                  Efron = ls.sum$Efron,
                  Peron = ls.sum$Peron)
  if(dir.exists(dir) == FALSE){dir.create(dir)}
  saveRDS(results, file = file.path(dir,"test-parallel.rds"))
}else if(identical(save, FALSE)){
   GS <- readRDS(file = file.path(dir,"test-parallel.rds"))
  
  test_that("comparison with the previous version", {
    expect_equalBT(ls.sum$Gehan, GS$Gehan)
    expect_equalBT(ls.sum$Peto, GS$Peto)
    expect_equalBT(ls.sum$Efron, GS$Efron)
    expect_equalBT(ls.sum$Peron, GS$Peron)
  })
}
