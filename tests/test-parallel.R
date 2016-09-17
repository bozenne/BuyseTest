library(BuyseTest)
library(testthat)
BuyseTest.options(trace = 0, keep.bootstrap = TRUE)

n.patients <- 100
n.bootstrap <- 100
save <- FALSE # TRUE to save results, FALSE to test, NULL to ignore
if(identical(save, TRUE)){
  dirSave <- paste0("Results-version",packageVersion("BuyseTest"))
  if(dir.exists(dirSave) == FALSE){dir.create(dirSave)}
}else{
  dirSave <- paste0("Results-version","1.0")
}
source(file.path("FCT","FCT_check.R"))

#### data ####
verboseContext("Check parallel boostrap")
set.seed(10)
dt.BT <- simulBT(n.patients)

#### parallel computation ####
ls.sum <- list()
for(method in c("Gehan","Peto","Efron","Peron")){ # method <- "Gehan"
  cat(method, " ")
  ls.sum[[method]] <- BuyseTest(data=dt.BT,endpoint="Y_TTE1",treatment="Treatment",
                                type="TTE",censoring="event1",threshold=0, cpus = "all",
                                n.bootstrap=n.bootstrap,method=method)
}
cat("\n")

##### export
if(identical(save, TRUE)){
  results <- list(data = dt.BT,
                  Gehan = ls.sum$Gehan,
                  Peto = ls.sum$Peto,
                  Efron = ls.sum$Efron,
                  Peron = ls.sum$Peron)
  saveRDS(results, file = file.path(dirSave,"test-parallel.rds"))
}else if(identical(save, FALSE)){
  cat("* Previous version \n")
  GS <- readRDS(file = file.path(dirSave,"test-parallel.rds"))
  
  test_that("comparison with the previous version", {
    expect_equalBT(ls.sum$Gehan, GS$Gehan)
    expect_equalBT(ls.sum$Peto, GS$Peto)
    expect_equalBT(ls.sum$Efron, GS$Efron)
    expect_equalBT(ls.sum$Peron, GS$Peron)
  })
}

ls.sum$Peron@delta_boot
GS$Peron@delta_boot




