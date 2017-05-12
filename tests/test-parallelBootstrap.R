### test-parallelBootstrap.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: maj 12 2017 (18:07) 
##           By: Brice
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: Check whether the parallel bootstrap can be run
## Note: it is not integrated in the testthat directory since the test was failing when I was putting it there
##
## path <- "~/GitHub/BuyseTest/tests/"
## setwd(path)
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
source(file.path("FCT","FCT_check.R")) # file containing additional function for performing the tests
verboseContext("Check parallel boostrap")

# {{{ load package
library(BuyseTest)
library(testthat)
source(file.path("FCT","FCT_check.R"))
# }}}

# {{{ parametrisation
BuyseTest.options(trace = 0, keep.bootstrap = TRUE)
n.patients <- 100
n.bootstrap <- 100

save <- TRUE # TRUE to save results, FALSE to test, NULL to ignore
if(identical(save, TRUE)){
  dirSave <- paste0("Results-version",utils::packageVersion("BuyseTest"))
  if(dir.exists(dirSave) == FALSE){dir.create(dirSave)}
}else{
  dirSave <- paste0("Results-version","1.1")
}
# }}}

# {{{ generate data
set.seed(10)
dt.BT <- simulBT(n.patients)
# }}}

# {{{ run and test parallel computation
ls.sum <- list()
for(method in c("Gehan","Peto","Efron","Peron")){ # method <- "Gehan"
  cat(method, " ")
  ls.sum[[method]] <- BuyseTest(data=dt.BT,endpoint="eventtime",treatment="Treatment",
                                type="TTE",censoring="status",threshold=0, cpus = "all",
                                n.bootstrap=n.bootstrap,method=method)
}
cat("\n")
# }}}

# {{{ compare results with previous version
if(identical(save, FALSE)){
    cat("* Previous version \n")
    GS <- readRDS(file = file.path(dirSave,"test-parallel.rds"))
  
    test_that("comparison with the previous version", {
        expect_equalBT(ls.sum$dt.BT, GS$dt.BT)
        expect_equalBT(ls.sum$Gehan, GS$Gehan)
        expect_equalBT(ls.sum$Peto, GS$Peto)
        expect_equalBT(ls.sum$Efron, GS$Efron)
        expect_equalBT(ls.sum$Peron, GS$Peron)
    })
}
# }}}

# {{{ export
if(identical(save, TRUE)){
    results <- list(data = dt.BT,
                    Gehan = ls.sum$Gehan,
                    Peto = ls.sum$Peto,
                    Efron = ls.sum$Efron,
                    Peron = ls.sum$Peron)
    saveRDS(results, file = file.path(dirSave,"test-parallel.rds"))
}
# }}}


#----------------------------------------------------------------------
### test-parallelBootstrap.R ends here
