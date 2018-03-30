# setwd("GitHub/BuyseTest/tests")
# setwd("tests")
# setwd("tests/testthat")
library(testthat)
library(BuyseTest)

#### additional packages
library(lava)
library(data.table)
library(survival)
source(file.path("FCT","FCT_check.R")) # file containing additional function for performing the tests

#### specifications
BuyseTest.option(trace = 0)
precision <- 10^{-7} # for expect_equal
save <- FALSE # TRUE to save results, FALSE to test, NULL to ignore
conv2df <- FALSE # should the data be in data.frame format or data.table format

if(identical(save, TRUE)){
  dirSave <- paste0("../Results-version",packageVersion("BuyseTest"))
  if(dir.exists(dirSave) == FALSE){dir.create(dirSave)}
}else{
  dirSave <- paste0("../Results-version","1.1")
}

#### run tests
test_check("BuyseTest")
