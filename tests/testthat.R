## * load packages
library(testthat)
library(BuyseTest)

library(lava)
library(data.table)
library(survival)
library(riskRegression)

## * run tests
test_check("BuyseTest")
