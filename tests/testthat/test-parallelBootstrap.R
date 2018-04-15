### test-parallelBootstrap.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: apr 16 2018 (00:08) 
##           By: Brice Ozenne
##     Update #: 10
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
context("Check parallel boostrap")

if(FALSE){
    library(BuyseTest)
    library(testthat)
}


## * settings
BuyseTest.options(trace = 0, n.permutation = 10)
n.patients <- 10

## * Simulate data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients,
                       n.C = n.patients,
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))

## * Bootstrap
method <- "Peron"
test_that("boostrap", {
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                    data = dt.sim, method = method)
    ## endpoint threshold pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf
    ## eventtime1     1e-12      100         48.6          23.86          0    27.54
    ## delta Delta CIinf.Delta CIsup.Delta n.permutation p.value 
    ## 0.247 0.247      -0.537       0.495            10     0.4 
})



#----------------------------------------------------------------------
### test-parallelBootstrap.R ends here
