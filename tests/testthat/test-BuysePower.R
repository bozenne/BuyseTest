### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: nov 12 2019 (11:10) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(testthat)
    library(BuyseTest)
    library(data.table)
}

context("Check BuysePower \n")

## * binary endpoint
test_that("binary endpoint", {
    e.powerBT <- powerBuyseTest(sim = simBuyseTest, sample.size = c(50,50), n.rep = 5,
                                formula = treatment ~ bin(toxicity),
                                method.inference = "u-statistic", trace = 0)
})

######################################################################
### test-BuysePower.R ends here
