### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: sep 13 2019 (09:32) 
##           By: Brice Ozenne
##     Update #: 7
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

    expect_true(all(e.powerBT@results[,.SD$netBenefit[1]==.SD$netBenefit[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$netBenefit.se[1]==.SD$netBenefit.se[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$netBenefit.lower[1]==.SD$netBenefit.lower[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$netBenefit.upper[1]==.SD$netBenefit.upper[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$netBenefit.p.value[1]==.SD$netBenefit.p.value[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$winRatio[1]==.SD$winRatio[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$winRatio.se[1]==.SD$winRatio.se[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$winRatio.lower[1]==.SD$winRatio.lower[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$winRatio.upper[1]==.SD$winRatio.upper[1],by="simulation"][[2]]))
    expect_true(all(e.powerBT@results[,.SD$winRatio.p.value[1]==.SD$winRatio.p.value[1],by="simulation"][[2]]))
})

######################################################################
### test-BuysePower.R ends here
