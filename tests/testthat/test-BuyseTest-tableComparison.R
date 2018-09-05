### test-BuyseTest-tableComparison.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:33) 
## Version: 
## Last-Updated: sep  5 2018 (09:47) 
##           By: Brice Ozenne
##     Update #: 20
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
}

context("Check tableComparison matches the summary of BuyseTest objects")

## * Settings
n.patients <- c(90,100)
BuyseTest.options(check = FALSE,
                  keep.comparison = TRUE,
                  method.inference = "none",
                  trace = 0)


## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))

dtRed.sim <- dt.sim[, .SD[1:50], by = "Treatment"]

## * test against tableComparison (no correction)
formula <- Treatment ~ tte(eventtime1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(eventtime1, 0.25, status1) + cont(score1, 0.5)
BT.mixed <- BuyseTest(formula,
                      data = dt.sim, method.tte = "Peron")

test_that("Full data", {
    test <- aggrTableComparison(table = BT.mixed@tableComparison,
                                correct.tte = BT.mixed@method.tte$correction)
        
    expect_equal(unname(tail(BT.mixed@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(BT.mixed@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])
})

## * test against tableComparison (correction)
formula <- Treatment ~ tte(eventtime1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(eventtime1, 0.25, status1) + cont(score1, 0.5)
BT.mixed <- BuyseTest(formula,
                      data = dt.sim, method.tte = "Peron corrected")

test_that("Full data", {
    test <- aggrTableComparison(table = BT.mixed@tableComparison,
                                correct.tte = BT.mixed@method.tte$correction)
        
    expect_equal(unname(tail(BT.mixed@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(BT.mixed@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])
})



### does not work because of the estimation of the survival
## test_that("First 50 patients in each arm", {
    ## BT.mixedRed <- BuyseTest(formula,
                             ## data = dtRed.sim, method.tte = "Peron")
    ## summary(BT.mixedRed, percentage = FALSE)

    ## BT.mixedRed@tableComparison[[1]][1:5]
    
    ## test <- tableComparison2Delta(BT.mixed@tableComparison,
                                  ## correct.tte = BT.mixed@method.tte$correction,
                                  ## maxData.T = which(dt.sim$Treatment==1)[50],
                                  ## maxData.C = which(dt.sim$Treatment==0)[50])
    ## expect_equal(BT.mixedRed@Delta.netChance,test[["Delta.netChance"]])
    ## expect_equal(BT.mixedRed@Delta.winRatio,test[["Delta.winRatio"]])
## })

##----------------------------------------------------------------------
### test-BuyseTest-tableComparison.R ends here
