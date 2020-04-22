### test-BuyseTest-engine.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2020 (09:46) 
## Version: 
## Last-Updated: apr 22 2020 (17:19) 
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

context("Check BuyseTest without strata")

## * Settings
n.patients <- c(60,65)
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1/(1:3), rates.Censoring.T = rep(1,3)))

## * Compare
test_that("TTE with decreasing thresholds",{
    iFormula <- treatment ~ tte(eventtime1,status1,5) + tte(eventtime1,status1,1) + tte(eventtime1,status1,0.5) + tte(eventtime1,status1,0)

    BuyseTest.options(engine = "GPC_cpp")
    e.BT1 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")

    BuyseTest.options(engine = "GPC2_cpp")
    e.BT2 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")

    expect_equal(confint(e.BT1)[,"estimate"],confint(e.BT2)[,"estimate"], tol = 1e-6)
    expect_equal(confint(e.BT1)[,"se"],confint(e.BT2)[,"se"], tol = 1e-2)
    ## expected difference because GPC2 compute the influence function over all pairs,
    ## while GPC only over pairs with informative scores.
    ## the difference is expected to be small though in large samples

    GS <- matrix(c(0, -0.05552322, -0.07409667, -0.12245246, 0, 0.10511706, 0.1165731, 0.12540316, 0, -0.25639241, -0.29494435, -0.35626146, 0, 0.14994324, 0.15426618, 0.12578738, 1, 0.59811202, 0.52654102, 0.3337042), 
                 nrow = 4, 
                 ncol = 5, 
                 dimnames = list(c("eventtime1_5", "eventtime1_1", "eventtime1_0.5", "eventtime1_1e-12"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    test <- confint(e.BT2)
    attr(test,"n.resampling") <- NULL
    expect_equal(GS,test, tol = 1e-6)
})

test_that("different TTE with decreasing thresholds",{
    iFormula <- treatment ~ tte(eventtime1,status1,1) + tte(eventtime2,status2,1) + tte(eventtime1,status1,0.25) + tte(eventtime2,status2,0)

    BuyseTest.options(engine = "GPC_cpp")
    e.BT1 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")
    BuyseTest.options(engine = "GPC2_cpp")
    e.BT2 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")
    expect_equal(confint(e.BT1)[,"estimate"],confint(e.BT2)[,"estimate"], tol = 1e-6)
    expect_equal(confint(e.BT1)[,"se"],confint(e.BT2)[,"se"], tol = 1e-2)
    ## expected difference because GPC2 compute the influence function over all pairs,
    ## while GPC only over pairs with informative scores.
    ## the difference is expected to be small though in large samples

    GS <- matrix(c(-0.05552322, -0.20777959, -0.22000332, -0.21987125, 0.10511706, 0.13520906, 0.13793648, 0.1365901, -0.25639241, -0.45247782, -0.46819889, -0.4659088, 0.14994324, 0.06601621, 0.06036837, 0.05772614, 0.59811202, 0.13567039, 0.12283402, 0.11939353), 
                 nrow = 4, 
                 ncol = 5, 
                 dimnames = list(c("eventtime1_1", "eventtime2_1", "eventtime1_0.25", "eventtime2_1e-12"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    test <- confint(e.BT2)
    attr(test,"n.resampling") <- NULL
    expect_equal(GS, test, tol = 1e-3)
})
##----------------------------------------------------------------------
### test-BuyseTest-engine.R ends here
