### test-BuyseTest-engine.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2020 (09:46) 
## Version: 
## Last-Updated: mar 23 2020 (10:04) 
##           By: Brice Ozenne
##     Update #: 3
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

GS <- matrix(c(0, -0.05552322, -0.07409667, -0.12245246, 0.00344144, 0.11142578, 0.12031357, 0.12441761, -0.00674501, -0.26794284, -0.30166002, -0.35454813, 0.00674501, 0.16204436, 0.16145411, 0.12385689, 1, 0.61899652, 0.53947556, 0.32988924), 
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

    GS <- matrix(c(-0.05552322, -0.20777959, -0.22000332, -0.21987125, 0.11142578, 0.1354407, 0.13797685, 0.13741368, -0.26794284, -0.45285508, -0.46826381, -0.46723574, 0.16204436, 0.06648861, 0.06045122, 0.05941651, 0.61899652, 0.13634053, 0.1229436, 0.12162459), 
                 nrow = 4, 
                 ncol = 5, 
                 dimnames = list(c("eventtime1_1", "eventtime2_1", "eventtime1_0.25", "eventtime2_1e-12"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    test <- confint(e.BT2)
    attr(test,"n.resampling") <- NULL
    expect_equal(GS, test, tol = 1e-6)
})
##----------------------------------------------------------------------
### test-BuyseTest-engine.R ends here
