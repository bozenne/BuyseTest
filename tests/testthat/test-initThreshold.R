### test-initThreshold.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 22 2017 (18:37) 
## Version: 
## Last-Updated: feb 20 2026 (14:39) 
##           By: Brice Ozenne
##     Update #: 31
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
context("Check the function initializing the thresholds")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

set.seed(10)
dt <- simBuyseTest(10)

## * binary outcomes
test_that("Do not accept threshold argument for binary outcomes", {
    expect_error(BuyseTest(treatment ~ bin(toxicity, threshold = 1),
                           data = dt,
                           method.inference = "none", trace = 0))
})

## * continuous outcomes

test_that("Reject non-decreasing thresholds", {

    expect_error(BuyseTest(treatment ~ cont(score, threshold = 1) + cont(score, threshold = 2),
                      data = dt,
                      method.inference = "none", trace = 0))

})

test_that("convert 0 to 1e-12 - threshold",{
    test <- BuyseTest(treatment ~ cont(score, threshold = 1) + cont(score),
                      data = dt,
                      method.inference = "none", trace = 0)

    expect_equal(as.double(test@threshold), c(1,1e-12))    
})

## * time to event outcomes

test_that("Reject non-decreasing thresholds", {

    expect_error(BuyseTest(treatment ~ cont(score, threshold = 1) + cont(score, threshold = 2),
                      data = dt,
                      method.inference = "none", trace = 0))

    expect_error(BuyseTest(treatment ~ tte(eventtime, status = status, threshold = 1) + tte(eventtime, censoring = status, threshold = 2),
                           data = dt,
                           method.inference = "none", trace = 0))
    
    
})

test_that("convert 0 to 1e-12 - threshold",{
    test <- BuyseTest(treatment ~ cont(score, threshold = 1) + cont(score),
                      data = dt,
                      method.inference = "none", trace = 0)

    expect_equal(as.double(test@threshold), c(1,1e-12))    
})

## * multiplicative threshold as log additive threshold
data(profil, package = "BuyseTest")
profil$duration_log <- log(profil$duration)
log_2 <- log(2)
profilR <- profil[profil$treatment %in% c("placebo","lowDose"),]
profilR$treatment <- droplevels(profilR$treatment)

test_that("multiplicative threshold as log additive threshold",{
    eMult <- BuyseTest(treatment ~ cont(duration, threshold = 2, operator = "*"), data = profilR, trace = 0, method.inference = "U-statistic")
    eLog <- BuyseTest(treatment ~ cont(duration_log, threshold = log_2), data = profilR, trace = 0, method.inference = "U-statistic")

    expect_equal(as.double(confint(eMult)-confint(eLog)),rep(0,6), tol =  1e-3)
})
##----------------------------------------------------------------------
### test-initThreshold.R ends here
