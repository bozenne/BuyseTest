### test-initThreshold.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 22 2017 (18:37) 
## Version: 
## Last-Updated: feb 27 2019 (22:34) 
##           By: Brice Ozenne
##     Update #: 24
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
test_that("Only accept NA or 1/2 for binary outcomes", {

    test <- BuyseTest(Treatment ~ bin(toxicity, threshold = 0.5) + bin(status),
                      data = dt,
                      method.inference = "none", trace = 0)
    expect_true(all(test@threshold==0.5))
    
    expect_error(BuyseTest(Treatment ~ bin(toxicity, threshold = 1),
                           data = dt,
                           method.inference = "none", trace = 0))
})

## * continuous outcomes

test_that("Reject non-decreasing thresholds", {

    expect_error(BuyseTest(Treatment ~ cont(score, threshold = 1) + cont(score, threshold = 2),
                      data = dt,
                      method.inference = "none", trace = 0))

})

test_that("convert 0 to 1e-12 - threshold",{
    test <- BuyseTest(Treatment ~ cont(score, threshold = 1) + cont(score),
                      data = dt,
                      method.inference = "none", trace = 0)

    expect_equal(test@threshold, c(1,1e-12))    
})

## * time to event outcomes

test_that("Reject non-decreasing thresholds", {

    expect_error(BuyseTest(Treatment ~ cont(score, threshold = 1) + cont(score, threshold = 2),
                      data = dt,
                      method.inference = "none", trace = 0))

    expect_error(BuyseTest(Treatment ~ tte(eventtime, censoring = status, threshold = 1) + tte(eventtime, censoring = status, threshold = 2),
                           data = dt,
                           method.inference = "none", trace = 0))
    
    
})

test_that("convert 0 to 1e-12 - threshold",{
    test <- BuyseTest(Treatment ~ cont(score, threshold = 1) + cont(score),
                      data = dt,
                      method.inference = "none", trace = 0)

    expect_equal(test@threshold, c(1,1e-12))    
})


##----------------------------------------------------------------------
### test-initThreshold.R ends here
