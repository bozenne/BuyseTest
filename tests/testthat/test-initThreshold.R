### test-initThreshold.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 22 2017 (18:37) 
## Version: 
## Last-Updated: apr 15 2018 (13:59) 
##           By: Brice Ozenne
##     Update #: 11
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
context("Check the function initializing the thresholds")

initThreshold <- BuyseTest:::initThreshold

## * binary outcomes
test_that("Only accept NA or 1/2 for binary outcomes", {

    initThreshold(threshold = NA, type = 1, D = 1, endpoint = c("Y1"))
    initThreshold(threshold = 1/2, type = 1, D = 1, endpoint = c("Y1"))
    expect_error(initThreshold(threshold = 0, type = 1, D = 1, endpoint = c("Y1")))

})

## * Continuous/TTE outcomes

test_that("Reject non-decreasing thresholds", {

    initThreshold(threshold = c(1,1/2), type = c(2,2), D = 2, endpoint = c("Y1","Y1"))
    expect_error(initThreshold(threshold = c(1,2), type = c(2,2), D = 2, endpoint = c("Y1","Y1")))
    expect_error(initThreshold(threshold = c(1,2), type = c(3,3), D = 1, endpoint = c("Y1","Y1")))
    
})

test_that("convert 0 to 1e-12 - threshold",{
  expect_equal(1, initThreshold(threshold=1, type=3, D=1, endpoint="time"))
  expect_equal(1e-12, initThreshold(threshold=0, type=3, D=1, endpoint="time"))
})

##----------------------------------------------------------------------
### test-initThreshold.R ends here
