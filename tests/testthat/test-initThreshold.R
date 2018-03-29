### test-initThreshold.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 22 2017 (18:37) 
## Version: 
## Last-Updated: dec 22 2017 (18:57) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
library(testthat)
verboseContext("Check the function initializing the thresholds")

initThreshold <- BuyseTest:::initThreshold

## * binary outcomes
test_that("Only accept NA or 1/2 for binary outcomes", {

    initThreshold(threshold = NA, type = 1, D = 1, endpoint = c("Y1"))
    initThreshold(threshold = 1/2, type = 1, D = 1, endpoint = c("Y1"))
    expect_error(initThreshold(threshold = 0, type = 1, D = 1, endpoint = c("Y1")))

})

## * Continuous/TTE outcomes

test_that("Reject non-decreasing thresholds", {

    expect_error(initThreshold(threshold = NA, type = 2, D = 1, endpoint = c("Y1")))
    expect_error(initThreshold(threshold = NA, type = 3, D = 1, endpoint = c("Y1")))
    
})

test_that("Reject non-decreasing thresholds", {

    initThreshold(threshold = c(1,1/2), type = c(2,2), D = 2, endpoint = c("Y1","Y1"))
    expect_error(initThreshold(threshold = 1:2, type = c(2,2), D = 2, endpoint = c("Y1","Y1")))
    
})

##----------------------------------------------------------------------
### test-initThreshold.R ends here
