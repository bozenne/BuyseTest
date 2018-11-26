### test-BuyseTest-operator.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  2 2018 (15:21) 
## Version: 
## Last-Updated: okt 30 2018 (16:21) 
##           By: Brice Ozenne
##     Update #: 18
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

context("Check that the option operator in BuyseTest is working correctly \n")


## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * one pair
test_that("check - 1 pair",{
    ## binary 
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),1)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
    ## getPairScore(BT)

    ## continuous
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ cont(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),1)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

})

##----------------------------------------------------------------------
### test-BuyseTest-operator.R ends here
