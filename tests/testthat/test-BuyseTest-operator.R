### test-BuyseTest-operator.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  2 2018 (15:21) 
## Version: 
## Last-Updated: apr  2 2018 (15:35) 
##           By: Brice Ozenne
##     Update #: 5
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

verboseContext("Check that the option operator in BuyseTest is working correctly \n")

## * settings
BuyseTest.option(n.permutation = 0, trace = 0)

## * settings
BuyseTest.option(n.permutation = 0, trace = 0, keepComparison = TRUE)

## * one pair
test_that("check - 1 pair",{
    ## binary 
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),1)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
    ## BT@tableComparison

    ## continuous
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ cont(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),1)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)

})

##----------------------------------------------------------------------
### test-BuyseTest-operator.R ends here
