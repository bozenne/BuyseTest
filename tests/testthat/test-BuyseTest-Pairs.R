### test-BuyseTest-Pairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 30 2018 (13:17) 
## Version: 
## Last-Updated: apr  1 2018 (20:48) 
##           By: Brice Ozenne
##     Update #: 37
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

context("Check BuyseTest on tractable examples")

## * settings
BuyseTest.option(n.permutation = 0, trace = 0, keepComparison = TRUE)

## * binary endpoint
## ** favorable
test_that("check favorable - Binary",{
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data=data)
    expect_equal(as.double(BT@count_favorable),1)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
    ## BT@tableComparison
  
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = rbind(data, data))
    expect_equal(as.double(BT@count_favorable),4)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
    ## BT@tableComparison
  
})

## ** unfavorable
test_that("check unfavorable - Binary",{
    data <- data.frame(toxicity1 = c(0,1), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data=data)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),1)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
  
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data=rbind(data,data))
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),4)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
})

## ** neutral
test_that("check neutral - Binary",{
  data <- data.frame(toxicity1 = c(1,1), Treatment = c(1,0))
  BT <- BuyseTest(Treatment ~ bin(toxicity1), data=data)
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),0)
  expect_equal(as.double(BT@count_neutral),1)
  expect_equal(as.double(BT@count_uninf),0)
  
  BT <- BuyseTest(Treatment ~ bin(toxicity1), data=rbind(data,data))
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),0)
  expect_equal(as.double(BT@count_neutral),4)
  expect_equal(as.double(BT@count_uninf),0)
})

## ** NA as uninformative
test_that("check NA - Binary",{
  data <- data.frame(toxicity1 = c(NA, 1, 1), Treatment = c(1, 1,0))
  BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data)
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),0)
  expect_equal(as.double(BT@count_neutral),1)
  expect_equal(as.double(BT@count_uninf),1)
})

## * continous endpoint
## ** favorable
test_that("check favorable - continous",{
    data <- data.frame(size = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(BT@count_favorable),1)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)

    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = rbind(data, data))
    expect_equal(as.double(BT@count_favorable),4)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
})

## ** unfavorable
test_that("check unfavorable - continous",{
    data <- data.frame(size = c(-1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),1)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)

    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = rbind(data, data))
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),4)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
})

## ** neutral
test_that("check neutral - continous",{
    ## 0 threshold
    data <- data.frame(size = c(1,1), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data=data)
    BT.bis <- BuyseTest(Treatment ~ continuous(size), data=data)
    expect_equal(BT.bis,BT)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),1)
    expect_equal(as.double(BT@count_uninf),0)

    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = rbind(data, data))
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),4)
    expect_equal(as.double(BT@count_uninf),0)

    ## non 0 threshold 
    data <- data.frame(size = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 2), data=data)
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),1)
    expect_equal(as.double(BT@count_uninf),0)

    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 2), data = rbind(data, data))
    expect_equal(as.double(BT@count_favorable),0)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),4)
    expect_equal(as.double(BT@count_uninf),0)
})
## ** NA as uninformative
test_that("check NA - continuous",{
  data <- data.frame(size = c(NA, 1, 1), Treatment = c(1, 1, 0))
  BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data = data)
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),0)
  expect_equal(as.double(BT@count_neutral),1)
  expect_equal(as.double(BT@count_uninf),1)

  BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data = rbind(data, data))
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),0)
  expect_equal(as.double(BT@count_neutral),4)
  expect_equal(as.double(BT@count_uninf),4)
})

## * Time to event endpoint

##----------------------------------------------------------------------
### test-BuyseTest-Pairs.R ends here
