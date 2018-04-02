### test-BuyseTest-Pairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 30 2018 (13:17) 
## Version: 
## Last-Updated: apr  2 2018 (15:21) 
##           By: Brice Ozenne
##     Update #: 41
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

## * one binary endpoint
## ** favorable
test_that("check favorable - Binary",{
    ## one pair
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data=data)
    expect_equal(as.double(BT@count_favorable),1)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
    ## BT@tableComparison

    ## seveal pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(BT@count_favorable),4)
    expect_equal(as.double(BT@count_unfavorable),0)
    expect_equal(as.double(BT@count_neutral),0)
    expect_equal(as.double(BT@count_uninf),0)
    ## BT@tableComparison

    ## strata pairs
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(BT@count_favorable),c(4,4))
    expect_equal(as.double(BT@count_unfavorable),c(0,0))
    expect_equal(as.double(BT@count_neutral),c(0,0))
    expect_equal(as.double(BT@count_uninf),c(0,0))
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

## * one continous endpoint
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

## * One time to event endpoint
## ** favorable
test_that("check favorable - TTE",{
    data <- data.frame(time = c(1,0), Treatment = c(1,0))
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


##----------------------------------------------------------------------
### test-BuyseTest-Pairs.R ends here
