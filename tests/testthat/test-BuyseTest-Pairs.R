### test-BuyseTest-Pairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 30 2018 (13:17) 
## Version: 
## Last-Updated: mar 30 2018 (13:44) 
##           By: Brice Ozenne
##     Update #: 7
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

verboseContext("Check BuyseTest on tractable examples")

## * settings
BuyseTest.option(n.permutation = 0)

## * binary endpoint 
test_that("check favorable - Binary",{
  data_Bin <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
  BT_Bin1 <- BuyseTest(data=data_Bin,endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"))
  expect_equal(as.double(BT_Bin1@count_favorable),1)
  expect_equal(as.double(BT_Bin1@count_unfavorable),0)
  expect_equal(as.double(BT_Bin1@count_neutral),0)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
  
  BT_Bin1 <- BuyseTest(data=rbind(data_Bin,data_Bin),endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"),
                       n.bootstrap=0)
  expect_equal(as.double(BT_Bin1@count_favorable),4)
  expect_equal(as.double(BT_Bin1@count_unfavorable),0)
  expect_equal(as.double(BT_Bin1@count_neutral),0)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
})

test_that("check unfavorable - Binary",{
  data_Bin <- data.frame(toxicity1 = c(0,1), Treatment = c(1,0))
  BT_Bin1 <- BuyseTest(data=data_Bin,endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"),
                       n.bootstrap=0)
  expect_equal(as.double(BT_Bin1@count_favorable),0)
  expect_equal(as.double(BT_Bin1@count_unfavorable),1)
  expect_equal(as.double(BT_Bin1@count_neutral),0)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
  
  BT_Bin1 <- BuyseTest(data=rbind(data_Bin,data_Bin),endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"),
                       n.bootstrap=0)
  expect_equal(as.double(BT_Bin1@count_favorable),0)
  expect_equal(as.double(BT_Bin1@count_unfavorable),4)
  expect_equal(as.double(BT_Bin1@count_neutral),0)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
})

test_that("check neutral - Binary",{
  data_Bin <- data.frame(toxicity1 = c(1,1), Treatment = c(1,0))
  BT_Bin1 <- BuyseTest(data=data_Bin,endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"),
                       n.bootstrap=0)
  expect_equal(as.double(BT_Bin1@count_favorable),0)
  expect_equal(as.double(BT_Bin1@count_unfavorable),0)
  expect_equal(as.double(BT_Bin1@count_neutral),1)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
  
  BT_Bin1 <- BuyseTest(data=rbind(data_Bin,data_Bin),endpoint=c("toxicity1"),
                       treatment="Treatment", type=c("bin"),
                       n.bootstrap=0)
  expect_equal(as.double(BT_Bin1@count_favorable),0)
  expect_equal(as.double(BT_Bin1@count_unfavorable),0)
  expect_equal(as.double(BT_Bin1@count_neutral),4)
  expect_equal(as.double(BT_Bin1@count_uninf),0)
})


##----------------------------------------------------------------------
### test-BuyseTest-Pairs.R ends here
