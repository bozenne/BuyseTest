### test-BuyseTest-Pairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 30 2018 (13:17) 
## Version: 
## Last-Updated: okt 30 2018 (16:33) 
##           By: Brice Ozenne
##     Update #: 151
##----------------------------------------------------------------------
## 
### Commentary: 
## Test BuyseTest on example where we have an explicit solution 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(testthat)
    library(BuyseTest)
}

context("Check BuyseTest on simple examples")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * one binary endpoint
## ** favorable
test_that("check favorable - 1 Binary",{
    
    ## one pair
    data <- data.frame(toxicity1 = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data=data)
    expect_equal(as.double(BT@count.favorable),1)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
    ## getPairScore(BT)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(BT@count.favorable),4)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
    ## getPairScore(BT)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(4,4))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
    ## getPairScore(BT)

})

## ** unfavorable
test_that("check unfavorable - 1 Binary",{
    ## one pair
    data <- data.frame(toxicity1 = c(0,1), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),1)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs
    data2 <- rbind(data,data)
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(4,4))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
    ## getPairScore(BT)

})

## ** neutral
test_that("check neutral - 1 Binary",{
    ## 1 pair
    data <- data.frame(toxicity1 = c(1,1), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),1)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),4)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(0,0))
    ## getPairScore(BT)
})

## ** NA as uninformative
test_that("check NA - 1 Binary",{
  data <- data.frame(toxicity1 = c(NA, 1, 1), Treatment = c(1, 1, 0))
  BT <- BuyseTest(Treatment ~ bin(toxicity1), data = data)
  expect_equal(as.double(BT@count.favorable),0)
  expect_equal(as.double(BT@count.unfavorable),0)
  expect_equal(as.double(BT@count.neutral),1)
  expect_equal(as.double(BT@count.uninf),1)

  ## with strata
  data2 <- rbind(cbind(data, strata = 0), cbind(data, strata = 1))
  BT <- BuyseTest(Treatment ~ bin(toxicity1) + strata, data = data2)
  expect_equal(as.double(BT@count.favorable),c(0,0))
  expect_equal(as.double(BT@count.unfavorable),c(0,0))
  expect_equal(as.double(BT@count.neutral),c(1,1))
  expect_equal(as.double(BT@count.uninf),c(1,1))
})

## * two binary endpoints
## ** unfavorable
test_that("check unfavorable - 2 Binary",{
    dt <- data.table(Treatment = 1:2,
                     toxicity1 = 0,
                     toxicity2 = c(1,0),
                     Id = 1:10)
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + bin(toxicity2), data = dt)
    ## getPairScore(BT)

    ## total pairs: 25
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,25))
    expect_equal(as.double(BT@count.neutral),c(25,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))

    expect_equal(as.double(BT@delta.netBenefit),c(0,-1))
    expect_equal(as.double(BT@Delta.netBenefit),c(0,-1))

    expect_equal(as.double(BT@delta.winRatio),c(NaN,0))
    expect_equal(as.double(BT@Delta.winRatio),c(NaN,0))
})

## ** mixed
test_that("check mixed - 2 Binary",{
    dt <- data.table(Treatment = c(1,1,2,2),
                     toxicity1 = 0,
                     toxicity2 = c(1,0,1,0))
    BT <- BuyseTest(Treatment ~ bin(toxicity1) + bin(toxicity2), data = dt,
                    keep.pairScore = TRUE)
    ## getPairScore(BT)

    expect_equal(as.double(BT@count.favorable),c(0,1))
    expect_equal(as.double(BT@count.unfavorable),c(0,1))
    expect_equal(as.double(BT@count.neutral),c(4,2))
    expect_equal(as.double(BT@count.uninf),c(0,0))

    expect_equal(as.double(BT@delta.netBenefit),c(0,0))
    expect_equal(as.double(BT@Delta.netBenefit),c(0,0))

    expect_equal(as.double(BT@delta.winRatio),c(NaN,1))
    expect_equal(as.double(BT@Delta.winRatio),c(NaN,1))
})

## * one continous endpoint
## ** favorable
test_that("check favorable - continous",{
    ## one pair
    data <- data.frame(size = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(BT@count.favorable),1)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(BT@count.favorable),4)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(4,4))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** unfavorable
test_that("check unfavorable - continous",{
    ## one pair
    data <- data.frame(size = c(-1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),1)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(4,4))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** neutral
test_that("check neutral - continous",{
    ## one pair, 0 threshold
    data <- data.frame(size = c(1,1), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data=data)
    BT.bis <- BuyseTest(Treatment ~ continuous(size), data=data)
    expect_equal(BT.bis,BT)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),1)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs, 0 threshold
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),4)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata, 0 threshold
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(0,0))

    ## 1 pair, non 0 threshold 
    data <- data.frame(size = c(1,0), Treatment = c(1,0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 2), data=data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),1)
    expect_equal(as.double(BT@count.uninf),0)

    ## several pairs, non 0 threshold
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 2), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),4)
    expect_equal(as.double(BT@count.uninf),0)

    ## with strata, non 0 threshold
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 2) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** NA as uninformative
test_that("check NA - continuous",{
    data <- data.frame(size = c(NA, 1, 1), Treatment = c(1, 1, 0))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),1)
    expect_equal(as.double(BT@count.uninf),1)

    ## more pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0), data = data2)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),4)
    expect_equal(as.double(BT@count.uninf),4)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ continuous(size, threshold = 0) + strata, data = data3)
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(4,4))
})

## * one time to event endpoint (no extrapolation needed or possible)
## ** favorable
test_that("check favorable - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(1,0), Treatment = c(1,0), status = 1)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data)
    expect_equal(as.double(BT@count.favorable),1)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data2)    
    expect_equal(as.double(BT@count.favorable),4)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold, strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status) + strata, data = data3)    
    expect_equal(as.double(BT@count.favorable),c(4,4))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))

    ## 1 threshold, strata
    BT <- BuyseTest(Treatment ~ tte(time, threshold = 1, censoring = status) + strata, data = data3)    
    expect_equal(as.double(BT@count.favorable),c(4,4))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** unfavorable
test_that("check unfavorable - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(0,1), Treatment = c(1,0), status = 1)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),1)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data2)    
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status) + strata, data = data3)    
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(4,4))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
    
    ## 1 threshold, strata
    BT <- BuyseTest(Treatment ~ tte(time, threshold = 1, censoring = status) + strata, data = data3)    
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(4,4))
    expect_equal(as.double(BT@count.neutral),c(0,0))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** neutral
test_that("check neutral - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(1,1), Treatment = c(1,0), status = 1)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data)
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),1)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data2)    
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),4)
    expect_equal(as.double(BT@count.uninf),0)

    ## 0 threshold strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status) + strata, data = data3)    
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(0,0))
    
    ## 1 threshold, strata
    data4 <- data3
    data4[data4$Treatment==1,"time"] <- 2
    BT <- BuyseTest(Treatment ~ tte(time, threshold = 3, censoring = status) + strata, data = data4)    
    expect_equal(as.double(BT@count.favorable),c(0,0))
    expect_equal(as.double(BT@count.unfavorable),c(0,0))
    expect_equal(as.double(BT@count.neutral),c(4,4))
    expect_equal(as.double(BT@count.uninf),c(0,0))
})

## ** NA as uninformative
test_that("check NA - time to event",{
    ## censored after the event in the other arm
    data <- data.frame(time = c(2,1), Treatment = c(1,0), status = c(0,1))
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data)
    expect_equal(as.double(BT@count.favorable),1)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)

    ## censored at the same time as the event in the other arm
    data <- data.frame(time = c(1,2,1), Treatment = c(1,1,0), status = c(0,1,1))
    BT <- BuyseTest(Treatment ~ tte(time, censoring = status), data = data)
    expect_equal(as.double(BT@count.favorable),2)
    expect_equal(as.double(BT@count.unfavorable),0)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
})

## * one time to event endpoint (extrapolation using Gehan/Peron)
## ** with 2 pairs
# one is censored
dt.2pairs <- data.table(id = 1:4,
                        time = c(10,20,12,32),
                        cens = c(0,1,1,1),
                        treat = c("control","6-MP","control","6-MP")) 


## *** Gehan
test_that("2 pairs - Gehan - no correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, censoring = cens), data = dt.2pairs,
                    method.tte="Gehan", correction.uninf = FALSE)
  
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),2)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),2)
 
    expect_equal(as.double(BT@delta.netBenefit),-1/2) ## 
    expect_equal(as.double(BT@Delta.netBenefit),-1/2) ##
    
    expect_equal(as.double(BT@delta.winRatio),0)
    expect_equal(as.double(BT@Delta.winRatio),0)
})

test_that("2 pairs - Gehan - correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, censoring = cens), data = dt.2pairs,
                    method.tte="Gehan", correction.uninf = TRUE)
  
    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
 
    expect_equal(as.double(BT@delta.netBenefit),-1) ## 
    expect_equal(as.double(BT@Delta.netBenefit),-1) ##
    
    expect_equal(as.double(BT@delta.winRatio),0)
    expect_equal(as.double(BT@Delta.winRatio),0)
})

## *** Peron
test_that("2 pairs - Peron - no correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, censoring = cens), data = dt.2pairs,
                    method.tte="Peron", correction.uninf = FALSE)
  
    ## different survival curve per groups (denoting S survival time and T group)
    ## P[T>=t|T=0,S>=10] = 1 (t=<12), 0 (t>12)
    ## P[T>=t|T=1,S>=10] = 1 (t=<20), 1/2 (20<t=<32), 0 (t>32)

    ## all pairScores (see getPairScore(BT))
    ## 10* vs 20 : unfavorable
    ## 10* vs 32 : unfavorable
    ## 12 vs 20 : unfavorable
    ## 12 vs 32 : unfavorable

    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
  
    expect_equal(as.double(BT@delta.netBenefit),-1)
    expect_equal(as.double(BT@Delta.netBenefit),-1)

    expect_equal(as.double(BT@delta.winRatio),0)
    expect_equal(as.double(BT@Delta.winRatio),0)
  
})

test_that("2 pairs - Peron - correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, censoring = cens), data = dt.2pairs,
                    method.tte="Peron", correction.uninf = TRUE)
  
    ## different survival curve per groups (denoting S survival time and T group)
    ## P[T>=t|T=0,S>=10] = 1 (t=<12), 0 (t>12)
    ## P[T>=t|T=1,S>=10] = 1 (t=<20), 1/2 (20<t=<32), 0 (t>32)

    ## all pairScores (see getPairScore(BT))
    ## 10* vs 20 : unfavorable
    ## 10* vs 32 : unfavorable
    ## 12 vs 20 : unfavorable
    ## 12 vs 32 : unfavorable

    expect_equal(as.double(BT@count.favorable),0)
    expect_equal(as.double(BT@count.unfavorable),4)
    expect_equal(as.double(BT@count.neutral),0)
    expect_equal(as.double(BT@count.uninf),0)
  
    expect_equal(as.double(BT@delta.netBenefit),-1)
    expect_equal(as.double(BT@Delta.netBenefit),-1)

    expect_equal(as.double(BT@delta.winRatio),0)
    expect_equal(as.double(BT@Delta.winRatio),0)
  
})

## ** Predictible events
## all events in the group 1 will happen at time 10
## all events in the group 0 will happen before time 2
## so only in favor of group 1
M.all <- rbind(c(time = 1, status = 1, trt = 0), 
               c(time = 1, status = 1, trt = 0), 
               c(time = 1, status = 1, trt = 0),
               c(time = 1, status = 0, trt = 0), ## obs 1 pair 1
               c(time = 2, status = 1, trt = 0),
               c(time = 2, status = 1, trt = 0),
               c(time = 1, status = 0, trt = 1),  ## obs 2 pair 1
               c(time = 10, status = 1, trt = 1))
df.all <- as.data.frame(M.all)

## plot(prodlim(Hist(time, status) ~ trt, data = df.all))

test_that("Peron - predictible events",{
    df.all$trt <- abs(df.all$trt)
    BT <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = status), data = df.all,
                    method.tte="Peron", correction.uninf = FALSE)
    expect_equal(as.double(BT@Delta.netBenefit),1)

    df.all$trt <- -abs(df.all$trt)
    BT <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = status), data = df.all,
                    method.tte="Peron", correction.uninf = FALSE)
    expect_equal(as.double(BT@Delta.netBenefit),-1)

})

##----------------------------------------------------------------------
### test-BuyseTest-Pairs.R ends here







