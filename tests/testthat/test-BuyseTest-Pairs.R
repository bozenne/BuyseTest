### test-BuyseTest-Pairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 30 2018 (13:17) 
## Version: 
## Last-Updated: apr  2 2020 (17:18) 
##           By: Brice Ozenne
##     Update #: 187
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
    library(data.table)
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
    data <- data.frame(toxicity1 = c(1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data=data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
    ## getPairScore(BT)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
    ## getPairScore(BT)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
    ## getPairScore(BT)
})

## ** unfavorable
test_that("check unfavorable - 1 Binary",{
    ## one pair
    data <- data.frame(toxicity1 = c(0,1),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs
    data2 <- rbind(data,data)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
    ## getPairScore(BT)
})

## ** neutral
test_that("check neutral - 1 Binary",{
    ## 1 pair
    data <- data.frame(toxicity1 = c(1,1),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ bin(toxicity1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ bin(toxicity1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
    ## getPairScore(BT)
})

## ** NA as uninformative
test_that("check NA - 1 Binary",{
  data <- data.frame(toxicity1 = c(NA, 1, 1),
                     treatment = c(1, 1, 0),
                     stringsAsFactors = FALSE)
  BT <- BuyseTest(treatment ~ bin(toxicity1), data = data)
  expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
  expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
  expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
  expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),1)

  ## with strata
  data2 <- rbind(cbind(data, strata = 0), cbind(data, strata = 1))
  BT <- BuyseTest(treatment ~ bin(toxicity1) + strata, data = data2)
  expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
  expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
  expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(1,1))
  expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(1,1))
})

## * two binary endpoints
## ** unfavorable
test_that("check unfavorable - 2 Binary",{
    dt <- data.table(treatment = 1:2,
                     toxicity1 = 0,
                     toxicity2 = c(1,0),
                     Id = 1:10)
    BT <- BuyseTest(treatment ~ bin(toxicity1) + bin(toxicity2), data = dt)
    
    ## total pairs: 25
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),c(0,25))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),c(25,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),c(0,0))

    expect_equal(as.double(coef(BT, statistic = "netBenefit")),c(0,-1))
    expect_equal(as.double(coef(BT, statistic = "winRatio")),c(NaN,0))
})

## ** mixed
test_that("check mixed - 2 Binary",{
    dt <- data.table(treatment = c(1,1,2,2),
                     toxicity1 = 0,
                     toxicity2 = c(1,0,1,0))
    BT <- BuyseTest(treatment ~ bin(toxicity1) + bin(toxicity2), data = dt,
                    keep.pairScore = TRUE)
    ## getPairScore(BT)
    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),c(0,1))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),c(0,1))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),c(4,2))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),c(0,0))

    expect_equal(as.double(coef(BT, statistic = "netBenefit")),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "winRatio")),c(NaN,1))
})

## * one continous endpoint
## ** favorable
test_that("check favorable - continous",{
    ## one pair
    data <- data.frame(size = c(1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** unfavorable
test_that("check unfavorable - continous",{
    ## one pair
    data <- data.frame(size = c(-1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** neutral
test_that("check neutral - continous",{
    ## one pair, 0 threshold
    data <- data.frame(size = c(1,1),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 0), data=data)
    BT.bis <- BuyseTest(treatment ~ continuous(size), data=data)
    expect_equal(BT.bis,BT)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs, 0 threshold
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata, 0 threshold
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 1) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))

    ## 1 pair, non 0 threshold 
    data <- data.frame(size = c(1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 2), data=data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## several pairs, non 0 threshold
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 2), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## with strata, non 0 threshold
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 2) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** NA as uninformative
test_that("check NA - continuous",{
    data <- data.frame(size = c(NA, 1, 1),
                       treatment = c(1, 1, 0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 0), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),1)

    ## more pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 0), data = data2)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),4)

    ## with strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ continuous(size, threshold = 0) + strata, data = data3)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(4,4))
})

## * one time to event endpoint (no extrapolation needed or possible)
## ** favorable
test_that("check favorable - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(1,0),
                       treatment = c(1,0),
                       status = 1,
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data2)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold, strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ tte(time, status = status) + strata, data = data3)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))

    ## 1 threshold, strata
    BT <- BuyseTest(treatment ~ tte(time, threshold = 1, status = status) + strata, data = data3)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** unfavorable
test_that("check unfavorable - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(0,1),
                       treatment = c(1,0),
                       status = 1,
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data2)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ tte(time, status = status) + strata, data = data3)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
    
    ## 1 threshold, strata
    BT <- BuyseTest(treatment ~ tte(time, threshold = 1, status = status) + strata, data = data3)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** neutral
test_that("check neutral - time to event",{
    ## 0 threshold, 1 pair
    data <- data.frame(time = c(1,1),
                       treatment = c(1,0),
                       status = 1,
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold, several pairs
    data2 <- rbind(data, data)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data2)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## 0 threshold strata
    data3 <- rbind(cbind(data2, strata = 0), cbind(data2, strata = 1))
    BT <- BuyseTest(treatment ~ tte(time, status = status) + strata, data = data3)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
    
    ## 1 threshold, strata
    data4 <- data3
    data4[data4$treatment==1,"time"] <- 2
    BT <- BuyseTest(treatment ~ tte(time, threshold = 3, status = status) + strata, data = data4)    
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE, stratified = TRUE)),c(0,0))
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE, stratified = TRUE)),c(4,4))
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE, stratified = TRUE)),c(0,0))
})

## ** NA as uninformative
test_that("check NA - time to event",{
    ## censored after the event in the other arm
    data <- data.frame(time = c(2,1),
                       treatment = c(1,0),
                       status = c(0,1),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ tte(time, status = status),
                    data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

    ## censored at the same time as the event in the other arm
    data <- data.frame(time = c(1,2,1),
                       treatment = c(1,1,0),
                       status = c(0,1,1),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ tte(time, status = status), data = data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),2)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
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
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, status = cens), data = dt.2pairs,
                    scoring.rule = "Gehan", correction.uninf = FALSE)

    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),2)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),2)
 
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),-1/2) ## 
    expect_equal(as.double(coef(BT, statistic = "winRatio")),c(0))
})

test_that("2 pairs - Gehan - correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, status = cens), data = dt.2pairs,
                    scoring.rule = "Gehan", correction.uninf = TRUE)
  
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
 
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),-1) ## 
    expect_equal(as.double(coef(BT, statistic = "winRatio")),0) ## 
})

## *** Peron
test_that("2 pairs - Peron - no correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, status = cens), data = dt.2pairs,
                    scoring.rule = "Peron", correction.uninf = FALSE)
  
    ## different survival curve per groups (denoting S survival time and T group)
    ## P[T>=t|T=0,S>=10] = 1 (t=<12), 0 (t>12)
    ## P[T>=t|T=1,S>=10] = 1 (t=<20), 1/2 (20<t=<32), 0 (t>32)

    ## all pairScores (see getPairScore(BT))
    ## 10* vs 20 : unfavorable
    ## 10* vs 32 : unfavorable
    ## 12 vs 20 : unfavorable
    ## 12 vs 32 : unfavorable

    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
  
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),-1)
    expect_equal(as.double(coef(BT, statistic = "winRatio")),0)
  
})

test_that("2 pairs - Peron - correction",{
    BT <- BuyseTest(treat ~ tte(time, threshold = 0, status = cens), data = dt.2pairs,
                    scoring.rule = "Peron", correction.uninf = TRUE)
  
    ## different survival curve per groups (denoting S survival time and T group)
    ## P[T>=t|T=0,S>=10] = 1 (t=<12), 0 (t>12)
    ## P[T>=t|T=1,S>=10] = 1 (t=<20), 1/2 (20<t=<32), 0 (t>32)

    ## all pairScores (see getPairScore(BT))
    ## 10* vs 20 : unfavorable
    ## 10* vs 32 : unfavorable
    ## 12 vs 20 : unfavorable
    ## 12 vs 32 : unfavorable

    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),4)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
  
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),-1)
    expect_equal(as.double(coef(BT, statistic = "winRatio")),0)
  
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
df.all <- as.data.frame(M.all, stringsAsFactors = FALSE)

## plot(prodlim(Hist(time, status) ~ trt, data = df.all))

test_that("Peron - predictible events",{
    df.all$trt <- abs(df.all$trt)
    BT <- BuyseTest(trt ~ tte(time, threshold = 0, status = status), data = df.all,
                    scoring.rule="Peron", correction.uninf = FALSE)
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),1)

    df.all$trt <- -abs(df.all$trt)
    BT <- BuyseTest(trt ~ tte(time, threshold = 0, status = status), data = df.all,
                    scoring.rule="Peron", correction.uninf = FALSE)
    expect_equal(as.double(coef(BT, statistic = "netBenefit")),-1)

})

##----------------------------------------------------------------------
### test-BuyseTest-Pairs.R ends here







