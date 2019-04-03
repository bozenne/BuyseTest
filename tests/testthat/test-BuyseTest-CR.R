### test-BuyseTest-CR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne, Eva Cantagallo
## Created: jul 12 2018 (16:58) 
## Version: 
## Last-Updated: feb 27 2019 (22:33) 
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
    library(data.table)
}

context("Check that BuyseTest with competing risks \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)


alphaE.X <- 2
alphaCR.X <- 1
alphaE.Y <- 3
alphaCR.Y <- 2
alpha.cens <- 1.5
n <- 1e2

## * Simulate CR data without censoring
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X), time2 = rexp(n, rate = alphaCR.X), group = "1"),
            data.frame(time1 = rexp(n, rate = alphaE.Y), time2 = rexp(n, rate = alphaCR.Y), group = "2"))
df$time <- pmin(df$time1,df$time2) ## first event
df$event <- (df$time2<df$time1)+1 ## type of event

## * Simulate CR data with censoring
set.seed(10)
df2 <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X), time2 = rexp(n, rate = alphaCR.X), time.cens = rexp(n, rate = alpha.cens), treatment = "1"),
            data.frame(time1 = rexp(n, rate = alphaE.Y), time2 = rexp(n, rate = alphaCR.Y), time.cens = rexp(n, rate = alpha.cens), treatment = "2"))
df2$time <- pmin(df2$time1, df2$time2, df2$time.cens) ## first status
df2$status <- ifelse(df2$time == df2$time1, 1, ifelse(df2$time == df2$time2, 2, 0)) ## type of event
df2$strata <- sample(c('a', 'b', 'c'), 2*n, replace = T)
df2$toxicity <- sample(0:1, 2*n, replace = T)
df2$strata2 <- sample(c('d', 'e', 'f'), 2*n, replace = T)


## * test
test_that("tte = 2 is equivalent to continuous with infty when cause=2", {
    e.BT <- BuyseTest(group ~ tte(time, censoring = event), data = df,
                      method.inference = "none", method.tte = "Gehan",
                      trace = 0)
    ## summary(e.BT)
    df$timeXX <- df$time
    df$timeXX[df$event==2] <- max(df$time)+1
    e.BT.bis <- BuyseTest(group ~ cont(timeXX), data = df,
                          method.inference = "none", trace = 0)
    ## summary(e.BT.bis)
    
    expect_equal(as.double(e.BT@Delta.netBenefit),
                 as.double(e.BT.bis@Delta.netBenefit))
    expect_equal(as.double(e.BT@delta.netBenefit),
                 as.double(e.BT.bis@delta.netBenefit))
    
})

test_that("BuyseTest package and Eva's R code give the same results with one endpoint and one stratum", {
  
  ## Net benefit computed with Eva's R code (see reproduce-results-CR.R)
  delta.R = 0.04377023
  
  ## Apply GPC with BuyseTest package
  BT = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df2)
  
  ## Test
  expect_equal(as.double(delta.R), as.double(BT@Delta.netBenefit))
  expect_equal(as.double(delta.R), as.double(BT@delta.netBenefit))
  
})

test_that("New package version gives the same results as previous one", {
  
  #### Net benefit computed with previous version 
  delta11 = 0.04377023 # one outcome, one stratum
  Delta11 = 0.04377023 # one outcome, one stratum
  delta13 = c(-0.06967276, 0.13925664, 0.13200801) # one outcome, 3 strata
  Delta13 = 0.0477878  # one outcome, 3 strata
  delta21 = c(0.04377023, 0.02028528) # 2 outcomes, one stratum
  Delta21 = c(0.04377023, 0.06405551) # 2 outcomes, one stratum
  delta23 = matrix(data = c(-0.06967276, 0.045860858, 0.13925664, 0.009916651, 0.13200801, 0.015133257), byrow = T, nrow = 3) # 2 outcomes, 3 strata
  Delta23 = c(0.04778780, 0.07473073) # 2 outcomes, 3 strata
  
  #### Apply GPC with new version of BuyseTest package
  ## One outcome, one stratum
  BT11.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df2)@delta.netBenefit
  BT11.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df2)@Delta.netBenefit

  ## One outcome, 3 strata
  BT13.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata, data = df2)@delta.netBenefit
  BT13.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata, data = df2)@Delta.netBenefit

  ## Two outcomes, one stratum
  BT21.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity), data = df2)@delta.netBenefit
  BT21.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity), data = df2)@Delta.netBenefit

  ## Two outcomes, one stratum
  BT23.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity) + strata, data = df2)@delta.netBenefit
  BT23.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity) + strata, data = df2)@Delta.netBenefit

  #### Tests
  expect_equal(delta11, as.double(BT11.d))
  expect_equal(Delta11, as.double(BT11.D))
  expect_equal(delta13, as.double(BT13.d))
  expect_equal(Delta13, as.double(BT13.D))
  expect_equal(delta21, as.double(BT21.d))
  expect_equal(Delta21, as.double(BT21.D))
  expect_equal(delta23[1,], as.double(BT23.d[1,]));expect_equal(delta23[2,], as.double(BT23.d[2,]));expect_equal(delta23[3,], as.double(BT23.d[3,]))
  expect_equal(Delta23, as.double(BT23.D))
  
})

test_that("Package give the same results when model.tte is (not) provided as an argument", {
  
  ## Create prodlim object to be inserted as an argument
  fit = prodlim(Hist(time, status) ~ treatment + strata + strata2, data = df2)
  
  ## Net benefit without passing model.tte
  B = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata + strata2, data = df2) 
  
  ## Net benefit with model.tte passed as an argument
  B.model = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata + strata2, data = df2, model.tte = fit) 
  
  ## Tests
  expect_equal(as.double(B@delta.netBenefit), as.double(B.model@delta.netBenefit))
  expect_equal(as.double(B@Delta.netBenefit), as.double(B.model@Delta.netBenefit))
  
})

test_that("When TTE endpoints are analyzed several times with different thresholds, the results are the same than those when analyzed once with lowest threshold", {

  ## Endpoint analyzed once with threshold = 0.5
  B1 = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df2) 
  
  ## Endpoint analyzed twice with threshold = c(1, 0.75, 0.5)
  B2 = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 1) + tte(time, censoring = status, threshold = 0.75) + tte(time, censoring = status, threshold = 0.5),
                 data = df2) 
  
  ## Tests
  expect_equal(as.double(B1@Delta.netBenefit), as.double(B2@Delta.netBenefit[3]))
  
})
######################################################################
### test-BuyseTest-CR.R ends here
