### test-BuyseTest-CR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne, Eva Cantagallo
## Created: jul 12 2018 (16:58) 
## Version: 
## Last-Updated: dec  3 2020 (13:41) 
##           By: Brice Ozenne
##     Update #: 53
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
library(prodlim)

context("Check that BuyseTest with competing risks \n")

## * Gehan scoring rule (Maintainer: Brice)
## ** settings
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
n.big <- 5e4
true.sHR <- c(0, 0.5, 1, 2) 
p.1C <- 0.55
v <- -0.30

## ** Simulate data
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X),
                       time2 = rexp(n, rate = alphaCR.X),
                       group = "1",
                       stringsAsFactors = FALSE),
            data.frame(time1 = rexp(n, rate = alphaE.Y),
                       time2 = rexp(n, rate = alphaCR.Y),
                       group = "2",
                       stringsAsFactors = FALSE))
df$time <- pmin(df$time1,df$time2) ## first event
df$event <- (df$time2<df$time1)+1 ## type of event

## ** test
test_that("tte = 2 is equivalent to continuous with infty when cause=2", {
    e.BT <- BuyseTest(group ~ tte(time, status = event), data = df,
                      method.inference = "none", scoring.rule = "Gehan",
                      trace = 0)
    ## summary(e.BT)
    df$timeXX <- df$time
    df$timeXX[df$event==2] <- max(df$time)+1
    e.BT.bis <- BuyseTest(group ~ cont(timeXX), data = df,
                          method.inference = "none", trace = 0)
    ## summary(e.BT.bis)
    
    expect_equal(as.double(coef(e.BT)),
                 as.double(coef(e.BT.bis)),
                 tol = 1e-6)
})


## * Peron scoring rule (Maintainer: Eva)
alphaE.X <- 2
alphaCR.X <- 1
alphaE.Y <- 3
alphaCR.Y <- 2
alpha.cens <- 1.5
n <- 1e2
n.big <- 7.5e4
sHR <- 0.5 # c(0, 0.5, 1, 3) 
true.Delta = 0.1519 # c(1/3, 0.1519, 0, -0.4012)

## ** Simulate CR data without censoring
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X),
                       time2 = rexp(n, rate = alphaCR.X),
                       group = "1",
                       stringsAsFactors = FALSE),
            data.frame(time1 = rexp(n, rate = alphaE.Y),
                       time2 = rexp(n, rate = alphaCR.Y),
                       group = "2",
                       stringsAsFactors = FALSE))
df$time <- pmin(df$time1,df$time2) ## first event
df$event <- (df$time2<df$time1)+1 ## type of event

## ** Simulate CR data with censoring
set.seed(10)
df2 <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X),
                        time2 = rexp(n, rate = alphaCR.X),
                        time.cens = rexp(n, rate = alpha.cens),
                        treatment = "1",
                        stringsAsFactors = FALSE),
            data.frame(time1 = rexp(n, rate = alphaE.Y),
                       time2 = rexp(n, rate = alphaCR.Y),
                       time.cens = rexp(n, rate = alpha.cens),
                       treatment = "2",
                       stringsAsFactors = FALSE))
df2$time <- pmin(df2$time1, df2$time2, df2$time.cens) ## first status
df2$status <- ifelse(df2$time == df2$time1, 1, ifelse(df2$time == df2$time2, 2, 0)) ## type of event
df2$strata <- sample(c('a', 'b', 'c'), 2*n, replace = T)
df2$toxicity <- sample(0:1, 2*n, replace = T)
df2$strata2 <- sample(c('d', 'e', 'f'), 2*n, replace = T)

test_that("BuyseTest package and Eva's R code give the same results with one endpoint and one stratum", {
  
  ## Net benefit computed with Eva's R code (see inst/Code/reproduce-results-CR.R)
  delta.R <- 0.04377023

  ## Apply GPC with BuyseTest package
  BT <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5), data = df2)
  
  ## Test
  expect_equal(as.double(delta.R), as.double(coef(BT, statistic = "netBenefit")), tol = 1e-5)
  ## BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5),
  ##           data = df2, method.inference = "u-statistic")
  
})

test_that("New package version gives the same results as previous one", {
  
  #### Net benefit computed with previous version 
  delta11 <- 0.04377023 # one outcome, one stratum
  Delta11 <- 0.04377023 # one outcome, one stratum
  delta13 <- c(-0.06967276, 0.13925664, 0.13200801) # one outcome, 3 strata
  Delta13 <- 0.0477878  # one outcome, 3 strata
  delta21 <- c(0.04377023, 0.02028528) # 2 outcomes, one stratum
  Delta21 <- c(0.04377023, 0.06405551) # 2 outcomes, one stratum
  delta23 <- matrix(data = c(-0.06967276, 0.045860858, 0.13925664, 0.009916651, 0.13200801, 0.015133257), byrow = T, nrow = 3) # 2 outcomes, 3 strata
  Delta23 <- c(0.04778780, 0.07473073) # 2 outcomes, 3 strata
  
  #### Apply GPC with new version of BuyseTest package
  ## One outcome, one stratum
  BT11.D <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5), data = df2)

  ## One outcome, 3 strata
  BT13.D <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + strata, data = df2)

  ## Two outcomes, one stratum
  BT21.D <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + bin(toxicity), data = df2)

  ## Two outcomes, 3 strata
  BT23.D <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + bin(toxicity) + strata, data = df2)

  #### Tests
  expect_equal(Delta11, as.double(coef(BT11.D, statistic = "netBenefit")), tol = 1e-5)
  ## expect_equal(Delta13, as.double(coef(BT13.D, statistic = "netBenefit")), tol = 1e-5)
  ## expect_equal(Delta21, as.double(coef(BT21.D, statistic = "netBenefit")), tol = 1e-5)
  ## expect_equal(Delta23, as.double(coef(BT23.D, statistic = "netBenefit")), tol = 1e-5)  
})

test_that("Package give the same results when model.tte is (not) provided as an argument", {
  
  ## Create prodlim object to be inserted as an argument
  fit = prodlim::prodlim(Hist(time, status) ~ treatment + strata, data = df2)
  #fit = prodlim(Hist(time, status) ~ treatment + strata + strata2, data = df2)
  
  ## Net benefit without passing model.tte
  B = BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + strata, data = df2) 
  #B = BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + strata + strata2, data = df2) 
  
  ## Net benefit with model.tte passed as an argument
  B.model = BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + strata, data = df2, model.tte = fit)
  #B.model = BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5) + strata + strata2, data = df2, model.tte = fit) 
  
  ## Tests
  expect_equal(coef(B), coef(B.model))
  
})

test_that("When TTE endpoints are analyzed several times with different thresholds, the results are the same than those when analyzed once with lowest threshold", {

  ## Endpoint analyzed once with threshold = 0.5
  B1 = BuyseTest(treatment ~ tte(time, status = status, threshold = 0.5), data = df2) 
  
  ## Endpoint analyzed twice with threshold = c(1, 0.75, 0.5)
  B2 = BuyseTest(treatment ~ tte(time, status = status, threshold = 1) + tte(time, status = status, threshold = 0.75) + tte(time, status = status, threshold = 0.5),
                 data = df2) 
  
  ## Tests
  expect_equal(as.double(coef(B1)), as.double(coef(B2)[3]))
  
})

## ** Simulate HR with proportional subdistribution hazard ratio
if(FALSE){ ## works but time consuming
test_that("The relationship between net benefit and subdistribution hazard ratio is verified", {
  
  set.seed(10)
  
  ## Compute true net benefit from true subdistribution hazard ratio
  b.1C <- v * log(1 - p.1C)
  b.1T <- b.1C * true.sHR
  p.1T <- 1 - exp(b.1T / v)
  
  true.Delta <- ((1 - true.sHR) / (1 + true.sHR)) * (1 - (1 - p.1T) * (1 - p.1C))
  
  ## Simulate big datasets with pre-specified subdistribution hazard ratio (based on Jeong and Fine, 2006)
  for (q in 1:length(true.sHR)) { ## q <- 1
    
    dd <- simCompetingRisks(n.T = n.big/2, n.C = n.big/2, p.1C = p.1C, v.1C = v,
                            v.1T = v, v.2C = v, v.2T = v, sHR = true.sHR[q])
    
    ## Compute net benefit with BuyseTest package
    B <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0),
                   data = dd,
                   scoring.rule = "Gehan",
                   keep.pairScore = FALSE,
                   method.inference = "none",
                   trace = FALSE) 
    
    ## Tests
    expect_equal(true.Delta[q], as.double(coef(B, statistic = "netBenefit")), tolerance = 1e-2) # tolerance because of the too small dataset
    
  }
})
}
######################################################################
### test-BuyseTest-CR.R ends here
