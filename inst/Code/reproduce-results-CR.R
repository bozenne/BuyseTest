### reproduce-results-CR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne, Eva Cantagallo
## Created: mar 28 2019 (16:01) 
## Version: 
## Last-Updated: 
##           By: 
##     Update #: 
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

setwd("D:\\Documents\\EORTC\\BENEFIT\\IT\\BuyseTest\\inst\\Code")
source("GPCV2.R")

## Simulate data
alphaE.X <- 2
alphaCR.X <- 1
alphaE.Y <- 3
alphaCR.Y <- 2
alpha.cens <- 1.5
n <- 1e2

## Simulate data
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X), time2 = rexp(n, rate = alphaCR.X), time.cens = rexp(n, rate = alpha.cens), treatment = "1"),
            data.frame(time1 = rexp(n, rate = alphaE.Y), time2 = rexp(n, rate = alphaCR.Y), time.cens = rexp(n, rate = alpha.cens), treatment = "2"))
df$time <- pmin(df$time1, df$time2, df$time.cens) ## first status
df$status <- ifelse(df$time == df$time1, 1, ifelse(df$time == df$time2, 2, 0)) ## type of event
df$strata <- sample(c('a', 'b', 'c'), 2*n, replace = T)
df$toxicity <- sample(0:1, 2*n, replace = T)

###############################  Apply GPC with Eva's code  ###############################
delta = GPCWithCR("time", "treatment", "status", threshold = 0.5, df, print = F)$delta.netChance


####  BuyseTest package version 1.7.4  ###############################
## One outcome, one stratum
BT11.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df)@delta.netBenefit
BT11.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df)@Delta.netBenefit

## One outcome, 3 strata
BT13.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata, data = df)@delta.netBenefit
BT13.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata, data = df)@Delta.netBenefit

## Two outcomes, one stratum
BT21.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity), data = df)@delta.netBenefit
BT21.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity), data = df)@Delta.netBenefit

## Two outcomes, one stratum
BT23.d = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity) + strata, data = df)@delta.netBenefit
BT23.D = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + bin(toxicity) + strata, data = df)@Delta.netBenefit

## Tests
expect_equal(as.double(BT11@Delta.netBenefit), as.double(BT@Delta.netBenefit))
