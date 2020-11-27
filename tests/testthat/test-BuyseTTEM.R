### test-BuyseTTEM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  2 2019 (11:54) 
## Version: 
## Last-Updated: nov 26 2020 (20:31) 
##           By: Brice Ozenne
##     Update #: 21
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
context("Check BuyseTTEM against predictCoxPL")

library(riskRegression)
library(prodlim)
library(survival)

## * Survival case
## ** no strata (other than treatment group)
## treatment arm with last observation NA
## control arm with last observation event
set.seed(10)
dt <- simBuyseTest(100)

e.BuyseTTEM <- BuyseTTEM(Hist(eventtime,status)~treatment, data = dt, iid = TRUE)
## ** with strata (other than treatment group)
## treatment arm with last observation NA
## control arm with last observation event

e.tempo <- coxph(Surv(eventtime,status)~1, data = dt, x = TRUE, y = TRUE)
iid.tempo <- predictCox(e.tempo,newdata = dt, times = 1,iid = TRUE)$survival.iid
test.riskRegression.version <- all(abs(apply(iid.tempo,2:3,mean, na.omit =TRUE))<1e-6) ## to be removed


test_that("iid with 1 strata variable", {
    e.prodlim <- prodlim(Hist(eventtime,status)~treatment, data = dt)
    e.survival <- coxph(Surv(eventtime,status)~strata(treatment), data = dt, x = TRUE, y = TRUE)

    iid.BT <- lava::iid(e.prodlim)
    iid.RR <- iidCox(e.survival, return.object = FALSE)

    iid.RR$IFsurvival[[1]]


    if(test.riskRegression.version){
    iid.RR$IFsurvival <- list(predictCox(e.survival, newdata = dt[dt$treatment == "C",][1], times = iid.RR$time[[1]], iid = TRUE)$survival.iid[,,1],
                              predictCox(e.survival, newdata = dt[dt$treatment == "T",][1], times = iid.RR$time[[2]], iid = TRUE)$survival.iid[,,1])

    expect_equal(unname(iid.BT$IFhazard[[1]]), unname(iid.RR$IFhazard[[1]]))
    expect_equal(unname(iid.BT$IFhazard[[2]]), unname(iid.RR$IFhazard[[2]]))

    expect_equal(unname(iid.BT$IFcumhazard[[1]]), unname(iid.RR$IFcumhazard[[1]]))
    expect_equal(unname(iid.BT$IFcumhazard[[2]]), unname(iid.RR$IFcumhazard[[2]]))

    expect_equal(unname(iid.BT$IFsurvival[[1]]), unname(iid.RR$IFsurvival[[1]]))
    expect_equal(unname(iid.BT$IFsurvival[[2]]), unname(iid.RR$IFsurvival[[2]]))
    }
})
## profvis::profvis(lava::iid(e.prodlim))

test_that("iid with 2 strata variables", {
    e.prodlim <- prodlim(Hist(eventtime,status)~treatment+toxicity, data = dt)
    e.survival <- coxph(Surv(eventtime,status)~strata(treatment)+strata(toxicity), data = dt, x = TRUE, y = TRUE)

    iid.BT <- lava::iid(e.prodlim)
    iid.RR <- iidCox(e.survival, return.object = FALSE)
    if(test.riskRegression.version){
        iid.RR$IFsurvival <- list(t(predictCox(e.survival, newdata = dt[dt$treatment == "C" & toxicity == 0,][1], times = iid.RR$time[[1]], iid = TRUE)$survival.iid[,,1]),
                                  t(predictCox(e.survival, newdata = dt[dt$treatment == "T" & toxicity == 0,][1], times = iid.RR$time[[2]], iid = TRUE)$survival.iid[,,1]))

        expect_equal(unname(iid.BT$IFhazard[[1]]), unname(iid.RR$IFhazard[[1]]))
        expect_equal(unname(iid.BT$IFhazard[[3]]), unname(iid.RR$IFhazard[[2]]))
        expect_equal(unname(iid.BT$IFhazard[[2]]), unname(iid.RR$IFhazard[[3]]))
        expect_equal(unname(iid.BT$IFhazard[[4]]), unname(iid.RR$IFhazard[[4]]))

        expect_equal(unname(iid.BT$IFcumhazard[[1]]), unname(iid.RR$IFcumhazard[[1]]))
        expect_equal(unname(iid.BT$IFcumhazard[[3]]), unname(iid.RR$IFcumhazard[[2]]))
        expect_equal(unname(iid.BT$IFcumhazard[[2]]), unname(iid.RR$IFcumhazard[[3]]))
        expect_equal(unname(iid.BT$IFcumhazard[[4]]), unname(iid.RR$IFcumhazard[[4]]))
    }
})


## * Competing risk case
## ** no strata (other than treatment group)
## treatment arm with last observation NA
## control arm with last observation event

## ** with strata (other than treatment group)
## treatment arm with last observation NA
## control arm with last observation event

######################################################################
### test-BuyseTTEM.R ends here