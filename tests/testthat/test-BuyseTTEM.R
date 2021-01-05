### test-BuyseTTEM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  2 2019 (11:54) 
## Version: 
## Last-Updated: jan  5 2021 (11:33) 
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
    library(data.table)
}
context("Check BuyseTTEM against predictCoxPL")

library(riskRegression)
library(prodlim)
library(survival)

if(FALSE){
## * Survival case
## treatment arm with last observation NA
## control arm with last observation event
set.seed(10)
dt <- simBuyseTest(100)
dt[,last:=eventtime==max(eventtime), by = c("treatment","toxicity")]
dt[treatment == "C" & last, status:=1]
dt[treatment == "T" & last, status:=0]
dt2 <- dt[status!=0]
## plot(prodlim(Hist(eventtime,status)~treatment+toxicity,data=dt))

seqTau <- c(0,unique(dt$eventtime),runif(100,0,max(dt$eventtime)),1e5)

## ** no strata (other than treatment group)
test_that("no strata, survival", {
    ## ## no censoring ## ##
    e.rr <- coxph(Surv(eventtime,status)~strata(treatment), data = dt2, x = TRUE, y = TRUE)
    e.bb <- BuyseTTEM(Hist(eventtime,status)~treatment, data = dt2, iid = TRUE, iid.surv = "exp") 

    ## Treatment
    test <- predict(e.bb, time = seqTau, treatment  = "T", iid = TRUE)
    GS <- predictCoxPL(e.rr, time = seqTau, newdata  = dt[treatment=="T",.SD[1]], iid = TRUE)

    expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
    expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)

    ## Control
    test <- predict(e.bb, time = seqTau, treatment  = "C", iid = TRUE)
    GS <- predictCoxPL(e.rr, time = seqTau, newdata  = dt[treatment=="C",.SD[1]], iid = TRUE)

    expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
    expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)

    ## ## censoring ## ##
    e.r <- coxph(Surv(eventtime,status)~strata(treatment), data = dt, x = TRUE, y = TRUE)
    e.b <- BuyseTTEM(Hist(eventtime,status)~treatment, data = dt, iid = TRUE, iid.surv = "exp")

    ## Treamtent
    test <- predict(e.b, time = seqTau, treatment  = "T", iid = TRUE)
    GS <- predictCoxPL(e.r, time = seqTau, newdata  = dt[treatment=="T",.SD[1]], iid = TRUE)

    expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
    expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)

    ## Control
    test <- predict(e.b, time = seqTau, treatment  = "C", iid = TRUE)
    GS <- predictCoxPL(e.r, time = seqTau, newdata  = dt[treatment=="C",.SD[1]], iid = TRUE)

    expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
    expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)
})

## ** with strata (other than treatment group)

test_that("strata, survival", {
    ## ## no censoring ## ##
    e.rr <- coxph(Surv(eventtime,status)~strata(treatment)+strata(toxicity), data = dt2, x = TRUE, y = TRUE)
    e.bb <- BuyseTTEM(Hist(eventtime,status)~treatment+toxicity, data = dt2, iid = TRUE, iid.surv = "exp",
                     treatment = "treatment")

    for(iStrata in c("yes","no")){ ## iStrata <- "yes"
        ## Treatment
        test <- predict(e.bb, time = seqTau, treatment  = "T", strata = iStrata, iid = TRUE)
        GS <- predictCoxPL(e.rr, time = seqTau, newdata  = dt[treatment=="T" & toxicity == iStrata,.SD[1]], iid = TRUE)

        expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
        expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)

        ## Control
        test <- predict(e.bb, time = seqTau, treatment  = "C", strata = iStrata, iid = TRUE)
        GS <- predictCoxPL(e.rr, time = seqTau, newdata  = dt[treatment=="C" & toxicity == iStrata,.SD[1]], iid = TRUE)

        expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
        expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)
    }

    ## ## censoring ## ##
    e.r <- coxph(Surv(eventtime,status)~strata(treatment)+strata(toxicity), data = dt, x = TRUE, y = TRUE)
    e.b <- BuyseTTEM(Hist(eventtime,status)~treatment+toxicity, data = dt, iid = TRUE, iid.surv = "exp",
                     treatment = "treatment")

    for(iStrata in c("yes","no")){
        ## Treatment
        test <- predict(e.b, time = seqTau, treatment  = "T", strata = iStrata, iid = TRUE)
        GS <- predictCoxPL(e.r, time = seqTau, newdata  = dt[treatment=="T" & toxicity == iStrata,.SD[1]], iid = TRUE)

        expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
        expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)

        ## Control
        test <- predict(e.b, time = seqTau, treatment  = "C", strata = iStrata, iid = TRUE)
        GS <- predictCoxPL(e.r, time = seqTau, newdata  = dt[treatment=="C" & toxicity == iStrata,.SD[1]], iid = TRUE)

        expect_equal(test$survival,GS$survival[1,], tol = 1e-6)
        expect_equal(test$survival.iid,GS$survival.iid[,,1], tol = 1e-6)
    }
})

## ** other
dataT <- data.table(time = 1:5,
                    treatment = "T",
                    status1 = c(1,0,1,1,1),
                    status2 = c(1,0,1,1,1),
                    status3 = c(1,1,1,1,1))
dataC <- data.table(time = c(1:5-0.1,5,5),
                    treatment = "C",
                    status1 = c(1,1,0,1,0,0,0),
                    status2 = c(1,1,0,1,0,1,1),
                    status3 = c(1,1,1,1,1,1,1))
data <- rbind(dataC, dataT)
seqThreshold <- c(1e-12,0.5,1)

e.TTEM <- BuyseTTEM(Hist(time, status1) ~ treatment, data = data, iid = FALSE, iid.surv = "prodlim")
expect_equal(predict(e.TTEM, treatment = "C", time = c(1,3,4,5)+1e-12)$survival,
             c(0.8571429, 0.7142857, 0.5357143, 0.5357143), tol = 1e-6)
             

## * Competing risk case
set.seed(10)
dt <- simBuyseTest(100, argsTTE = list(CR=TRUE))
setkeyv(dt, "eventtime")
dt[,last:=eventtime==max(eventtime), by = c("treatment","toxicity")]
dt[treatment == "C" & last, status:=1]
dt[treatment == "T" & last, status:=0]
## plot(prodlim(Hist(eventtime,status)~treatment+toxicity,data=dt))
dt2 <- dt[status!=0]

seqTau <- c(0,unique(dt$eventtime),runif(100,0,max(dt$eventtime)),1e5)

## ** no strata (other than treatment group)

test_that("no strata, competing risks", {
    ## ## no censoring ## ##
    e.rr <- CSC(Hist(eventtime,status)~strata(treatment), data = dt2)
    e.bb <- BuyseTTEM(Hist(eventtime,status)~treatment, data = dt2, iid = TRUE, iid.surv = "prodlim")

    for(iCause in 1:2){
        ## Treatment
        test <- predict(e.bb, time = seqTau, treatment  = "T", iid = TRUE, cause = iCause)
        GS <- predict(e.rr, time = seqTau, newdata  = dt[treatment=="T",.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

        expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
        expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)

        ## Control
        test <- predict(e.bb, time = seqTau, treatment  = "C", iid = TRUE, cause = iCause)
        GS <- predict(e.rr, time = seqTau, newdata  = dt[treatment=="C",.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

        expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
        expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)
    }
    
    ## ## censoring ## ##
    e.r <- CSC(Hist(eventtime,status)~strata(treatment), data = dt)
    e.b <- BuyseTTEM(Hist(eventtime,status)~treatment, data = dt, iid = TRUE, iid.surv = "prodlim")

    for(iCause in 1:2){
        ## Treatment
        test <- predict(e.b, time = seqTau, treatment  = "T", iid = TRUE, cause = iCause)
        GS <- predict(e.r, time = seqTau, newdata  = dt[treatment=="T",.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

        expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
        expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)

        ## Control
        test <- predict(e.b, time = seqTau, treatment  = "C", iid = TRUE, cause = iCause)
        GS <- predict(e.r, time = seqTau, newdata  = dt[treatment=="C",.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

        expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
        expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)
    }
})

## ** with strata (other than treatment group)

test_that("with strata, competing risks", {
    ## ## no censoring ## ##
    e.rr <- CSC(Hist(eventtime,status)~strata(treatment)+strata(toxicity), data = dt2)
    e.bb <- BuyseTTEM(Hist(eventtime,status)~treatment+toxicity, data = dt2, iid = TRUE, iid.surv = "prodlim", treatment = "treatment")

    for(iCause in 1:2){
        for(iStrata in c("yes","no")){
            ## Treatment
            test <- predict(e.bb, time = seqTau, treatment  = "T", strata = iStrata, iid = TRUE, cause = iCause)
            GS <- predict(e.rr, time = seqTau, newdata  = dt[treatment=="T" & toxicity == iStrata,.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

            expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
            expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)

            ## Control
            test <- predict(e.bb, time = seqTau, treatment  = "C", strata = iStrata, iid = TRUE, cause = iCause)
            GS <- predict(e.rr, time = seqTau, newdata  = dt[treatment=="C" & toxicity == iStrata,.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

            expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
            expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)
        }
    }
    
    ## ## censoring ## ##
    e.r <- CSC(Hist(eventtime,status)~strata(treatment)+strata(toxicity), data = dt)
    e.b <- BuyseTTEM(Hist(eventtime,status)~treatment+toxicity, data = dt, iid = TRUE, iid.surv = "prodlim", treatment = "treatment")

    for(iCause in 1:2){
        for(iStrata in c("yes","no")){
            ## Treatment
            test <- predict(e.b, time = seqTau, treatment  = "T", strata = iStrata, iid = TRUE, cause = iCause)
            GS <- predict(e.r, time = seqTau, newdata  = dt[treatment=="T" & toxicity == iStrata,.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

        expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
        expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)

        ## Control
        test <- predict(e.b, time = seqTau, treatment  = "C", strata = iStrata, iid = TRUE, cause = iCause)
        GS <- predict(e.r, time = seqTau, newdata  = dt[treatment=="C" & toxicity == iStrata,.SD[1]], iid = TRUE, prodlim = TRUE, cause = iCause)

            expect_equal(test$cif,GS$absRisk[1,], tol = 1e-6)
            expect_equal(test$cif.iid,GS$absRisk.iid[,,1], tol = 1e-6)
        }
    }
})
}
######################################################################
### test-BuyseTTEM.R ends here
