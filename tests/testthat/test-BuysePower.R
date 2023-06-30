### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: jun 30 2023 (13:06) 
##           By: Brice Ozenne
##     Update #: 51
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

context("Check BuysePower \n")

BuyseTest.options(order.Hprojection = 1,
                  pool.strata = "Buyse")

## * 1 binary endpoint
if(FALSE){ ## to save computation time for CRAN check
    test_that("1 binary endpoint", {
        seqN <- c(10,20,30,40,50)
        nrep <- 5
        formula <- treatment ~ bin(toxicity)

        ## automatic
        e.bin <- powerBuyseTest(sim = simBuyseTest,
                                sample.size = seqN,
                                n.rep = nrep,
                                formula = formula,
                                method.inference = "u-statistic", trace = 0,
                                seed = 5)
        vec.seed <- e.bin@seed

        ## manual
        GS <- NULL
        for(iRep in 1:nrep){
            set.seed(vec.seed[iRep])
            d <- simBuyseTest(max(seqN), name.cluster = NULL)
            d[, id := 1:.N, by = "treatment"]
            iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
                data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", trace = 0)))
            })
            GS <- rbind(GS,do.call(rbind,iLs))
        }

        GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
        test <- model.tables(e.bin)[, .SD,.SDcols = names(GS.S)]
        expect_equivalent(unlist(GS.S),unlist(test), tol = 1e-6)

        GS.bis <- data.frame("n.T" = c(10, 20, 30, 40, 50), 
                             "mean.estimate" = c(0.06, 0.04, 0.06666667, 0.06, 0.016), 
                             "sd.estimate" = c(0.23021729, 0.19811613, 0.0781736, 0.10547512, 0.11081516), 
                             "mean.se" = c(0.21583773, 0.1537722, 0.12812014, 0.11087229, 0.09903152), 
                             "rejection.rate" = c(0, 0, 0, 0.2, 0.2))
        expect_equivalent(GS.bis, as.data.frame(test), tol = 1e-6)
    })
}

## * 1 tte endpoint
## ** Gehan
test_that("1 tte endpoint - Gehan", {    
    seqN <- c(10,30,50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status)

    ## automatic
    e.tte <- powerBuyseTest(sim = simBuyseTest,
                            sample.size = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Gehan",
                            seed = 5)
    vec.seed <- e.tte@seed

    ## manual
    GS <- NULL
    for(iRep in 1:nrep){
        set.seed(vec.seed[iRep])
        d <- simBuyseTest(max(seqN), name.cluster = NULL)
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
            data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Gehan", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    test <- model.tables(e.tte)[,.SD,.SDcols= names(GS.S)]
    expect_equivalent(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(10, 30, 50), 
                         "mean.estimate" = c(-0.05, -0.00822222, -0.01248), 
                         "sd.estimate" = c(0.20639767, 0.18556986, 0.16888822), 
                         "mean.se" = c(0.21809705, 0.12183075, 0.0950564), 
                         "rejection.rate" = c(0, 0.2, 0.2))
    expect_equivalent(GS.bis, as.data.frame(test), tol = 1e-6)
})

## ** Peron
test_that("1 tte endpoint - Peron", {    
    seqN <- c(50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status)

    ## automatic
    e.tte <- powerBuyseTest(sim = simBuyseTest,
                            sample.size = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Peron",
                            seed = 5)
    vec.seed <- e.tte@seed

    ## manual
    GS <- NULL
    for(iRep in 1:nrep){
        set.seed(vec.seed[iRep])
        d <- simBuyseTest(max(seqN), name.cluster = NULL)
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
            data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    test <- model.tables(e.tte)[,.SD,.SDcols = names(GS.S)]
    expect_equivalent(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(50), 
                         "mean.estimate" = c(-0.02962616), 
                         "sd.estimate" = c(0.20130054), 
                         "mean.se" = c(0.12733956), 
                         "rejection.rate" = c(0.2))
    expect_equivalent(GS.bis, as.data.frame(test), tol = 1e-6)
})



## * Multiple endpoints
test_that("Multiple endpoints", {    
    seqN <- c(10,50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status, threshold = 0.25) + bin(toxicity) + tte(eventtime, status = status, threshold = 0)

    ## automatic
    e.tte <- powerBuyseTest(sim = simBuyseTest,
                            sample.size = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Peron",
                            seed = 5)
    vec.seed <- e.tte@seed

    ## manual
    GS <- NULL
    for(iRep in 1:nrep){
        set.seed(vec.seed[iRep])
        d <- simBuyseTest(max(seqN), name.cluster = NULL)
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 1
            iBT <- BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
            iCI <- confint(iBT)
            data.table(n.T = seqN[iN], endpoint = rownames(iCI), iCI)
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[endpoint == "eventtime", .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    
    test <- model.tables(e.tte)[,.SD,.SDcols = names(GS.S)]
    expect_equivalent(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(10, 50), 
                         "mean.estimate" = c(-0.01050476, -0.04230385), 
                         "sd.estimate" = c(0.18982451, 0.19697267), 
                         "mean.se" = c(0.28263784, 0.12597809), 
                         "rejection.rate" = c(0, 0.2))
    expect_equivalent(GS.bis, as.data.frame(test), tol = 1e-3)
})




######################################################################
### test-BuysePower.R ends here
