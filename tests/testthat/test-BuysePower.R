### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: Dec 21 2021 (17:43) 
##           By: Brice Ozenne
##     Update #: 43
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

BuyseTest.options(order.Hprojection = 1)

## * 1 binary endpoint
if(FALSE){ ## to save computation time for CRAN check
    test_that("1 binary endpoint", {
        seqN <- c(10,20,30,40,50)
        nrep <- 5
        formula <- treatment ~ bin(toxicity)

        ## automatic
        e.bin <- powerBuyseTest(sim = simBuyseTest,
                                sample.sizeT = seqN,
                                sample.sizeC = seqN,
                                n.rep = nrep,
                                formula = formula,
                                method.inference = "u-statistic", trace = 0,
                                seed = 10)

        ## manual
        set.seed(10)
        GS <- NULL
        for(iRep in 1:nrep){
            d <- simBuyseTest(max(seqN))
            d[, id := 1:.N, by = "treatment"]
            iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
                data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", trace = 0)))
            })
            GS <- rbind(GS,do.call(rbind,iLs))
        }

        GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
        test <- summary(e.bin, print = FALSE)[, .SD,.SDcols = names(GS.S)]
        expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

        GS.bis <- data.frame("n.T" = c(10, 20, 30, 40, 50), 
                             "mean.estimate" = c(0.1, 0.05, 0.02666667, 0.045, 0.024), 
                             "sd.estimate" = c(0.31622777, 0.24748737, 0.22656861, 0.20186629, 0.16816658), 
                             "mean.se" = c(0.2077989, 0.15264937, 0.12585806, 0.10933957, 0.09831528), 
                             "rejection.rate" = c(0.2, 0.2, 0.4, 0.4, 0.4))
        expect_equal(GS.bis, as.data.frame(test), tol = 1e-6)
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
                            sample.sizeT = seqN,
                            sample.sizeC = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Gehan",
                            seed = 10)

    ## manual
    set.seed(10)
    GS <- NULL
    for(iRep in 1:nrep){
        d <- simBuyseTest(max(seqN))
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
            data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Gehan", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    test <- summary(e.tte, print = FALSE)[,.SD,.SDcols= names(GS.S)]
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(10, 30, 50), 
                         "mean.estimate" = c(0.016, 0.03066667, -0.01008), 
                         "sd.estimate" = c(0.22266567, 0.0879134, 0.1003021), 
                         "mean.se" = c(0.20307347, 0.12681144, 0.09572955), 
                         "rejection.rate" = c(0, 0, 0))
    expect_equal(GS.bis, as.data.frame(test), tol = 1e-6)
})

## ** Peron
test_that("1 tte endpoint - Peron", {    
    seqN <- c(50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status)

    ## automatic
    e.tte <- powerBuyseTest(sim = simBuyseTest,
                            sample.sizeT = seqN,
                            sample.sizeC = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Peron",
                            seed = 10)

    ## manual
    set.seed(10)
    GS <- NULL
    for(iRep in 1:nrep){
        d <- simBuyseTest(max(seqN))
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 2
            data.table(n.T = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    test <- summary(e.tte, print = FALSE)[,.SD,.SDcols = names(GS.S)]
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(50), 
                         "mean.estimate" = c(-0.02685844), 
                         "sd.estimate" = c(0.12158465), 
                         "mean.se" = c(0.12782063), 
                         "rejection.rate" = c(0))
    expect_equal(GS.bis, as.data.frame(test), tol = 1e-6)
})



## * Multiple endpoints
test_that("Multiple endpoints", {    
    seqN <- c(10,50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status, threshold = 0.25) + bin(toxicity) + tte(eventtime, status = status, threshold = 0)

    ## automatic
    e.tte <- powerBuyseTest(sim = simBuyseTest,
                            sample.sizeT = seqN,
                            sample.sizeC = seqN,
                            n.rep = nrep,
                            formula = formula,
                            method.inference = "u-statistic", trace = 0,
                            scoring.rule = "Peron",
                            seed = 10)

    ## manual
    set.seed(10)
    GS <- NULL
    for(iRep in 1:nrep){
        d <- simBuyseTest(max(seqN))
        d[, id := 1:.N, by = "treatment"]
        iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 1
            iBT <- BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
            iCI <- confint(iBT)
            data.table(n.T = seqN[iN], endpoint = rownames(iCI), iCI)
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[endpoint == "eventtime", .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    
    test <- summary(e.tte, print = FALSE)[,.SD,.SDcols = names(GS.S)]
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(10, 50), 
                         "mean.estimate" = c(0.00469841, -0.01343508), 
                         "sd.estimate" = c(0.1225553, 0.10901633), 
                         "mean.se" = c(0.3019124, 0.12739454), 
                         "rejection.rate" = c(0, 0))
    expect_equal(GS.bis, as.data.frame(test), tol = 1e-3)
})




######################################################################
### test-BuysePower.R ends here
