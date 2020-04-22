### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: apr 22 2020 (17:17) 
##           By: Brice Ozenne
##     Update #: 35
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
                             "mean.estimate" = c(-0.1, -0.07, -0.04, 0.005, 0.004), 
                             "sd.estimate" = c(0.3082207, 0.20493902, 0.13207742, 0.09905806, 0.08876936), 
                             "mean.se" = c(0.21097018, 0.15363469, 0.12747414, 0.11086051, 0.09920359), 
                             "rejection.rate" = c(0, 0.2, 0, 0, 0))
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
                         "mean.estimate" = c(0.13, 0.08177778, 0.04248), 
                         "sd.estimate" = c(0.31819805, 0.07265936, 0.09814088), 
                         "mean.se" = c(0.18568528, 0.11933458, 0.09174785), 
                         "rejection.rate" = c(0.2, 0, 0.2))
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
                         "mean.estimate" = c(0.04281066), 
                         "sd.estimate" = c(0.12750526), 
                         "mean.se" = c(0.12987798), 
                         "rejection.rate" = c(0))
    expect_equal(GS.bis, as.data.frame(test), tol = 1e-6)
})



## * Multiple endpoints
test_that("Multiple endpoints", {    
    seqN <- c(10,50)
    nrep <- 5
    formula <- treatment ~ tte(eventtime, status = status, threshold = 0.25) + bind(toxicity) + tte(eventtime, status = status, threshold = 0)

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

    GS.S <- GS[endpoint == "eventtime_1e-12", .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    
    test <- summary(e.tte, print = FALSE)[,.SD,.SDcols = names(GS.S)]
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- data.frame("n.T" = c(10, 50), 
                         "mean.estimate" = c(-0.01278659, 0.02538207), 
                         "sd.estimate" = c(0.24334861, 0.13532166), 
                         "mean.se" = c(0.28827423, 0.12892174), 
                         "rejection.rate" = c(0, 0))
    expect_equal(GS.bis, as.data.frame(test), tol = 1e-3)

    ## seqN <- c(10,25,50)
    ## nrep <- 5
    ## formula <- treatment ~ tte(eventtime, status = status, threshold = 0.25) + bind(toxicity) + tte(eventtime, status = status, threshold = 0)

    ## ## automatic
    ## iCorrection <- 2
    
    ## e.tte <- powerBuyseTest(sim = simBuyseTest,
    ##                         sample.sizeT = seqN,
    ##                         sample.sizeC = seqN,
    ##                         n.rep = nrep,
    ##                         formula = formula,
    ##                         method.inference = "u-statistic", trace = 0,
    ##                         scoring.rule = "Gehan", correction.uninf = iCorrection,
    ##                         seed = 10)

    ## ## manual
    ## set.seed(10)
    ## GS <- NULL
    ## for(iRep in 1:nrep){
    ##     d <- simBuyseTest(max(seqN))
    ##     d[, id := 1:.N, by = "treatment"]
    ##     iLs <-  lapply(1:length(seqN), function(iN){ ## iN <- 1
    ##         iBT <- BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Gehan", trace = 0, correction.uninf = iCorrection)
    ##         iCI <- confint(iBT)
    ##         data.table(n.T = seqN[iN], n.C = seqN[iN], endpoint = rownames(iCI), iCI)
    ##     })
    ##     GS <- rbind(GS,do.call(rbind,iLs))
    ## }

    ## GS.S <- GS[endpoint == "eventtime_1e-12", .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection.rate" =  mean(p.value <= 0.05)), by = "n.T"]
    ## test <- summary(e.tte, print = FALSE)[,.SD,.SDcols = names(GS.S)]
    ## expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    
})




######################################################################
### test-BuysePower.R ends here
