### test-BuysePower.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2019 (18:24) 
## Version: 
## Last-Updated: mar 23 2020 (13:45) 
##           By: Brice Ozenne
##     Update #: 18
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
            data.table(n = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection (TRUE)" =  mean(p.value <= 0.05)), by = "n"]
    test <- summary(e.bin, print = FALSE)[["netBenefit"]][, .SD,.SDcols = c("n.T","mean.estimate", "sd.estimate", "mean.se", "rejection (TRUE)")]
    setnames(test, old = "n.T", new = "n")
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- c("n1" = 10, "n2" = 20, "n3" = 30, "n4" = 40, "n5" = 50, "mean.estimate1" = -0.1, "mean.estimate2" = -0.07, "mean.estimate3" = -0.04, "mean.estimate4" = 0.005, "mean.estimate5" = 0.004, "sd.estimate1" = 0.3082207, "sd.estimate2" = 0.20493902, "sd.estimate3" = 0.13207742, "sd.estimate4" = 0.09905806, "sd.estimate5" = 0.08876936, "mean.se1" = 0.21097018, "mean.se2" = 0.15363469, "mean.se3" = 0.12747414, "mean.se4" = 0.11086051, "mean.se5" = 0.09920359, "rejection (TRUE)1" = 0, "rejection (TRUE)2" = 0.2, "rejection (TRUE)3" = 0, "rejection (TRUE)4" = 0, "rejection (TRUE)5" = 0)
    expect_equal(GS.bis, unlist(test), tol = 1e-6)
})

## * 1 tte endpoint
## ** Gehan
test_that("1 tte endpoint - Gehan", {    
    seqN <- c(10,20,30,40,50)
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
            data.table(n = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Gehan", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection (TRUE)" =  mean(p.value <= 0.05)), by = "n"]
    test <- summary(e.tte, print = FALSE)[["netBenefit"]][, .SD,.SDcols = c("n.T","mean.estimate", "sd.estimate", "mean.se", "rejection (TRUE)")]
    setnames(test, old = "n.T", new = "n")
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- c("n1" = 10, "n2" = 20, "n3" = 30, "n4" = 40, "n5" = 50, "mean.estimate1" = 0.13, "mean.estimate2" = 0.057, "mean.estimate3" = 0.08177778, "mean.estimate4" = 0.079375, "mean.estimate5" = 0.04248, "sd.estimate1" = 0.31819805, "sd.estimate2" = 0.06836026, "sd.estimate3" = 0.07265936, "sd.estimate4" = 0.07469731, "sd.estimate5" = 0.09814088, "mean.se1" = 0.18568528, "mean.se2" = 0.14375078, "mean.se3" = 0.11933458, "mean.se4" = 0.10166738, "mean.se5" = 0.09174785, "rejection (TRUE)1" = 0.2, "rejection (TRUE)2" = 0, "rejection (TRUE)3" = 0, "rejection (TRUE)4" = 0.2, "rejection (TRUE)5" = 0.2)
    expect_equal(GS.bis, unlist(test), tol = 1e-6)
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
            data.table(n = seqN[iN], confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)))
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection (TRUE)" =  mean(p.value <= 0.05)), by = "n"]
    test <- summary(e.tte, print = FALSE)[["netBenefit"]][, .SD,.SDcols = c("n.T","mean.estimate", "sd.estimate", "mean.se", "rejection (TRUE)")]
    setnames(test, old = "n.T", new = "n")
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- c("n" = 50, "mean.estimate" = 0.04281066, "sd.estimate" = 0.12750526, "mean.se" = 0.12883854, "rejection (TRUE)" = 0)
    expect_equal(GS.bis, unlist(test), tol = 1e-6)
})



## * Multiple endpoints
test_that("Multiple endpoints", {    
    seqN <- c(50)
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
            iConfint <- confint(BuyseTest(formula, data = d[id <= seqN[iN]], method.inference = "u-statistic", scoring.rule = "Peron", trace = 0))
            data.table(n = seqN[iN], iConfint[NROW(iConfint),,drop=FALSE])
        })
        GS <- rbind(GS,do.call(rbind,iLs))
    }

    GS.S <- GS[, .(mean.estimate = mean(estimate), sd.estimate = sd(estimate), mean.se =  mean(se), "rejection (TRUE)" =  mean(p.value <= 0.05)), by = "n"]
    test <- summary(e.tte, print = FALSE)[["netBenefit"]][, .SD,.SDcols = c("n.T","mean.estimate", "sd.estimate", "mean.se", "rejection (TRUE)")]
    setnames(test, old = "n.T", new = "n")
    expect_equal(unlist(GS.S),unlist(test), tol = 1e-6)

    GS.bis <- c("n" = 50, "mean.estimate" = 0.02538207, "sd.estimate" = 0.13532166, "mean.se" = 0.12851884, "rejection (TRUE)" = 0)
    expect_equal(GS.bis, unlist(test), tol = 1e-4)
})




######################################################################
### test-BuysePower.R ends here
