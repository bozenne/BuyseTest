### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: sep 15 2023 (16:33) 
##           By: Brice Ozenne
##     Update #: 209
#----------------------------------------------------------------------
## 
### Commentary: Check 
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
context("Check resampling")

if(FALSE){
    library(BuyseTest)
    library(testthat)
    library(data.table)
}

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  keep.survival = FALSE,
                  order.Hprojection = 1,
                  add.1.presample = FALSE,
                  pool.strata = "Buyse",
                  trace = 0)
n.patients <- 100
method <- "Peron"

## * Simulate data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients,
                       n.C = n.patients,
                       argsBin = list(p.T = list(c(0.5,0.5),c(0.25,0.75))),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(scale.T = 1:3, scale.censoring.T = rep(1,3)),
                       n.strata = 3)

## * Permutation
test_that("permutation", {
    BT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1),
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "studentized permutation", n.resampling = 5)

    ## ** summary (two.sided)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
    ## outSummaryPerc
    ##       Generalized pairwise comparisons with 2 prioritized endpoints

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : permutation test with 20 samples 
    ##                     p-value computed using the studentized permutation distribution 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > right-censored pairs: probabilistic score based on the survival curves
    ## > neutral pairs   : ignored at lower priority endpoints
    ## > results
    ##   endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta   Delta p.value 
    ## eventtime1         1   100.00        20.87          21.91      57.22        0 -0.0104 -0.0104    0.85 
    ##  toxicity1       0.5    57.22        10.92          17.62      28.68        0 -0.0670 -0.0774    0.60 
    
    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"])) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"]))),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"])) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"]))))
    expect_equal(outSummaryStud$p.value, p.value)

    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"])),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"])))
    expect_equal(outSummaryPerc$p.value, p.value)
       
    ## ** summary (greater)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "greater", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "greater", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) <= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) <= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$p.value, p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] <= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] <= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$p.value, p.value)

    ## ** summary (less)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "less", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "less", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) >= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) >= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$p.value, p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] >= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] >= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$p.value, p.value)
    
    ## ** check permutation
    vec.seed <- BT.perm@seed
    for(iResample in 1:2){ ## iResample <- 1
        set.seed(vec.seed[iResample])
        dt.perm <- copy(dt.sim)

        dt.perm[, treatment := treatment[sample.int(.N, size = .N, replace = FALSE)] ]
        ## expect_equal(table(dt.perm$treatment), table(dt.sim$treatment))

        iBT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1),
                              data = dt.perm, scoring.rule = method,
                              method.inference = "u-statistic")

        expect_equal(as.double(iBT.perm@Delta[,"netBenefit"]),
                     as.double(BT.perm@DeltaResampling[iResample,,"netBenefit"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.perm@Delta[,"winRatio"]),
                     as.double(BT.perm@DeltaResampling[iResample,,"winRatio"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.perm@covariance),
                     as.double(BT.perm@covarianceResampling[iResample,,]),
                     tol = 1e-6)
    }

})


## * Stratified permutation
test_that("stratified permutation", {
    BT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "studentized permutation", n.resampling = 5)

    ## ** summary (two.sided)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
    ##       Generalized pairwise comparisons with 2 prioritized endpoints and 3 strata

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : permutation test with 5 samples 
    ##                     confidence intervals/p-values computed using the quantiles of the empirical distribution 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > right-censored pairs: probabilistic score based on the survival curves
    ## > neutral pairs   : ignored at lower priority endpoints
    ## > uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints
    ## > results
    ##   endpoint threshold strata total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta   Delta p.value 
    ## eventtime1         1 global   100.00        20.77          19.57      55.25     4.40  0.0120   0.012     0.6 
    ##                           0    33.56         9.38           7.61      16.58     0.00  0.0527                 
    ##                           1    34.74         4.50           4.58      24.41     1.25 -0.0023                 
    ##                           2    31.70         6.90           7.38      14.27     3.15 -0.0153                 
    ##  toxicity1       0.5 global    59.66        11.53          18.18      29.95     0.00 -0.0665 -0.0545     0.6 
    ##                           0    16.58         2.96           5.19       8.43     0.00 -0.0662                 
    ##                           1    25.66         4.64           8.63      12.39     0.00 -0.1151                 
    ##                           2    17.42         3.93           4.36       9.13     0.00 -0.0138                 

    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"])) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"]))),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"])) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"]))))
    expect_equal(outSummaryStud[outSummaryStud$strata == "global","p.value"], p.value)

    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"])),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"])))
    expect_equal(outSummaryPerc[outSummaryPerc$strata=="global","p.value"], p.value)

    ## ** summary (greater)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "greater", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "greater", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) <= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) <= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud[outSummaryStud$strata=="global","p.value"], p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] <= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] <= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc[outSummaryPerc$strata=="global","p.value"], p.value)

    ## ** summary (less)
    outSummaryPerc <- suppressWarnings(model.tables(BT.perm, alternative = "less", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(model.tables(BT.perm, alternative = "less", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) >= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) >= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud[outSummaryStud$strata=="global","p.value"], p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] >= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] >= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc[outSummaryPerc$strata=="global","p.value"], p.value)
    
    ## ** check permutation
    vec.seed <- BT.perm@seed
    for(iResample in 1:2){ ## iResample <- 1
        set.seed(vec.seed[iResample])
        dt.perm <- copy(dt.sim)
        dt.perm[, treatment := treatment[sample.int(.N, size = .N, replace = FALSE)]]

        iBT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                              data = dt.perm, scoring.rule = method,
                              method.inference = "u-statistic")

        expect_equal(as.double(iBT.perm@Delta[,"netBenefit"]),
                     as.double(BT.perm@DeltaResampling[iResample,,"netBenefit"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.perm@Delta[,"winRatio"]),
                     as.double(BT.perm@DeltaResampling[iResample,,"winRatio"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.perm@covariance),
                     as.double(BT.perm@covarianceResampling[iResample,,]),
                     tol = 1e-6)
    }

    for(strataVar in c("strata","toxicity1")){
        BT.perm2 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                              data = dt.sim, scoring.rule = method, seed = 10, 
                              method.inference = "studentized permutation", strata.resampling = strataVar, n.resampling = 5)

        for(iResample in 1:2){ ## iResample <- 1 
            set.seed(vec.seed[iResample])
            dt.perm2 <- copy(dt.sim)
            setkeyv(dt.perm2, cols = strataVar)
            dt.perm2[, treatment := treatment[sample.int(.N, size = .N, replace = FALSE)], by = strataVar]

            iBT.perm2 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                                   data = dt.perm2, scoring.rule = method,
                                   method.inference = "u-statistic")

            expect_equal(as.double(iBT.perm2@Delta[,"netBenefit"]),
                         as.double(BT.perm2@DeltaResampling[iResample,,"netBenefit"]),
                         tol = 1e-6)
            expect_equal(as.double(iBT.perm2@Delta[,"winRatio"]),
                         as.double(BT.perm2@DeltaResampling[iResample,,"winRatio"]),
                         tol = 1e-6)
            expect_equal(as.double(iBT.perm2@covariance),
                         as.double(BT.perm2@covarianceResampling[iResample,,]),
                         tol = 1e-6)
        }
    }
})

## * Bootstrap
test_that("Bootstrap", {

    BT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1)  + bin(toxicity1),
                         data = dt.sim, scoring.rule = method, seed = 10,
                         method.inference = "studentized bootstrap", n.resampling = 10)
    
    ## ** summary (two.sided)
    ## warnings because too few bootstrap samples
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
    ## summary(BT.boot)
    ##     Generalized pairwise comparisons with 2 prioritized endpoints

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : bootstrap resampling with 10 samples 
    ##                     CI computed using the studentized method; p-value by test inversion 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > right-censored pairs: probabilistic score based on the survival curves
    ## > neutral pairs   : ignored at lower priority endpoints
    ## > results
    ##   endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta   Delta p.value   CI [2.5 ; 97.5]
    ## eventtime1         1   100.00        20.87          21.91      57.22        0 -0.0104 -0.0104     0.9  [-0.1762;0.1284]
    ##  toxicity1       0.5    57.22        10.92          17.62      28.68        0 -0.0670 -0.0774     0.5  [-0.2693;0.0596]

    CI <- t(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.025, 0.975)))
    expect_equal(as.double(unlist(outSummaryPerc[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- t(apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = c(0.025, 0.975)))
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,1],
                BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,2])
    expect_equal(as.double(unlist(outSummaryStud[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** greater
    ## warnings because too few bootstrap samples
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "greater", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "greater", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.05)), Inf)
    expect_equal(as.double(unlist(outSummaryPerc[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.05)
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz, Inf)
    expect_equal(as.double(unlist(outSummaryStud[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    ## ** lower
    ## warnings because too few bootstrap samples
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "less", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "less", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(-Inf, apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.95)))
    expect_equal(as.double(unlist(outSummaryPerc[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.95)
    CI <- cbind(-Inf, BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz)
    
    expect_equal(as.double(unlist(outSummaryStud[,c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** check bootstrap
    vec.seed <- BT.boot@seed

    for(iResample in 1:2){ ## iResample <- 1
        set.seed(vec.seed[iResample])
        dt.boot <- dt.sim[sample.int(.N, size = .N, replace = TRUE)]

        iBT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1),
                              data = dt.boot, scoring.rule = method,
                              method.inference = "u-statistic")

        expect_equal(as.double(iBT.boot@Delta[,"netBenefit"]),
                     as.double(BT.boot@DeltaResampling[iResample,,"netBenefit"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.boot@Delta[,"winRatio"]),
                     as.double(BT.boot@DeltaResampling[iResample,,"winRatio"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.boot@covariance),
                     as.double(BT.boot@covarianceResampling[iResample,,]),
                     tol = 1e-6)

    }
})


## * Stratified bootstrap
test_that("stratified bootstrap", {
    BT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1)  + bin(toxicity1) + strata,
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "studentized bootstrap", n.resampling = 10)
    
    ## ** summary (two.sided)
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
    ## summary(BT.boot)
    ##       Generalized pairwise comparisons with 2 prioritized endpoints and 3 strata

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : bootstrap resampling with 10 samples 
    ##                     CI computed using the studentized method; p-value by test inversion 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > right-censored pairs: probabilistic score based on the survival curves
    ## > neutral pairs   : ignored at lower priority endpoints
    ## > uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints
    ## > results
    ##   endpoint threshold strata total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta   Delta p.value   CI [2.5 ; 97.5]
    ## eventtime1         1 global   100.00        20.77          19.57      55.25     4.40  0.0120   0.012     0.6  [-0.1158;0.1104]
    ##                           0    33.56         9.38           7.61      16.58     0.00  0.0527                                  
    ##                           1    34.74         4.50           4.58      24.41     1.25 -0.0023                                  
    ##                           2    31.70         6.90           7.38      14.27     3.15 -0.0153                                  
    ##  toxicity1       0.5 global    59.66        11.53          18.18      29.95     0.00 -0.0665 -0.0545     0.5  [-0.2217;0.0803]
    ##                           0    16.58         2.96           5.19       8.43     0.00 -0.0662                                  
    ##                           1    25.66         4.64           8.63      12.39     0.00 -0.1151                                  
    ##                           2    17.42         3.93           4.36       9.13     0.00 -0.0138                                  
    CI <- t(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.025, 0.975)))
    expect_equal(as.double(unlist(outSummaryPerc[outSummaryPerc$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- t(apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = c(0.025, 0.975)))
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,1],
                BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,2])
    expect_equal(as.double(unlist(outSummaryStud[outSummaryStud$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** greater
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "greater", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "greater", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.05)), Inf)
    expect_equal(as.double(unlist(outSummaryPerc[outSummaryPerc$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.05)
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz, Inf)
    
    expect_equal(as.double(unlist(outSummaryStud[outSummaryStud$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    ## ** lower
    outSummaryPerc <- suppressMessages(model.tables(BT.boot, alternative = "less", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressMessages(model.tables(BT.boot, alternative = "less", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(-Inf, apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.95)))
    expect_equal(as.double(unlist(outSummaryPerc[outSummaryPerc$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(sweep(BT.boot@DeltaResampling[,,"netBenefit"], MARGIN = 2, FUN = "-", STATS = coef(BT.boot))/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.95)
    CI <- cbind(-Inf, BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz)
    
    expect_equal(as.double(unlist(outSummaryStud[outSummaryStud$strata=="global",c("lower.ci","upper.ci")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** check bootstrap
    vec.seed <- BT.boot@seed

    for(iResample in 1:2){ ## iResample <- 1
        set.seed(vec.seed[iResample])
        dt.boot <- dt.sim[sample.int(.N, size = .N, replace = TRUE)]

        iBT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                              data = dt.boot, scoring.rule = method,
                              method.inference = "u-statistic")

        expect_equal(as.double(iBT.boot@Delta[,"netBenefit"]),
                     as.double(BT.boot@DeltaResampling[iResample,,"netBenefit"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.boot@Delta[,"winRatio"]),
                     as.double(BT.boot@DeltaResampling[iResample,,"winRatio"]),
                     tol = 1e-6)
        expect_equal(as.double(iBT.boot@covariance),
                     as.double(BT.boot@covarianceResampling[iResample,,]),
                     tol = 1e-6)

    }

    for(strataVar in c("treatment","strata")){
        BT.boot2 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1)  + bin(toxicity1) + strata,
                              data = dt.sim, scoring.rule = method, seed = 10, 
                              method.inference = "studentized bootstrap", strata.resampling = strataVar, n.resampling = 10)

        for(iResample in 1:2){ ## iResample <- 1
        set.seed(vec.seed[iResample])
            dt.boot2 <- copy(dt.sim)
            setkeyv(dt.boot2, cols = strataVar)
            dt.boot2 <- dt.boot2[,.SD[sample.int(.N, size = .N, replace = TRUE)], by = strataVar]

            iBT.boot2 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + bin(toxicity1) + strata,
                                   data = dt.boot2, scoring.rule = method,
                                   method.inference = "u-statistic")

            expect_equal(as.double(iBT.boot2@Delta[,"netBenefit"]),
                         as.double(BT.boot2@DeltaResampling[iResample,,"netBenefit"]),
                         tol = 1e-6)
            expect_equal(as.double(iBT.boot2@Delta[,"winRatio"]),
                         as.double(BT.boot2@DeltaResampling[iResample,,"winRatio"]),
                         tol = 1e-6)
            expect_equal(as.double(iBT.boot2@covariance),
                         as.double(BT.boot2@covarianceResampling[iResample,,]),
                         tol = 1e-6)

        }
    }
})


## * t-test example
## ** data
set.seed(10)
df <- rbind(data.frame(Group = "T",
                       score = rnorm(25, mean = 0),
                       stringsAsFactors = FALSE),
            data.frame(Group = "C",
                       score = rnorm(25, mean = 0.75),
                       stringsAsFactors = FALSE)
            )

## ** BT
e.perm <- BuyseTest(Group ~ cont(score),
                    data = df,
                    method.inference = "studentized permutation", n.resampling = 200,
                    trace = 0)
e.boot <- BuyseTest(Group ~ cont(score),
                    data = df,
                    method.inference = "studentized bootstrap", n.resampling = 200,
                    trace = 0)
e.ustat <- BuyseTest(Group ~ cont(score),
                     data = df,
                     method.inference = "u-statistic",
                     trace = 0)

## ** confint (two.sided)
test_that("compare with t-test (two.sided)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "two.sided")

    ls.res <- list(perm = suppressWarnings(confint(e.perm, alternative = "two.sided")),
                   percboot = suppressMessages(confint(e.boot, alternative = "two.sided", method.ci.resampling = "percentile")),
                   gausboot = suppressMessages(confint(e.boot, alternative = "two.sided", method.ci.resampling = "gaussian", transformation = FALSE)),
                   gausboot.trans = suppressMessages(confint(e.boot, alternative = "two.sided", method.ci.resampling = "gaussian", transformation = TRUE)),
                   studboot = suppressMessages(confint(e.boot, alternative = "two.sided", method.ci.resampling = "studentized", transformation = FALSE)),
                   studboot.trans = suppressMessages(confint(e.boot, alternative = "two.sided", method.ci.resampling = "studentized", transformation = TRUE)),
                   ustat = confint(e.ustat, alternative = "two.sided", transformation = FALSE),
                   ustat.trans = confint(e.ustat, alternative = "two.sided", transformation = TRUE)
                   )
    M.res <- do.call(rbind,ls.res)
    rownames(M.res) <- names(ls.res)

    ## same estimates for all
    expect_true(all(abs(diff(M.res[,"estimate"]))<1e-6))
    ## same variance for percentile and gaussian bootstrap
    expect_true(all(abs(diff(M.res[c("percboot","gausboot","gausboot.trans"),"se"])<1e-6)))
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(all(abs(diff(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"])<1e-6)))
    ## lower.ci smaller than upper ci
    expect_true(all(is.na(M.res[1,c("lower.ci","upper.ci")])))
    expect_true(all(M.res[-1,"lower.ci"]<M.res[-1,"upper.ci"]))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- data.frame("estimate" = c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112), 
                     "se" = c(0.14543358, 0.15109607, 0.15109607, 0.15109607, 0.14543358, 0.14543358, 0.14543358, 0.14543358), 
                     "lower.ci" = c(NA, -0.65632966, -0.70734285, -0.6657388, -0.70907164, -0.6244038, -0.69624457, -0.65276625), 
                     "upper.ci" = c(NA, -0.02091063, -0.11505715, -0.07093912, -0.07297107, -0.00322604, -0.12615543, -0.09372943), 
                     "null" = c(0, 0, 0, 0, 0, 0, 0, 0), 
                     "p.value" = c(0.01, 0.03, 0.00649967, 0.01925831, 0, 0.04, 0.00469266, 0.01252311))

    expect_equivalent(GS, M.res, tol = 1e-4)
})

## ** confint (greater)
test_that("compare with t-test (greater)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "greater")

    ls.res <- list(perm = suppressWarnings(confint(e.perm, alternative = "greater")),
                   percboot = suppressMessages(confint(e.boot, alternative = "greater", method.ci.resampling = "percentile")),
                   gausboot = suppressMessages(confint(e.boot, alternative = "greater", method.ci.resampling = "gaussian", transformation = FALSE)),
                   gausboot.trans = suppressMessages(confint(e.boot, alternative = "greater", method.ci.resampling = "gaussian", transformation = TRUE)),
                   studboot = suppressMessages(confint(e.boot, alternative = "greater", method.ci.resampling = "studentized", transformation = FALSE)),
                   studboot.trans = suppressMessages(confint(e.boot, alternative = "greater", method.ci.resampling = "studentized", transformation = TRUE)),
                   ustat = confint(e.ustat, alternative = "greater", transformation = FALSE),
                   ustat.trans = confint(e.ustat, alternative = "greater", transformation = TRUE)
                   )
    M.res <- do.call(rbind,ls.res)
    rownames(M.res) <- names(ls.res)

    
    ## same estimates for all
    expect_true(all(abs(diff(M.res[,"estimate"]))<1e-6))
    ## same variance for percentile and gaussian bootstrap
    expect_true(all(abs(diff(M.res[c("percboot","gausboot","gausboot.trans"),"se"])<1e-6)))
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(all(abs(diff(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"])<1e-6)))
    ## lower.ci smaller than upper ci
    expect_true(all(is.na(M.res[1,c("lower.ci","upper.ci")])))
    expect_true(all(M.res[-1,"lower.ci"]<M.res[-1,"upper.ci"]))
    ## upper ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"upper.ci"]==1))
    expect_true(all(is.infinite(M.res[-grep("trans|perm",rownames(M.res)),"upper.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- data.frame("estimate" = c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112), 
                     "se" = c(0.14543358, 0.15109607, 0.15109607, 0.15109607, 0.14543358, 0.14543358, 0.14543358, 0.14543358), 
                     "lower.ci" = c(NA, -0.62211538, -0.65973091, -0.6316809, -0.66549951, -0.60288984, -0.65041694, -0.61996641), 
                     "upper.ci" = c(NA, Inf, Inf, 1, Inf, 1, Inf, 1), 
                     "null" = c(0, 0, 0, 0, 0, 0, 0, 0), 
                     "p.value" = c(0.995, 0.985, 0.99675016, 0.99037084, 1, 0.98, 0.99765367, 0.99373845))
    expect_equivalent(GS, M.res, tol = 1e-4)
})

## ** confint (less)
test_that("compare with t-test (less)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "less")

    ls.res <- list(perm = suppressWarnings(confint(e.perm, alternative = "less")),
                   percboot = suppressMessages(confint(e.boot, alternative = "less", method.ci.resampling = "percentile")),
                   gausboot = suppressMessages(confint(e.boot, alternative = "less", method.ci.resampling = "gaussian", transformation = FALSE)),
                   gausboot.trans = suppressMessages(confint(e.boot, alternative = "less", method.ci.resampling = "gaussian", transformation = TRUE)),
                   studboot = suppressMessages(confint(e.boot, alternative = "less", method.ci.resampling = "studentized", transformation = FALSE)),
                   studboot.trans = suppressMessages(confint(e.boot, alternative = "less", method.ci.resampling = "studentized", transformation = TRUE)),
                   ustat = confint(e.ustat, alternative = "less", transformation = FALSE),
                   ustat.trans = confint(e.ustat, alternative = "less", transformation = TRUE)
                   )
    M.res <- do.call(rbind,ls.res)
    rownames(M.res) <- names(ls.res)

    
    ## same estimates for all
    expect_true(all(abs(diff(M.res[,"estimate"]))<1e-6))
    ## same variance for percentile and gaussian bootstrap
    expect_true(all(abs(diff(M.res[c("percboot","gausboot","gausboot.trans"),"se"])<1e-6)))
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(all(abs(diff(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"])<1e-6)))
    ## lower.ci smaller than upper ci
    expect_true(all(is.na(M.res[1,c("lower.ci","upper.ci")])))
    expect_true(all(M.res[-1,"lower.ci"]<M.res[-1,"upper.ci"]))
    ## lower ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"lower.ci"]==-1))
    expect_true(all(is.infinite(M.res[-grep("trans|perm",rownames(M.res)),"lower.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- data.frame("estimate" = c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112), 
           "se" = c(0.14543358, 0.15109607, 0.15109607, 0.15109607, 0.14543358, 0.14543358, 0.14543358, 0.14543358), 
           "lower.ci" = c(NA, -Inf, -Inf, -1, -Inf, -1, -Inf, -1), 
           "upper.ci" = c(NA, -0.15025183, -0.16266909, -0.1291752, -0.18137423, -0.14017782, -0.17198306, -0.14806218), 
           "null" = c(0, 0, 0, 0, 0, 0, 0, 0), 
           "p.value" = c(0.005, 0.015, 0.00324984, 0.00962916, 0, 0.02, 0.00234633, 0.00626155))
    expect_equivalent(GS, M.res, tol = 1e-4)
})

