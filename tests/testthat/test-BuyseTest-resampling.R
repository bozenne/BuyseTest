### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: Mar 13 2023 (09:32) 
##           By: Brice Ozenne
##     Update #: 190
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
                  add.1.pperm = FALSE,
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
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
    ## summary(BT.perm)
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
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"])),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"])))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)
       
    ## ** summary (greater)
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "greater", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "greater", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) <= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) <= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] <= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] <= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "less", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "less", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) >= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) >= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] >= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] >= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
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
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
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
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(abs(BT.perm@Delta[1,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,1,"netBenefit"])),
                 mean(abs(BT.perm@Delta[2,"netBenefit"]) <= abs(BT.perm@DeltaResampling[,2,"netBenefit"])))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "greater", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "greater", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) <= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) <= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] <= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] <= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummaryPerc <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "less", method.ci.resampling = "percentile"))
    outSummaryStud <- suppressWarnings(summary(BT.perm, print = FALSE, alternative = "less", method.ci.resampling = "studentized"))
    
    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"]/sqrt(BT.perm@covariance[1,"netBenefit"]) >= BT.perm@DeltaResampling[,1,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,1,"netBenefit"])),
                 mean(BT.perm@Delta[2,"netBenefit"]/sqrt(BT.perm@covariance[2,"netBenefit"]) >= BT.perm@DeltaResampling[,2,"netBenefit"]/sqrt(BT.perm@covarianceResampling[,2,"netBenefit"])))
    expect_equal(outSummaryStud$table[outSummaryStud$table$strata=="global","p.value"],
                 p.value)

    p.value <- c(mean(BT.perm@Delta[1,"netBenefit"] >= BT.perm@DeltaResampling[,1,"netBenefit"]),
                 mean(BT.perm@Delta[2,"netBenefit"] >= BT.perm@DeltaResampling[,2,"netBenefit"]))
    expect_equal(outSummaryPerc$table[outSummaryPerc$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1 
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

        set.seed(10)
        for(iResample in 1:2){ ## iResample <- 1 
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
    outSummaryPerc <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE))
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
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- t(apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = c(0.025, 0.975)))
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,1],
                BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,2])
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** greater
    ## warnings because too few bootstrap samples
    outSummaryPerc <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "greater", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "greater", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.05)), Inf)
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.05)
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz, Inf)
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    ## ** lower
    ## warnings because too few bootstrap samples
    outSummaryPerc <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "less", method.ci.resampling = "percentile", transform = FALSE))
    outSummaryStud <- suppressWarnings(summary(BT.boot, print = FALSE, alternative = "less", method.ci.resampling = "studentized", transform = FALSE))

    CI <- cbind(-Inf, apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.95)))
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.95)
    CI <- cbind(-Inf, BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz)
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** check bootstrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
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
    outSummaryPerc <- summary(BT.boot, print = FALSE, alternative = "two.sided", method.ci.resampling = "percentile", transform = FALSE)
    outSummaryStud <- summary(BT.boot, print = FALSE, alternative = "two.sided", method.ci.resampling = "studentized", transform = FALSE)
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
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- t(apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = c(0.025, 0.975)))
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,1],
                BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz[,2])
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** greater
    outSummaryPerc <- summary(BT.boot, print = FALSE, alternative = "greater", method.ci.resampling = "percentile", transform = FALSE)
    outSummaryStud <- summary(BT.boot, print = FALSE, alternative = "greater", method.ci.resampling = "studentized", transform = FALSE)

    CI <- cbind(apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.05)), Inf)
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.05)
    CI <- cbind(BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz, Inf)
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    ## ** lower
    outSummaryPerc <- summary(BT.boot, print = FALSE, alternative = "less", method.ci.resampling = "percentile", transform = FALSE)
    outSummaryStud <- summary(BT.boot, print = FALSE, alternative = "less", method.ci.resampling = "studentized", transform = FALSE)

    CI <- cbind(-Inf, apply(BT.boot@DeltaResampling[,,"netBenefit"], 2, quantile, probs = c(0.95)))
    expect_equal(as.double(unlist(outSummaryPerc$table[outSummaryPerc$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)

    qz <- apply(apply(BT.boot@DeltaResampling[,,"netBenefit"],2,scale,scale=FALSE)/sqrt(BT.boot@covarianceResampling[,,"netBenefit"]), 2, quantile, probs = 0.95)
    CI <- cbind(-Inf, BT.boot@Delta[,"netBenefit"] + sqrt(BT.boot@covariance[,"netBenefit"]) * qz)
    
    expect_equal(as.double(unlist(outSummaryStud$table[outSummaryStud$table$strata=="global",c("CIinf.Delta","CIsup.Delta")])),
                 as.double(CI), tol = 1e-6)
    
    ## ** check bootstrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
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

        set.seed(10)
        for(iResample in 1:2){ ## iResample <- 1
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
                   percboot = confint(e.boot, alternative = "two.sided", method.ci.resampling = "percentile"),
                   gausboot = confint(e.boot, alternative = "two.sided", method.ci.resampling = "gaussian", transformation = FALSE),
                   gausboot.trans = confint(e.boot, alternative = "two.sided", method.ci.resampling = "gaussian", transformation = TRUE),
                   studboot = confint(e.boot, alternative = "two.sided", method.ci.resampling = "studentized", transformation = FALSE),
                   studboot.trans = confint(e.boot, alternative = "two.sided", method.ci.resampling = "studentized", transformation = TRUE),
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
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, 0.145434, 0.147342, 0.147342, 0.147342, 0.145434, 0.145434, 0.145434, 0.145434, -0.657594, -0.665138, -0.699984, -0.66126, -0.745069, -0.632558, -0.696245, -0.652766, -0.095572, -0.067126, -0.122416, -0.078896, -0.126607, -0.057924, -0.126155, -0.093729, 0.01, 0, 0.005258, 0.01672, 0, 0.01, 0.004693, 0.012523), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    ## expect_equal(GS, M.res, tol = 1e-4)
})

## ** confint (greater)
test_that("compare with t-test (greater)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "greater")

    ls.res <- list(perm = suppressWarnings(confint(e.perm, alternative = "greater")),
                   percboot = confint(e.boot, alternative = "greater", method.ci.resampling = "percentile"),
                   gausboot = confint(e.boot, alternative = "greater", method.ci.resampling = "gaussian", transformation = FALSE),
                   gausboot.trans = confint(e.boot, alternative = "greater", method.ci.resampling = "gaussian", transformation = TRUE),
                   studboot = confint(e.boot, alternative = "greater", method.ci.resampling = "studentized", transformation = FALSE),
                   studboot.trans = confint(e.boot, alternative = "greater", method.ci.resampling = "studentized", transformation = TRUE),
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
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## upper ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"upper.ci"]==1))
    expect_true(all(is.infinite(M.res[-grep("trans|perm",rownames(M.res)),"upper.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, 0.145434, 0.147342, 0.147342, 0.147342, 0.145434, 0.145434, 0.145434, 0.145434, -0.636778, -0.622851, -0.653555, -0.627628, -0.685803, -0.606192, -0.650417, -0.619966, 1, Inf, Inf, 1, Inf, 1, Inf, 1, 1, 1, 0.997371, 0.99164, 0.995, 0.995, 0.997654, 0.993738), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 )
    ## expect_equal(GS, M.res, tol = 1e-5)
})

## ** confint (less)
test_that("compare with t-test (less)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "less")

    ls.res <- list(perm = suppressWarnings(confint(e.perm, alternative = "less")),
                   percboot = confint(e.boot, alternative = "less", method.ci.resampling = "percentile"),
                   gausboot = confint(e.boot, alternative = "less", method.ci.resampling = "gaussian", transformation = FALSE),
                   gausboot.trans = confint(e.boot, alternative = "less", method.ci.resampling = "gaussian", transformation = TRUE),
                   studboot = confint(e.boot, alternative = "less", method.ci.resampling = "studentized", transformation = FALSE),
                   studboot.trans = confint(e.boot, alternative = "less", method.ci.resampling = "studentized", transformation = TRUE),
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
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## lower ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"lower.ci"]==-1))
    expect_true(all(is.infinite(M.res[-grep("trans|perm",rownames(M.res)),"lower.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(-0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, -0.4112, 0.145434, 0.147342, 0.147342, 0.147342, 0.145434, 0.145434, 0.145434, 0.145434, -1, -Inf, -Inf, -1, -Inf, -1, -Inf, -1, -0.138781, -0.157787, -0.168845, -0.135772, -0.19426, -0.144137, -0.171983, -0.148062, 0, 0, 0.002629, 0.00836, 0, 0.005, 0.002346, 0.006262), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 )
    ## expect_equal(GS, M.res, tol = 1e-5)
})

