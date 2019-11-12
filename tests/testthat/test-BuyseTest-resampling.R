### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: nov 12 2019 (11:05) 
##           By: Brice Ozenne
##     Update #: 136
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

student <- FALSE ## need to fix asymptotics for several outcome to use studentized bootstrap

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  keep.survival = FALSE,
                  order.Hprojection = 1,
                  trace = 0)
n.patients <- 100
method <- "Peron"

## * Simulate data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients,
                       n.C = n.patients,
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censoring.T = rep(1,3)),
                       n.strata = 3)

## * Permutation
test_that("permutation", {
    BT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "permutation", n.resampling = 20)

    ## ** summary (two.sided)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "two.sided")
    ## summary(BT.perm)
    ##       Hierarchical generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : permutation test with 20 samples 
    ##                     confidence intervals/p-values computed using the quantiles of the empirical distribution 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > censored pairs  : use Kaplan Meier survival curves to compute the score
    ## > uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints (if any)
    ## > results
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf   delta  Delta p.value 
    ## eventtime1     1e-12 global 100.00     49.15       48.28    0.00  2.57  0.0087 0.0087       1 
    ##                           0  33.56     18.07       15.50    0.00  0.00  0.0766                
    ##                           1  34.74     17.24       17.03    0.00  0.47  0.0062                
    ##                           2  31.70     13.84       15.76    0.00  2.10 -0.0605                
    ##  toxicity1       0.5 global   2.57      0.61        0.51    1.46  0.00  0.0010 0.0097       1 
    ##                           0   0.00      0.00        0.00    0.00  0.00  0.0000                
    ##                           1   0.47      0.11        0.12    0.24  0.00 -0.0002                
    ##                           2   2.10      0.50        0.39    1.22  0.00  0.0034                
    
    p.value <- c(mean(abs(BT.perm@Delta.netBenefit[1]) <= abs(BT.perm@DeltaResampling.netBenefit[,1])),
                 mean(abs(BT.perm@Delta.netBenefit[2]) <= abs(BT.perm@DeltaResampling.netBenefit[,2])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
       
    ## ** summary (greater)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] <= BT.perm@DeltaResampling.netBenefit[,1]),
                 mean(BT.perm@Delta.netBenefit[2] <= BT.perm@DeltaResampling.netBenefit[,2]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] > BT.perm@DeltaResampling.netBenefit[,1]),
                 mean(BT.perm@Delta.netBenefit[2] > BT.perm@DeltaResampling.netBenefit[,2]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.perm <- copy(dt.sim)

        dt.perm[, treatment := treatment[sample.int(.N, size = .N, replace = FALSE)] ]
        expect_equal(table(dt.perm$treatment), table(dt.sim$treatment))

        iBT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                              data = dt.perm, scoring.rule = method,
                              method.inference = "none")

        expect_equal(as.double(iBT.perm@delta.netBenefit),
                     as.double(BT.perm@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.perm@delta.winRatio),
                     as.double(BT.perm@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.perm@Delta.netBenefit),
                     as.double(BT.perm@DeltaResampling.netBenefit[iResample,]))
        expect_equal(as.double(iBT.perm@Delta.winRatio),
                     as.double(BT.perm@DeltaResampling.winRatio[iResample,]))
    }

})


## * Stratified permutation
test_that(paste0("stratified permutation"), {
    BT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "permutation", strata.resampling = "strata", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "two.sided")

    p.value <- c(mean(abs(BT.perm@Delta.netBenefit[1]) <= abs(BT.perm@DeltaResampling.netBenefit[,1])),
                 mean(abs(BT.perm@Delta.netBenefit[2]) <= abs(BT.perm@DeltaResampling.netBenefit[,2])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] <= BT.perm@DeltaResampling.netBenefit[,1]),
                 mean(BT.perm@Delta.netBenefit[2] <= BT.perm@DeltaResampling.netBenefit[,2]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] >= BT.perm@DeltaResampling.netBenefit[,1]),
                 mean(BT.perm@Delta.netBenefit[2] >= BT.perm@DeltaResampling.netBenefit[,2]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1 
        dt.perm <- copy(dt.sim)
        setkeyv(dt.perm, cols = "strata")
        
        dt.perm[, treatment := treatment[sample.int(.N, size = .N, replace = FALSE)], by = "strata"]

        iBT.perm <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                              data = dt.perm, scoring.rule = method,
                              method.inference = "none")

        expect_equal(as.double(iBT.perm@delta.netBenefit),
                     as.double(BT.perm@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.perm@delta.winRatio),
                     as.double(BT.perm@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.perm@Delta.netBenefit),
                     as.double(BT.perm@DeltaResampling.netBenefit[iResample,]))
        expect_equal(as.double(iBT.perm@Delta.winRatio),
                     as.double(BT.perm@DeltaResampling.winRatio[iResample,]))
    }
})

## * Bootstrap
test_that("Bootstrap", {

    BT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0)  + bin(toxicity1) + strata,
                         data = dt.sim, scoring.rule = method, seed = 10, 
                         method.inference = "bootstrap", n.resampling = 20)

    if(student){
        BT.bootT <- suppressWarnings(BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0)  + bin(toxicity1) + strata,
                                               data = dt.sim, scoring.rule = method, seed = 10,
                                               method.inference = "studentized bootstrap", n.resampling = 20))
    }
    
    ## same point estimate with or without computation of the variance
    expect_true(sum(dim(BT.boot@covarianceResampling))==0)
    if(student){
        expect_equal(BT.boot@Delta.netBenefit,BT.bootT@Delta.netBenefit)
        expect_equal(BT.boot@Delta.winRatio,BT.bootT@Delta.winRatio)
        expect_true(sum(dim(BT.bootT@covarianceResampling))!=0)
    }
    
    
    ## ** summary (two.sided)
    outSummary <- summary(BT.boot, print = FALSE, alternative = "two.sided")
    ## summary(BT.boot)
    ##       Hierarchical generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ## > statistic       : net benefit (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > confidence level: 0.95 
    ## > inference       : bootstrap resampling with 20 samples 
    ##                     confidence intervals/p-values computed using the quantiles of the empirical distribution 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > censored pairs  : use Kaplan Meier survival curves to compute the score
    ## > uninformative pairs: no contribution at the current endpoint, analyzed at later endpoints (if any)
    ## > results
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf   delta  Delta  CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12 global 100.00     49.15       48.28    0.00  2.57  0.0087 0.0087   [-0.132;0.193]    0.95 
    ##                           0  33.56     18.07       15.50    0.00  0.00  0.0766                                 
    ##                           1  34.74     17.24       17.03    0.00  0.47  0.0062                                 
    ##                           2  31.70     13.84       15.76    0.00  2.10 -0.0605                                 
    ##  toxicity1       0.5 global   2.57      0.61        0.51    1.46  0.00  0.0010 0.0097 [-0.1415;0.1925]    0.95 
    ##                           0   0.00      0.00        0.00    0.00  0.00  0.0000                                 
    ##                           1   0.47      0.11        0.12    0.24  0.00 -0.0002                                 
    ##                           2   2.10      0.50        0.39    1.22  0.00  0.0034

    ## CI
    if(student){
        perBoot.confint <- confint(BT.bootT, method.ci.resampling = "percentile", alternative = "two.sided")
        perBoot.manual <- t(apply(BT.bootT@DeltaResampling.netBenefit, 2, quantile, probs = c(0.025, 0.975)))

        expect_equal(unname(perBoot.confint[,c("lower.ci","upper.ci")]),unname(perBoot.manual), tol = 1e-6)

        gausBoot.confint <- confint(BT.bootT, method.ci.resampling = "gaussian", alternative = "two.sided", transform = FALSE)
        expect_equal(unname(gausBoot.confint[,"se"]), unname(apply(BT.bootT@DeltaResampling.netBenefit, 2, sd)), tol = 1e-6)
        expect_equal(unname(gausBoot.confint[,"lower.ci"]), unname(gausBoot.confint[,"estimate"] + qnorm(0.025) * gausBoot.confint[,"se"]), tol = 1e-6)
        expect_equal(unname(gausBoot.confint[,"upper.ci"]), unname(gausBoot.confint[,"estimate"] + qnorm(0.975) * gausBoot.confint[,"se"]), tol = 1e-6)
    }
    
    ## ** greater
    if(student){
    perBoot.confint <- confint(BT.bootT, method.ci.resampling = "percentile", alternative = "greater")
    perBoot.manual <- cbind(apply(BT.bootT@DeltaResampling.netBenefit, 2, quantile, probs = 0.05),Inf)

    expect_equal(unname(perBoot.confint[,c("lower.ci","upper.ci")]),unname(perBoot.manual), tol = 1e-6)

    gausBoot.confint <- confint(BT.bootT, method.ci.resampling = "gaussian", alternative = "greater", transform = FALSE)
    expect_equal(unname(gausBoot.confint[,"se"]), unname(apply(BT.bootT@DeltaResampling.netBenefit, 2, sd)), tol = 1e-6)
    expect_equal(unname(gausBoot.confint[,"lower.ci"]), unname(gausBoot.confint[,"estimate"] + qnorm(0.05) * gausBoot.confint[,"se"]), tol = 1e-6)
    }
    
    ## ** lower
    if(student){
        perBoot.confint <- confint(BT.bootT, method.ci.resampling = "percentile", alternative = "less")
        perBoot.manual <- cbind(-Inf,apply(BT.bootT@DeltaResampling.netBenefit, 2, quantile, probs = 0.95))

        expect_equal(unname(perBoot.confint[,c("lower.ci","upper.ci")]),unname(perBoot.manual), tol = 1e-6)

        gausBoot.confint <- confint(BT.bootT, method.ci.resampling = "gaussian", alternative = "less", transform = FALSE)
        expect_equal(unname(gausBoot.confint[,"se"]), unname(apply(BT.bootT@DeltaResampling.netBenefit, 2, sd)), tol = 1e-6)
        expect_equal(unname(gausBoot.confint[,"upper.ci"]), unname(gausBoot.confint[,"estimate"] + qnorm(0.95) * gausBoot.confint[,"se"]), tol = 1e-6)
    }
    
    ## ** check bootstrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.boot <- dt.sim[sample.int(.N, size = .N, replace = TRUE)]
        
        iBT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                              data = dt.boot, scoring.rule = method, method.inference = "none")

        expect_equal(as.double(iBT.boot@delta.netBenefit),
                     as.double(BT.boot@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.boot@delta.winRatio),
                     as.double(BT.boot@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.boot@Delta.netBenefit),
                     as.double(BT.boot@DeltaResampling.netBenefit[iResample,]))
        expect_equal(as.double(iBT.boot@Delta.winRatio),
                     as.double(BT.boot@DeltaResampling.winRatio[iResample,]))

        if(student){
            iBT.bootT <- suppressWarnings(BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                                                    data = dt.boot, scoring.rule = method,                                               
                                                    method.inference = "u-statistic"))
            expect_equal(as.double(iBT.bootT@delta.netBenefit),
                         as.double(BT.bootT@deltaResampling.netBenefit[,,iResample]))
            expect_equal(as.double(iBT.bootT@delta.winRatio),
                         as.double(BT.bootT@deltaResampling.winRatio[,,iResample]))
            expect_equal(as.double(iBT.bootT@Delta.netBenefit),
                         as.double(BT.bootT@DeltaResampling.netBenefit[iResample,]))
            expect_equal(as.double(iBT.bootT@Delta.winRatio),
                         as.double(BT.bootT@DeltaResampling.winRatio[iResample,]))

            expect_equal(as.double(iBT.bootT@covariance),
                         as.double(BT.bootT@covarianceResampling[iResample,,]))
        }
    }
})


## * Stratified bootstrap
for(iStrataVar in c("treatment","strata")){ ## iStrataVar <- "treatment"
    test_that(paste0("Stratified bootstrap (",iStrataVar,")"), {
        BT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0)  + bin(toxicity1) + strata,
                             data = dt.sim, scoring.rule = method, seed = 10, 
                             method.inference = "bootstrap", strata.resampling = iStrataVar, n.resampling = 20)
        if(student){
            BT.bootT <- suppressWarnings(BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0)  + bin(toxicity1) + strata,
                                                   data = dt.sim, scoring.rule = method, seed = 10, 
                                                   method.inference = "studentized bootstrap", strata.resampling = iStrataVar, n.resampling = 20))
        }
    
        ## same point estimate with or without computation of the variance
        if(student){
            expect_equal(BT.boot@Delta.winRatio,BT.bootT@Delta.winRatio)
            expect_equal(BT.boot@Delta.netBenefit,BT.bootT@Delta.netBenefit)
            expect_true(sum(dim(BT.boot@covarianceResampling))==0)
            expect_true(sum(dim(BT.bootT@covarianceResampling))!=0)
        }
    
    ## ** check bootstrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.boot <- copy(dt.sim)
        setkeyv(dt.boot, cols = iStrataVar)
        dt.boot <- dt.boot[, .SD[sample.int(.N, size = .N, replace = TRUE)], by = iStrataVar]

        iBT.boot <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                              data = dt.boot, scoring.rule = method, method.inference = "none")

        expect_equal(as.double(iBT.boot@delta.netBenefit),
                     as.double(BT.boot@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.boot@delta.winRatio),
                     as.double(BT.boot@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.boot@Delta.netBenefit),
                     as.double(BT.boot@DeltaResampling.netBenefit[iResample,]))
        expect_equal(as.double(iBT.boot@Delta.winRatio),
                     as.double(BT.boot@DeltaResampling.winRatio[iResample,]))

        if(student){
            iBT.bootT <- suppressWarnings(BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0) + bin(toxicity1) + strata,
                                              data = dt.boot, scoring.rule = method,
                                              method.inference = "u-statistic"))

            expect_equal(as.double(iBT.bootT@delta.netBenefit),
                         as.double(BT.boot@deltaResampling.netBenefit[,,iResample]))
            expect_equal(as.double(iBT.bootT@delta.winRatio),
                         as.double(BT.boot@deltaResampling.winRatio[,,iResample]))
            expect_equal(as.double(iBT.bootT@Delta.netBenefit),
                         as.double(BT.boot@DeltaResampling.netBenefit[iResample,]))
            expect_equal(as.double(iBT.bootT@Delta.winRatio),
                         as.double(BT.boot@DeltaResampling.winRatio[iResample,]))

            expect_equal(as.double(iBT.bootT@covariance),
                         as.double(BT.bootT@covarianceResampling[iResample,,]))
        }
    }
    })
}


## * t-test example
## ** data
set.seed(10)
df <- rbind(data.frame(Group = "T", score = rnorm(25, mean = 0)),
            data.frame(Group = "C", score = rnorm(25, mean = 0.75))
            )

## ** BT
e.perm <- BuyseTest(Group ~ cont(score),
                    data = df,
                    method.inference = "permutation", n.resampling = 200,
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

    ls.res <- list(perm = confint(e.perm, alternative = "two.sided"),
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
    expect_true(length(unique(round(M.res[,"estimate"],6)))==1)
    ## same variance for percentile and gaussian bootstrap
    expect_true(length(unique(round(M.res[c("percboot","gausboot","gausboot.trans"),"se"],6)))==1)
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(length(unique(round(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"],6)))==1)
    ## lower.ci smaller than upper ci
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112,
                   0.176403, 0.153859, 0.153859, 0.153859, 0.145434, 0.145434, 0.145434, 0.145434,
                   0.08928, 0.135485, 0.109643, 0.056216, 0.088518, 0.144998, 0.126155, 0.093729,
                   0.74576, 0.698294, 0.712757, 0.673889, 0.92824, 0.612199, 0.696245, 0.652766,
                   0.015, 0.005, 0.007527, 0.024473, 0.02, 0.005, 0.004693, 0.012523), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),
                                 c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 )
    ## expect_equal(GS, M.res, tol = 1e-5)
})

## ** confint (greater)
test_that("compare with t-test (greater)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "greater")

    ls.res <- list(perm = confint(e.perm, alternative = "greater"),
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
    expect_true(length(unique(round(M.res[,"estimate"],6)))==1)
    ## same variance for percentile and gaussian bootstrap
    expect_true(length(unique(round(M.res[c("percboot","gausboot","gausboot.trans"),"se"],6)))==1)
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(length(unique(round(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"],6)))==1)
    ## lower.ci smaller than upper ci
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## upper ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"upper.ci"]==1))
    expect_true(all(is.infinite(M.res[-grep("trans",rownames(M.res)),"upper.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112,
                   0.176403, 0.153859, 0.153859, 0.153859, 0.145434, 0.145434, 0.145434, 0.145434,
                   0.13072, 0.179298, 0.158125, 0.116957, 0.132322, 0.19367, 0.171983, 0.148062,
                   Inf, Inf, Inf, 1, Inf, 1, Inf, 1,
                   0.005, 0.005, 0.003763, 0.012236, 0.01, 0, 0.002346, 0.006262), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),
                                 c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    ## expect_equal(GS, M.res, tol = 1e-5)
})

## ** confint (less)
test_that("compare with t-test (less)", {
    ## just to see what to expect
    res.tt <- t.test(y = df[df$Group=="T","score"], x = df[df$Group=="C","score"], alternative = "less")

    ls.res <- list(perm = confint(e.perm, alternative = "less"),
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
    expect_true(length(unique(round(M.res[,"estimate"],6)))==1)
    ## same variance for percentile and gaussian bootstrap
    expect_true(length(unique(round(M.res[c("percboot","gausboot","gausboot.trans"),"se"],6)))==1)
    ## same variance for studentized bootstrap and asymptotic 
    expect_true(length(unique(round(M.res[c("studboot","studboot.trans","ustat","ustat.trans"),"se"],6)))==1)
    ## lower.ci smaller than upper ci
    expect_true(all(M.res[,"lower.ci"]<M.res[,"upper.ci"]))
    ## upper ci
    expect_true(all(M.res[grep("trans",rownames(M.res)),"lower.ci"]==-1))
    expect_true(all(is.infinite(M.res[-grep("trans",rownames(M.res)),"lower.ci"])))
    ## check values
    ## butils::object2script(M.res, digits = 6)
    GS <- matrix(c(0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112, 0.4112,
                   0.176403, 0.153859, 0.153859, 0.153859, 0.145434, 0.145434, 0.145434, 0.145434,
                   -Inf, -Inf, -Inf, -1, -Inf, -1, -Inf, -1,
                   0.7136, 0.663419, 0.664275, 0.639078, 0.788808, 0.581402, 0.650417, 0.619966,
                   0.995, 0.995, 0.996237, 0.987764, 0.99, 0.995, 0.997654, 0.993738), 
                 nrow = 8, 
                 ncol = 5, 
                 dimnames = list(c("perm", "percboot", "gausboot", "gausboot.trans", "studboot", "studboot.trans", "ustat", "ustat.trans"),
                                 c("estimate", "se", "lower.ci", "upper.ci", "p.value")) 
                 ) 
    ## expect_equal(GS, M.res, tol = 1e-5)
})

