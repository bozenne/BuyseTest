### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: okt 30 2018 (16:19) 
##           By: Brice Ozenne
##     Update #: 93
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
}


## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,                  
                  trace = 0)
n.patients <- 100
method <- "Peron"

## * Simulate data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients,
                       n.C = n.patients,
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)),
                       n.strata = 3)

## * Permutation
test_that("permutation", {
    BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                         data = dt.sim, method.tte = method, seed = 10, 
                         method.inference = "permutation", n.resampling = 20)

    ## ** summary (two.sided)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "two.sided")
    ## summary(BT.perm, alternative = "two.sided")
    ##             Generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ##  > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
    ##  > null hypothesis : Delta == 0 
    ##  > permutation test: 50 samples, confidence level 0.95 
    ##  > treatment groups: 0 (control) vs. 1 (treatment) 
    ##  > censored pairs  : imputation using Kaplan Meier stratified by treatment group

    ##  > results
    ##    endpoint threshold strata  total favorable unfavorable neutral uninf   delta  Delta  CI [2.5 ; 97.5] p.value 
    ##  eventtime1     1e-12 global 100.00     48.83       47.77    0.00  3.40  0.0106 0.0106 [-0.1698;0.2138]    0.92 
    ##                            a  33.56     18.07       15.50    0.00  0.00  0.0257     NA                          
    ##                            b  34.74     16.92       17.03    0.00  0.79 -0.0011     NA                          
    ##                            c  31.70     13.84       15.24    0.00  2.62 -0.0140     NA                          
    ##   toxicity1       0.5 global   3.40      0.61        0.86    1.94  0.00 -0.0025 0.0082 [-0.1696;0.2127]    0.94 
    ##                            a   0.00      0.00        0.00    0.00  0.00  0.0000     NA                          
    ##                            b   0.79      0.11        0.24    0.43  0.00 -0.0013     NA                          
    ##                            c   2.62      0.50        0.61    1.51  0.00 -0.0012     NA                          
    ## NOTE: confidence intervals computed under the null hypothesis
    
    p.value <- c(mean(abs(BT.perm@Delta.netBenefit[1]) < abs(BT.perm@DeltaResampling.netBenefit[1,])),
                 mean(abs(BT.perm@Delta.netBenefit[2]) < abs(BT.perm@DeltaResampling.netBenefit[2,])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] < BT.perm@DeltaResampling.netBenefit[1,]),
                 mean(BT.perm@Delta.netBenefit[2] < BT.perm@DeltaResampling.netBenefit[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] > BT.perm@DeltaResampling.netBenefit[1,]),
                 mean(BT.perm@Delta.netBenefit[2] > BT.perm@DeltaResampling.netBenefit[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){
        dt.perm <- copy(dt.sim)
        setkeyv(dt.perm, cols = c("strata","Treatment"))

        dt.perm[, Treatment := Treatment[sample.int(.N, size = .N, replace = FALSE)] ]
        expect_equal(table(dt.perm$Treatment), table(dt.sim$Treatment))

        iBT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                              data = dt.perm, method.tte = method,
                              method.inference = "none")

        expect_equal(as.double(iBT.perm@delta.netBenefit),
                     as.double(BT.perm@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.perm@delta.winRatio),
                     as.double(BT.perm@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.perm@Delta.netBenefit),
                     as.double(BT.perm@DeltaResampling.netBenefit[,iResample]))
        expect_equal(as.double(iBT.perm@Delta.winRatio),
                     as.double(BT.perm@DeltaResampling.winRatio[,iResample]))
    }

    ## ** with more samples
    if(FALSE){
        BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                        data = dt.sim, method.tte = "Peron", trace = 3,
                        cpus = 4, 
                        method.inference = "permutation", n.resampling = 1000)

        ## summary(BT)       
        ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
        ## > null hypothesis : Delta == 0 
        ## > permutation test: 1000 samples, confidence level 0.95 
        ## > treatment groups: 0 (control) vs. 1 (treatment) 
        ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

        ## > results
        ##   endpoint threshold total favorable unfavorable neutral uninf  delta  Delta CI [2.5 ; 97.5] p.value 
        ## eventtime1     1e-12   100     48.57       51.43       0     0 -0.029 -0.029  [-0.231;0.164]   0.789 
       
    }
})


## * Stratified permutation
test_that("stratified permutation", {
    BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                    data = dt.sim, method.tte = method, seed = 10, 
                    method.inference = "stratified permutation", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "two.sided")
 ##       endpoint threshold strata  total favorable unfavorable neutral uninf  delta Delta CI [2.5 ; 97.5] p.value 
 ## eventtime1     1e-12 global 100.00     50.87       46.92       0   2.2  0.039 0.039   [-0.353;0.31]     0.9 
 ##                           a  31.50     11.81       19.69       0   0.0 -0.079                               
 ##                           b  35.43     20.79       14.64       0   0.0  0.062                               
 ##                           c  33.07     18.27       12.60       0   2.2  0.057                               
 ##  toxicity1       0.5 global   2.20      2.20        0.00       0   0.0  0.022 0.062  [-0.324;0.351]     0.9 
 ##                           a   0.00      0.00        0.00       0   0.0  0.000                               
 ##                           b   0.00      0.00        0.00       0   0.0  0.000                               
 ##                           c   2.20      2.20        0.00       0   0.0  0.022

    p.value <- c(mean(abs(BT.perm@Delta.netBenefit[1]) < abs(BT.perm@DeltaResampling.netBenefit[1,])),
                 mean(abs(BT.perm@Delta.netBenefit[2]) < abs(BT.perm@DeltaResampling.netBenefit[2,])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] < BT.perm@DeltaResampling.netBenefit[1,]),
                 mean(BT.perm@Delta.netBenefit[2] < BT.perm@DeltaResampling.netBenefit[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

       ## ** summary (less)
    outSummary <- summary(BT.perm, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT.perm@Delta.netBenefit[1] > BT.perm@DeltaResampling.netBenefit[1,]),
                 mean(BT.perm@Delta.netBenefit[2] > BT.perm@DeltaResampling.netBenefit[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1 
        dt.perm <- copy(dt.sim)
        setkeyv(dt.perm, cols = c("strata","Treatment"))
        
        dt.perm[, Treatment := Treatment[sample.int(.N, size = .N, replace = FALSE)], by = "strata"]

        iBT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.perm, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(iBT.perm@delta.netBenefit),
                     as.double(BT.perm@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.perm@delta.winRatio),
                     as.double(BT.perm@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.perm@Delta.netBenefit),
                     as.double(BT.perm@DeltaResampling.netBenefit[,iResample]))
        expect_equal(as.double(iBT.perm@Delta.winRatio),
                     as.double(BT.perm@DeltaResampling.winRatio[,iResample]))
    }

    ## ** with more samples
    if(FALSE){
        BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                        data = dt.sim, method.tte = "Peron", trace = 3,
                        cpus = 4, 
                        method.inference = "stratified permutation", n.resampling = 1000)

        ## summary(BT)       
        ##                Generalized pairwise comparison with 1 prioritized endpoint

        ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
        ## > null hypothesis : Delta == 0 
        ## > permutation test: 1000 samples, confidence level 0.95 
        ## > treatment groups: 0 (control) vs. 1 (treatment) 
        ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

        ## > results
        ##   endpoint threshold total favorable unfavorable neutral uninf  delta  Delta CI [2.5 ; 97.5] p.value 
        ## eventtime1     1e-12   100     48.57       51.43       0     0 -0.029 -0.029   [-0.24;0.163]   0.773
        
       
    }

})

## * Bootstrap
test_that("Bootstrap", {

    BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1)  + bin(toxicity1) + strata,
                         data = dt.sim, method.tte = method, seed = 10, 
                         method.inference = "bootstrap", n.resampling = 20)
    tol.boot <- 1.1/NCOL(BT.boot@DeltaResampling.netBenefit)  ## ok to have a difference of 1 unit

    ## method <- "Peron"
    ## BT.boot@deltaResampling.netBenefit
    ## ** summary (two.sided)
    outSummary <- summary(BT.boot, print = FALSE, alternative = "two.sided")
    ## summary(BT.boot)
    ##          Generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > bootstrap resampling: 50 samples, confidence level 0.95 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group

    ## > results
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf   delta  Delta  CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12 global 100.00     48.83       47.77    0.00  3.40  0.0106 0.0106 [-0.1829;0.2106]    0.68 
    ##                           a  33.56     18.07       15.50    0.00  0.00  0.0257     NA                          
    ##                           b  34.74     16.92       17.03    0.00  0.79 -0.0011     NA                          
    ##                           c  31.70     13.84       15.24    0.00  2.62 -0.0140     NA                          
    ##  toxicity1       0.5 global   3.40      0.61        0.86    1.94  0.00 -0.0025 0.0082  [-0.1857;0.213]    0.66 
    ##                           a   0.00      0.00        0.00    0.00  0.00  0.0000     NA                          
    ##                           b   0.79      0.11        0.24    0.43  0.00 -0.0013     NA                          
    ##                           c   2.62      0.50        0.61    1.51  0.00 -0.0012     NA                          

    ## ** check bootsrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.boot <- copy(dt.sim)
        setkeyv(dt.boot, cols = c("strata","Treatment"))
        dt.boot <- dt.boot[sample.int(.N, size = .N, replace = TRUE)]
        
        ## BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
        iBT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                              data = dt.boot, method.tte = method,
                              method.inference = "none")

        expect_equal(as.double(iBT.boot@delta.netBenefit),
                     as.double(BT.boot@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(iBT.boot@delta.winRatio),
                     as.double(BT.boot@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(iBT.boot@Delta.netBenefit),
                     as.double(BT.boot@DeltaResampling.netBenefit[,iResample]))
        expect_equal(as.double(iBT.boot@Delta.winRatio),
                     as.double(BT.boot@DeltaResampling.winRatio[,iResample]))
    }
})


## * Stratified bootstrap
test_that("Stratified bootstrap", {
    ## BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1)  + bin(toxicity1) + strata,
                    data = dt.sim, method.tte = method, seed = 10, 
                    method.inference = "stratified bootstrap", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT, print = FALSE, alternative = "two.sided")

    ## summary(BT)
    ##            Generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > stratified bootstrap resampling: 10 samples, confidence level 0.95 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

    ## > results
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf  delta Delta CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12 global 100.00     48.83       47.77    0.00  3.40  0.011 0.011  [-0.133;0.148]     0.4 
    ##                           a  33.56     18.07       15.50    0.00  0.00  0.026                               
    ##                           b  34.74     16.92       17.03    0.00  0.79 -0.001                               
    ##                           c  31.70     13.84       15.24    0.00  2.62 -0.014                               
    ##  toxicity1       0.5 global   3.40      0.61        0.86    1.94  0.00 -0.002 0.008  [-0.135;0.166]     0.4 
    ##                           a   0.00      0.00        0.00    0.00  0.00  0.000                               
    ##                           b   0.79      0.11        0.24    0.43  0.00 -0.001                               
    ##                           c   2.62      0.50        0.61    1.51  0.00 -0.001     

    
    ## ** check bootsrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.boot <- copy(dt.sim)
        setkeyv(dt.boot, cols = c("strata","Treatment"))
        dt.boot <- dt.boot[, .SD[sample.int(.N, size = .N, replace = TRUE)], by = "strata"]
        
        ## BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
        BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.boot, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(BT.boot@delta.netBenefit),
                     as.double(BT@deltaResampling.netBenefit[,,iResample]))
        expect_equal(as.double(BT.boot@delta.winRatio),
                     as.double(BT@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(BT.boot@Delta.netBenefit),
                     as.double(BT@DeltaResampling.netBenefit[,iResample]))
        expect_equal(as.double(BT.boot@Delta.winRatio),
                     as.double(BT@DeltaResampling.winRatio[,iResample]))
    }
})
