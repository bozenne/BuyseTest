### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: maj 26 2018 (17:20) 
##           By: Brice Ozenne
##     Update #: 67
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
BuyseTest.options(check = FALSE,
                  keep.comparison = TRUE,                  
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
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                    data = dt.sim, method.tte = method, seed = 10, 
                    method.inference = "permutation", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT, print = FALSE, alternative = "two.sided")
    ## summary(BT, alternative = "two.sided")
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf  delta Delta CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12 global 100.00     48.83       47.77    0.00  3.40  0.011 0.011   [-0.099;0.33]       1 
    ##                           a  33.56     18.07       15.50    0.00  0.00  0.026                               
    ##                           b  34.74     16.92       17.03    0.00  0.79 -0.001                               
    ##                           c  31.70     13.84       15.24    0.00  2.62 -0.014                               
    ##  toxicity1       0.5 global   3.40      0.61        0.86    1.94  0.00 -0.002 0.008  [-0.104;0.335]       1 
    ##                           a   0.00      0.00        0.00    0.00  0.00  0.000                               
    ##                           b   0.79      0.11        0.24    0.43  0.00 -0.001                               
    ##                           c   2.62      0.50        0.61    1.51  0.00 -0.001      
    

    p.value <- c(mean(abs(BT@Delta.netChance[1]) < abs(BT@DeltaResampling.netChance[1,])),
                 mean(abs(BT@Delta.netChance[2]) < abs(BT@DeltaResampling.netChance[2,])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT@Delta.netChance[1] < BT@DeltaResampling.netChance[1,]),
                 mean(BT@Delta.netChance[2] < BT@DeltaResampling.netChance[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT@Delta.netChance[1] > BT@DeltaResampling.netChance[1,]),
                 mean(BT@Delta.netChance[2] > BT@DeltaResampling.netChance[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){
        dt.perm <- copy(dt.sim)
        dt.perm[, Treatment := Treatment[sample.int(.N, size = .N, replace = FALSE)] ]
        expect_equal(table(dt.perm$Treatment), table(dt.sim$Treatment))

        BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.perm, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(BT.perm@delta.netChance),
                     as.double(BT@deltaResampling.netChance[,,iResample]))
        expect_equal(as.double(BT.perm@delta.winRatio),
                     as.double(BT@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(BT.perm@Delta.netChance),
                     as.double(BT@DeltaResampling.netChance[,iResample]))
        expect_equal(as.double(BT.perm@Delta.winRatio),
                     as.double(BT@DeltaResampling.winRatio[,iResample]))
    }

    ## ** with more samples
    if(FALSE){
        BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
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
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                    data = dt.sim, method.tte = method, seed = 10, 
                    method.inference = "stratified permutation", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT, print = FALSE, alternative = "two.sided")
 ##       endpoint threshold strata  total favorable unfavorable neutral uninf  delta Delta CI [2.5 ; 97.5] p.value 
 ## eventtime1     1e-12 global 100.00     50.87       46.92       0   2.2  0.039 0.039   [-0.353;0.31]     0.9 
 ##                           a  31.50     11.81       19.69       0   0.0 -0.079                               
 ##                           b  35.43     20.79       14.64       0   0.0  0.062                               
 ##                           c  33.07     18.27       12.60       0   2.2  0.057                               
 ##  toxicity1       0.5 global   2.20      2.20        0.00       0   0.0  0.022 0.062  [-0.324;0.351]     0.9 
 ##                           a   0.00      0.00        0.00       0   0.0  0.000                               
 ##                           b   0.00      0.00        0.00       0   0.0  0.000                               
 ##                           c   2.20      2.20        0.00       0   0.0  0.022

    p.value <- c(mean(abs(BT@Delta.netChance[1]) < abs(BT@DeltaResampling.netChance[1,])),
                 mean(abs(BT@Delta.netChance[2]) < abs(BT@DeltaResampling.netChance[2,])))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT@Delta.netChance[1] < BT@DeltaResampling.netChance[1,]),
                 mean(BT@Delta.netChance[2] < BT@DeltaResampling.netChance[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

       ## ** summary (less)
    outSummary <- summary(BT, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT@Delta.netChance[1] > BT@DeltaResampling.netChance[1,]),
                 mean(BT@Delta.netChance[2] > BT@DeltaResampling.netChance[2,]))
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check permutation
    set.seed(10)
    for(iResample in 1:2){
        dt.perm <- copy(dt.sim)
        dt.perm[, Treatment := Treatment[sample.int(.N, size = .N, replace = FALSE)], by = "strata"]

        BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.perm, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(BT.perm@delta.netChance),
                     as.double(BT@deltaResampling.netChance[,,iResample]))
        expect_equal(as.double(BT.perm@delta.winRatio),
                     as.double(BT@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(BT.perm@Delta.netChance),
                     as.double(BT@DeltaResampling.netChance[,iResample]))
        expect_equal(as.double(BT.perm@Delta.winRatio),
                     as.double(BT@DeltaResampling.winRatio[,iResample]))
    }

    ## ** with more samples
    if(FALSE){
        BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
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
    ## BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1)  + bin(toxicity1) + strata,
                    data = dt.sim, method.tte = method, seed = 10, 
                    method.inference = "bootstrap", n.resampling = 10)

    ## ** summary (two.sided)
    outSummary <- summary(BT, print = FALSE, alternative = "two.sided")
    ## summary(BT)
    ##            Generalized pairwise comparison with 2 prioritized endpoints and 3 strata

    ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
    ## > null hypothesis : Delta == 0 
    ## > bootstrap resampling: 10 samples, confidence level 0.95 
    ## > treatment groups: 0 (control) vs. 1 (treatment) 
    ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

    ## > results
    ##   endpoint threshold strata  total favorable unfavorable neutral uninf  delta Delta CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12 global 100.00     48.83       47.77    0.00  3.40  0.011 0.011  [-0.043;0.265]     0.9 
    ##                           a  33.56     18.07       15.50    0.00  0.00  0.026                               
    ##                           b  34.74     16.92       17.03    0.00  0.79 -0.001                               
    ##                           c  31.70     13.84       15.24    0.00  2.62 -0.014                               
    ##  toxicity1       0.5 global   3.40      0.61        0.86    1.94  0.00 -0.002 0.008  [-0.045;0.263]     0.9 
    ##                           a   0.00      0.00        0.00    0.00  0.00  0.000                               
    ##                           b   0.79      0.11        0.24    0.43  0.00 -0.001                               
    ##                           c   2.62      0.50        0.61    1.51  0.00 -0.001   
    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] > 0 ), ## >0 because >0 punctual estimate
                 mean(BT@DeltaResampling.netChance[2,] > 0 ) ## >0 because >0 punctual estimate
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] > 0 ), 
                 mean(BT@DeltaResampling.netChance[2,] > 0 ) 
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] < 0 ), 
                 mean(BT@DeltaResampling.netChance[2,] < 0 ) 
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check bootsrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1 
        dt.boot <- dt.sim[sample.int(.N, size = .N, replace = TRUE)]
        
        ## BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
        BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.boot, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(BT.boot@delta.netChance),
                     as.double(BT@deltaResampling.netChance[,,iResample]))
        expect_equal(as.double(BT.boot@delta.winRatio),
                     as.double(BT@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(BT.boot@Delta.netChance),
                     as.double(BT@DeltaResampling.netChance[,iResample]))
        expect_equal(as.double(BT.boot@Delta.winRatio),
                     as.double(BT@DeltaResampling.winRatio[,iResample]))
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

    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] > 0 ), ## >0 because >0 punctual estimate
                 mean(BT@DeltaResampling.netChance[2,] > 0 ) ## >0 because >0 punctual estimate
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (greater)
    outSummary <- summary(BT, print = FALSE, alternative = "greater")
    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] > 0 ), 
                 mean(BT@DeltaResampling.netChance[2,] > 0 ) 
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)

    ## ** summary (less)
    outSummary <- summary(BT, print = FALSE, alternative = "less")
    
    p.value <- c(mean(BT@DeltaResampling.netChance[1,] < 0 ), 
                 mean(BT@DeltaResampling.netChance[2,] < 0 ) 
                 )
    expect_equal(outSummary$table[outSummary$table$strata=="global","p.value"],
                 p.value)
    
    ## ** check bootsrap
    set.seed(10)
    for(iResample in 1:2){ ## iResample <- 1
        dt.boot <- dt.sim[, .SD[sample.int(.N, size = .N, replace = TRUE)], by = "strata"]
        
        ## BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
        BT.boot <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1) + bin(toxicity1) + strata,
                             data = dt.boot, method.tte = method,
                             method.inference = "none")

        expect_equal(as.double(BT.boot@delta.netChance),
                     as.double(BT@deltaResampling.netChance[,,iResample]))
        expect_equal(as.double(BT.boot@delta.winRatio),
                     as.double(BT@deltaResampling.winRatio[,,iResample]))
        expect_equal(as.double(BT.boot@Delta.netChance),
                     as.double(BT@DeltaResampling.netChance[,iResample]))
        expect_equal(as.double(BT.boot@Delta.winRatio),
                     as.double(BT@DeltaResampling.winRatio[,iResample]))
    }
})
