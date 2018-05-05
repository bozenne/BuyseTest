### test-BuyseTest-resampling.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:34) 
## Version: 
## last-updated: maj  5 2018 (22:59) 
##           By: Brice Ozenne
##     Update #: 22
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
n.patients <- 10

## * Simulate data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients,
                       n.C = n.patients,
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))

## * Bootstrap
method <- "Peron"
test_that("permutation", {
    BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                    data = dt.sim, method.tte = method, seed = 10,
                    method.inference = "permutation", n.resampling = 10)

    ## set.seed(10)
    ## dt.perm <- copy(dt.sim)
    ## indexT <- sample.int(NROW(dt.sim), size = sum(dt.sim$Treatment==1), replace = FALSE)
    ## dt.perm[indexT, Treatment := 1]
    ## dt.perm[-indexT, Treatment := 0]
    ## BT.perm <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                         ## data = dt.perm, method.tte = method, n.resampling = 0)
    ## summary(BT)
    ## new
    ## endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
    ## eventtime1     1e-12   100      48.6       23.86       0 27.54 0.247 0.247    [-0.3;0.251]     0.4

    ## old
    ## endpoint threshold pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf
    ## eventtime1     1e-12      100         48.6          23.86          0    27.54
    ## delta Delta CIinf.Delta CIsup.Delta n.permutation p.value 
    ## 0.247 0.247      -0.537       0.495            10     0.4 

    if(FALSE){
        BT <- BuyseTest(Treatment ~ tte(eventtime1, 0, status1),
                    data = dt.sim, method.tte = "Peron", trace = 3,
                    cpus = 4, 
                    method.inference = "permutation", n.resampling = 1000)

        ## new
        ##         Generalized pairwise comparison with 1 prioritized endpoint

        ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
        ## > null hypothesis : Delta == 0 
        ## > permutation test: 1000 samples, confidence level 0.95 
        ## > groups          : 0 (control) vs. 1 (treatment) 
        ## > results
        ## endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
        ## eventtime1     1e-12   100      48.6       23.86       0 27.54 0.247 0.247  [-0.321;0.769]   0.348 

        ## old
        ## Generalized pairwise comparison with 1 prioritized endpoint

        ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
        ## > null hypothesis : Delta == 0 
        ## > permutation test: 1000 samples, confidence level 0.95 
        ## > groups          : 0(control) vs. 1(treatment) 
        ## > results
        ## endpoint threshold total favorable unfavorable neutral uninf delta Delta
        ## eventtime1     1e-12   100      48.6       23.86       0 27.54 0.247 0.247
        ## CI [2.5 ; 97.5] p.value 
        ## [-0.284;0.765]   0.377 
    }
})

#----------------------------------------------------------------------
### test-BuyseTest-resampling.R ends here
