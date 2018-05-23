### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: maj 23 2018 (13:50) 
##           By: Brice Ozenne
##     Update #: 31
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
}

context("Check that bugs that have been reported are fixed \n")


## * settings
BuyseTest.options(check = FALSE,
                  keep.comparison = TRUE,
                  method.inference = "none",
                  trace = 0)

## * Joris: jeudi 5 avril 2018 à 14:57
dt.sim <- data.table(
    ttt=c(rep(0,3),rep(1,3)),
    timeOS = c(10,20,30,15,20,35),
    eventOS = c(1,1,0,0,1,1),
    Mgrade.tox = -c(1,2,3,2,4,2)
)

test_that("number of pairs - argument neutral.as.uninf", {
    BT.T <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = TRUE)
    BTS.T <- as.data.table(summary(BT.T, print = FALSE)$table)
    
    BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = FALSE)
    BTS.F <- as.data.table(summary(BT.F, print = FALSE)$table)

    expect_equal(BTS.T[1,],BTS.F[1,])
    expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", pc.total],
                 BTS.T[endpoint == "timeOS" & strata == "global", pc.neutral+pc.uninf])
    expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", pc.total],
                 BTS.F[endpoint == "timeOS" & strata == "global", pc.uninf])

    GS <- data.table("endpoint" = c("timeOS", "timeOS", "Mgrade.tox", "Mgrade.tox"), 
                     "threshold" = c(1e-12, 1e-12, 1e-12, 1e-12), 
                     "strata" = c("global", "1", "global", "1"), 
                     "pc.total" = c(100.00000, 100.00000,  33.33333,  33.33333), 
                     "pc.favorable" = c(50.00000, 50.00000, 16.66667, 16.66667), 
                     "pc.unfavorable" = c(16.66667, 16.66667, 11.11111, 11.11111), 
                     "pc.neutral" = c(11.111111, 11.111111,  5.555556,  5.555556), 
                     "pc.uninf" = c(22.22222, 22.22222,  0.00000,  0.00000), 
                     "delta" = c(0.33333333, 0.33333333, 0.05555556, 0.05555556), 
                     "Delta" = c(0.3333333, NA, 0.3888889, NA), 
                     "CIinf.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                     "CIsup.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                     "p.value" = as.numeric(c(NA, NA, NA, NA)),
                     "n.resampling" = as.numeric(c(NA, NA, NA, NA)))
    ##    butils::object2script(BTS.T)

    attr(BTS.T,"index") <- NULL
    expect_equal(BTS.T, GS, tol = 1e-6)
    ## class(BTS.T[["n.resampling"]])
    ## class(GS[["n.resampling"]])
    
})

## * Joris: jeudi 5 avril 2018 à 14:57
if(FALSE){
    system.time(
        BT_Gehan <- BuyseTest(treatment = "trt",
                              strata = "celltype",
                              type = "TTE",
                              endpoint = "time",
                              threshold = 0,
                              censoring = "status", 
                              data=veteran, method.tte = "Gehan",
                              n.resampling = 1e4)
    )
    1
    confint(BT_Gehan)
    dt.sim <- data.table(
        ttt = c(0,0,0,0,0,1,1,1,1,1),
        y = c(1,2,3,4,5,2,3,4,5,6),
        end = c(1,1,1,1,1,1,1,1,1,1),
        time = c(10,11,12,13,14,11,12,13,14,15)
    )

    BT <- BuyseTest(ttt ~ cont(y,threshold=4) + TTE(time,censoring=end) + cont(y,threshold=3),
                    data = dt.sim, method.inference = "permutation", method.tte = "Peron", n.resampling = 1000)
    summary(BT, statistic = "netChance")
 ##        Generalized pairwise comparison with 3 prioritized endpoints

 ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
 ## > null hypothesis : Delta == 0 
 ## > permutation test: 1000 samples, confidence level 0.95 
 ## > treatment groups: 0 (control) vs. 1 (treatment) 
 ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

 ## > results
 ## endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
 ##        y         4   100        12           0      88     0  0.12  0.12   [-0.04;0.281]   0.184 
 ##     time     1e-12    88        48          24      16     0  0.24  0.36    [-0.32;1.04]   0.297 
 ##        y         3    16         0           0      16     0  0.00  0.36    [-0.32;1.04]   0.297 

  xx <- summary(BT,statistic="winRatio")
## Generalized pairwise comparison with 3 prioritized endpoints

##  > statistic       : win ratio (delta: endpoint specific, Delta: global) 
##  > null hypothesis : Delta == 1 
##  > permutation test: 1000 samples, confidence level 0.95 
##  > treatment groups: 0 (control) vs. 1 (treatment) 
##  > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

##  > results
##  endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
##         y         4   100        12           0      88     0   Inf   Inf                         
##      time     1e-12    88        48          24      16     0     2   2.5    [1.65;8.167]   0.118 
##         y         3    16         0           0      16     0         2.5    [1.65;8.167]   0.118
}

## ** table Comparison
## $y_4
##    strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1       1       6       1                   1                   1         0           0       1             0
## 2       1       6       2                   1                   2         0           0       1             0
## 3       1       6       3                   1                   3         0           0       1             0
## 4       1       6       4                   1                   4         0           0       1             0
## 5       1       6       5                   1                   5         0           0       1             0
## 6       1       7       1                   2                   1         0           0       1             0
## 7       1       7       2                   2                   2         0           0       1             0
## 8       1       7       3                   2                   3         0           0       1             0
## 9       1       7       4                   2                   4         0           0       1             0
## 10      1       7       5                   2                   5         0           0       1             0
## 11      1       8       1                   3                   1         0           0       1             0
## 12      1       8       2                   3                   2         0           0       1             0
## 13      1       8       3                   3                   3         0           0       1             0
## 14      1       8       4                   3                   4         0           0       1             0
## 15      1       8       5                   3                   5         0           0       1             0
## 16      1       9       1                   4                   1         1           0       0             0
## 17      1       9       2                   4                   2         0           0       1             0
## 18      1       9       3                   4                   3         0           0       1             0
## 19      1       9       4                   4                   4         0           0       1             0
## 20      1       9       5                   4                   5         0           0       1             0
## 21      1      10       1                   5                   1         1           0       0             0
## 22      1      10       2                   5                   2         1           0       0             0
## 23      1      10       3                   5                   3         0           0       1             0
## 24      1      10       4                   5                   4         0           0       1             0
## 25      1      10       5                   5                   5         0           0       1             0

## $`time_1e-12`
##    strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1       1       6       1                   1                   1      1.00        0.00       0          0.00
## 2       1       6       2                   1                   2      1.00        0.00       0          0.00
## 3       1       6       3                   1                   3      0.75        0.00       0          0.25
## 4       1       6       4                   1                   4      0.50        0.25       0          0.25
## 5       1       6       5                   1                   5      0.25        0.50       0          0.25
## 6       1       7       1                   2                   1      1.00        0.00       0          0.00
## 7       1       7       2                   2                   2      1.00        0.00       0          0.00
## 8       1       7       3                   2                   3      0.00        0.00       1          0.00
## 9       1       7       4                   2                   4      0.00        1.00       0          0.00
## 10      1       7       5                   2                   5      0.00        1.00       0          0.00
## 11      1       8       1                   3                   1      1.00        0.00       0          0.00
## 12      1       8       2                   3                   2      1.00        0.00       0          0.00
## 13      1       8       3                   3                   3      1.00        0.00       0          0.00
## 14      1       8       4                   3                   4      0.00        0.00       1          0.00
## 15      1       8       5                   3                   5      0.00        1.00       0          0.00
## 16      1       9       2                   4                   2      1.00        0.00       0          0.00
## 17      1       9       3                   4                   3      1.00        0.00       0          0.00
## 18      1       9       4                   4                   4      1.00        0.00       0          0.00
## 19      1       9       5                   4                   5      0.00        0.00       1          0.00
## 20      1      10       3                   5                   3      1.00        0.00       0          0.00
## 21      1      10       4                   5                   4      1.00        0.00       0          0.00
## 22      1      10       5                   5                   5      1.00        0.00       0          0.00

## $y_3
##   strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1      1       7       3                   2                   3         0        0.00    1.00             0
## 2      1       8       4                   3                   4         0        0.00    1.00             0
## 3      1       9       5                   4                   5         0        0.00    1.00             0
## 4      1       6       3                   1                   3         0        0.00    0.25             0
## 5      1       6       4                   1                   4         0        0.00    0.25             0
## 6      1       6       5                   1                   5         0        0.25    0.00             0


## * running time (hidden)
## tps <- system.time(
##     BT.T <- .BuyseTest(alternative = "two.sided",
##                        censoring = NULL,
##                        correctionTTE = FALSE,
##                        cpus = 1,
##                        data = dt.sim,
##                        endpoint = c("timeOS","Mgrade.tox"),
##                        index.survivalM1 = numeric(0),
##                        keep.comparison = FALSE,
##                        method.tte = 0,
##                        method.inference = "none",
##                        neutral.as.uninf = TRUE,
##                        operator = c(">0",">0"),
##                        seed = 10,
##                        strata = NULL,
##                        threshold = c(1e-12,1e-12),
##                        threshold.TTEM1 = numeric(0),
##                        trace = 0,
##                        treatment = "ttt",
##                        type = c(1,1),
##                        Wscheme = matrix(nrow = 0, ncol = 0)
##                        )
## )    
## tps  

######################################################################
### test-BuyseTest-previousBug.R ends here


