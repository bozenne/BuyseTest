### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: apr 30 2018 (15:22) 
##           By: Brice Ozenne
##     Update #: 16
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
    BTS.T <- as.data.table(summary(BT.T, show = FALSE)$table)
    
    BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = FALSE)
    BTS.F <- as.data.table(summary(BT.F, show = FALSE)$table)

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
                     "CIinf.Delta" = c(NA, NA, NA, NA), 
                     "CIsup.Delta" = c(NA, NA, NA, NA), 
                     "p.value" = c(NA, NA, NA, NA), 
                     "n.resampling" = c(NA, NA, NA, NA))
    ##    butils::object2script(BTS.T)
    expect_equal(BTS.T, GS, tol = 1e-6)
})


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
##                        method = 0,
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


