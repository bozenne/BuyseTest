### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: apr 17 2018 (16:56) 
##           By: Brice Ozenne
##     Update #: 9
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
BuyseTest.options(n.permutation = 0, trace = 0, keep.comparison = TRUE)

## * Joris: jeudi 5 avril 2018 à 14:57
dt.sim <- data.table(
    ttt=c(rep(0,3),rep(1,3)),
    timeOS = c(10,20,30,15,20,35),
    eventOS= c(1,1,0,0,1,1),
    moins.grade.tox = -c(1,2,3,2,4,2)
)

test_that("number of pairs - argument neutral.as.uninf", {
    BT.T <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(moins.grade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = TRUE)
    BTS.T <- as.data.table(summary(BT.T, show = FALSE))
    
    BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(moins.grade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = FALSE)
    BTS.F <- as.data.table(summary(BT.F, show = FALSE))

    expect_equal(BTS.T[1,],BTS.F[1,])
    expect_equal(BTS.T[2,pc.total],BTS.T[1,pc.neutral+pc.uninf])
    expect_equal(BTS.F[2,pc.total],BTS.F[1,pc.uninf])   
})

######################################################################
### test-BuyseTest-previousBug.R ends here
