### test-BuyseTest-correctionTTE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 30 2018 (23:45) 
## Version: 
## Last-Updated: maj  5 2018 (22:57) 
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

context("Check that the argument correctTTE in BuyseTest is working correctly \n")

## * settings
BuyseTest.options(check = FALSE,
                  keep.comparison = TRUE,
                  method.inference = "none",
                  trace = 0)
## * 1 endpoint
df <- data.frame("survie" = c(2.1, 4.1, 6.1, 8.1, 4, 6, 8, 10),
                 "event" = c(1, 1, 1, 0, 1, 0, 0, 1),
                 "group" = c(0, 0, 0, 0, 1, 1, 1, 1),
                 "score" = 1)

## ** Gehan
test_that("1 endpoint - Gehan", {
    Gehan <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                       data = df, 
                       method.tte = "Gehan")

    expect_equal(as.double(Gehan@count.favorable), c(9,0))
    expect_equal(as.double(Gehan@count.unfavorable), c(2,0))
    expect_equal(as.double(Gehan@count.neutral), c(1,5))
    expect_equal(as.double(Gehan@count.uninf), c(4,0))

    GehanC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Gehan", correctionTTE = TRUE)

    factor <- 16/12 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(GehanC@count.favorable), c(9*factor,0))
    expect_equal(as.double(GehanC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(GehanC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(GehanC@count.uninf), c(0,0))
})

## ** Peron
test_that("1 endpoint - Peron", {
    Peron <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                       data = df, 
                       method.tte = "Peron")

    expect_equal(as.double(Peron@count.favorable), c(10,0))
    expect_equal(as.double(Peron@count.unfavorable), c(2,0))
    expect_equal(as.double(Peron@count.neutral), c(1,4))
    expect_equal(as.double(Peron@count.uninf), c(3,0))

    PeronC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Peron", correctionTTE = TRUE)

    factor <- 16/13 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(PeronC@count.favorable), c(10*factor,0))
    expect_equal(as.double(PeronC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(PeronC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(PeronC@count.uninf), c(0,0))
})

##----------------------------------------------------------------------
### test-BuyseTest-correctionTTE.R ends here
