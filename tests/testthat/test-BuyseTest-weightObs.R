### test-BuyseTest-weightObs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Mar 20 2023 (12:17) 
## Version: 
## Last-Updated: jul 17 2023 (14:24) 
##           By: Brice Ozenne
##     Update #: 6
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
    library(data.table)
    library(survival)
}

context("Check calculations with weighted observations")

## * same weight for all
veteran$weight <- 2
veteran2 <- rbind(veteran,veteran)

test_that("same weight for all",{
    e.wBT <- BuyseTest(trt ~ bin(prior) + cont(time) + celltype, data = veteran, trace = FALSE,
                       method.inference = "u-statistic", weightObs = "weight")
    e.BT2 <- BuyseTest(trt ~ bin(prior) + cont(time) + celltype, data = veteran2, trace = FALSE,
                       method.inference = "u-statistic")

    expect_equal(unlist(confint(e.wBT)), unlist(confint(e.BT2)), tol = 1e-6)
    expect_equal(unlist(confint(e.wBT, strata = TRUE)), unlist(confint(e.BT2, strata = TRUE)), tol = 1e-6)

    ## note when using time to event endpoint, there might not be equality 
    ## because the weights are only used for the GPC, not when fitting the survival model
})

## * different weight per individual
veteran$weight <- 1
veteran[1:10,"weight"] <- 2
veteran[11:20,"weight"] <- 3
veteran[21:30,"weight"] <- 4
veteran2 <- rbind(veteran,veteran[1:30,],veteran[11:30,],veteran[21:30,])


test_that("individual specific weight",{
    e.wBT <- BuyseTest(trt ~ bin(prior) + cont(time) + celltype, data = veteran, trace = FALSE,
                       method.inference = "u-statistic", weightObs = "weight")
    e.BT2 <- BuyseTest(trt ~ bin(prior) + cont(time) + celltype, data = veteran2, trace = FALSE,
                       method.inference = "u-statistic")

    expect_equal(unlist(confint(e.wBT)), unlist(confint(e.BT2)), tol = 1e-6)
    expect_equal(unlist(confint(e.wBT, strata = TRUE)), unlist(confint(e.BT2, strata = TRUE)), tol = 1e-6)

    ## note when using time to event endpoint, there might not be equality 
    ## because the weights are only used for the GPC, not when fitting the survival model
})
##----------------------------------------------------------------------
### test-BuyseTest-weightObs.R ends here
