### test-BuyseTest-multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 28 2021 (11:01) 
## Version: 
## Last-Updated: May  1 2023 (10:07) 
##           By: Brice Ozenne
##     Update #: 8
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
}

context("Check BuyseTest when performing multiple comparisons")

## * Settings
BuyseTest.options(check = TRUE,
                  pool.strata = "Buyse",
                  trace = 0)

## * Case 1: multiple univariate analyses

## ** no strata
n <- 100
set.seed(10)
d <- simBuyseTest(n)

test_that("univariate result - no strata",{
    ## check correctly extracted univariate from the multivariate analysis 
    eG.BT <- BuyseTest(treatment ~ bin(toxicity) + cont(score,1) + cont(eventtime), data = d,
                       hierarchical = FALSE, method.inference = "u-statistic")
    eG1.BT <- BuyseTest(treatment ~ bin(toxicity), data = d, method.inference = "u-statistic")
    eG2.BT <- BuyseTest(treatment ~ cont(score,1), data = d, method.inference = "u-statistic")
    eG3.BT <- BuyseTest(treatment ~ cont(eventtime), data = d, method.inference = "u-statistic")

    test <- confint(eG.BT, statistic =  "netBenefit", cumulative = FALSE)
    GS <- rbind(confint(eG1.BT, statistic =  "netBenefit", cumulative = FALSE),
                confint(eG2.BT, statistic =  "netBenefit", cumulative = FALSE),
                confint(eG3.BT, statistic =  "netBenefit", cumulative = FALSE))
    expect_equal(as.double(unlist(test)),as.double(unlist(GS)), tol = 1e-6)

    test <- confint(eG.BT, statistic =  "winRatio", cumulative = FALSE)
    GS <- rbind(confint(eG1.BT, statistic =  "winRatio", cumulative = FALSE),
                confint(eG2.BT, statistic =  "winRatio", cumulative = FALSE),
                confint(eG3.BT, statistic =  "winRatio", cumulative = FALSE))
    expect_equal(as.double(unlist(test)),as.double(unlist(GS)), tol = 1e-6)
    
    xx <- BuyseMultComp(eG.BT, cumulative = FALSE, endpoint = 1:3)    
})

## ** strata
n <- 100
set.seed(10)
dS <- simBuyseTest(n, n.strata = 3)

test_that("Separate BuyseTest vs. single BuyseTest for obtaining univariate result",{
    eG.BT <- BuyseTest(treatment ~ bin(toxicity) + cont(score,1) + cont(eventtime) + strata, data = dS,
                       hierarchical = FALSE, method.inference = "u-statistic")
    eG1.BT <- BuyseTest(treatment ~ bin(toxicity) + strata, data = dS, method.inference = "u-statistic")
    eG2.BT <- BuyseTest(treatment ~ cont(score,1) + strata, data = dS, method.inference = "u-statistic")
    eG3.BT <- BuyseTest(treatment ~ cont(eventtime) + strata, data = dS, method.inference = "u-statistic")

    test <- confint(eG.BT, statistic =  "netBenefit", cumulative = FALSE)
    GS <- rbind(confint(eG1.BT, statistic =  "netBenefit", cumulative = FALSE),
                confint(eG2.BT, statistic =  "netBenefit", cumulative = FALSE),
                confint(eG3.BT, statistic =  "netBenefit", cumulative = FALSE))
    expect_equal(as.double(unlist(test)),as.double(unlist(GS)), tol = 1e-6)

    test <- confint(eG.BT, statistic =  "winRatio", cumulative = FALSE)
    GS <- rbind(confint(eG1.BT, statistic =  "winRatio", cumulative = FALSE),
                confint(eG2.BT, statistic =  "winRatio", cumulative = FALSE),
                confint(eG3.BT, statistic =  "winRatio", cumulative = FALSE))
    expect_equal(as.double(unlist(test)),as.double(unlist(GS)), tol = 1e-6)

    xx <- BuyseMultComp(eG.BT, cumulative = FALSE, endpoint = 1:3)    
})

##----------------------------------------------------------------------
### test-BuyseTest-multcomp.R ends here
