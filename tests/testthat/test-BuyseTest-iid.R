### test-BuyseTest-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  8 2019 (11:54) 
## Version: 
## Last-Updated: mar 28 2019 (14:31) 
##           By: Brice Ozenne
##     Update #: 32
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

context("Check correct computation of the variance \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  order.Hprojection = 1,
                  trace = 0)

## * Binary case
## ** no strata
## equal number in each group
d <- data.table(id = 1:4, group = c("C","C","T","T"), toxicity = c(1,0,1,0))

test_that("iid: binary and no strata (balanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity),
                      data = d, 
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(1/16,1/16,-1/16, 1/4, 4) )

    expect_equal(iid(e.BT)[,"favorable"], c(-1/8,1/8,1/8,-1/8))
    expect_equal(iid(e.BT)[,"unfavorable"], c(1/8,-1/8,-1/8,1/8))
})

## unequal number in each group
d.bis <- data.table(id = 1:4, group = c("C","T","T","T"), toxicity = c(1,1,1,0))

test_that("iid: binary and no strata (unbalanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity),
                      data = d.bis, 
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0,2/27,0) )

    expect_equal(iid(e.BT)[,"favorable"], c(0,0,0,0))
    expect_equal(iid(e.BT)[,"unfavorable"], c(0,-1/9,-1/9,2/9))
})


## ** strata
d2 <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2),
            cbind(d, strata = 3))

test_that("iid: binary with strata (balanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                      data = d2, 
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(1/16,1/16,-1/16,1/4,4)/3 )

    expect_equal(iid(e.BT)[,"favorable"], rep(c(-1/8,1/8,1/8,-1/8),3)/3)
    expect_equal(iid(e.BT)[,"unfavorable"], rep(c(1/8,-1/8,-1/8,1/8),3)/3)
})

d2.bis <- rbind(cbind(d.bis, strata = 1),
                cbind(d.bis, strata = 2),
                cbind(d.bis, strata = 3))

test_that("iid: binary and no strata (unbalanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                      data = d2.bis, 
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0,2/27,0)/3 )

    expect_equal(iid(e.BT)[,"favorable"], rep(c(0,0,0,0),3)/3)
    expect_equal(iid(e.BT)[,"unfavorable"], rep(c(0,-1/9,-1/9,2/9),3)/3)
})

## * Two endpoints
## ** no strata
set.seed(10)
d <- simBuyseTest(50)

test_that("iid: two endpoints (no strata)", {
    ## different endpoints
    e.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["toxicity_0.5",]),
                 c(0.002499994, 0.002499994, -0.002492006, 0.009984000, 0.160256410), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score_1",]),
                 c(0.003049562, 0.003202234, -0.002925978, 0.012103750, 0.077787351), tol = 1e-6 )

    ## favorable unfavorable   covariance
    ## toxicity 0.002499994 0.002499994 -0.002492006
    ## score    0.008041562 0.008194234 -0.002925978

    ## same endpoint
    e.BT <- BuyseTest(Treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["score_2",]),
                 c(0.00036192, 0.000613760, -0.000166400, 0.00130848, 0.13086775), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score_1",]),
                 c(0.00190759, 0.002360218, -0.001708275, 0.007684358, 0.088893573), tol = 1e-6 )

    ## same endpoint tte
    e.BT <- suppressWarnings(BuyseTest(Treatment ~  tte(eventtime, threshold = 2, censoring = status) + tte(eventtime, threshold = 1, censoring = status),
                                       data = d,
                                       method.inference = "u-statistic"))
    e2.BT <- suppressWarnings(BuyseTest(Treatment ~  tte(eventtime, threshold = 2, censoring = status) + tte(eventtime, threshold = 1, censoring = status),
                                        data = d, keep.pairScore = TRUE,
                                        method.inference = "u-statistic-bebu"))
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["eventtime_2",]),
                 c(0.0000000000, 0.0004251505, 0.0000000000, 0.0004251505, 0.0000000000), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["eventtime_1",]),
                 c(0.0008748297, 0.0015243380, -0.0009593909, 0.0043179495, 0.0171570309), tol = 1e-6 )
})


## ** strata
d2 <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2))
d2$score1 <- d2$score
d2[strata == 1, score1 := 1]
d2$score2 <- d2$score
d2[strata == 2, score2 := 1]

test_that("iid: two endpoints (strata)", {

    e0.BT <- BuyseTest(Treatment ~ cont(score, threshold = 1) + strata,
                       data = d2,
                       method.inference = "u-statistic")

    e.BT <- BuyseTest(Treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2,
                      method.inference = "u-statistic")

    e2.BT <- BuyseTest(Treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(as.double(e0.BT@covariance), as.double(e.BT@covariance[2,]))
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e0.BT@covariance), c(0.0009537952, 0.001180109, -0.0008541376, 0.0038421792, 0.0444467864), tol = 1e-6 )
})

######################################################################
### test-BuyseTest-iid.R ends here
