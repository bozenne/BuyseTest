### test-BuyseTest-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  8 2019 (11:54) 
## Version: 
## Last-Updated: jan  8 2019 (15:08) 
##           By: Brice Ozenne
##     Update #: 10
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
                      method.inference = "asymptotic")
    e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d, 
                       method.inference = "asymptotic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(1/16,1/16,-1/16) )

    expect_equal(iid(e.BT)[,"favorable"], c(-1/8,1/8,1/8,-1/8))
    expect_equal(iid(e.BT)[,"unfavorable"], c(1/8,-1/8,-1/8,1/8))
})

## unequal number in each group
d.bis <- data.table(id = 1:4, group = c("C","T","T","T"), toxicity = c(1,1,1,0))

test_that("iid: binary and no strata (unbalanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity),
                      data = d.bis, 
                      method.inference = "asymptotic")
    e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d.bis, 
                       method.inference = "asymptotic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0) )

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
                      method.inference = "asymptotic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2, 
                       method.inference = "asymptotic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(1/16,1/16,-1/16)/3 )

    expect_equal(iid(e.BT)[,"favorable"], rep(c(-1/8,1/8,1/8,-1/8),3)/3)
    expect_equal(iid(e.BT)[,"unfavorable"], rep(c(1/8,-1/8,-1/8,1/8),3)/3)
})

d2.bis <- rbind(cbind(d.bis, strata = 1),
                cbind(d.bis, strata = 2),
                cbind(d.bis, strata = 3))

test_that("iid: binary and no strata (unbalanced groups)", {
    e.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                      data = d2.bis, 
                      method.inference = "asymptotic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2.bis, 
                       method.inference = "asymptotic-bebu")
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0)/3 )

    expect_equal(iid(e.BT)[,"favorable"], rep(c(0,0,0,0),3)/3)
    expect_equal(iid(e.BT)[,"unfavorable"], rep(c(0,-1/9,-1/9,2/9),3)/3)
})

######################################################################
### test-BuyseTest-iid.R ends here

## * Two endpoints
set.seed(10)
d <- simBuyseTest(50)


test_that("iid: two endpoints", {
    e.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d,
                      method.inference = "asymptotic")
    e2.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                       data = d,
                       method.inference = "asymptotic-bebu")
   
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["toxicity",]), c(0.002499994, 0.002499994, -0.002492006), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score",]), c(0.008041562, 0.008194234, -0.002925978), tol = 1e-6 )

})


