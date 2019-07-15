### test-BuyseTest-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  8 2019 (11:54) 
## Version: 
## Last-Updated: jul 15 2019 (23:23) 
##           By: Brice Ozenne
##     Update #: 48
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
                  keep.pairScore = FALSE,
                  trace = 0)

## * Binary case
## ** no strata
## equal number in each group
d <- data.table(id = 1:4, group = c("C","C","T","T"), toxicity = c(1,0,1,0))

test_that("iid: binary and no strata (balanced groups)", {
    ## first order
    BuyseTest.options(order.Hprojection = 1)
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

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(group ~ bin(toxicity),
                      data = d, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(group ~ bin(toxicity),
                     data = d, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e1.BT@tablePairScore[[1]]$index.pair,1:4) ## correct ordering of the pairs
    expect_equal(e.BT@covariance, e1.BT@covariance) ## assumes vs. does not assumes binary score when computing second order terms
    expect_equal(e.BT@covariance, e2.BT@covariance) 
    expect_equal(as.double(e.BT@covariance), c(5/64,5/64,-3/64, 1/4, 4) )
})

## unequal number in each group
d.bis <- data.table(id = 1:4, group = c("C","T","T","T"), toxicity = c(1,1,1,0))

test_that("iid: binary and no strata (unbalanced groups)", {
    ## first order
    BuyseTest.options(order.Hprojection = 1)
    suppressWarnings(e.BT <- BuyseTest(group ~ bin(toxicity),
                                       data = d.bis, 
                                       method.inference = "u-statistic"))
    suppressWarnings(e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu"))
    
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0,2/27,0) )

    expect_equal(iid(e.BT)[,"favorable"], c(0,0,0,0))
    expect_equal(iid(e.BT)[,"unfavorable"], c(0,-1/9,-1/9,2/9))

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    suppressWarnings(e.BT <- BuyseTest(group ~ bin(toxicity),
                      data = d.bis, keep.pairScore = FALSE,
                      method.inference = "u-statistic"))
    suppressWarnings(e1.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic"))
    suppressWarnings(e2.BT <- BuyseTest(group ~ bin(toxicity),
                       data = d.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu"))
    
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair,1:3) ## correct ordering of the pairs
    expect_equal(e.BT@covariance, e1.BT@covariance) ## assumes vs. does not assumes binary score when computing second order terms
    expect_equal(e.BT@covariance, e2.BT@covariance) 
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0,2/27,0) )
})


## ** strata
d2 <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2),
            cbind(d, strata = 3))

test_that("iid: binary with strata (balanced groups)", {
    ## first order
    BuyseTest.options(order.Hprojection = 1)
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

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                      data = d2, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair,1:12) ## correct ordering of the pairs
    expect_equal(e.BT@covariance, e1.BT@covariance) ## assumes vs. does not assumes binary score when computing second order terms
    expect_equal(e.BT@covariance, e2.BT@covariance)
    
    expect_equal(as.double(e.BT@covariance), c(5/64,5/64,-3/64, 1/4, 4)/3 )
})

d2.bis <- rbind(cbind(d.bis, strata = 1),
                cbind(d.bis, strata = 2),
                cbind(d.bis, strata = 3))

test_that("iid: binary and no strata (unbalanced groups)", {
    ## first order
    BuyseTest.options(order.Hprojection = 1)
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

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                      data = d2.bis, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic")
    e2.BT <- BuyseTest(group ~ bin(toxicity) + strata,
                       data = d2.bis, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair,1:9) ## correct ordering of the pairs
    expect_equal(e.BT@covariance, e1.BT@covariance) ## assumes vs. does not assumes binary score when computing second order terms
    expect_equal(e.BT@covariance, e2.BT@covariance)
    
    expect_equal(as.double(e.BT@covariance), c(0,2/27,0,2/27,0)/3 )
})

## * Two endpoints
## ** no strata
set.seed(10)
d <- simBuyseTest(50)

BuyseTest.options(order.Hprojection = 1)
test_that("iid: two endpoints (no strata - first order)", {
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
    e.BT <- BuyseTest(Treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                                       data = d, scoring.rule = "Gehan",
                                       method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                                        data = d, keep.pairScore = TRUE, scoring.rule = "Gehan",
                                        method.inference = "u-statistic-bebu")
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["eventtime_1",]),
                 c(0.0001065024, 0.0001116864, -0.0000153344, 0.0002488576, 0.4454834743), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["eventtime_1e-12",]),
                 c(0.00194112, 0.00143968, -0.00028544, 0.00395168, 0.38625867), tol = 1e-6 )
})

BuyseTest.options(order.Hprojection = 2)
test_that("iid: two endpoints (no strata - second order)", {
    ## different endpoints
    e.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  bin(toxicity) + cont(score, threshold = 1),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e1.BT@tablePairScore[[1]]$index.pair, 1:2500)
    expect_equal(e1.BT@tablePairScore[[1]][e1.BT@tablePairScore[[2]]$index.pair,.(index.C,index.T)],
                 e1.BT@tablePairScore[[2]][,.(index.C,index.T)])
    expect_equal(e.BT@covariance, e1.BT@covariance)
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["toxicity_0.5",]),
                 c(0.002524914, 0.002524914, -0.002467086, 0.009984000, 0.160256410), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score_1",]),
                 c(0.003079997, 0.003232467, -0.002921262, 0.012154988, 0.078117626), tol = 1e-6 )

    ## same endpoint
    e.BT <- BuyseTest(Treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(Treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e.BT@covariance, e1.BT@covariance)
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair, 1:2500)
    expect_equal(e1.BT@tablePairScore[[1]][e1.BT@tablePairScore[[2]]$index.pair,.(index.C,index.T)],
                 e1.BT@tablePairScore[[2]][,.(index.C,index.T)])
    expect_equal(as.double(e.BT@covariance["score_2",]),
                 c(0.000374400, 0.0006309248, -0.000164736, 0.001334797, 0.13361289), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score_1",]),
                 c(0.001935052, 0.002390279, -0.001695749, 0.007716830, 0.089279985), tol = 1e-6 )

    ## same endpoint tte
    e.BT <- BuyseTest(Treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                      data = d, scoring.rule = "Gehan", keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(Treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                      data = d, scoring.rule = "Gehan", keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                       data = d, keep.pairScore = TRUE, scoring.rule = "Gehan",
                       method.inference = "u-statistic-bebu")

    expect_equal(e.BT@covariance, e1.BT@covariance)
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair, 1:2500)
    expect_equal(e1.BT@tablePairScore[[1]][e1.BT@tablePairScore[[2]]$index.pair,.(index.C,index.T)],
                 e1.BT@tablePairScore[[2]][,.(index.C,index.T)])
    expect_equal(as.double(e.BT@covariance["eventtime_1",]),
                 c(1.126726e-04, 1.183647e-04, -1.522106e-05,  2.614794e-04, 4.680545e-01), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["eventtime_1e-12",]),
                 c(0.0019582080,  0.0014531264, -0.0002877952,  0.0039869248,  0.3897334933), tol = 1e-6 )
})

## ** strata
d2 <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2))
d2$score1 <- d2$score
d2[strata == 1, score1 := 1]
d2$score2 <- d2$score
d2[strata == 2, score2 := 1]

test_that("iid: two endpoints (strata)", {

    ## first order
    BuyseTest.options(order.Hprojection = 1)
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
    expect_equal(as.double(e0.BT@covariance),
                 c(0.0009537952, 0.001180109, -0.0008541376, 0.0038421792, 0.0444467864), tol = 1e-6 )

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(Treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(Treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(Treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(e.BT@covariance, e1.BT@covariance)
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair, 1:5000)
    expect_equal(e1.BT@tablePairScore[[1]][e1.BT@tablePairScore[[2]]$index.pair,.(index.C,index.T)],
                 e1.BT@tablePairScore[[2]][,.(index.C,index.T)])
    expect_equal(as.double(e.BT@covariance["score2_1",]),
                 c(0.0009675260, 0.0011951397, -0.0008478746, 0.0038584150, 0.0446399926), tol = 1e-6 )

})



## * iid Peron
if(FALSE){
    library(data.table)
    library(prodlim)
    butils.base:::sourcePackage("BuyseTest", c.code = TRUE)
    dtAllC <- data.table(time = 1:10,
                         status = 1,
                         group = "C")
    dtAllT <- data.table(time = (1:10) + 0.5,
                         status = 1,
                         group = "T")
    dtAll <-  rbind(dtAllC, dtAllT)
    e.tte <- prodlim(Hist(time,status)~group, data = dtAll)

    BuyseTest.options(keep.survival = TRUE)
    e.BT <- BuyseTest(group ~ tte(time,status,threshold = 0.1),
                      data = data.frame(time = c(2.1,2.3),
                                        status = c(0,1),
                                        group = c("T","C")),
                      method.inference = "none",
                      model.tte = e.tte,
                      scoring.rule = "Peron",
                      keep.pairScore = TRUE)
    getPairScore(e.BT)
    scoreFavorable <- getSurvival(e.BT)$survTimeC[[1]][[1]][7]/getSurvival(e.BT)$survTimeT[[1]][[1]][5]
    scoreUnfavorable <- 1 - getSurvival(e.BT)$survTimeC[[1]][[1]][5]/getSurvival(e.BT)$survTimeT[[1]][[1]][5]
}

######################################################################
### test-BuyseTest-iid.R ends here
