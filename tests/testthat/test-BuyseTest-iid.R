### test-BuyseTest-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  8 2019 (11:54) 
## Version: 
## Last-Updated: nov 12 2019 (18:25) 
##           By: Brice Ozenne
##     Update #: 76
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
var2 <- function(x){var(x)*(length(x)-1)/length(x)}
cov2 <- function(x,y){cov(x,y)*(length(x)-1)/length(x)}

## * Settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = FALSE,
                  trace = 0)

## * iid average
## ** 1 binary variable
## *** no strata
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


## *** strata
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

## ** 1 continuous variable

test_that("Manual calculation of second order H projection (no strata)",{
    n <- 5
    set.seed(10)
    dt <- simBuyseTest(n)

    BuyseTest.options(order.Hprojection = 1)
    e.BT_c1 <- BuyseTest(treatment ~ cont(score),
                         data = dt, trace = 0, 
                         method.inference = "u-statistic")
    BuyseTest.options(order.Hprojection = 2)
    e.BT_c2 <- BuyseTest(treatment ~ cont(score),
                         data = dt, trace = 0, 
                         method.inference = "u-statistic")
    e.BT_c3 <- BuyseTest(treatment ~ cont(score),
                         data = dt, trace = 0, keep.pairScore = TRUE, 
                         method.inference = "u-statistic")

    ## manual calculation
    dt.pair <- getPairScore(e.BT_c3)[,.(index.C,index.T,favorable,unfavorable)]
    dt.pair[, H1C.favorable := mean(favorable), by = index.C]
    dt.pair[, H1T.favorable := mean(favorable), by = index.T]
    dt.pair[, H1C.favorable := H1C.favorable - mean(favorable)]
    dt.pair[, H1T.favorable := H1T.favorable - mean(favorable)]

    dt.pair[, H1C.unfavorable := mean(unfavorable), by = index.C]
    dt.pair[, H1T.unfavorable := mean(unfavorable), by = index.T]
    dt.pair[, H1C.unfavorable := H1C.unfavorable - mean(unfavorable)]
    dt.pair[, H1T.unfavorable := H1T.unfavorable - mean(unfavorable)]

    dt.pair[, H2.favorable := favorable - H1C.favorable - H1T.favorable]
    dt.pair[, H2.unfavorable := unfavorable - H1C.unfavorable - H1T.unfavorable]

    ## check H1
    expect_true(all(abs(dt.pair[!duplicated(index.C),.(H1C.favorable/.N,H1C.unfavorable/.N)]-iid(e.BT_c3)[1:n,])<1e-6))
    expect_true(all(abs(dt.pair[!duplicated(index.T),.(H1T.favorable/.N,H1T.unfavorable/.N)]-iid(e.BT_c3)[(n+1):(2*n),])<1e-6))
        
    ## check H2
    manual <- dt.pair[,.(favorable = var2(H2.favorable)/.N,
                         unfavorable = var2(H2.unfavorable)/.N,
                         covariance = cov2(H2.favorable,H2.unfavorable)/.N)]
    expect_equal(as.double(manual), as.double((e.BT_c2@covariance - e.BT_c1@covariance)[1,c("favorable","unfavorable","covariance")]))
    expect_equal(as.double(manual), c(0.002304, 0.002304, -0.002304), tol = 1e-5)   
})

test_that("Manual calculation of second order H projection (strata)",{
    n <- 5
    set.seed(10)
    dt <- simBuyseTest(n)
    dtS <- rbind(cbind(S = 1, dt), cbind(S = 2, dt), cbind(S = 3, dt))
    
    BuyseTest.options(order.Hprojection = 1)
    e.BT_c1 <- BuyseTest(treatment ~ cont(score) + S,
                         data = dtS, trace = 0, 
                         method.inference = "u-statistic")
    BuyseTest.options(order.Hprojection = 2)
    e.BT_c2 <- BuyseTest(treatment ~ cont(score) + S,
                         data = dtS, trace = 0, 
                         method.inference = "u-statistic")
    e.BT_c3 <- BuyseTest(treatment ~ cont(score) + S,
                         data = dtS, trace = 0, keep.pairScore = TRUE, 
                         method.inference = "u-statistic")

    ## manual calculation
    dt.pair <- getPairScore(e.BT_c3)[,.(strata,index.C,index.T,favorable,unfavorable)]
    dt.pair[, H1C.favorable := mean(favorable), by = index.C]
    dt.pair[, H1T.favorable := mean(favorable), by = index.T]
    dt.pair[, H1C.favorable := H1C.favorable - mean(favorable), by = "strata"]
    dt.pair[, H1T.favorable := H1T.favorable - mean(favorable), by = "strata"]

    dt.pair[, H1C.unfavorable := mean(unfavorable), by = index.C]
    dt.pair[, H1T.unfavorable := mean(unfavorable), by = index.T]
    dt.pair[, H1C.unfavorable := H1C.unfavorable - mean(unfavorable), by = "strata"]
    dt.pair[, H1T.unfavorable := H1T.unfavorable - mean(unfavorable), by = "strata"]

    dt.pair[, H2.favorable := favorable - H1C.favorable - H1T.favorable]
    dt.pair[, H2.unfavorable := unfavorable - H1C.unfavorable - H1T.unfavorable]

    ## check H1
    expect_true(all(abs(dt.pair[!duplicated(index.C),.(H1C.favorable/.N,H1C.unfavorable/.N)]-iid(e.BT_c3)[which(dtS$treatment=="C"),])<1e-6))
    expect_true(all(abs(dt.pair[!duplicated(index.T),.(H1T.favorable/.N,H1T.unfavorable/.N)]-iid(e.BT_c3)[which(dtS$treatment=="T"),])<1e-6))
        
    ## check H2
    manual <- dt.pair[,.(favorable = var2(H2.favorable)/.N,
                         unfavorable = var2(H2.unfavorable)/.N,
                         covariance = cov2(H2.favorable,H2.unfavorable)/.N)]
    expect_equal(as.double(manual), as.double((e.BT_c2@covariance - e.BT_c1@covariance)[1,c("favorable","unfavorable","covariance")]))
    expect_equal(as.double(manual), c(0.000768,  0.000768, -0.000768), tol = 1e-5)
})

## ** 1 TTE variable
test_that("iid: TTE and no strata",{
    BuyseTest.options(order.Hprojection = 1)

    n <- 5
    set.seed(10)
    dt <- simBuyseTest(n)
    ## dt

    e.BT_tte1 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           scoring.rule = "Gehan",
                           method.inference = "u-statistic")
    e.BT_tte2 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1e5) + tte(eventtime, status, threshold = 1-1e-5),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           scoring.rule = "Gehan",
                           method.inference = "u-statistic")
    e.BT_tte3 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1) + tte(eventtime, status, threshold = 1-1e-5),
                           data = dt, trace = 0,
                           keep.pairScore = TRUE,
                           scoring.rule = "Gehan",
                           method.inference = "u-statistic")

    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte2)[2,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte3)[1,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte3)[1,],
                 confint(e.BT_tte3)[2,],
                 tol = 1e-6)
})

## ** Two endpoints
## *** no strata
set.seed(10)
d <- simBuyseTest(50)

BuyseTest.options(order.Hprojection = 1)
test_that("iid: two endpoints (no strata - first order)", {
    ## different endpoints
    e.BT <- BuyseTest(treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  bin(toxicity) + cont(score, threshold = 1),
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
    e.BT <- BuyseTest(treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                       data = d, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e.BT@covariance["score_2",]),
                 c(0.00036192, 0.000613760, -0.000166400, 0.00130848, 0.13086775), tol = 1e-6 )
    expect_equal(as.double(e.BT@covariance["score_1",]),
                 c(0.00190759, 0.002360218, -0.001708275, 0.007684358, 0.088893573), tol = 1e-6 )

    ## same endpoint tte
    e.BT <- BuyseTest(treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                                       data = d, scoring.rule = "Gehan",
                                       method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
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
    e.BT <- BuyseTest(treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(treatment ~  bin(toxicity) + cont(score, threshold = 1),
                      data = d, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  bin(toxicity) + cont(score, threshold = 1),
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
    e.BT <- BuyseTest(treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d, keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
                      data = d, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  cont(score, threshold = 2) + cont(score, threshold = 1),
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
    e.BT <- BuyseTest(treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                      data = d, scoring.rule = "Gehan", keep.pairScore = FALSE,
                      method.inference = "u-statistic")
    e1.BT <- BuyseTest(treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
                      data = d, scoring.rule = "Gehan", keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~  tte(eventtime, threshold = 1, censoring = status) + tte(eventtime, threshold = 0, censoring = status),
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

## *** strata
d2 <- rbind(cbind(d, strata = 1),
            cbind(d, strata = 2))
d2$score1 <- d2$score
d2[strata == 1, score1 := 1]
d2$score2 <- d2$score
d2[strata == 2, score2 := 1]

test_that("iid: two endpoints (strata)", {

    ## first order
    BuyseTest.options(order.Hprojection = 1)
    e0.BT <- BuyseTest(treatment ~ cont(score, threshold = 1) + strata,
                       data = d2,
                       method.inference = "u-statistic")
    e.BT <- BuyseTest(treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2,
                      method.inference = "u-statistic")

    e2.BT <- BuyseTest(treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu")

    expect_equal(as.double(e0.BT@covariance), as.double(e.BT@covariance[2,]))
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(as.double(e0.BT@covariance),
                 c(0.0009537952, 0.001180109, -0.0008541376, 0.0038421792, 0.0444467864), tol = 1e-6 )

    ## second order
    BuyseTest.options(order.Hprojection = 2)
    e.BT <- BuyseTest(treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2,
                      method.inference = "u-statistic") ## neglect some terms
    e1.BT <- BuyseTest(treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                      data = d2, keep.pairScore = TRUE,
                      method.inference = "u-statistic")
    e2.BT <- BuyseTest(treatment ~ cont(score1, threshold = 1) + cont(score2, threshold = 1) + strata,
                       data = d2, keep.pairScore = TRUE,
                       method.inference = "u-statistic-bebu") ## neglect some terms

    ## expect_equal(e.BT@covariance, e1.BT@covariance, tol = 1e-4) ## imperfect match
    expect_equal(e.BT@covariance, e2.BT@covariance)
    expect_equal(e1.BT@tablePairScore[[1]]$index.pair, 1:5000)
    expect_equal(e1.BT@tablePairScore[[1]][e1.BT@tablePairScore[[2]]$index.pair,.(index.C,index.T)],
                 e1.BT@tablePairScore[[2]][,.(index.C,index.T)])
    expect_equal(as.double(e1.BT@covariance["score2_1",]),
                 c(0.0009675260, 0.0011951397, -0.0008478746, 0.0038584150, 0.0446399926), tol = 1e-6 )

})





## * iid Peron

## ** 1 TTE variable
test_that("iid with nuisance parameters: 1 TTE",{
    BuyseTest.options(order.Hprojection = 1)

    n <- 5
    set.seed(10)
    dt <- simBuyseTest(n)
    dt$X0 <- 0
    dt$treatment2 <- 1-as.numeric(dt$treatment)
    ## dt <- data.table("treatment" = c("C", "C", "C", "C", "C", "T", "T", "T", "T", "T"), 
    ##                  "toxicity" = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 0), 
    ##                  "score" = c( 0.54361539, -0.70762484, -0.36944577, -1.32197565,  1.28059746,  0.01874617, -0.18425254, -1.37133055, -0.59916772,  0.29454513), 
    ##                  "eventtime" = c(1.8252132, 2.9489056, 0.7213402, 0.6322603, 0.2212117, 0.1453481, 0.4855601, 0.2547505, 1.0340368, 0.3579324), 
    ##                  "status" = c(0, 1, 0, 1, 0, 0, 0, 0, 0, 1))

    e.BT_tte1 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           method.inference = "u-statistic")
    e.BT_tte1.bis <- BuyseTest(treatment2 ~ tte(eventtime, status, threshold = 1),
                               data = dt, trace = 0, 
                               keep.pairScore = TRUE,
                               method.inference = "u-statistic")
        
    e.BT_tte2 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1e5) + tte(eventtime, status, threshold = 1-1e-5),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           method.inference = "u-statistic")

    e.BT_tte3 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1) + tte(eventtime, status, threshold = 1-1e-5),
                           data = dt, trace = 0,
                           keep.pairScore = TRUE,
                           method.inference = "u-statistic")

    e.BT_tte4 <- BuyseTest(treatment ~ bin(X0) + tte(eventtime, status, threshold = 1),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           method.inference = "u-statistic")

    e.BT_tte5 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 2) + bin(X0) + tte(eventtime, status, threshold = 1),
                           data = dt, trace = 0, 
                           keep.pairScore = TRUE,
                           method.inference = "u-statistic")
    e.BT_tte5.bis <- BuyseTest(treatment2 ~ tte(eventtime, status, threshold = 2) + bin(X0) + tte(eventtime, status, threshold = 1),
                               data = dt, trace = 0, 
                               keep.pairScore = TRUE,
                               method.inference = "u-statistic")

    ## unchanged when switching around the group
    expect_equal(e.BT_tte1@count.unfavorable,e.BT_tte1.bis@count.favorable)
    expect_equal(e.BT_tte1@covariance["unfavorable"],e.BT_tte1.bis@covariance["favorable"])

    expect_equal(e.BT_tte5@count.unfavorable,e.BT_tte5.bis@count.favorable)
    expect_equal(e.BT_tte5@covariance["unfavorable"],e.BT_tte5.bis@covariance["favorable"])

    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte2)[2,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte3)[1,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte2)[2,],
                 confint(e.BT_tte3)[2,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte3)[1,],
                 confint(e.BT_tte3)[2,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte4)[2,],
                 tol = 1e-6)
    expect_equal(confint(e.BT_tte1)[1,],
                 confint(e.BT_tte5)[3,],
                 tol = 1e-6)
    
})


## ** 1 TTE variable and 1 binary
test_that("iid with nuisance parameters: 1 TTE + 1 binary",{
    BuyseTest.options(order.Hprojection = 1)
    
    n <- 5
    set.seed(10)
    dt <- simBuyseTest(n)
    dt$bin0 <- 0
    dt$binF <- 1-as.numeric(dt$treatment)
    
    ## dt
    e.BT_ttebin <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 1) + bin(bin),
                             data = dt, trace = 0, 
                             keep.pairScore = TRUE,
                             method.inference = "u-statistic")

    e.BT_tte0 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0),
                             data = dt, trace = 0, 
                             keep.pairScore = TRUE,
                             method.inference = "u-statistic")
    e.BT_tte0@covariance
    e.BT_ttebin@covariance
    e.BT_ttebin@iidNuisance$unfavorable
}
######################################################################
### test-BuyseTest-iid.R ends here
