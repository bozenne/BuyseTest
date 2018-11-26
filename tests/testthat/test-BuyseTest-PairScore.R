### test-BuyseTest-tableComparison.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:33) 
## Version: 
## Last-Updated: okt 30 2018 (16:21) 
##           By: Brice Ozenne
##     Update #: 46
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

context("Check tableComparison matches the summary of BuyseTest objects")

## * Settings
n.patients <- c(90,100)
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  keep.survival = TRUE,
                  method.inference = "none",
                  trace = 0)


## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))
dt.sim[eventtime1 >= 1, status1 := 0]
dt.sim[, time1 := eventtime1]
dt.sim[eventtime1 >= 1, time1 := 1]


## * test against tableComparison (no correction)
formula <- Treatment ~ tte(time1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(time1, 0.25, status1) + cont(score1, 0.5)
test_that("Full data - no correction", {

    BT.mixed <- BuyseTest(formula, data = dt.sim, method.tte = "Peron", correction.uninf = FALSE)

    expect_equal(as.double(BT.mixed@n.pairs),
                 prod(table(dt.sim$Treatment)))

    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorable*weight),
                                                                  unfavorable = sum(unfavorable*weight),
                                                                  neutral = sum(neutral*weight),
                                                                  uninf = sum(uninf*weight))]
        manualScore <- rbind(manualScore,iScore)
    }
    expect_equal(as.double(manualScore$favorable),
                 as.double(BT.mixed@count.favorable))
    expect_equal(as.double(manualScore$unfavorable),
                 as.double(BT.mixed@count.unfavorable))
    expect_equal(as.double(manualScore$neutral),
                 as.double(BT.mixed@count.neutral))
    expect_equal(as.double(manualScore$uninf),
                 as.double(BT.mixed@count.uninf))

    expect_equal(as.double(cumsum(BT.mixed@count.favorable-BT.mixed@count.unfavorable)/BT.mixed@n.pairs),
                 as.double(BT.mixed@Delta.netBenefit))
    expect_equal(as.double(cumsum(BT.mixed@count.favorable)/cumsum(BT.mixed@count.unfavorable)),
                 as.double(BT.mixed@Delta.winRatio))
})

## * test against tableComparison (correction)
formula <- Treatment ~ tte(time1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(time1, 0.25, status1) + cont(score1, 0.5)

test_that("Full data", {

    BT.mixed <- BuyseTest(formula,
                          data = dt.sim, method.tte = "Peron", correction.uninf = TRUE)

    expect_equal(as.double(BT.mixed@n.pairs),
                 prod(table(dt.sim$Treatment)))

    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorableC),
                                                                  unfavorable = sum(unfavorableC),
                                                                  neutral = sum(neutralC))]
        manualScore <- rbind(manualScore,iScore)
    }

    expect_equal(unname(BT.mixed@Delta.netBenefit),manualScore[,cumsum(favorable-unfavorable)]/BT.mixed@n.pairs)
    expect_equal(unname(BT.mixed@Delta.winRatio),manualScore[,cumsum(favorable)/cumsum(unfavorable)])
})



##----------------------------------------------------------------------
### test-BuyseTest-tableComparison.R ends here
