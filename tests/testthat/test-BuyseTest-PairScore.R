### test-BuyseTest-tableComparison.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:33) 
## Version: 
## Last-Updated: okt 14 2021 (19:04) 
##           By: Brice Ozenne
##     Update #: 69
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
                       argsBin = list(p.T = list(c(0.5,0.5),c(0.25,0.75))),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(scale.T = 1:3, scale.Censoring.T = rep(1,3)))
dt.sim[eventtime1 >= 1, status1 := 0]
dt.sim[, time1 := eventtime1]
dt.sim[eventtime1 >= 1, time1 := 1]


## * test against tableComparison (no correction)
formula <- treatment ~ tte(time1, status1, threshold = 0.5) + cont(score1, 1) + bin(toxicity1) + tte(time1, status1, threshold = 0.25) + cont(score1, 0.5)
test_that("Full data - no correction", {

    BT.mixed <- BuyseTest(formula, data = dt.sim, scoring.rule = "Peron", correction.uninf = FALSE)

    expect_equal(as.double(BT.mixed@n.pairs),
                 prod(table(dt.sim$treatment)))

    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){ ## iEndpoint <- 1
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorable*weight),
                                                                  unfavorable = sum(unfavorable*weight),
                                                                  neutral = sum(neutral*weight),
                                                                  uninf = sum(uninf*weight))]
        manualScore <- rbind(manualScore,iScore)
    }

    ## check tablePairScore
    expect_equal(as.double(manualScore$favorable),
                 as.double(coef(BT.mixed, statistic = "count.favorable", cumulative = FALSE)))
    expect_equal(as.double(manualScore$unfavorable),
                 as.double(coef(BT.mixed, statistic = "count.unfavorable", cumulative = FALSE)))
    expect_equal(as.double(manualScore$neutral),
                 as.double(coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE)))
    expect_equal(as.double(manualScore$uninf),
                 as.double(coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE)))

    expect_equal(as.double(cumsum(BT.mixed@count.favorable-BT.mixed@count.unfavorable)/BT.mixed@n.pairs),
                 as.double(coef(BT.mixed, statistic = "netBenefit")))
    expect_equal(as.double(cumsum(BT.mixed@count.favorable)/cumsum(BT.mixed@count.unfavorable)),
                 as.double(coef(BT.mixed, statistic = "winRatio")))

    ## check number of pairs
    D <- length(BT.mixed@endpoint)
    vec.pair <- (coef(BT.mixed, statistic = "count.favorable", cumulative = FALSE) + coef(BT.mixed, statistic = "count.unfavorable", cumulative = FALSE) + coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE) + coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE))
    vec.RP <- (coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE) + coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE))
    expect_equal(as.double(vec.RP[-D]),as.double(vec.pair[-1]))
    
})

## * test against tableComparison (correction)
formula <- treatment ~ tte(time1, status1, threshold = 0.5) + cont(score1, 1) + bin(toxicity1) + tte(time1, status1, threshold = 0.25) + cont(score1, 0.5)

test_that("Full data", {
    BT.mixed <- suppressWarnings(BuyseTest(formula, data = dt.sim, scoring.rule = "Peron", correction.uninf = TRUE))

    expect_equal(as.double(BT.mixed@n.pairs),
                 prod(table(dt.sim$treatment)))

    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){ ## iEndpoint <- 1
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorableC),
                                                                  unfavorable = sum(unfavorableC),
                                                                  neutral = sum(neutralC))]
        manualScore <- rbind(manualScore,iScore)
    }

    ## check tablePairScore
    expect_equal(unname(coef(BT.mixed, statistic = "netBenefit")),manualScore[,cumsum(favorable-unfavorable)]/BT.mixed@n.pairs)
    expect_equal(unname(coef(BT.mixed, statistic = "winRatio")),manualScore[,cumsum(favorable)/cumsum(unfavorable)])

    ## check number of pairs
    D <- length(BT.mixed@endpoint)
    vec.pair <- (coef(BT.mixed, statistic = "count.favorable", cumulative = FALSE) + coef(BT.mixed, statistic = "count.unfavorable", cumulative = FALSE) + coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE) + coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE))
    vec.RP <- (coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE) + coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE))
    expect_equal(as.double(vec.RP[-D]),as.double(vec.pair[-1]))
})



##----------------------------------------------------------------------
### test-BuyseTest-tableComparison.R ends here
