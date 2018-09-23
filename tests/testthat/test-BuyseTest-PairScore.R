### test-BuyseTest-tableComparison.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:33) 
## Version: 
## Last-Updated: sep 23 2018 (11:41) 
##           By: Brice Ozenne
##     Update #: 31
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
BuyseTest.options(check = FALSE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)


## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1:3, rates.Censor = rep(1,3)))

dtRed.sim <- dt.sim[, .SD[1:50], by = "Treatment"]

## * test against tableComparison (no correction)
formula <- Treatment ~ tte(eventtime1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(eventtime1, 0.25, status1) + cont(score1, 0.5)
test_that("Full data - no correction", {

    BT.mixed <- BuyseTest(formula, data = dt.sim, method.tte = "Peron", correction.uninf.tte = FALSE)

    summary(BT.mixed, percentage = FALSE)
    getPairScore(BT.mixed, endpoint = 1)
    
    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorable*weight),
                                                                  unfavorable = sum(unfavorable*weight),
                                                                  neutral = sum(neutral*weight),
                                                                  uninformative = sum(uninformative*weight))]
        manualScore <- rbind(manualScore,iScore)
    }
    expect_equal(as.double(manualScore$favorable),
                 as.double(BT.mixed@count.favorable))
    expect_equal(as.double(manualScore$unfavorable),
                 as.double(BT.mixed@count.unfavorable))
    expect_equal(as.double(manualScore$neutral),
                 as.double(BT.mixed@count.neutral))
    expect_equal(as.double(manualScore$uninformative),
                 as.double(BT.mixed@count.uninf))

    summary(BT.mixed, percentage = FALSE)
    3133.39     +2361.10+ 3918.03+ 90.95
    
    expect_equal(as.double(cumsum(BT.mixed@count.favorable-BT.mixed@count.unfavorable)/BT.mixed@n.pairs),
                 as.double(BT.mixed@Delta.netChance))
    expect_equal(as.double(cumsum(BT.mixed@count.favorable)/cumsum(BT.mixed@count.unfavorable)),
                 as.double(BT.mixed@Delta.winRatio))
})

## * test against tableComparison (correction)
if(FALSE){
formula <- Treatment ~ tte(eventtime1, 0.5, status1) + cont(score1, 1) + bin(toxicity1) + tte(eventtime1, 0.25, status1) + cont(score1, 0.5)
BT.mixed <- BuyseTest(formula,
                      data = dt.sim, method.tte = "Peron corrected")

test_that("Full data", {

    manualScore <- NULL
    for(iEndpoint in 1:length(BT.mixed@endpoint)){
        iScore <- getPairScore(BT.mixed, endpoint = iEndpoint)[,.(favorable = sum(favorable*weight),
                                                                  unfavorable = sum(unfavorable*weight),
                                                                  neutral = sum(neutral*weight),
                                                                  uninformative = sum(uninformative*weight))]
        manualScore <- rbind(manualScore,iScore)
    }

    expect_equal(unname(tail(BT.mixed@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(BT.mixed@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])
})
}


### does not work because of the estimation of the survival
## test_that("First 50 patients in each arm", {
    ## BT.mixedRed <- BuyseTest(formula,
                             ## data = dtRed.sim, method.tte = "Peron")
    ## summary(BT.mixedRed, percentage = FALSE)

    ## BT.mixedRed@tableComparison[[1]][1:5]
    
    ## test <- tableComparison2Delta(BT.mixed@tableComparison,
                                  ## correct.tte = BT.mixed@method.tte$correction,
                                  ## maxData.T = which(dt.sim$Treatment==1)[50],
                                  ## maxData.C = which(dt.sim$Treatment==0)[50])
    ## expect_equal(BT.mixedRed@Delta.netChance,test[["Delta.netChance"]])
    ## expect_equal(BT.mixedRed@Delta.winRatio,test[["Delta.winRatio"]])
## })

##----------------------------------------------------------------------
### test-BuyseTest-tableComparison.R ends here