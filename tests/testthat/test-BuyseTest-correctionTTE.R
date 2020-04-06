### test-BuyseTest-correctionTTE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 30 2018 (23:45) 
## Version: 
## Last-Updated: apr  6 2020 (11:45) 
##           By: Brice Ozenne
##     Update #: 129
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

context("Check that scoring.rule = corrected  in BuyseTest is working correctly \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)
## * 1 endpoint
df <- data.frame("survie" = c(2.1, 4.1, 6.1, 8.1, 4, 6, 8, 10),
                 "event" = c(1, 1, 1, 0, 1, 0, 0, 1),
                 "group" = c(0, 0, 0, 0, 1, 1, 1, 1),
                 "score" = 1,
                 stringsAsFactors = FALSE)

## ** Gehan
test_that("1 TTE endpoint - Gehan (no correction)", {
    Gehan <- BuyseTest(group ~ tte(survie, status = event, threshold = 1) + cont(score),
                       data = df,
                       scoring.rule = "Gehan", correction.uninf = FALSE)

    expect_equal(as.double(coef(Gehan, statistic = "count.favorable", cumulative = FALSE)), c(9,0))
    expect_equal(as.double(coef(Gehan, statistic = "count.unfavorable", cumulative = FALSE)), c(2,0))
    expect_equal(as.double(coef(Gehan, statistic = "count.neutral", cumulative = FALSE)), c(1,5))
    expect_equal(as.double(coef(Gehan, statistic = "count.uninf", cumulative = FALSE)), c(4,0))

    iScore <- copy(getPairScore(Gehan))
    iScore[[1]][,c("endpoint","n") := list(1,NROW(iScore[[1]]))]
    iScore[[2]][,c("endpoint","n") := list(2,NROW(iScore[[1]]))]
    iiScore <- do.call(rbind, iScore)
    iiScoreS <- iiScore[,.(n = n[1], favorable = sum(favorable), unfavorable = sum(unfavorable)),by = "endpoint"]
        
    expect_equal(as.double(coef(Gehan, statistic = "netBenefit"))
,iiScoreS[,cumsum(favorable-unfavorable)/n])
    expect_equal(as.double(coef(Gehan, statistic = "winRatio")),iiScoreS[,cumsum(favorable)/cumsum(unfavorable)])
})

test_that("1 TTE endpoint - Gehan (correction at the pair level)", {
    ## survival first    
    GehanC <- BuyseTest(group ~ tte(survie, status = event, threshold = 1) + cont(score),
                        data = df, 
                        scoring.rule = "Gehan", correction.uninf = TRUE)
    ## getPairScore(GehanC)
    ## summary(GehanC)
    
    expect_equal(as.double(coef(GehanC, statistic = "count.favorable", cumulative = FALSE)), c(12,0))
    expect_equal(as.double(coef(GehanC, statistic = "count.unfavorable", cumulative = FALSE)), c(2+2/3,0))
    expect_equal(as.double(coef(GehanC, statistic = "count.neutral", cumulative = FALSE)), c(1 + 1/3,1 + 1/3))
    expect_equal(as.double(coef(GehanC, statistic = "count.uninf", cumulative = FALSE)), c(0,0))

    iScore <- copy(getPairScore(GehanC, endpoint = 1, strata = 1))
    expect_equal(iScore[, favorable + uninf * sum(favorable)/sum(favorable + unfavorable + neutral)], iScore[["favorableC"]])
    expect_equal(iScore[, unfavorable + uninf * sum(unfavorable)/sum(favorable + unfavorable + neutral)], iScore[["unfavorableC"]])
    expect_equal(iScore[, neutral + uninf * sum(neutral)/sum(favorable + unfavorable + neutral)], iScore[["neutralC"]])
    
    expect_equal(as.double(coef(GehanC, statistic = "netBenefit")[1]),iScore[,sum(favorableC-unfavorableC)/.N])
    expect_equal(as.double(coef(GehanC, statistic = "winRatio")[1]),iScore[,sum(favorableC)/sum(unfavorableC)])

})

test_that("1 TTE endpoint - Gehan (correction IPCW)", {
    ## survival first    
    GehanC <- BuyseTest(group ~ tte(survie, status = event, threshold = 1) + cont(score),
                        data = df, 
                        scoring.rule = "Gehan", correction.uninf = 2)

    factor <- 16/12 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(coef(GehanC, statistic = "count.favorable", cumulative = FALSE)), c(9*factor,0))
    expect_equal(as.double(coef(GehanC, statistic = "count.unfavorable", cumulative = FALSE)), c(2*factor,0))
    expect_equal(as.double(coef(GehanC, statistic = "count.neutral", cumulative = FALSE)), c(1*factor,1*factor))
    expect_equal(as.double(coef(GehanC, statistic = "count.uninf", cumulative = FALSE)), c(0,0))

    iScore <- copy(getPairScore(GehanC, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable),
                         factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(coef(GehanC, statistic = "netBenefit"))[1], iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(coef(GehanC, statistic = "winRatio"))[1], iScoreS[,sum(favorable*factor)/sum(unfavorable*factor)])
})

## ** Peron
test_that("1 TTE endpoint - Peron (no correction)", {
    Peron <- BuyseTest(group ~ tte(survie, status = event, threshold = 1) + cont(score),
                       data = df, 
                       scoring.rule = "Peron", correction.uninf = FALSE)

    ## summary(Peron, percentage = FALSE)
    expect_equal(as.double(coef(Peron, statistic = "count.favorable", cumulative = FALSE)), c(10,0))
    expect_equal(as.double(coef(Peron, statistic = "count.unfavorable", cumulative = FALSE)), c(2,0))
    expect_equal(as.double(coef(Peron, statistic = "count.neutral", cumulative = FALSE)), c(1,4))
    expect_equal(as.double(coef(Peron, statistic = "count.uninf", cumulative = FALSE)), c(3,0))

    iScore <- copy(getPairScore(Peron, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable))]

    expect_equal(as.double(coef(Peron, statistic = "netBenefit")[1]),iScoreS[,sum(favorable-unfavorable)/n])
    expect_equal(as.double(coef(Peron, statistic = "winRatio")[1]),iScoreS[,sum(favorable)/cumsum(unfavorable)])
})

    
test_that("1 TTE endpoint - Peron (IPCW)", {
    ## survival first
    PeronC <- BuyseTest(group ~ tte(survie, status = event, threshold = 1) + cont(score),
                        data = df, 
                        scoring.rule = "Peron", correction.uninf = 2)

    ## summary(PeronC, percentage = FALSE)
    factor <- as.double(PeronC@n.pairs/(PeronC@n.pairs-3)) ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(coef(PeronC, statistic = "count.favorable", cumulative = FALSE)), c(10*factor,0))
    expect_equal(as.double(coef(PeronC, statistic = "count.unfavorable", cumulative = FALSE)), c(2*factor,0))
    expect_equal(as.double(coef(PeronC, statistic = "count.neutral", cumulative = FALSE)), c(1*factor,1*factor))
    expect_equal(as.double(coef(PeronC, statistic = "count.uninf", cumulative = FALSE)), c(0,0))

    iScore <- copy(getPairScore(PeronC, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable), 
                         factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(coef(PeronC, statistic = "netBenefit")[1]),iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(coef(PeronC, statistic = "winRatio")[1]),iScoreS[,sum(favorable*factor)/cumsum(unfavorable*factor)])
})

## * 2 endpoints
## ** categorical variables

## simulate data
n <- 10
set.seed(10)
dt <- data.table(trt = c(rep(0,n/2),rep(1,n/2)),
                 Y1 = 1,
                 C = (1:n) %in% sample.int(n, n/2, replace = FALSE))
dt[C==FALSE, Y1c := Y1]
dt[,Y2 := as.numeric(NA)]
dt[trt==0, Y2 := as.numeric(C)]
dt[trt==1, Y2 := as.numeric(1-C)]

## table(dt$trt,dt$Y1,dt$Y2)

test_that("2 endpoints - IPW induces bias when censoring is correlated with 2nd endpoint", {
    BT.all <- BuyseTest(trt ~ cont(Y1, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 0, method.inference = "none")
    expect_equal(as.double(coef(BT.all, statistic = "count.favorable", cumulative = FALSE)), c(0,4))
    expect_equal(as.double(coef(BT.all, statistic = "count.unfavorable", cumulative = FALSE)), c(0,4))
    expect_equal(as.double(coef(BT.all, statistic = "count.neutral", cumulative = FALSE)), c(25,17))
    expect_equal(as.double(coef(BT.all, statistic = "count.uninf", cumulative = FALSE)), c(0,0))

    BT.uninf <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                          correction.uninf = 0, method.inference = "none")
    expect_equal(as.double(coef(BT.uninf, statistic = "count.favorable", cumulative = FALSE)), c(0,4))
    expect_equal(as.double(coef(BT.uninf, statistic = "count.unfavorable", cumulative = FALSE)), c(0,4))
    expect_equal(as.double(coef(BT.uninf, statistic = "count.neutral", cumulative = FALSE)), c(4,17))
    expect_equal(as.double(coef(BT.uninf, statistic = "count.uninf", cumulative = FALSE)), c(21,0))

    BT.ipw <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 2, method.inference = "none")
    expect_equal(as.double(coef(BT.ipw, statistic = "count.favorable", cumulative = FALSE)), c(0,25))
    expect_equal(as.double(coef(BT.ipw, statistic = "count.unfavorable", cumulative = FALSE)), c(0,0))
    expect_equal(as.double(coef(BT.ipw, statistic = "count.neutral", cumulative = FALSE)), c(25,0))
    expect_equal(as.double(coef(BT.ipw, statistic = "count.uninf", cumulative = FALSE)), c(0,0))

    BT.esp <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 1, method.inference = "none", keep.pairScore = TRUE)
    expect_equal(as.double(coef(BT.esp, statistic = "count.favorable", cumulative = FALSE)), as.double(coef(BT.all, statistic = "count.favorable", cumulative = FALSE)))
    expect_equal(as.double(coef(BT.esp, statistic = "count.unfavorable", cumulative = FALSE)), as.double(coef(BT.all, statistic = "count.unfavorable", cumulative = FALSE)))
    expect_equal(as.double(coef(BT.esp, statistic = "count.neutral", cumulative = FALSE)), as.double(coef(BT.all, statistic = "count.neutral", cumulative = FALSE)))
    expect_equal(as.double(coef(BT.esp, statistic = "count.uninf", cumulative = FALSE)), as.double(coef(BT.all, statistic = "count.uninf", cumulative = FALSE)))
})

## ** time to event variables
set.seed(10)
dt.sim <- simBuyseTest(n.T = 20,
                       n.C = 20)
test_that("2 TTE endpoints - check consistency across threshold when using corrections", {
    e.BT2 <- BuyseTest(treatment ~ tte(eventtime,status,1)  + tte(eventtime,status,0.75) + tte(eventtime,status,0.5),
                       correction.uninf = 1, data = dt.sim)
    e.BT1 <- BuyseTest(treatment ~ tte(eventtime,status,0.5),
                       correction.uninf = 1, data = dt.sim)

    expect_equal(sum(coef(e.BT2, statistic = "count.favorable", cumulative = FALSE)), as.double(coef(e.BT1, statistic = "count.favorable", cumulative = FALSE)[1]))
    expect_equal(sum(coef(e.BT2, statistic = "count.unfavorable", cumulative = FALSE)), as.double(coef(e.BT1, statistic = "count.unfavorable", cumulative = FALSE)[1]))
    expect_equal(as.double(coef(e.BT2, statistic = "count.neutral", cumulative = FALSE)[3]), as.double(coef(e.BT1, statistic = "count.neutral", cumulative = FALSE)[1]))
    expect_equal(as.double(coef(e.BT2, statistic = "count.uninf", cumulative = FALSE)[3]), as.double(coef(e.BT1, statistic = "count.uninf", cumulative = FALSE)[1]))

    e.BT2 <- BuyseTest(treatment ~ tte(eventtime,status,1)  + tte(eventtime,status,0.75) + tte(eventtime,status,0.5),
                       correction.uninf = 2, data = dt.sim)
    e.BT1 <- BuyseTest(treatment ~ tte(eventtime,status,0.5),
                       correction.uninf = 2, data = dt.sim)

    expect_equal(sum(coef(e.BT2, statistic = "count.favorable", cumulative = FALSE)), as.double(coef(e.BT1, statistic = "count.favorable", cumulative = FALSE)[1]))
    expect_equal(sum(coef(e.BT2, statistic = "count.unfavorable", cumulative = FALSE)), as.double(coef(e.BT1, statistic = "count.unfavorable", cumulative = FALSE)[1]))
    expect_equal(as.double(coef(e.BT2, statistic = "count.neutral", cumulative = FALSE)[3]), as.double(coef(e.BT1, statistic = "count.neutral", cumulative = FALSE)[1]))
    expect_equal(as.double(coef(e.BT2, statistic = "count.uninf", cumulative = FALSE)[3]), as.double(coef(e.BT1, statistic = "count.uninf", cumulative = FALSE)[1]))
})

test_that("TTE,cont,TTE endpoints - check consistency across threshold when using corrections", {
    e.BT2 <- BuyseTest(treatment ~ tte(eventtime,status,1)  + cont(score,1) + tte(eventtime,status,0.5),
                       correction.uninf = 1, data = dt.sim)
    expect_equal(as.double(coef(e.BT2, statistic = "count.neutral", cumulative = FALSE)[2]),
                 as.double(coef(e.BT2, statistic = "count.favorable", cumulative = FALSE)[3] + coef(e.BT2, statistic = "count.unfavorable", cumulative = FALSE)[3] + coef(e.BT2, statistic = "count.neutral", cumulative = FALSE)[3]))
})
##----------------------------------------------------------------------
### test-BuyseTest-correctionTTE.R ends here
