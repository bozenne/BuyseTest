### test-BuyseTest-correctionTTE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 30 2018 (23:45) 
## Version: 
## Last-Updated: okt 30 2018 (17:22) 
##           By: Brice Ozenne
##     Update #: 97
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

context("Check that method.tte = corrected  in BuyseTest is working correctly \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)
## * 1 endpoint
df <- data.frame("survie" = c(2.1, 4.1, 6.1, 8.1, 4, 6, 8, 10),
                 "event" = c(1, 1, 1, 0, 1, 0, 0, 1),
                 "group" = c(0, 0, 0, 0, 1, 1, 1, 1),
                 "score" = 1)

## ** Gehan
test_that("1 TTE endpoint - Gehan (no correction)", {
    Gehan <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                       data = df,
                       method.tte = "Gehan", correction.uninf = FALSE)

    expect_equal(as.double(Gehan@count.favorable), c(9,0))
    expect_equal(as.double(Gehan@count.unfavorable), c(2,0))
    expect_equal(as.double(Gehan@count.neutral), c(1,5))
    expect_equal(as.double(Gehan@count.uninf), c(4,0))

    iScore <- copy(getPairScore(Gehan))
    iScore[[1]][,c("endpoint","n") := list(1,NROW(iScore[[1]]))]
    iScore[[2]][,c("endpoint","n") := list(2,NROW(iScore[[1]]))]
    iiScore <- do.call(rbind, iScore)
    iiScoreS <- iiScore[,.(n = n[1], favorable = sum(favorable), unfavorable = sum(unfavorable)),by = "endpoint"]
        
    expect_equal(as.double(Gehan@Delta.netBenefit),iiScoreS[,cumsum(favorable-unfavorable)/n])
    expect_equal(as.double(Gehan@Delta.winRatio),iiScoreS[,cumsum(favorable)/cumsum(unfavorable)])
})

test_that("1 TTE endpoint - Gehan (correction at the pair level)", {
    ## survival first    
    GehanC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Gehan", correction.uninf = TRUE)
    ## getPairScore(GehanC)
    
    expect_equal(as.double(GehanC@count.favorable), c(12,0))
    expect_equal(as.double(GehanC@count.unfavorable), c(2+2/3,0))
    expect_equal(as.double(GehanC@count.neutral), c(1 + 1/3,1 + 1/3))
    expect_equal(as.double(GehanC@count.uninf), c(0,0))

    iScore <- copy(getPairScore(GehanC, endpoint = 1, strata = 1))
    expect_equal(iScore[, favorable + uninf * sum(favorable)/sum(favorable + unfavorable + neutral)], iScore[["favorableC"]])
    expect_equal(iScore[, unfavorable + uninf * sum(unfavorable)/sum(favorable + unfavorable + neutral)], iScore[["unfavorableC"]])
    expect_equal(iScore[, neutral + uninf * sum(neutral)/sum(favorable + unfavorable + neutral)], iScore[["neutralC"]])
    
    expect_equal(as.double(GehanC@Delta.netBenefit)[1],iScore[,sum(favorableC-unfavorableC)/.N])
    expect_equal(as.double(GehanC@Delta.winRatio)[1],iScore[,sum(favorableC)/sum(unfavorableC)])

})

test_that("1 TTE endpoint - Gehan (correction IPCW)", {
    ## survival first    
    GehanC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Gehan", correction.uninf = 2)

    factor <- 16/12 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(GehanC@count.favorable), c(9*factor,0))
    expect_equal(as.double(GehanC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(GehanC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(GehanC@count.uninf), c(0,0))

    iScore <- copy(getPairScore(GehanC, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable),
                         factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(GehanC@Delta.netBenefit[1]),iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(GehanC@Delta.winRatio[1]),iScoreS[,sum(favorable*factor)/sum(unfavorable*factor)])

    ## survival second
    GehanC2 <- BuyseTest(group ~  cont(score) + tte(survie, censoring = event, threshold = 1),
                         data = df, 
                         method.tte = "Gehan", correction.uninf = 2)
    
    expect_equal(GehanC@count.favorable[1], GehanC2@count.favorable[2])
    expect_equal(GehanC@count.unfavorable[1], GehanC2@count.unfavorable[2])
    expect_equal(GehanC@count.neutral[1], GehanC@count.neutral[2])
    expect_equal(GehanC@count.uninf[1], GehanC2@count.uninf[2])

    iScore <- copy(getPairScore(GehanC2, endpoint = 2, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable),
                           factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(GehanC2@Delta.netBenefit[2]),iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(GehanC2@Delta.winRatio[2]),iScoreS[,sum(favorable*factor)/cumsum(unfavorable*factor)])
})

## ** Peron
test_that("1 TTE endpoint - Peron (no correction)", {

    Peron <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                       data = df, 
                       method.tte = "Peron", correction.uninf = FALSE)

    ## summary(Peron, percentage = FALSE)
    expect_equal(as.double(Peron@count.favorable), c(10,0))
    expect_equal(as.double(Peron@count.unfavorable), c(2,0))
    expect_equal(as.double(Peron@count.neutral), c(1,4))
    expect_equal(as.double(Peron@count.uninf), c(3,0))

    iScore <- copy(getPairScore(Peron, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable))]

    expect_equal(as.double(Peron@Delta.netBenefit[1]),iScoreS[,sum(favorable-unfavorable)/n])
    expect_equal(as.double(Peron@Delta.winRatio[1]),iScoreS[,sum(favorable)/cumsum(unfavorable)])
})

    
test_that("1 TTE endpoint - Peron (IPCW)", {
    ## survival first
    PeronC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Peron", correction.uninf = 2)

    ## summary(PeronC, percentage = FALSE)
    factor <- PeronC@n.pairs/(PeronC@n.pairs-3) ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(PeronC@count.favorable), c(10*factor,0))
    expect_equal(as.double(PeronC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(PeronC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(PeronC@count.uninf), c(0,0))

    iScore <- copy(getPairScore(PeronC, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable), 
                         factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(PeronC@Delta.netBenefit[1]),iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(PeronC@Delta.winRatio[1]),iScoreS[,sum(favorable*factor)/cumsum(unfavorable*factor)])

    ## survival second
    PeronC2 <- BuyseTest(group ~  cont(score) + tte(survie, censoring = event, threshold = 1),
                         data = df, 
                         method.tte = "Peron", correction.uninf = 2)
    
    expect_equal(PeronC@count.favorable[1], PeronC2@count.favorable[2])
    expect_equal(PeronC@count.unfavorable[1], PeronC2@count.unfavorable[2])
    expect_equal(PeronC@count.neutral[1], PeronC@count.neutral[2])
    expect_equal(PeronC@count.uninf[1], PeronC2@count.uninf[2])

    iScore <- copy(getPairScore(PeronC2, endpoint = 1, strata = 1))
    iScoreS <- iScore[,.(n = .N, favorable = sum(favorable), unfavorable = sum(unfavorable), 
                         factor = .N/sum(favorable+unfavorable+neutral))]

    expect_equal(as.double(PeronC2@Delta.netBenefit[1]),iScoreS[,sum(favorable*factor-unfavorable*factor)/n])
    expect_equal(as.double(PeronC2@Delta.winRatio[1]),iScoreS[,sum(favorable*factor)/cumsum(unfavorable*factor)])
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

test_that("2 TTE endpoint - IPW induces bias when censoring is correlated with 2nd endpoint", {
    BT.all <- BuyseTest(trt ~ cont(Y1, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 0, method.inference = "none")
    expect_equal(as.double(BT.all@count.favorable), c(0,4))
    expect_equal(as.double(BT.all@count.unfavorable), c(0,4))
    expect_equal(as.double(BT.all@count.neutral), c(25,17))
    expect_equal(as.double(BT.all@count.uninf), c(0,0))

    BT.uninf <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                          correction.uninf = 0, method.inference = "none")
    expect_equal(as.double(BT.uninf@count.favorable), c(0,4))
    expect_equal(as.double(BT.uninf@count.unfavorable), c(0,4))
    expect_equal(as.double(BT.uninf@count.neutral), c(4,17))
    expect_equal(as.double(BT.uninf@count.uninf), c(21,0))

    BT.ipw <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 2, method.inference = "none")
    expect_equal(as.double(BT.ipw@count.favorable), c(0,25))
    expect_equal(as.double(BT.ipw@count.unfavorable), c(0,0))
    expect_equal(as.double(BT.ipw@count.neutral), c(25,0))
    expect_equal(as.double(BT.ipw@count.uninf), c(0,0))

    BT.esp <- BuyseTest(trt ~ cont(Y1c, threshold = 1) + bin(Y2), data = dt,
                        correction.uninf = 1, method.inference = "none", keep.pairScore = TRUE)
    expect_equal(as.double(BT.esp@count.favorable), as.double(BT.all@count.favorable))
    expect_equal(as.double(BT.esp@count.unfavorable), as.double(BT.all@count.unfavorable))
    expect_equal(as.double(BT.esp@count.neutral), as.double(BT.all@count.neutral))
    expect_equal(as.double(BT.esp@count.uninf), as.double(BT.all@count.uninf))

})

##----------------------------------------------------------------------
### test-BuyseTest-correctionTTE.R ends here
