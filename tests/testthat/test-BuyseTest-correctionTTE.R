### test-BuyseTest-correctionTTE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 30 2018 (23:45) 
## Version: 
## Last-Updated: sep  7 2018 (13:39) 
##           By: Brice Ozenne
##     Update #: 53
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
BuyseTest.options(check = FALSE,
                  keep.individualScore = TRUE,
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
                       method.tte = "Gehan")

    expect_equal(as.double(Gehan@count.favorable), c(9,0))
    expect_equal(as.double(Gehan@count.unfavorable), c(2,0))
    expect_equal(as.double(Gehan@count.neutral), c(1,5))
    expect_equal(as.double(Gehan@count.uninf), c(4,0))

    test <- aggrTableIndividualScore(table = Gehan@tableIndividualScore,
                                correct.tte = Gehan@method.tte$correction)
        
    expect_equal(unname(tail(Gehan@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(Gehan@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])

})

test_that("1 TTE endpoint - Gehan (correction IPCW)", {
    ## survival first    
    GehanC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Gehan IPCW")

    factor <- 16/12 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(GehanC@count.favorable), c(9*factor,0))
    expect_equal(as.double(GehanC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(GehanC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(GehanC@count.uninf), c(0,0))

    test <- aggrTableIndividualScore(table = GehanC@tableIndividualScore,
                                correct.tte = GehanC@method.tte$correction)
        
    expect_equal(unname(tail(GehanC@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(GehanC@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])

    ## survival second
    GehanC2 <- BuyseTest(group ~  cont(score) + tte(survie, censoring = event, threshold = 1),
                         data = df, 
                         method.tte = "Gehan IPCW")
    expect_equal(GehanC@count.favorable[1], GehanC2@count.favorable[2])
    expect_equal(GehanC@count.unfavorable[1], GehanC2@count.unfavorable[2])
    expect_equal(GehanC@count.neutral[1], GehanC@count.neutral[2])
    expect_equal(GehanC@count.uninf[1], GehanC2@count.uninf[2])

    test <- aggrTableIndividualScore(table = GehanC2@tableIndividualScore,
                                correct.tte = GehanC2@method.tte$correction)
        
    expect_equal(unname(tail(GehanC2@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(GehanC2@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])

})

## ** Peron
test_that("1 TTE endpoint - Peron (no correction)", {
    Peron <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                       data = df, 
                       method.tte = "Peron")

    expect_equal(as.double(Peron@count.favorable), c(10,0))
    expect_equal(as.double(Peron@count.unfavorable), c(2,0))
    expect_equal(as.double(Peron@count.neutral), c(1,4))
    expect_equal(as.double(Peron@count.uninf), c(3,0))

    test <- aggrTableIndividualScore(table = Peron@tableIndividualScore,
                                correct.tte = Peron@method.tte$correction)
        
    expect_equal(unname(tail(Peron@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(Peron@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])

})

    
test_that("1 TTE endpoint - Peron (IPCW)", {
    ## survival first
    PeronC <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 1) + cont(score),
                        data = df, 
                        method.tte = "Peron IPCW")

    factor <- 16/13 ## n.pairs/(n.pairs-n.uninf)
    expect_equal(as.double(PeronC@count.favorable), c(10*factor,0))
    expect_equal(as.double(PeronC@count.unfavorable), c(2*factor,0))
    expect_equal(as.double(PeronC@count.neutral), c(1*factor,1*factor))
    expect_equal(as.double(PeronC@count.uninf), c(0,0))

    test <- aggrTableIndividualScore(table = PeronC@tableIndividualScore,
                                correct.tte = PeronC@method.tte$correction)
        
    expect_equal(unname(tail(PeronC@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(PeronC@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])

    ## survival second
    PeronC2 <- BuyseTest(group ~  cont(score) + tte(survie, censoring = event, threshold = 1),
                         data = df, 
                         method.tte = "Peron IPCW")
    expect_equal(PeronC@count.favorable[1], PeronC2@count.favorable[2])
    expect_equal(PeronC@count.unfavorable[1], PeronC2@count.unfavorable[2])
    expect_equal(PeronC@count.neutral[1], PeronC@count.neutral[2])
    expect_equal(PeronC@count.uninf[1], PeronC2@count.uninf[2])

    test <- aggrTableIndividualScore(table = PeronC2@tableIndividualScore,
                                correct.tte = PeronC2@method.tte$correction)
        
    expect_equal(unname(tail(PeronC2@Delta.netChance,1)),test[,mean(favorable-unfavorable)])
    expect_equal(unname(tail(PeronC2@Delta.winRatio,1)),test[,sum(favorable)/sum(unfavorable)])
})


##----------------------------------------------------------------------
### test-BuyseTest-correctionTTE.R ends here
