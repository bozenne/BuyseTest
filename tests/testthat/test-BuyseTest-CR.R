### test-BuyseTest-CR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 12 2018 (16:58) 
## Version: 
## Last-Updated: okt 30 2018 (16:21) 
##           By: Brice Ozenne
##     Update #: 9
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

context("Check that BuyseTest with competing risks \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)


alphaE.X <- 2
alphaCR.X <- 1
alphaE.Y <- 3
alphaCR.Y <- 2
n <- 1e2

## * Simulate data
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X), time2 = rexp(n, rate = alphaCR.X), group = "1"),
            data.frame(time1 = rexp(n, rate = alphaE.Y), time2 = rexp(n, rate = alphaCR.Y), group = "2"))
df$time <- pmin(df$time1,df$time2) ## first event
df$event <- (df$time2<df$time1)+1 ## type of event

## * test
test_that("tte = 2 is equivalent to continuous with infty when cause=2", {
    e.BT <- BuyseTest(group ~ tte(time, censoring = event), data = df,
                      method.inference = "none", method.tte = "Gehan",
                      trace = 0)
    ## summary(e.BT)
    df$timeXX <- df$time
    df$timeXX[df$event==2] <- max(df$time)+1
    e.BT.bis <- BuyseTest(group ~ cont(timeXX), data = df,
                          method.inference = "none", trace = 0)
    ## summary(e.BT.bis)
    
    expect_equal(as.double(e.BT@Delta.netBenefit),
                 as.double(e.BT.bis@Delta.netBenefit))
    expect_equal(as.double(e.BT@delta.netBenefit),
                 as.double(e.BT.bis@delta.netBenefit))
    
})
######################################################################
### test-BuyseTest-CR.R ends here
