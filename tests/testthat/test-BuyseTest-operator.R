### test-BuyseTest-operator.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  2 2018 (15:21) 
## Version: 
## Last-Updated: jan  8 2021 (15:43) 
##           By: Brice Ozenne
##     Update #: 36
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

context("Check that the option operator in BuyseTest is working correctly \n")


## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * one pair
test_that("check - 1 pair",{
    ## binary 
    data <- data.frame(toxicity1 = c(1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ bin(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)
    ## getPairScore(BT)

    ## continuous
    data <- data.frame(toxicity1 = c(1,0),
                       treatment = c(1,0),
                       stringsAsFactors = FALSE)
    BT <- BuyseTest(treatment ~ cont(toxicity1, operator = "<0"), data=data)
    expect_equal(as.double(coef(BT, statistic = "count.favorable", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.unfavorable", cumulative = FALSE)),1)
    expect_equal(as.double(coef(BT, statistic = "count.neutral", cumulative = FALSE)),0)
    expect_equal(as.double(coef(BT, statistic = "count.uninf", cumulative = FALSE)),0)

})

## * many pairs
data(veteran, package = "survival")
test_that("check - many pairs",{
    ## original - twice same tte
    BT0 <- BuyseTest(trt ~ tte(time,status,threshold=5)+tte(time,status,threshold=0),
                     data=veteran, method.inference = "u-statistic")
    ## switched - twice same tte
    BuyseTest.options(engine = "GPC_cpp")
    BT <- BuyseTest(trt ~ tte(time,status,threshold=5,operator="<0")+tte(time,status,threshold=0,operator="<0"),
                    data=veteran, method.inference = "u-statistic")
    BuyseTest.options(engine = "GPC2_cpp")
    BT2 <- BuyseTest(trt ~ tte(time,status,threshold=5,operator="<0")+tte(time,status,threshold=0,operator="<0"),
                     data=veteran, method.inference = "u-statistic")

    expect_equal(confint(BT0)[,"estimate"], -confint(BT)[,"estimate"], tol = 1e-6)
    expect_equal(confint(BT0)[,c("se","p.value")], confint(BT)[,c("se","p.value")], tol = 1e-6)
    expect_equal(confint(BT), confint(BT2))
    expect_equal(confint(BT), confint(BT2))

    ## original -  tte + bin
    BT0 <- BuyseTest(trt ~ tte(time,status,threshold=5)+bin(prior),
                     data=veteran, method.inference = "u-statistic")

    ## switched -  tte + bin
    BuyseTest.options(engine = "GPC_cpp")
    BT <- BuyseTest(trt ~ tte(time,status,threshold=5,operator="<0")+bin(prior,operator="<0"),
                    data=veteran, method.inference = "u-statistic")
    BuyseTest.options(engine = "GPC2_cpp")
    BT2 <- BuyseTest(trt ~ tte(time,status,threshold=5,operator="<0")+bin(prior,operator="<0"),
                     data=veteran, method.inference = "u-statistic")

    expect_equal(confint(BT0)[,"estimate"], -confint(BT)[,"estimate"], tol = 1e-6)
    expect_equal(confint(BT0)[,c("se","p.value")], confint(BT)[,c("se","p.value")], tol = 1e-2)
    expect_equal(confint(BT), confint(BT2), tol = 1e-2)
    expect_equal(confint(BT), confint(BT2), tol = 1e-2)
    ## expected difference in se because GPC2 compute the influence function over all pairs,
    ## while GPC only over pairs with informative scores.
    ## the difference is expected to be small though in large samples

})

##----------------------------------------------------------------------
### test-BuyseTest-operator.R ends here
