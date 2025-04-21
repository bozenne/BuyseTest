### test-auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:55) 
## Version: 
## Last-Updated: Apr 21 2025 (12:09) 
##           By: Brice Ozenne
##     Update #: 62
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
    library(pROC)
    library(cvAUC)
}

context("Check auc calculation vs. cvAUC")

## * Settings
BuyseTest.options(check = TRUE,
                  method.inference = "none",
                  trace = 0)


## * Compare AUC and CI
n <- 200
set.seed(10)
X <- rnorm(n)
dt <- data.table(Y = as.factor(rbinom(n, size = 1, prob = 1/(1+exp(1/2-X)))),
                 X = X,
                 fold = unlist(lapply(1:10,function(iL){rep(iL,n/10)})))

## boxplot(X~Y, data = dt)
## ** no CV
test_that("AUC - BuyseTest vs pROC",{
    ## example from BuyseTest
    test <- BuyseTest::auc(labels = dt$Y, predictions = dt$X, direction = ">", transformation = FALSE, pooling = "mean")
    test2 <- BuyseTest::auc(labels = as.character(dt$Y), predictions = dt$X, direction = ">", transformation = FALSE, pooling = "mean")
    ## GS2 <- cvAUC::cvAUC(predictions = dt$X, labels = dt$Y)
    GS <- pROC::roc(dt$Y, -dt$X, ci = TRUE, direction = ">", quiet = TRUE)
    GS2 <- cvAUC::ci.cvAUC(predictions = dt$X, labels = dt$Y)

    expect_equal(GS$auc[1], as.double(test[test$fold == "global","estimate"]), tol = 1e-6)
    expect_equal(GS$auc[1], as.double(test2[test$fold == "global","estimate"]), tol = 1e-6)
    expect_equal(GS$ci[1], as.double(test[test$fold == "global","lower"]), tol = 1e-3)
    expect_equal(GS$ci[3], as.double(test[test$fold == "global","upper"]), tol = 1e-3)

    expect_equal(GS2$cvAUC[1], as.double(test[test$fold == "global","estimate"]), tol = 1e-6)
    expect_equal(GS2$se, as.double(test[test$fold == "global","se"]), tol = 1e-3)
    expect_equal(GS2$ci[1], as.double(test[test$fold == "global","lower"]), tol = 1e-3)
    expect_equal(GS2$ci[2], as.double(test[test$fold == "global","upper"]), tol = 1e-3)

    ## butils::object2script(test, digit = 6)
    expect_equal(test$estimate, c(0.7054435), tol = 1e-6)
    expect_equal(test$se, c(0.03639812), tol = 1e-6)

    ## example from pROC
    data(aSAH, package = "pROC")
    roc.s100b <- pROC::roc(aSAH$outcome, aSAH$s100b, quiet = TRUE)
    e.BTauc <- BuyseTest::auc(labels = aSAH$outcome, prediction = aSAH$s100b, transformation = TRUE, pooling = "mean")
    expect_equal(e.BTauc[1,"estimate"],as.numeric(pROC::auc(roc.s100b)), tol = 1e-6)
})

## ** with CV
test_that("AUC after CV - BuyseTest vs cvAUC",{

    ## cvAUC incorrectly estimate standard error on unequally sized folds
    if(FALSE){
        dtS <- rbind(cbind(dt[,.(Y,X)], fold3 = "1", fold2 = "1", id = 1:NROW(dt)),
                     cbind(dt[,.(Y,X)], fold3 = "2", fold2 = "2", id = 1:NROW(dt)),
                     cbind(dt[,.(Y,X)], fold3 = "3", fold2 = "2", id = 1:NROW(dt)))
    
        testS3 <- cvAUC::ci.cvAUC(predictions = dtS$X, labels = dtS$Y, folds = dtS$fold3)
        testS2 <- cvAUC::ci.cvAUC(predictions = dtS$X, labels = dtS$Y, folds = dtS$fold2)
        GS <- cvAUC::ci.cvAUC(predictions = dt$X, labels = dt$Y)

        testBT2 <- auc(labels = dtS$Y, prediction = dtS$X,
                       fold = dtS$fold2, observation = 1:NROW(dtS), pooling = "mean")

        testBT3 <- auc(labels = dtS$Y, prediction = dtS$X,
                       fold = dtS$fold3, observation = 1:NROW(dtS), pooling = "mean")

        ## testS3$cvAUC - GS$cvAUC
        ## testS2$cvAUC - GS$cvAUC
        ## testS3$se - GS$se/sqrt(3)
        testS2$se - testS3$se ## this does not seems right
        testS2$se - testBT2$se[3]
        weight <- c(1/2,1/2)
        sum(weight * testBT2$estimate[1:2]) - testBT2$estimate[3]
        sqrt(sum(weight^2 * testBT2$se[1:2]^2)) - testBT2$se[3]
    }


    ## GS0 <- cvAUC::cvAUC(predictions = dt$X, labels = dt$Y, folds = dt$fold) ## gives wrong results as it ignores fold argument
    GS1 <- cvAUC::ci.cvAUC(predictions = dt$X, labels = dt$Y, folds = dt$fold)

    ## remove warning for uncertainty estimation (P-value/confidence intervals will not be valid with only one observation.)
    test <- suppressWarnings(BuyseTest::auc(labels = dt$Y, prediction = dt$X,
                                            fold = dt$fold, observation = 1:NROW(dt), pooling = "mean"))
    
    expect_equal(test[test$fold=="global", "estimate"], GS1$cvAUC, tol = 1e-6)
    expect_equal(test[test$fold=="global", "estimate"], 0.703265, tol = 1e-6)
    
    ## expect_equal(test[test$fold=="global", "se"], GS1$se, tol = 1e-6) ## differs (cf issue in cvAUC above)
    expect_equal(test[test$fold=="global", "se"], 0.0357925, tol = 1e-6)

    ## another example
    ## remove warning for uncertainty estimation (P-value/confidence intervals will not be valid with only one observation.)
    test.auc <- suppressWarnings(BuyseTest::auc(label = c(0,0,0,0,1,1,1,1),
                                                prediction = c(0.3,0.4,0.5,0.6,0.55,0.5,0.6,0.4,0.7,0.8),
                                                observation = c(1,2,3,4,7,1,2,6,7,8),
                                                fold = c(1,1,1,1,1,2,2,2,2,2),
                                                direction = ">"))
    test.auc
    
})



test_that("AUC after CV - ties",{
    y  <- c(1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1)
    fit <- c(0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789)
    fold <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10)
    observation <- c(62, 77, 50, 52, 64, 58, 74, 2, 26, 58, 2, 52, 12, 28, 17, 52, 68, 69, 14, 65, 22, 58, 26, 68, 9, 13, 32, 8, 21, 30, 20, 55, 36, 30, 22, 9, 31, 18, 65, 11)

    test <- suppressWarnings(BuyseTest::auc(labels = y, predictions = fit, fold = fold, observation = observation,
                                            add.halfNeutral = TRUE))
    expect_true(all(abs(test$estimate-0.5)<1e-6))
})

######################################################################
### test-auc.R ends here
