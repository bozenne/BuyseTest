### test-auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:55) 
## Version: 
## Last-Updated: aug 24 2021 (15:04) 
##           By: Brice Ozenne
##     Update #: 40
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
    test <- auc(labels = dt$Y, predictions = dt$X, direction = ">", transformation = FALSE)
    test <- auc(labels = as.character(dt$Y), predictions = dt$X, direction = ">", transformation = FALSE)
    test2 <- cvAUC::cvAUC(predictions = dt$X,
                   labels = dt$Y)
    test3 <- cvAUC::ci.cvAUC(predictions = dt$X,
                      labels = dt$Y)
    GS <- pROC::roc(dt$Y, -dt$X, ci = TRUE, direction = ">")

    expect_equal(GS$auc[1], as.double(test[test$fold == "global","estimate"]), tol = 1e-6)
    expect_equal(GS$ci[1], as.double(test[test$fold == "global","lower"]), tol = 1e-3)
    expect_equal(GS$ci[3], as.double(test[test$fold == "global","upper"]), tol = 1e-3)

    expect_equal(GS$auc[1], as.double(test2$cvAUC), tol = 1e-6)
    expect_equal(GS$auc[1], as.double(test3$cvAUC), tol = 1e-6)
    expect_equal(GS$ci[1], as.double(test3$ci[1]), tol = 1e-2)
    expect_equal(GS$ci[3], as.double(test3$ci[2]), tol = 1e-2)

    ## butils::object2script(test, digit = 6)
    expect_equal(test$estimate, c(0.705443), tol = 1e-6)
    expect_equal(test$se, c(0.036287), tol = 1e-6)

    ## example from pROC
    data(aSAH, package = "pROC")
    roc.s100b <- pROC::roc(aSAH$outcome, aSAH$s100b)
    e.BTauc <- BuyseTest::auc(labels = aSAH$outcome, prediction = aSAH$s100b, transformation = FALSE)
    expect_equal(e.BTauc[1,"estimate"],as.numeric(pROC::auc(roc.s100b)), tol = 1e-6)

})

## ** with CV
test_that("AUC after CV - BuyseTest vs cvAUC",{
    dt$fold0 <- c(rep(1,100),rep(2,100))
    
    test0 <- auc(labels = dt$Y, prediction = dt$X,
                 fold = dt$fold0, observation = 1:NROW(dt))
    GS0 <- cvAUC::ci.cvAUC(predictions = dt$X,
                           labels = dt$Y,
                           folds = dt$fold0)

    expect_equal(test0[test0$fold=="global", "estimate"],
                 GS0$cvAUC, tol = 1e-6)
    expect_equal(test0[test0$fold=="global", "se"],
                 GS0$se, tol = 1e-6)

    ## e.glm <- glm(Y~X, data = dt, family = binomial(link="logit"))
    ## e.glmiid <- attr(riskRegression::predictRisk(e.glm, newdata = dt, iid = TRUE), "iid")
    ## dim(e.glmiid)

    ## test1 <- auc(labels = dt$Y, prediction = dt$X,
    ##              fold = dt$fold0, observation = 1:NROW(dt))
    
})


test.auc <- auc(label = c(0,0,0,0,1,1,1,1),
                prediction = c(0.3,0.4,0.5,0.6,0.55,0.5,0.6,0.4,0.7,0.8),
                observation = c(1,2,3,4,7,1,2,6,7,8),
                fold = c(1,1,1,1,1,2,2,2,2,2),
                direction = ">")

## * example
y  <- c(1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1)
fit <- c(0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789, 0.52631579, 0.52631579, 0.52631579, 0.52631579, 0.53947368, 0.53947368, 0.53947368, 0.53947368, 0.51315789, 0.51315789, 0.51315789, 0.51315789)
fold <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10)
observation <- c(62, 77, 50, 52, 64, 58, 74, 2, 26, 58, 2, 52, 12, 28, 17, 52, 68, 69, 14, 65, 22, 58, 26, 68, 9, 13, 32, 8, 21, 30, 20, 55, 36, 30, 22, 9, 31, 18, 65, 11)

BuyseTest::auc(labels = y, predictions = fit, fold = fold, observation = observation,
               add.halfNeutral = TRUE)


######################################################################
### test-auc.R ends here
