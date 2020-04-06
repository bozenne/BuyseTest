### test-auc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2019 (16:55) 
## Version: 
## Last-Updated: apr  6 2020 (11:44) 
##           By: Brice Ozenne
##     Update #: 29
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

context("Check auc calculation vs. cvAUC")
## library(pROC)
## library(cvAUC)

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
    test <- auc(labels = dt$Y, predictions = dt$X, direction = ">")
    test <- auc(labels = as.character(dt$Y), predictions = dt$X, direction = ">")
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
    expect_equal(test$estimate, c(0.705443, 0.705443), tol = 1e-6)
    expect_equal(test$se, c(0.036287, 0.036287), tol = 1e-6)
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

######################################################################
### test-auc.R ends here
