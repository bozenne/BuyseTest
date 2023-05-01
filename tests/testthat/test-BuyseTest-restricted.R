### test-BuyseTest-restricted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 21 2021 (17:50) 
## Version: 
## Last-Updated: May  1 2023 (10:08) 
##           By: Brice Ozenne
##     Update #: 17
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
    library(survival)
}

context("Check restricted net benefit \n")

## * Setting
BuyseTest.options(pool.strata = "Buyse")

## * survival case (Peron)
n.obs <- 100
tau <- c(0.1,0.5,1,1.5) ## restriction

set.seed(10)
df.data <- as.data.frame(simBuyseTest(n.obs))
## tapply(df.data$eventtime,df.data$treatment,max)
##        C        T 
## 1.464086 1.668044 

test_that("BuyseTest - restriction",{
    ## no threshold
    test <- list(BuyseTest(treatment ~ tte(eventtime,status,restriction = 0.1), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 0.5), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 1), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 1.5), data = df.data, trace = FALSE)) ## after last event in one group

    df.data$eventtime.0.1 <- pmin(df.data$eventtime,0.1)
    df.data$status.0.1 <- ifelse(df.data$eventtime > 0.1, 1, df.data$status)
    df.data$eventtime.0.5 <- pmin(df.data$eventtime,0.5)
    df.data$status.0.5 <- ifelse(df.data$eventtime > 0.5, 1, df.data$status)
    df.data$eventtime.1 <- pmin(df.data$eventtime,1)
    df.data$status.1 <- ifelse(df.data$eventtime > 1, 1, df.data$status)
    df.data$eventtime.1.5 <- pmin(df.data$eventtime,1.5)
    df.data$status.1.5 <- ifelse(df.data$eventtime > 1.5, 1, df.data$status)

    GS <- list(BuyseTest(treatment ~ tte(eventtime.0.1,status.0.1), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.0.5,status.0.5), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.1,status.1), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.1.5,status.1.5), data = df.data, trace = FALSE))
    
    expect_equal(as.double(sapply(GS,coef)),as.double(sapply(test,coef)), tol = 1e-6)
    ## difference in se as the "GS" may modify how the last step is computed ...

    ## GSGS <- rnetbenefit(endpoint = "eventtime", treatment = "treatment", censoring = "status", threshold = .000000001, epsilon = 0.1, data = df.data)
    ## coef(GS[[1]])-GSGS$rdelta

    ## with threshold
    test <- list(BuyseTest(treatment ~ tte(eventtime,status,restriction = 0.1, threshold = 0.1), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 0.5, threshold = 0.1), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 1, threshold = 0.1), data = df.data, trace = FALSE),
                 BuyseTest(treatment ~ tte(eventtime,status,restriction = 1.5, threshold = 0.1), data = df.data, trace = FALSE)) ## after last event in one group

    GS <- list(BuyseTest(treatment ~ tte(eventtime.0.1,status.0.1, threshold = 0.1), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.0.5,status.0.5, threshold = 0.1), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.1,status.1, threshold = 0.1), data = df.data, trace = FALSE),
               BuyseTest(treatment ~ tte(eventtime.1.5,status.1.5, threshold = 0.1), data = df.data, trace = FALSE))
    
    expect_equal(as.double(sapply(GS,coef)),as.double(sapply(test,coef)), tol = 1e-6)

    ## min(df.data$eventtime)
    ## [1] 0.002049558
    test <- BuyseTest(treatment ~ tte(eventtime,status,restriction = 0.000001, threshold = 0.1), data = df.data, trace = FALSE)
    expect_equal(as.double(coef(test)),0, tol = 1e-6)
    ## GSGS <- rnetbenefit(endpoint = "eventtime", treatment = "treatment", censoring = "status", threshold = 0.1, epsilon = 0.1, data = df.data)
    ## coef(GS[[1]])-GSGS$rdelta
    ## GSGS <- rnetbenefit(endpoint = "eventtime", treatment = "treatment", censoring = "status", threshold = 0.1, epsilon = 0.5, data = df.data)
    ## [1] -0.1661543
    ## coef(GS[[2]])-GSGS$rdelta
    ## coef(test[[2]])-GSGS$rdelta
})


##----------------------------------------------------------------------
### test-BuyseTest-restricted.R ends here
