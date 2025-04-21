### test-BuyseTest-standardization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2025 (15:42) 
## Version: 
## Last-Updated: Apr 21 2025 (14:12) 
##           By: Brice Ozenne
##     Update #: 59
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
    library(prodlim)
}

context("Check standardization of the net benefit across strata \n")

## * Setting
BuyseTest.options(order.Hprojection = 1,
                  method.inference = "u-statistic",
                  trace = 0)

## * Binary outcome and binary covariate

## ** simulate data
seed <- 1
n <- 50
tau <- 0

set.seed(seed)
A   <- rbinom(n, size = 1, prob = 0.75)           
X   <- sample(0:1, size = n, replace = TRUE)      
df <- data.frame(index = 1:n,
                 Y = 0.5 * A + X + rnorm(n),
                 arm = A, X = X)

## all possible pairs of strata
df00 <- df[df$X==0,]
df11 <- df[df$X==1,]
df01 <- rbind(df[df$X==0 & df$arm==0,],df[df$X==1 & df$arm==1,])
df10 <- rbind(df[df$X==0 & df$arm==1,],df[df$X==1 & df$arm==0,])

## ** standardization (first order)
test_that("BuyseTest - standarization (first order, binary outcome)",{

    BT.std <- BuyseTest(arm ~ cont(Y, threshold = tau) + X, trace = FALSE, data = df, pool.strata = "standardization")
    expect_equal(as.double(confint(BT.std)), c(0.05811117, 0.19991436, -0.3229847, 0.42299107, 0, 0.77179675), tol = 1e-5)
    expect_equal(as.double(BT.std@covariance[,"netBenefit"]), 0.03996575, tol = 1e-5)
    expect_equal(mean(BT.std@iidAverage$favorable[,1]), 0, tol = 1e-5)
    expect_equal(mean(BT.std@iidAverage$unfavorable[,1]), 0, tol = 1e-5)

    ## similar to weighting observations except that the iid gets wrong
    df$w <- as.numeric((table(df$arm)[as.character(df$arm)]*table(df$X)[as.character(df$X)])/(NROW(df)*table(interaction(df$arm,df$X))[interaction(df$arm,df$X)]))
    BT.std2 <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = df, weightObs = df$w)
    expect_equal(as.double(coef(BT.std)), as.double(coef(BT.std2)), tol =  1e-5)
    BT.std2@covariance[,"netBenefit"]
    ## [1] 0.04562422

    ## retrieve standardization result via strata-specific analyses
    BT.strat00 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df00)
    BT.strat11 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df11)
    BT.strat01 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df01)
    BT.strat10 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df10)

    BT.res <- rbind(data.frame(X = 0, Xstar = 0, Delta = coef(BT.strat00), nC.X = sum(df00$arm==0), nE.Xstar = sum(df00$arm==1)),
                    data.frame(X = 1, Xstar = 1, Delta = coef(BT.strat11), nC.X = sum(df11$arm==0), nE.Xstar = sum(df11$arm==1)),
                    data.frame(X = 0, Xstar = 1, Delta = coef(BT.strat01), nC.X = sum(df01$arm==0), nE.Xstar = sum(df01$arm==1)),
                    data.frame(X = 1, Xstar = 0, Delta = coef(BT.strat10), nC.X = sum(df10$arm==0), nE.Xstar = sum(df10$arm==1)))
    BT.res$n.X <- as.numeric(table(df$X)[as.character(BT.res$X)])
    BT.res$n.Xstar <- as.numeric(table(df$X)[as.character(BT.res$Xstar)])
    BT.res
    ##   X Xstar       Delta nC.X nE.Xstar n.X n.Xstar
    ## 1 0     0 -0.08571429    7       15  22      22
    ## 2 1     1  0.12121212    6       22  28      28
    ## 3 0     1  0.41558442    7       22  22      28
    ## 4 1     0 -0.26666667    6       15  28      22

    Myiid <- matrix(0, nrow = n, ncol = 4)
    Myiid[df00$index,1] <- getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)]
    Myiid[df11$index,2] <- getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)]
    Myiid[df01$index,3] <- getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)]
    Myiid[df10$index,4] <- getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)]

    Myiid.norm <- sweep(Myiid, MARGIN = 1, FUN = "/", STATS = table(df$arm)[as.character(1-df$arm)])


    ## *** version 1 (joint normal)
    BT.res$w.unadjusted <- BT.res$nC.X * BT.res$nE.Xstar / prod(table(df$arm))
    BT.res$w.standardize <- BT.res$n.X * BT.res$n.Xstar / NROW(df$arm)^2
    expect_equal(as.double(coef(BT.std)), as.double(BT.res$Delta %*% BT.res$w.standardize), tol = 1e-5)
    MyiidSTD.norm <- sweep(Myiid.norm, MARGIN = 2, FUN = "*", STATS = BT.res$w.standardize/BT.res$w.unadjusted)
    expect_equal(as.double(coef(BT.std)), as.double(mean(rowSums(MyiidSTD.norm))), tol = 1e-5)
    expect_equal(as.double(BT.std@covariance[,"netBenefit"]),
                 as.double(crossprod(as.vector(rowSums(MyiidSTD.norm) - mean(rowSums(MyiidSTD.norm)))/table(df$arm)[as.character(df$arm)])), tol = 1e-5)
    

    ## *** version 2 (compact)
    myiid <- rep(0, n)
    myiid[df00$index] <- myiid[df00$index] + getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)] * BT.res$w.standardize[1]/BT.res$w.unadjusted[1]
    myiid[df11$index] <- myiid[df11$index] + getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)] * BT.res$w.standardize[2]/BT.res$w.unadjusted[2]
    myiid[df01$index] <- myiid[df01$index] + getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)] * BT.res$w.standardize[3]/BT.res$w.unadjusted[3]
    myiid[df10$index] <- myiid[df10$index] + getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)] * BT.res$w.standardize[4]/BT.res$w.unadjusted[4]

    myiid <- myiid/table(df$arm)[as.character(1-df$arm)]
    expect_equal(as.double(coef(BT.std)), as.double(mean(myiid)), tol = 1e-5)
    expect_equal(as.double(BT.std@covariance[,"netBenefit"]),
                 as.double(crossprod((myiid - mean(myiid))/table(df$arm)[as.character(df$arm)]) ), tol = 1e-5)
    
    ## ** retrive marginal estimate after stratification
    if(FALSE){ ## sanity check
        ## *** marginal estimate
        BT.raw <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = df)
        confint(BT.raw)
        ##    estimate        se  lower.ci  upper.ci null   p.value
        ## Y 0.0977131 0.1908652 -0.272599 0.4428102    0 0.6109743
        mean(getIid(BT.raw, center = FALSE, scale = FALSE))
        ## [1] 0.0977131
        crossprod(getIid(BT.raw))
        ##            [,1]
        ## [1,] 0.03642954
        crossprod(((getIid(BT.raw, center = FALSE, scale = FALSE) * table(df$arm)[as.character(1-df$arm)]) / table(df$arm)[as.character(1-df$arm)] - coef(BT.raw))/table(df$arm)[as.character(df$arm)])
        ## 0.03642954
        c(BT.raw@covariance[,"netBenefit"],sqrt(BT.raw@covariance[,"netBenefit"]))
        ## [1] 0.03642954 0.19086523

        ## *** marginalization version 1 (joint normal)
        BT.res$Delta %*% BT.res$w.unadjusted
        ## [1] 0.0977131
        range(getIid(BT.raw, center = FALSE, scale = FALSE)-rowSums(Myiid.norm))
        ## [1] -3.330669e-16  1.110223e-16
        mean(rowSums(Myiid.norm))
        ## [1] 0.0977131
        crossprod(as.vector(rowSums(Myiid.norm) - mean(rowSums(Myiid.norm)))/table(df$arm)[as.character(df$arm)])
        ##            [,1]
        ## [1,] 0.03642954

        ## *** marginalization version 2 (compact)
        myiid <- rep(0, n)
        myiid[df00$index] <- myiid[df00$index] + getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)] 
        myiid[df11$index] <- myiid[df11$index] + getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)] 
        myiid[df01$index] <- myiid[df01$index] + getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)] 
        myiid[df10$index] <- myiid[df10$index] + getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)] 
        myiid <- myiid / table(df$arm)[as.character(1-df$arm)]

        mean(myiid)
        ## [1] 0.0977131
        crossprod((myiid - mean(myiid))/table(df$arm)[as.character(df$arm)]) 
        ##            [,1]
        ## [1,] 0.03642954
    }
})

## ** standardization (second order)
BuyseTest.options(order.Hprojection = 2)

test_that("BuyseTest - standarization (second order, binary outcome)",{
    ## with balanced dataset standardization equals marginal analysis
    dfB <- do.call(rbind,by(df, interaction(df$arm,df$X), function(iDF){iDF[1:min(table(df$arm,df$X)),,drop=FALSE]}))
    rownames(dfB) <- NULL

    BTB.raw <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = dfB)
    BTB.std <- BuyseTest(arm ~ cont(Y, threshold = tau) + X, trace = FALSE, data = dfB, pool.strata = "standardization")
    expect_equal(as.double(confint(BTB.raw)),as.double(confint(BTB.std)))

})

BuyseTest.options(order.Hprojection = 1)

## * check Peron scoring rule in a perfectly balanced dataset
n <- 50

set.seed(10)
dt.strata <- simBuyseTest(n,
                          n.strata = 2,
                          names.strata = "stratum")

dtB.strata <- do.call(rbind,by(dt.strata, interaction(dt.strata$treatment,dt.strata$stratum), function(iDT){iDT[1:min(table(dt.strata$treatment,dt.strata$stratum))]}))
rownames(dtB.strata) <- NULL


test_that("BuyseTest - standarization (tte outcome)",{
    BTB.raw <- BuyseTest(treatment ~ tte(eventtime, status), trace = FALSE, data = dtB.strata,
                         model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
    ##             estimate        se   lower.ci  upper.ci null   p.value
    ## eventtime 0.02881503 0.1417496 -0.2441965 0.2975942    0 0.8390032

    BTB.std <- BuyseTest(treatment ~ tte(eventtime, status) + stratum, trace = FALSE, data = dtB.strata,
                         pool.strata = "standardization", 
                         model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
    expect_equal(as.double(confint(BTB.raw)), as.double(confint(BTB.std)), tol = 1e-5)
    expect_equal(as.double(confint(BTB.raw)), c(0.028815, 0.1417496, -0.2441965, 0.2975942, 0, 0.8390032), tol = 1e-5)
})


##----------------------------------------------------------------------
### test-BuyseTest-standardization.R ends here
