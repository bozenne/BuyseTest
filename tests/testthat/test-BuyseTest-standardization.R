### test-BuyseTest-standardization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2025 (15:42) 
## Version: 
## Last-Updated: Jul  7 2025 (12:30) 
##           By: Brice Ozenne
##     Update #: 109
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
n <- 100
tau <- 0

set.seed(seed)
A   <- rbinom(n, size = 1, prob = 0.75)           
X   <- sample(0:1, size = n, replace = TRUE)      
df <- data.frame(id = 1:n,
                 Y = 0.5 * A + X + rnorm(n),
                 arm = A, X = X)

## ** standardization (first order)
test_that("BuyseTest - standarization (first order, binary outcome)",{

    ## 1 - marginal
    ## BT.unstrat         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau), trace = F)
    ## BT.unstrat2        <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau), trace = TRUE, method.inference = "bootstrap")
    ## rbind(confint(BT.unstrat, transformation = FALSE),
    ##       confint(BT.unstrat2, method.ci = "gaussian"))
    ##     estimate        se    lower.ci  upper.ci null   p.value
    ## Y  0.1811263 0.1188404 -0.05179662 0.4140493    0 0.1274803
    ## Y1 0.1811263 0.1220947 -0.06875553 0.4096230    0 0.1543344

    ## 2- standardized
    BT.std         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = FALSE)
    expect_equal(mean(getIid(BT.std)),0, tol = 1e-6)
    expect_equivalent(confint(BT.std, transformation = FALSE),
                      data.frame("estimate" = c(0.23088881),"se" = c(0.10551968),"lower.ci" = c(0.02407403),"upper.ci" = c(0.43770359),"null" = c(0),"p.value" = c(0.02866149)),,
                      tol = 1e-6)
    expect_equivalent(coef(BT.std, strata = TRUE), c("0" = 0.03658537, "1.0" = -0.30406504, "0.1" = 0.75, "1" = 0.49583333))
    expect_error(confint(BT.std, strata = TRUE)) ## cannot compute strata specific uncertainty since the same individual appear in sevearl strata

    BT.std2         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = FALSE,
                                 method.inference = "bootstrap", strata.resampling = c("arm","X"), seed = 10)
    expect_equivalent(confint(BT.std2, method.ci = "gaussian", transformation = FALSE),
                      data.frame("estimate" = c(0.23088881),"se" = c(0.10915709),"lower.ci" = c(0.01694484),"upper.ci" = c(0.44483278),"null" = c(0),"p.value" = c(0.03441312)),
                      tol = 1e-6)

    ## 3- stratified to stratified (by hand)
    df00 <- df[df$X==0,]
    df10 <- rbind(df[df$X==0 & df$arm==1,],df[df$X==1 & df$arm==0,])
    df01 <- rbind(df[df$X==0 & df$arm==0,],df[df$X==1 & df$arm==1,])
    df11 <- df[df$X==1,]

    BT.strat00 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df00)
    BT.strat10 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df10)
    BT.strat01 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df01)
    BT.strat11 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df11)

    grid.strata <- expand.grid(X = sort(unique(df$X)), Xstar = sort(unique(df$X)))
    grid.strata$X.X <- paste(grid.strata$X,grid.strata$Xstar,sep=".")
    df.extended <- do.call(rbind,apply(grid.strata, MARGIN = 1, function(iRow){
        cbind(rbind(df[df$arm == 0 & df$X == iRow[1],],
                    df[df$arm == 1 & df$X == iRow[2],]),
              X.X = paste(iRow[1],iRow[2],sep="."))
    }))
    grid.strata$nC.X <- tapply(df.extended$arm==0,df.extended$X.X,sum)[grid.strata$X.X]
    grid.strata$nE.Xstar <- tapply(df.extended$arm==1,df.extended$X.X,sum)[grid.strata$X.X]
    grid.strata$n.X <- as.numeric(table(df$X)[as.character(grid.strata$X)])
    grid.strata$n.Xstar <- as.numeric(table(df$X)[as.character(grid.strata$Xstar)])
    grid.strata$w.unadjusted <- grid.strata$nC.X * grid.strata$nE.Xstar / prod(table(df$arm))
    grid.strata$w.standardize <- grid.strata$n.X * grid.strata$n.Xstar / NROW(df$arm)^2

    ## direct
    Myiid <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
    Myiid[df00$id,"0.0"] <- getIid(BT.strat00)
    Myiid[df10$id,"1.0"] <- getIid(BT.strat10)
    Myiid[df01$id,"0.1"] <- getIid(BT.strat01)
    Myiid[df11$id,"1.1"] <- getIid(BT.strat11)
    MyWiid <- sweep(Myiid, FUN = "*", MARGIN = 2, STATS = grid.strata$w.standardize)
    
    ## mean(rowSums(MyWiid))
    expect_equivalent(BT.std@covariance[,"netBenefit"], crossprod(rowSums(MyWiid)), tol = 1e-6)
    
    ## from the partial win/loss sums
    Myiid2 <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
    Myiid2[df00$id,"0.0"] <- getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)]
    Myiid2[df10$id,"1.0"] <- getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)]
    Myiid2[df01$id,"0.1"] <- getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)]
    Myiid2[df11$id,"1.1"] <- getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)]

    ## MyWiid2 <- sweep(Myiid2, FUN = "*", MARGIN = 2, STATS = grid.strata$w.standardize)
    ## MyWiid2.norm <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
    ## MyWiid2.norm[df00$id,"0.0"] <- ((MyWiid2[df00$id,"0.0"] / table(df00$arm)[as.character(1-df00$arm)]) - coef(BT.strat00)*BT.std@weightStrata[1])/table(df00$arm)[as.character(df00$arm)]
    ## MyWiid2.norm[df10$id,"1.0"] <- ((MyWiid2[df10$id,"1.0"] / table(df10$arm)[as.character(1-df10$arm)]) - coef(BT.strat10)*BT.std@weightStrata[2])/table(df10$arm)[as.character(df10$arm)]
    ## MyWiid2.norm[df01$id,"0.1"] <- ((MyWiid2[df01$id,"0.1"] / table(df01$arm)[as.character(1-df01$arm)]) - coef(BT.strat01)*BT.std@weightStrata[3])/table(df01$arm)[as.character(df01$arm)]
    ## MyWiid2.norm[df11$id,"1.1"] <- ((MyWiid2[df11$id,"1.1"] / table(df11$arm)[as.character(1-df11$arm)]) - coef(BT.strat11)*BT.std@weightStrata[4])/table(df11$arm)[as.character(df11$arm)]
    ## mean(rowSums(MyWiid2.norm))
    ## crossprod(rowSums(MyWiid2.norm))

    ## this is what is implemented in the C++ function
    MyWiid3 <- rowSums(sweep(Myiid2, FUN = "*", MARGIN = 2, STATS = grid.strata$w.standardize/(grid.strata$w.unadjusted*prod(table(df$arm)))))
    ## grid.strata$w.unadjusted*prod(table(df$arm)) - c(prod(table(df00$arm)),prod(table(df10$arm)),prod(table(df01$arm)),prod(table(df11$arm)))
    MyWiid3[df00$id] <- MyWiid3[df00$id] - coef(BT.strat00)*BT.std@weightStrata[1]/table(df00$arm)[as.character(df00$arm)]
    MyWiid3[df10$id] <- MyWiid3[df10$id] - coef(BT.strat10)*BT.std@weightStrata[2]/table(df10$arm)[as.character(df10$arm)]
    MyWiid3[df01$id] <- MyWiid3[df01$id] - coef(BT.strat01)*BT.std@weightStrata[3]/table(df01$arm)[as.character(df01$arm)]
    MyWiid3[df11$id] <- MyWiid3[df11$id] - coef(BT.strat11)*BT.std@weightStrata[4]/table(df11$arm)[as.character(df11$arm)]
    expect_equal(mean(MyWiid3), 0, tol = 1e-6)
    expect_equivalent(BT.std@covariance[,"netBenefit"], crossprod(MyWiid3), tol = 1e-6)
})

## ** standardization (second order)
BuyseTest.options(order.Hprojection = 2)

test_that("BuyseTest - standarization (second order, binary outcome, balanced)",{
    ## with balanced dataset standardization equals marginal analysis
    dfB <- do.call(rbind,by(df, interaction(df$arm,df$X), function(iDF){iDF[1:min(table(df$arm,df$X)),,drop=FALSE]}))
    dfB$id <- 1:NROW(dfB)
    rownames(dfB) <- NULL

    BTB.marginal         <- BuyseTest(data = dfB, arm ~ cont(Y, threshold = tau), trace = FALSE)
    BTB.marginal2         <- BuyseTest(data = dfB, arm ~ cont(Y, threshold = tau), trace = F, keep.pairScore = TRUE)
    BTB.std         <- BuyseTest(data = dfB, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = F, keep.pairScore = FALSE)
    BTB.std2         <- BuyseTest(data = dfB, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = F, keep.pairScore = TRUE)
    expect_equal(coef(BTB.marginal), coef(BTB.std), tol = 1e-5)
    expect_equivalent(coef(BTB.std), mean(coef(BTB.std, strata = TRUE)), tol = 1e-5)

    expect_true(confint(BTB.marginal)$se > confint(BTB.marginal, order.Hprojection = 1)$se) ## confint(BTB.marginal)$se^2 - confint(BTB.marginal, order.Hprojection = 1)$se^2
    expect_true(confint(BTB.std)$se > confint(BTB.std, order.Hprojection = 1)$se) ## confint(BTB.std)$se^2 - confint(BTB.std, order.Hprojection = 1)$se^2
    expect_equivalent(confint(BTB.std), confint(BTB.std2), tol = 1e-6)

    ## by hand (marginal)
    sigmaMB_h1fav <- tapply(getIid(BTB.marginal, statistic = "favorable"), dfB$arm, crossprod)
    termM1 <- (coef(BTB.marginal, statistic = "favorable") - coef(BTB.marginal, statistic = "favorable")^2)/prod(table(dfB$arm)) - sum(sigmaMB_h1fav/rev(table(dfB$arm)))
    sigmaMB_h1unfav <- tapply(getIid(BTB.marginal, statistic = "unfavorable"), dfB$arm, crossprod)
    termM2 <- (coef(BTB.marginal, statistic = "unfavorable") - coef(BTB.marginal, statistic = "unfavorable")^2)/prod(table(dfB$arm)) - sum(sigmaMB_h1unfav/rev(table(dfB$arm)))
    sigmaMB_h1cross <- tapply(getIid(BTB.marginal, statistic = "favorable")*getIid(BTB.marginal, statistic = "unfavorable"), dfB$arm, sum)
    termM3 <- 0 - coef(BTB.marginal, statistic = "favorable")*coef(BTB.marginal, statistic = "unfavorable")/prod(table(dfB$arm)) - sum(sigmaMB_h1cross/rev(table(dfB$arm)))
    expect_equivalent(termM1+termM2-2*termM3, confint(BTB.marginal, order.Hprojection = 2)$se^2 - confint(BTB.marginal, order.Hprojection = 1)$se^2, tol = 1e-6)

    mytableM <- getPairScore(BTB.marginal2)
    HprojM_fav <- getIid(BTB.marginal2, statistic = "favorable", scale = FALSE)
    estM_fav <- coef(BTB.marginal2, statistic = "favorable")
    mytableM[, hproj.0 := HprojM_fav[index.0]]
    mytableM[, hproj.1 := HprojM_fav[index.1]]
    expect_equivalent(termM1, mytableM[, sum((favorable - hproj.0 - hproj.1 - estM_fav)^2)/.N^2], tol = 1e-6)

    ## by hand (standardized)
    sigmaB_h1fav <- tapply(getIid(BTB.std, statistic = "favorable"), dfB$arm, crossprod)
    term1 <- (coef(BTB.std, statistic = "favorable") - coef(BTB.std, statistic = "favorable")^2)/prod(table(dfB$arm)) - sum(sigmaB_h1fav/rev(table(dfB$arm)))
    sigmaB_h1unfav <- tapply(getIid(BTB.std, statistic = "unfavorable"), dfB$arm, crossprod)
    term2 <- (coef(BTB.std, statistic = "unfavorable") - coef(BTB.std, statistic = "unfavorable")^2)/prod(table(dfB$arm)) - sum(sigmaB_h1unfav/rev(table(dfB$arm)))
    sigmaB_h1cross <- tapply(getIid(BTB.std, statistic = "favorable")*getIid(BTB.std, statistic = "unfavorable"), dfB$arm, sum)
    term3 <- 0 - coef(BTB.std, statistic = "favorable")*coef(BTB.std, statistic = "unfavorable")/prod(table(dfB$arm)) - sum(sigmaB_h1cross/rev(table(dfB$arm)))
    expect_equivalent(term1+term2-2*term3, confint(BTB.std, order.Hprojection = 2)$se^2 - confint(BTB.std, order.Hprojection = 1)$se^2, tol = 1e-6)

    mytable <- getPairScore(BTB.std2)
    mytable[, weightStrata := stats::setNames(BTB.std2@weightStrata,BTB.std2@level.strata)[strata]]
    mytable[, m.X := length(unique(index.0)), by = "strata"]
    mytable[, n.X := length(unique(index.1)), by = "strata"]
    mytable[, m := length(attr(BTB.std2@level.treatment,"indexC"))]
    mytable[, n := length(attr(BTB.std2@level.treatment,"indexT"))]
    mytable[, Wfavorable := favorable * weightStrata/((m.X*n.X)/(n*m))]
    ## mean(mytable$Wfavorable)
    ## coef(BTB.std2, statistic = "favorable")

    Hproj_fav <- getIid(BTB.std2, statistic = "favorable", scale = FALSE)
    mytable[, hproj.0 := Hproj_fav[index.0]]
    mytable[, hproj.1 := Hproj_fav[index.1]]
    mytable[, Tfavorable := coef(BTB.std2, statistic = "favorable")]

    expect_equivalent(term1, mytable[, sum((Wfavorable - hproj.0 - hproj.1 - Tfavorable)^2)/.N^2], tol = 1e-6)
})    
    
test_that("BuyseTest - standarization (second order, binary outcome, unbalanced)",{

    ## 2- standardized
    BTUB.std         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = FALSE)
    BTUB.std2         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", keep.pairScore = TRUE, trace = FALSE)

    expect_equivalent(confint(BTUB.std), confint(BTUB.std2), tol = 1e-6)
    expect_equivalent(confint(BTUB.std), data.frame("estimate" = c(0.23088881),"se" = c(0.10675249),"lower.ci" = c(0.01411401),"upper.ci" = c(0.42693406),"null" = c(0),"p.value" = c(0.03705691)), tol = 1e-6)

    ## by hand (standardized)
    mytableUB <- getPairScore(BTUB.std)
    mytableUB[, weightStrata := stats::setNames(BTUB.std@weightStrata,BTUB.std@level.strata)[strata]]
    mytableUB[, m.X := length(unique(index.0)), by = "strata"]
    mytableUB[, n.X := length(unique(index.1)), by = "strata"]
    mytableUB[, m := length(attr(BTUB.std@level.treatment,"indexC"))]
    mytableUB[, n := length(attr(BTUB.std@level.treatment,"indexT"))]
    mytableUB[, Wfavorable := favorable * weightStrata/((m.X*n.X)/(n*m))]
    ## mean(mytableUB$Wfavorable)
    ## coef(BTUB.std, statistic = "favorable")
    Hproj_fav <- getIid(BTUB.std2, statistic = "favorable", scale = FALSE)
    mytableUB[, hproj.0 := Hproj_fav[index.0]]
    mytableUB[, hproj.1 := Hproj_fav[index.1]]
    mytableUB[, Tfavorable := coef(BTUB.std2, statistic = "favorable")]
 
    term1 <- confint(BTUB.std, statistic = "favorable", order.Hprojection = 2)$se^2 - confint(BTUB.std, statistic = "favorable", order.Hprojection = 1)$se^2
    expect_equal(term1, mytableUB[, sum((Wfavorable - hproj.0 - hproj.1 - Tfavorable)^2)/.N^2], tol = 1e-6)

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
    BT.raw <- BuyseTest(treatment ~ tte(eventtime, status), trace = FALSE, data = dt.strata)
    ## confint(BT.raw)
    ##             estimate        se   lower.ci  upper.ci null   p.value
    ## eventtime -0.008226969 0.1286523 -0.2546679 0.2392174    0 0.9490145
    BT.std <- BuyseTest(treatment ~ tte(eventtime, status) + stratum, trace = FALSE, data = dt.strata,
                        pool.strata = "standardization")
    ##               estimate        se   lower.ci  upper.ci null   p.value
    ## eventtime -0.009499187 0.1207078 -0.2412526 0.2232793    0 0.9372784
    expect_equivalent(confint(BT.std), data.frame("estimate" = c(-0.00949919),"se" = c(0.12070778),"lower.ci" = c(-0.24125263),"upper.ci" = c(0.22327925),"null" = c(0),"p.value" = c(0.9372784)), tol = 1e-5)

    BTB.raw <- BuyseTest(treatment ~ tte(eventtime, status), trace = FALSE, data = dtB.strata,
                         model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
    ## confint(BTB.raw)
    ##             estimate        se   lower.ci  upper.ci null   p.value
    ## eventtime 0.02881503 0.1417496 -0.2441965 0.2975942    0 0.8390032

    BTB.std <- BuyseTest(treatment ~ tte(eventtime, status) + stratum, trace = FALSE, data = dtB.strata,
                         pool.strata = "standardization", 
                         model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
    ##             estimate        se   lower.ci  upper.ci null   p.value
    ## eventtime 0.02881503 0.1325542 -0.2271614 0.2810671    0 0.8280037
    expect_equivalent(confint(BTB.std), data.frame("estimate" = c(0.02881503),"se" = c(0.13255421),"lower.ci" = c(-0.22716139),"upper.ci" = c(0.28106715),"null" = c(0),"p.value" = c(0.82800367)), tol = 1e-5)
})

##----------------------------------------------------------------------
### test-BuyseTest-standardization.R ends here
