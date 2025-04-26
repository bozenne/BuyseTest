### test-BuyseTest-standardization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2025 (15:42) 
## Version: 
## Last-Updated: apr 25 2025 (16:38) 
##           By: Brice Ozenne
##     Update #: 61
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
## test_that("BuyseTest - standarization (first order, binary outcome)",{

##     BT.std <- BuyseTest(arm ~ cont(Y, threshold = tau) + X, trace = FALSE, data = df, pool.strata = "standardization")
##     expect_equal(as.double(confint(BT.std)), c(0.05811117, 0.19991436, -0.3229847, 0.42299107, 0, 0.77179675), tol = 1e-5)
##     expect_equal(as.double(BT.std@covariance[,"netBenefit"]), 0.03996575, tol = 1e-5)
##     expect_equal(mean(BT.std@iidAverage$favorable[,1]), 0, tol = 1e-5)
##     expect_equal(mean(BT.std@iidAverage$unfavorable[,1]), 0, tol = 1e-5)

##     ## similar to weighting observations except that the iid gets wrong
##     df$w <- as.numeric((table(df$arm)[as.character(df$arm)]*table(df$X)[as.character(df$X)])/(NROW(df)*table(interaction(df$arm,df$X))[interaction(df$arm,df$X)]))
##     BT.std2 <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = df, weightObs = df$w)
##     expect_equal(as.double(coef(BT.std)), as.double(coef(BT.std2)), tol =  1e-5)
##     BT.std2@covariance[,"netBenefit"]
##     ## [1] 0.04562422

##     ## retrieve standardization result via strata-specific analyses
##     BT.strat00 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df00)
##     BT.strat11 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df11)
##     BT.strat01 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df01)
##     BT.strat10 <- BuyseTest(arm ~ cont(Y, threshold = tau) , method.inference = "U-statistic", trace = F, data = df10)

##     BT.res <- rbind(data.frame(X = 0, Xstar = 0, Delta = coef(BT.strat00), nC.X = sum(df00$arm==0), nE.Xstar = sum(df00$arm==1)),
##                     data.frame(X = 1, Xstar = 1, Delta = coef(BT.strat11), nC.X = sum(df11$arm==0), nE.Xstar = sum(df11$arm==1)),
##                     data.frame(X = 0, Xstar = 1, Delta = coef(BT.strat01), nC.X = sum(df01$arm==0), nE.Xstar = sum(df01$arm==1)),
##                     data.frame(X = 1, Xstar = 0, Delta = coef(BT.strat10), nC.X = sum(df10$arm==0), nE.Xstar = sum(df10$arm==1)))
##     BT.res$n.X <- as.numeric(table(df$X)[as.character(BT.res$X)])
##     BT.res$n.Xstar <- as.numeric(table(df$X)[as.character(BT.res$Xstar)])
##     BT.res
##     ##   X Xstar       Delta nC.X nE.Xstar n.X n.Xstar
##     ## 1 0     0 -0.08571429    7       15  22      22
##     ## 2 1     1  0.12121212    6       22  28      28
##     ## 3 0     1  0.41558442    7       22  22      28
##     ## 4 1     0 -0.26666667    6       15  28      22

##     Myiid <- matrix(0, nrow = n, ncol = 4)
##     Myiid[df00$index,1] <- getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)]
##     Myiid[df11$index,2] <- getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)]
##     Myiid[df01$index,3] <- getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)]
##     Myiid[df10$index,4] <- getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)]

##     Myiid.norm <- sweep(Myiid, MARGIN = 1, FUN = "/", STATS = table(df$arm)[as.character(1-df$arm)])


##     ## *** version 1 (joint normal)
##     BT.res$w.unadjusted <- BT.res$nC.X * BT.res$nE.Xstar / prod(table(df$arm))
##     BT.res$w.standardize <- BT.res$n.X * BT.res$n.Xstar / NROW(df$arm)^2
##     expect_equal(as.double(coef(BT.std)), as.double(BT.res$Delta %*% BT.res$w.standardize), tol = 1e-5)
##     MyiidSTD.norm <- sweep(Myiid.norm, MARGIN = 2, FUN = "*", STATS = BT.res$w.standardize/BT.res$w.unadjusted)
##     expect_equal(as.double(coef(BT.std)), as.double(mean(rowSums(MyiidSTD.norm))), tol = 1e-5)
##     expect_equal(as.double(BT.std@covariance[,"netBenefit"]),
##                  as.double(crossprod(as.vector(rowSums(MyiidSTD.norm) - mean(rowSums(MyiidSTD.norm)))/table(df$arm)[as.character(df$arm)])), tol = 1e-5)
    

##     ## *** version 2 (compact)
##     myiid <- rep(0, n)
##     myiid[df00$index] <- myiid[df00$index] + getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)] * BT.res$w.standardize[1]/BT.res$w.unadjusted[1]
##     myiid[df11$index] <- myiid[df11$index] + getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)] * BT.res$w.standardize[2]/BT.res$w.unadjusted[2]
##     myiid[df01$index] <- myiid[df01$index] + getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)] * BT.res$w.standardize[3]/BT.res$w.unadjusted[3]
##     myiid[df10$index] <- myiid[df10$index] + getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)] * BT.res$w.standardize[4]/BT.res$w.unadjusted[4]

##     myiid <- myiid/table(df$arm)[as.character(1-df$arm)]
##     expect_equal(as.double(coef(BT.std)), as.double(mean(myiid)), tol = 1e-5)
##     expect_equal(as.double(BT.std@covariance[,"netBenefit"]),
##                  as.double(crossprod((myiid - mean(myiid))/table(df$arm)[as.character(df$arm)]) ), tol = 1e-5)
    
##     ## ** retrive marginal estimate after stratification
##     if(FALSE){ ## sanity check
##         ## *** marginal estimate
##         BT.raw <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = df)
##         confint(BT.raw)
##         ##    estimate        se  lower.ci  upper.ci null   p.value
##         ## Y 0.0977131 0.1908652 -0.272599 0.4428102    0 0.6109743
##         mean(getIid(BT.raw, center = FALSE, scale = FALSE))
##         ## [1] 0.0977131
##         crossprod(getIid(BT.raw))
##         ##            [,1]
##         ## [1,] 0.03642954
##         crossprod(((getIid(BT.raw, center = FALSE, scale = FALSE) * table(df$arm)[as.character(1-df$arm)]) / table(df$arm)[as.character(1-df$arm)] - coef(BT.raw))/table(df$arm)[as.character(df$arm)])
##         ## 0.03642954
##         c(BT.raw@covariance[,"netBenefit"],sqrt(BT.raw@covariance[,"netBenefit"]))
##         ## [1] 0.03642954 0.19086523

##         ## *** marginalization version 1 (joint normal)
##         BT.res$Delta %*% BT.res$w.unadjusted
##         ## [1] 0.0977131
##         range(getIid(BT.raw, center = FALSE, scale = FALSE)-rowSums(Myiid.norm))
##         ## [1] -3.330669e-16  1.110223e-16
##         mean(rowSums(Myiid.norm))
##         ## [1] 0.0977131
##         crossprod(as.vector(rowSums(Myiid.norm) - mean(rowSums(Myiid.norm)))/table(df$arm)[as.character(df$arm)])
##         ##            [,1]
##         ## [1,] 0.03642954

##         ## *** marginalization version 2 (compact)
##         myiid <- rep(0, n)
##         myiid[df00$index] <- myiid[df00$index] + getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)] 
##         myiid[df11$index] <- myiid[df11$index] + getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)] 
##         myiid[df01$index] <- myiid[df01$index] + getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)] 
##         myiid[df10$index] <- myiid[df10$index] + getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)] 
##         myiid <- myiid / table(df$arm)[as.character(1-df$arm)]

##         mean(myiid)
##         ## [1] 0.0977131
##         crossprod((myiid - mean(myiid))/table(df$arm)[as.character(df$arm)]) 
##         ##            [,1]
##         ## [1,] 0.03642954
##     }
## })


## ** standardization (second order)
BuyseTest.options(order.Hprojection = 2)

## test_that("BuyseTest - standarization (second order, binary outcome)",{
##     ## with balanced dataset standardization equals marginal analysis
##     dfB <- do.call(rbind,by(df, interaction(df$arm,df$X), function(iDF){iDF[1:min(table(df$arm,df$X)),,drop=FALSE]}))
##     rownames(dfB) <- NULL

##     BTB.raw <- BuyseTest(arm ~ cont(Y, threshold = tau), trace = FALSE, data = dfB)
##     BTB.std <- BuyseTest(arm ~ cont(Y, threshold = tau) + X, trace = FALSE, data = dfB, pool.strata = "standardization")
##     expect_equal(as.double(confint(BTB.raw)),as.double(confint(BTB.std)))

## })

BuyseTest.options(order.Hprojection = 1)

## * check Peron scoring rule in a perfectly balanced dataset
n <- 50

set.seed(10)
dt.strata <- simBuyseTest(n,
                          n.strata = 2,
                          names.strata = "stratum")

dtB.strata <- do.call(rbind,by(dt.strata, interaction(dt.strata$treatment,dt.strata$stratum), function(iDT){iDT[1:min(table(dt.strata$treatment,dt.strata$stratum))]}))
rownames(dtB.strata) <- NULL


## test_that("BuyseTest - standarization (tte outcome)",{
##     BTB.raw <- BuyseTest(treatment ~ tte(eventtime, status), trace = FALSE, data = dtB.strata,
##                          model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
##     ##             estimate        se   lower.ci  upper.ci null   p.value
##     ## eventtime 0.02881503 0.1417496 -0.2441965 0.2975942    0 0.8390032

##     BTB.std <- BuyseTest(treatment ~ tte(eventtime, status) + stratum, trace = FALSE, data = dtB.strata,
##                          pool.strata = "standardization", 
##                          model.tte = prodlim(Hist(eventtime, status) ~ treatment, data = dtB.strata))
##     expect_equal(as.double(confint(BTB.raw)), as.double(confint(BTB.std)), tol = 1e-5)
##     expect_equal(as.double(confint(BTB.raw)), c(0.028815, 0.1417496, -0.2441965, 0.2975942, 0, 0.8390032), tol = 1e-5)
## })


## * debug
## grid.prob <- expand.grid(A=0:1, X=1:2)
## grid.prob$prob <- c(0.05,0.3,0.3,NA)
## grid.prob[NROW(grid.prob),"prob"] <- 1-sum(grid.prob$prob,na.rm =TRUE)

## set.seed(1)
## n <- 500
## eps <- rnorm(n)
## index.AX <- sample(1:NROW(grid.prob), size = n, replace = TRUE, prob = grid.prob$prob)
## A   <- grid.prob[index.AX,"A"]           # treatment allocation
## X   <- grid.prob[index.AX,"X"]   # value of stratification factor
## Y   <- .5 * A + X + eps                     # outcome
  
## df <- data.frame(id = 1:length(Y), 'Y' = Y, 'arm' = A, 'X' = X)
  
## ## 1 - marginal
## BT.unstrat         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau), trace = F)
## BT.unstrat2        <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau), trace = TRUE, method.inference = "bootstrap")
## rbind(confint(BT.unstrat, transformation = FALSE),
##       confint(BT.unstrat2, method.ci = "gaussian"))
## ##       estimate         se   lower.ci  upper.ci null   p.value
## ## Y  0.006467689 0.05539344 -0.1021015 0.1150368    0 0.9070509
## ## Y1 0.006467689 0.05702373 -0.1052614 0.1180356    0 0.9099811

## ## 2- standardized
## BT.std         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = F)
## BT.std2         <- BuyseTest(data = df, arm ~ cont(Y, threshold = tau) + X, pool.strata = "standardization", trace = TRUE,
##                              method.inference = "bootstrap", strata.resampling = c("arm","X"))
## rbind(confint(BT.std, transformation = FALSE),
##       confint(BT.std2, method.ci = "gaussian"))
## ##     estimate         se   lower.ci  upper.ci null      p.value
## ## Y  0.1939914 0.05303433 0.09004605 0.2979368    0 0.0002543434
## ## Y1 0.1939914 0.05337804 0.08730577 0.2962777    0 0.0004085061

## ## 3- stratified
## df00 <- df[df$X==1,]
## df10 <- rbind(df[df$X==1 & df$arm==1,],df[df$X==2 & df$arm==0,])
## df01 <- rbind(df[df$X==1 & df$arm==0,],df[df$X==2 & df$arm==1,])
## df11 <- df[df$X==2,]

## grid.strata <- expand.grid(X = sort(unique(df$X)), Xstar = sort(unique(df$X)))
## grid.strata$X.X <- paste(grid.strata$X,grid.strata$Xstar,sep=".")
## df.extended <- do.call(rbind,apply(grid.strata, MARGIN = 1, function(iRow){
##     cbind(rbind(df[df$arm == 0 & df$X == iRow[1],],
##                 df[df$arm == 1 & df$X == iRow[2],]),
##           X.X = paste(iRow[1],iRow[2],sep="."))
## }))
## grid.strata$nC.X <- tapply(df.extended$arm==0,df.extended$X.X,sum)[grid.strata$X.X]
## grid.strata$nE.Xstar <- tapply(df.extended$arm==1,df.extended$X.X,sum)[grid.strata$X.X]
## grid.strata$n.X <- as.numeric(table(df$X)[as.character(grid.strata$X)])
## grid.strata$n.Xstar <- as.numeric(table(df$X)[as.character(grid.strata$Xstar)])
## grid.strata$w.unadjusted <- grid.strata$nC.X * grid.strata$nE.Xstar / prod(table(df$arm))
## grid.strata$w.standardize <- grid.strata$n.X * grid.strata$n.Xstar / NROW(df$arm)^2


## BT.strata  <- BuyseTest(data = df.extended, arm ~ cont(Y, threshold = tau) + X.X, trace = F)

## BT.strat00 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df00)
## BT.strat10 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df10)
## BT.strat01 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df01)
## BT.strat11 <- BuyseTest(arm ~ cont(Y, threshold = tau), method.inference = "U-statistic", trace = F, data = df11)


## BT.strata2 <- do.call(rbind,lapply(1:1000, function(iSeed){ ## iSeed <- 1
##     set.seed(iSeed)
##     iDF <- do.call(rbind,by(df,interaction(df$X,df$arm), function(iDF){iDF[sample(NROW(iDF), size = NROW(iDF), replace = TRUE),]}))
##     iDF.extended <- do.call(rbind,apply(grid.strata, MARGIN = 1, function(iRow){
##         cbind(rbind(iDF[iDF$arm == 0 & iDF$X == iRow[1],],
##                     iDF[iDF$arm == 1 & iDF$X == iRow[2],]),
##               X.X = paste(iRow[1],iRow[2],sep="."))
##     }))
##     iBT <- BuyseTest(data = iDF.extended, arm ~ cont(Y, threshold = tau) + X.X, method.inference = "none", trace = F)
##     c(global = coef(iBT), setNames(coef(iBT, strata = TRUE), names(coef(iBT, strata = TRUE))))
## }))
## rbind(confint(BT.strata),confint(BT.strata, strata = TRUE))
## ##         estimate         se   lower.ci   upper.ci null      p.value
## ## Y      0.1623210 0.02631187  0.1103524  0.2134043    0 1.359732e-09
## ## Y.1.1  0.3389333 0.07521176  0.1842154  0.4772609    0 3.282531e-05
## ## Y.1.2  0.7600000 0.04412518  0.6592417  0.8339468    0 0.000000e+00
## ## Y.2.1 -0.2311820 0.04611988 -0.3193607 -0.1390343    0 1.351153e-06
## ## Y.2.2  0.3083950 0.04304014  0.2217990  0.3901645    0 2.056422e-11
## ## rbind(confint(BT.strat00),confint(BT.strat01), confint(BT.strat10),confint(BT.strat11))

## Myiid <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
## Myiid[df00$id,"1.1"] <- getIid(BT.strat00, scale = FALSE, center = FALSE) * table(df00$arm)[as.character(1-df00$arm)]
## Myiid[df10$id,"2.1"] <- getIid(BT.strat10, scale = FALSE, center = FALSE) * table(df10$arm)[as.character(1-df10$arm)]
## Myiid[df01$id,"1.2"] <- getIid(BT.strat01, scale = FALSE, center = FALSE) * table(df01$arm)[as.character(1-df01$arm)]
## Myiid[df11$id,"2.2"] <- getIid(BT.strat11, scale = FALSE, center = FALSE) * table(df11$arm)[as.character(1-df11$arm)]

## ## Myiid.merge <- rowSums(sweep(Myiid, MARGIN = 2, FUN = "*", STATS = grid.strata$w.standardize/grid.strata$w.unadjusted))
## ## sum(Myiid.merge[df$arm==0])/prod(table(df$arm))
## ## sum(Myiid.merge[df$arm==1])/prod(table(df$arm))
## ## sum(BT.std@weightStrata*coef(BT.strata, strata = TRUE))

## MyWiid <- sweep(Myiid, FUN = "*", MARGIN = 2, STATS = BT.std@weightStrata)
## MyWiid.norm <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
## MyWiid.norm[df00$id,"1.1"] <- ((MyWiid[df00$id,"1.1"] / table(df00$arm)[as.character(1-df00$arm)]) - coef(BT.strat00)*BT.std@weightStrata[1])/table(df00$arm)[as.character(df00$arm)]
## MyWiid.norm[df10$id,"2.1"] <- ((MyWiid[df10$id,"2.1"] / table(df10$arm)[as.character(1-df10$arm)]) - coef(BT.strat10)*BT.std@weightStrata[2])/table(df10$arm)[as.character(df10$arm)]
## MyWiid.norm[df01$id,"1.2"] <- ((MyWiid[df01$id,"1.2"] / table(df01$arm)[as.character(1-df01$arm)]) - coef(BT.strat01)*BT.std@weightStrata[3])/table(df01$arm)[as.character(df01$arm)]
## MyWiid.norm[df11$id,"2.2"] <- ((MyWiid[df11$id,"2.2"] / table(df11$arm)[as.character(1-df11$arm)]) - coef(BT.strat11)*BT.std@weightStrata[4])/table(df11$arm)[as.character(df11$arm)]

## mean(rowSums(MyWiid.norm))
## sqrt(crossprod(rowSums(MyWiid.norm)))

## MyWiid2.norm <- matrix(0, nrow = NROW(df), ncol = 4, dimnames = list(NULL,unique(df.extended$X.X)))
## MyWiid2.norm[df00$id,"1.1"] <- MyWiid[df00$id,"1.1"] / prod(table(df00$arm)) - coef(BT.strat00)*BT.std@weightStrata[1]/table(df00$arm)[as.character(df00$arm)]
## MyWiid2.norm[df10$id,"2.1"] <- MyWiid[df10$id,"2.1"] / prod(table(df10$arm)) - coef(BT.strat10)*BT.std@weightStrata[2]/table(df10$arm)[as.character(df10$arm)]
## MyWiid2.norm[df01$id,"1.2"] <- MyWiid[df01$id,"1.2"] / prod(table(df01$arm)) - coef(BT.strat01)*BT.std@weightStrata[3]/table(df01$arm)[as.character(df01$arm)]
## MyWiid2.norm[df11$id,"2.2"] <- MyWiid[df11$id,"2.2"] / prod(table(df11$arm)) - coef(BT.strat11)*BT.std@weightStrata[4]/table(df11$arm)[as.character(df11$arm)]
## mean(rowSums(MyWiid2.norm))
## sqrt(crossprod(rowSums(MyWiid2.norm)))

## MyWiid3 <- rowSums(sweep(Myiid, FUN = "*", MARGIN = 2, STATS = grid.strata$w.standardize/(grid.strata$w.unadjusted*prod(table(df$arm)))))
## ## grid.strata$w.unadjusted*prod(table(df$arm))
## ## prod(table(df00$arm))
## MyWiid3[df00$id] <- MyWiid3[df00$id] - coef(BT.strat00)*BT.std@weightStrata[1]/table(df00$arm)[as.character(df00$arm)]
## MyWiid3[df10$id] <- MyWiid3[df10$id] - coef(BT.strat10)*BT.std@weightStrata[2]/table(df10$arm)[as.character(df10$arm)]
## MyWiid3[df01$id] <- MyWiid3[df01$id] - coef(BT.strat01)*BT.std@weightStrata[3]/table(df01$arm)[as.character(df01$arm)]
## MyWiid3[df11$id] <- MyWiid3[df11$id] - coef(BT.strat11)*BT.std@weightStrata[4]/table(df11$arm)[as.character(df11$arm)]
## mean(MyWiid3)
## sqrt(crossprod(MyWiid3))

## MyWiid3


##----------------------------------------------------------------------
### test-BuyseTest-standardization.R ends here
