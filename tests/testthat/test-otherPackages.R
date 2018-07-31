### test-otherPackages.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 27 2018 (17:10) 
## Version: 
## Last-Updated: maj 27 2018 (20:42) 
##           By: Brice Ozenne
##     Update #: 16
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

context("Comparison with other softwares")


## * other packages:
## ** winRatioAnalysis
## Author:
## David A. Schoenfeld

## Description:
## Fits a model to data separately for each treatment group and then calculates the win-
## Ratio as a function of follow-up time.

## Reference:
## Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a composite outcome based
## on prioritized components. Biostatistics. 2015 Sep 8;17(1):178-87.

## library(winRatioAnalysis)

## winRatio(Surv(tdeath,cdeath),
         ## treatmentVariable='rx',
         ## treatmentCodes = c(1, 0),
         ## data=dat,
         ## secondSurvivalObject=Surv(tprog0,tprog1,type='interval'),
         ## pssmIntervals=1,
         ## method = "pssm",
         ## plotPoints =3,
         ## integrationIntervals=1)

## ** WWR
## Description:
## Calculate the (weighted) win loss statistics including the win ratio, win difference and win product
## and their variances, with which the p-values are also calculated.

## Author:
## Xiaodong Luo [aut, cre] et al.

## References:
## Pocock S.J., Ariti C.A., Collier T. J. and Wang D. 2012. The win ratio: a new approach to the anal-
## ysis of composite endpoints in clinical trials based on clinical priorities.  European Heart Journal,
## 33, 176-182.
## Luo X., Tian H., Mohanty S. and Tsai W.-Y. 2015.  An alternative approach to confidence interval
## estimation for the win ratio statistic. Biometrics, 71, 139-145.
## Bebu I. and Lachin J.M. 2016.   Large sample inference for a win ratio analysis of a composite
## outcome based on prioritized components. Biostatistics, 17, 178-187.
## Luo X., Qiu J., Bai S. and Tian H. 2017.  Weighted win loss approach for analyzing prioritized
## outcomes. Statistics in Medicine, <doi: 10.1002/sim.7284>

## library(WWR)

## set.seed(10)
## dt.sim <- rbind(data.table(time1 = runif(10,0,5),
                           ## status1 = 1,
                           ## time2 = 20,
                           ## status2 = 1,
                           ## treatment = 1),
                ## data.table(time1 = runif(10,0,5),
                           ## status1 = 1,
                           ## time2 = 0,
                           ## status2 = 1,
                           ## treatment = 0)
                ## )

## wtest <- WWR::winratio(y1 = dt.sim$time1,
                       ## y2 = dt.sim$time2,
                       ## d1 = dt.sim$status1,
                       ## d2 = dt.sim$status2,
                       ## z = dt.sim$treatment)

## summary(wtest)
## str(wtest)

## BT <- BuyseTest(treatment ~ tte(time1, censoring = status1) + tte(time2, censoring = status2),
                ## data = dt.sim,
                ## method.inference = "none")

## summary(BT, statistic = "winRatio", percentage = FALSE)
## 1



##  > results
##    endpoint threshold total favorable unfavorable neutral uninf    delta    Delta
## eventtime1     1e-12   400       287         113       0     0 2.539823 2.539823
## eventtime2     1e-12     0         0           0       0     0       NA 2.539823

## ** WLreg
## Description:
## Use various regression models for the analysis of win loss endpoints adjusting for non-binary and multivariate covariates.

## Author:
## Xiaodong Luo

## References:
## Pocock S.J., Ariti C.A., Collier T. J. and Wang D. 2012. The win ratio: a new approach to the anal-
## ysis of composite endpoints in clinical trials based on clinical priorities.  European Heart Journal,
## 33, 176-182.
## Luo X., Tian H., Mohanty S. and Tsai W.-Y. 2015.  An alternative approach to confidence interval
## estimation for the win ratio statistic. Biometrics, 71, 139-145.
## Luo X., Qiu J., Bai S. and Tian H. 2017.  Weighted win loss approach for analyzing prioritized
## outcomes. Statistics in Medicine, to appear.

## library(WLreg)

## set.seed(10)
## dt.sim <- rbind(data.table(time1 = runif(10,0,5),
##                            status1 = 1,
##                            time2 = 20,
##                            status2 = 1,
##                            treatment = 1),
##                 data.table(time1 = runif(10,0,5),
##                            status1 = 1,
##                            time2 = 0,
##                            status2 = 1,
##                            treatment = 0)
##                 )

## aa<-winreg(dt.sim$time1,
##            dt.sim$time2,
##            dt.sim$status1,
##            dt.sim$status2,
##            dt.sim$treatment)
## aa
##----------------------------------------------------------------------
### test-otherPackages.R ends here
