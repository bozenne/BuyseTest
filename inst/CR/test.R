### test.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  9 2018 (16:52) 
## Version: 
## Last-Updated: jul 10 2018 (18:23) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(riskRegression)

## * survival case
n <- 1e4
set.seed(10)
d <- sampleData(n, outcome = "survival")
d[,event1 := 1]

## ** no censoring
e.BT <- BuyseTest(X1 ~ tte(eventtime, censoring = event1), data = d, method.inference = "none")
e.BT

grid <- expand.grid(d[X1==0,eventtime],d[X1==1,eventtime])
mean(grid[,2]>grid[,1])-mean(grid[,2]<grid[,1])

d[,mean(time<=1), by = "X1"]

## ** censoring
e.BT <- BuyseTest(X1 ~ tte(time, censoring = event), data = d, method.inference = "none",
                  method.tte = "Gehan", keep.comparison = TRUE)
e.BT
e.BT@tableComparison[[1]][,mean(favorable)-mean(unfavorable)]

e.BTC <- BuyseTest(X1 ~ tte(time, censoring = event), data = d, method.inference = "none",
                   method.tte = "Gehan corrected", keep.comparison = TRUE)
e.BTC
e.BTC@tableComparison[[1]][,(mean(favorable)-mean(unfavorable))/(1-mean(uninformative))] ## inverse probability weighting


######################################################################
### test.R ends here
