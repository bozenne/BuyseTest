### test-BuyseTest-CR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne, Eva Cantagallo
## Created: jul 12 2018 (16:58) 
## Version: 
## Last-Updated: feb 27 2019 (22:33) 
##           By: Brice Ozenne
##     Update #: 11
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

context("Check that BuyseTest with competing risks \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0, keep.survival = T)


alpha.X <- 2
alpha.Y <- 3
alpha.cens <- 1.5
n <- 1e2

## * Simulate data with censoring
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alpha.X), time.cens = rexp(n, rate = alpha.cens), treatment = "1"),
             data.frame(time1 = rexp(n, rate = alpha.Y), time.cens = rexp(n, rate = alpha.cens), treatment = "2"))
df$time <- pmin(df$time1, df$time.cens) ## first status
df$status <- as.numeric((df$time1<df$time.cens)) ## type of event 
df$strata <- sample(c('a', 'b', 'c'), 2*n, replace = T)
df$toxicity <- sample(0:1, 2*n, replace = T)
df$strata2 <- sample(c('d', 'e', 'f'), 2*n, replace = T)

## * test
test_that("Package give the same results when model.tte is (not) provided as an argument", {
  
  ## Create prodlim object to be inserted as an argument
  fit = prodlim(Hist(time, status) ~ treatment +  strata + strata2, data = df)
  
  ## Net benefit without passing model.tte
  B = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata + strata2, data = df,
                method.inference = "none") 
  
  ## Net benefit with model.tte passed as an argument
  B.model = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5) + strata + strata2, data = df, model.tte = fit,
                      method.inference = "none") 
  
  
  # # treatment strata strata2 ..strata..
  # # 1          0      a       d          1
  # # 2          1      a       d          1
  # # 3          0      b       d          2
  # # 4          1      b       d          2
  # # 5          0      c       d          3
  # # 6          1      c       d          3
  # # 7          0      a       e          4
  # # 8          1      a       e          4
  # # 9          0      b       e          5
  # # 10         1      b       e          5
  # # 11         0      c       e          6
  # # 12         1      c       e          6
  # # 13         0      a       f          7
  # # 14         1      a       f          7
  # # 15         0      b       f          8
  # # 16         1      b       f          8
  # # 17         0      c       f          9
  # # 18         1      c       f          9
  # 
  # summary(B.model)
  # summary(B)
  # getSurvival(B.model)
  # 
  # identical(getSurvival(B),getSurvival(B.model))
  # 
  # getSurvival(B)[[1]][[1]][[1]] - getSurvival(B.model)[[1]][[1]][[1]]
  # getSurvival(B)[[2]][[1]][[1]] - getSurvival(B.model)[[2]][[1]][[1]]
  # getSurvival(B)[[2]][[1]][[1]] - getSurvival(B.model)[[2]][[1]][[1]]
  # 
  # getSurvival(B)[["lastSurv"]]
  # getSurvival(B.model)[["lastSurv"]]
  # names(getSurvival(B))
  
  ## Tests
  expect_equal(as.double(B@delta.netBenefit), as.double(B.model@delta.netBenefit))
  expect_equal(as.double(B@Delta.netBenefit), as.double(B.model@Delta.netBenefit))
  
})

test_that("When TTE endpoints are analyzed several times with different thresholds, the results are the same than those when analyzed once with lowest threshold", {
  
  ## Endpoint analyzed once with threshold = 0.5
  B1 = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 0.5), data = df) 
  
  ## Endpoint analyzed twice with threshold = c(1, 0.75, 0.5)
  B2 = BuyseTest(treatment ~ tte(time, censoring = status, threshold = 1) + tte(time, censoring = status, threshold = 0.75) + tte(time, censoring = status, threshold = 0.5), data = df) 
  
  ## Tests
  expect_equal(as.double(B1@Delta.netBenefit), as.double(B2@Delta.netBenefit[3]))
  
})
######################################################################
### test-BuyseTest-CR.R ends here
