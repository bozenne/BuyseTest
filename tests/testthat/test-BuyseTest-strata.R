### test-BuyseTest-strata.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan  5 2023 (11:45) 
## Version: 
## Last-Updated: Apr 21 2025 (12:23) 
##           By: Brice Ozenne
##     Update #: 65
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(BuyseTest)
    library(data.table)

    library(testthat)
}

context("Check stratification in BuyseTest")

## * setting
BuyseTest.options(pool.strata = "Buyse")

## * simulate data
n <- 100

set.seed(10)
dt.strata <- simBuyseTest(n,
                          n.strata = 3,
                          argsBin = list(name = "Y_1"),
                          argsCont = list(name = "Y_2"),
                          argsTTE = NULL,
                          names.strata = "stratum")
dt.strata$Y_1 <- as.numeric(dt.strata$Y_1)-1
dt.strata$stratum <- as.numeric(dt.strata$stratum)
dt.strata <- dt.strata[order(dt.strata$stratum),]

all.stat <- c("count.favorable","count.unfavorable","count.neutral","count.uninf","netBenefit","winRatio","favorable","unfavorable","neutral","uninf")

## * No strata
test_that("no strata (equivalence glm)",{
    ## GS.glm <- WINS::win.stat(data = dt.strata, summary.print = FALSE,
    ##                      ep_type = c("bin"),
    ##                      arm.name = c("T","C"), tau = 0, priority = 1,
    ##                      alpha = 0.05, digit = 3, censoring_adjust = "No",
    ##                      weight = "unstratified", pvalue = "two-sided")
    ## c(GS.glm$Win_statistic$Win_Ratio[1],GS.glm$Win_statistic$Net_Benefit[1])

    e.glmLOGIT <- glm(Y_1 ~ treatment, data = dt.strata, family = binomial(link="logit"))
    e.glmID <- glm(Y_1 ~ treatment, data = dt.strata, family = binomial(link="identity"))
    GS.glm <- c(WR = as.double(exp(coef(e.glmLOGIT)["treatmentT"])),
                NB = as.double(coef(e.glmID)["treatmentT"])
                )

    ## ** check equivalence glm
    e.BT <- BuyseTest(treatment ~ bin(Y_1), data = dt.strata, trace = FALSE,
                      method.inference = "u-statistic")
    ## summary(e.BT)
    test <- c(WR = as.double(coef(e.BT, statistic = "winRatio")),
              NB = as.double(coef(e.BT, statistic = "netBenefit"))
              )

    expect_equal(test, GS.glm, tol = 1e-5)

})

test_that("no strata (check coef)",{

    ## ** one endpoint
    e.BT <- BuyseTest(treatment ~ bin(Y_1), data = dt.strata, trace = FALSE,
                      method.inference = "bootstrap", n.resampling = 10)

    ## GS.coef <- sapply(all.stat, function(iStat){unname(coef(e.BT, statistic = iStat))})
    GS.coef <- c("count.favorable" = 2448, "count.unfavorable" = 2548, "count.neutral" = 5004, "count.uninf" = 0,
                 "netBenefit" = -0.01, "winRatio" = 0.96075353, "favorable" = 0.2448, "unfavorable" = 0.2548, "neutral" = 0.5004, "uninf" = 0)

    for(iStat in all.stat){ ## iStat <- all.stat[1]
        expect_equal(as.double(coef(e.BT, statistic = iStat)),
                     as.double(GS.coef[iStat]), tol = 1e-5)        
        if(iStat %in% c("count.neutral","neutral","count.uninf","uninf") == FALSE){
            expect_equal(coef(e.BT, statistic = iStat),
                         coef(e.BT, statistic = iStat, strata = FALSE, cumulative = TRUE), tol = 1e-5)
            expect_equal(as.double(coef(e.BT, statistic = iStat)),
                         as.double(coef(e.BT, statistic = iStat, strata = TRUE, cumulative = TRUE)), tol = 1e-5)
        }
        expect_equal(coef(e.BT, statistic = iStat),
                     coef(e.BT, statistic = iStat, strata = FALSE, cumulative = FALSE), tol = 1e-5)
        expect_equal(as.double(coef(e.BT, statistic = iStat)),
                     as.double(coef(e.BT, statistic = iStat, strata = TRUE, cumulative = FALSE)), tol = 1e-5)

        if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf") == FALSE){
            if(iStat %in% c("neutral","uninf") == FALSE){
                expect_equal(coef(e.BT, statistic = iStat, resampling = TRUE),
                             coef(e.BT, statistic = iStat, strata = FALSE, cumulative = TRUE, resampling = TRUE), tol = 1e-5)
                expect_equal(as.double(coef(e.BT, statistic = iStat, resampling = TRUE)),
                             as.double(coef(e.BT, statistic = iStat, strata = TRUE, cumulative = TRUE, resampling = TRUE)), tol = 1e-5)
            }
            expect_equal(coef(e.BT, statistic = iStat, resampling = TRUE),
                         coef(e.BT, statistic = iStat, strata = FALSE, cumulative = FALSE, resampling = TRUE), tol = 1e-5)
            expect_equal(as.double(coef(e.BT, statistic = iStat, resampling = TRUE)),
                         as.double(coef(e.BT, statistic = iStat, strata = TRUE, cumulative = FALSE, resampling = TRUE)), tol = 1e-5)
        }
    }

    ## ** two endpoints
    e.BT2 <- BuyseTest(treatment ~ bin(Y_1) + cont(Y_2), data = dt.strata, trace = FALSE,
                      method.inference = "bootstrap", n.resampling = 10)
    
    test <- c(WR = as.double(coef(e.BT2, statistic = "winRatio", endpoint = "Y_1")),
              NB = as.double(coef(e.BT2, statistic = "netBenefit", endpoint = "Y_1"))
              )

    expect_equal(test, c("WR" = 0.96075353, "NB" = -0.01), tol = 1e-5)

    for(iStat in all.stat){ ## iStat <- "count.favorable"
        expect_equal(as.double(coef(e.BT2, statistic = iStat, endpoint = "Y_1")),
                     as.double(GS.coef[iStat]), tol = 1e-5)
        if(iStat %in% c("count.neutral","neutral","count.uninf","uninf") == FALSE){
            expect_equal(coef(e.BT2, statistic = iStat),
                         coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = TRUE), tol = 1e-5)
            expect_equal(as.double(coef(e.BT2, statistic = iStat)),
                         as.double(coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = TRUE)), tol = 1e-5)
        }

        test1 <- coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = FALSE)
        test2 <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = FALSE)
        if(iStat %in% c("count.neutral","neutral","count.uninf","uninf")){
            expect_equal(coef(e.BT2, statistic = iStat), test1, tol = 1e-5)
            expect_equal(unname(coef(e.BT2, statistic = iStat)), as.double(test2), tol = 1e-5)
        }else if(iStat != "winRatio"){
            expect_equal(coef(e.BT2, statistic = iStat), cumsum(test1), tol = 1e-5)
            expect_equal(coef(e.BT2, statistic = iStat), cumsum(test2), tol = 1e-5)
        }

        if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf") == FALSE){
            if(iStat %in% c("neutral","uninf") == FALSE){
                expect_equal(coef(e.BT2, statistic = iStat, resampling = TRUE),
                             coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = TRUE, resampling = TRUE), tol = 1e-5)
                expect_equal(as.double(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                             as.double(coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = TRUE, resampling = TRUE)), tol = 1e-5)
            }

            test1 <- coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = FALSE, resampling = TRUE)
            test2 <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = FALSE, resampling = TRUE)
            if(iStat %in% c("neutral","uninf")){
                expect_equal(coef(e.BT2, statistic = iStat, resampling = TRUE), test1, tol = 1e-5)
                expect_equal(as.double(coef(e.BT2, statistic = iStat, resampling = TRUE)), as.double(test2), tol = 1e-5)
            }else if(iStat != "winRatio"){        
                expect_equal(unname(coef(e.BT2, statistic = iStat, resampling = TRUE)), BuyseTest:::.rowCumSum_cpp(test1), tol = 1e-5)
                expect_equal(unname(coef(e.BT2, statistic = iStat, resampling = TRUE)), unname(cbind(test2[,"Y_1"], test2[,"Y_1"]+test2[,"Y_2"])), tol = 1e-5)
            }
        }
    }
})

## * Strata (historical weight)
test_that("strata (historical weight)",{

    ## GS <- WINS::win.stat(data = dt.strata, summary.print = FALSE,
    ##                      ep_type = c("bin"),
    ##                      arm.name = c("T","C"), tau = 0, priority = 1,
    ##                      alpha = 0.05, digit = 3, censoring_adjust = "No",
    ##                      weight = "equal", pvalue = "two-sided")
    ## c(GS$Win_statistic$Win_Ratio[1],GS$Win_statistic$Net_Benefit[1])
    GS <- c("WR" = 0.99196326, "NB" = -0.00209895)

    ## ** check equivalence WINS
    e.BT <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, trace = FALSE,
                      method.inference = "u-statistic")
    ls.eBT <- by(dt.strata, INDICES = dt.strata$stratum, FUN = function(iData){
        BuyseTest(treatment ~ bin(Y_1), data = iData, method.inference = "u-statistic", trace = FALSE)
    })

    ## *** point estimate
    test <- c(WR = as.double(coef(e.BT, statistic = "winRatio")),
              NB = as.double(coef(e.BT, statistic = "netBenefit"))
              )
    expect_equal(test, GS, tol = 1e-5)

    count.favorable <- sapply(ls.eBT, function(iBT){iBT@count.favorable})
    count.unfavorable <- sapply(ls.eBT, function(iBT){iBT@count.unfavorable})
    n.pairs <- sapply(ls.eBT, function(iBT){iBT@n.pairs})
    weight <- n.pairs/sum(n.pairs)

    expect_equal(unname(test["NB"]), weighted.mean(sapply(ls.eBT, coef, statistic = "netBenefit"), weight))
    expect_equal(unname(test["WR"]), sum(count.favorable)/sum(count.unfavorable))

    ## *** iid
    iid.via.ls <- do.call(cbind,lapply(c("favorable","unfavorable","netBenefit"), function(iStat){do.call(c,lapply(ls.eBT, getIid, statistic = iStat))}))
    colnames(iid.via.ls) <- c("favorable","unfavorable","netBenefit")
    iid.direct <- cbind(netBenefit = rowSums(getIid(e.BT, strata = TRUE, statistic = "netBenefit")),
                        winRatio = rowSums(getIid(e.BT, strata = TRUE, statistic = "winRatio")))
    Delta.direct <- cbind(favorable = coef(e.BT, strata = TRUE, statistic = "favorable"),
                          unfavorable = coef(e.BT, strata = TRUE, statistic = "unfavorable"))[dt.strata$stratum,]
    
    expect_equal(as.double(iid.direct[,"netBenefit"]),
                 as.double(iid.via.ls[,"netBenefit"]),
                 tol = 1e-7)

    expect_equal(as.double(iid.direct[,"winRatio"]),
                 as.double(iid.via.ls[,"favorable"]/Delta.direct[,2] - iid.via.ls[,"unfavorable"]*Delta.direct[,1]/Delta.direct[,2]^2),
                 tol = 1e-7)

    expect_equal(c(sum(getIid(e.BT, statistic = "netBenefit")^2), sum(getIid(e.BT, statistic = "winRatio")^2)),
                 unname(e.BT@covariance[,c("netBenefit","winRatio")]),
                 tol = 1e-7)

    ## *** stratified iid
    iid.directNB <- getIid(e.BT, strata = TRUE, statistic = "netBenefit")
    iid.directWR <- getIid(e.BT, strata = TRUE, statistic = "winRatio")

    for(iStrata in 1:3){ ## iStrata <- 1
        expect_equal(unname(getIid(ls.eBT[[iStrata]], statistic = "netBenefit")),
                     unname(iid.directNB[dt.strata$stratum==iStrata,iStrata]),
                     tol = 1e-7)
        expect_equal(unname(getIid(ls.eBT[[iStrata]], statistic = "winRatio")),
                     unname(iid.directWR[dt.strata$stratum==iStrata,iStrata]),
                     tol = 1e-7)
    }
  
    ## *** confint
    test <- confint(e.BT, strata = TRUE)
    GS <- do.call(rbind,lapply(ls.eBT,confint))

    expect_equal(as.double(unlist(test)), as.double(unlist(GS)), tol = 1e-6)
})

test_that("strata (check coef)",{

  ## ** one endpoint
  e.BT <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, trace = FALSE,
                    method.inference = "bootstrap", n.resampling = 10)

  ## GS.coef <- sapply(all.stat, function(iStat){unname(coef(e.BT, statistic = iStat))})
  GS.coef <- c("count.favorable" = 864, "count.unfavorable" = 871, "count.neutral" = 1600, "count.uninf" = 0,
               "netBenefit" = -0.00209895, "winRatio" = 0.99196326, "favorable" = 0.25907046, "unfavorable" = 0.26116942, "neutral" = 0.47976012, "uninf" = 0)

  for(iStat in all.stat){ ## iStat <- "favorable"
      expect_equal(as.double(coef(e.BT, statistic = iStat)),
                   as.double(GS.coef[iStat]), tol = 1e-5)        
      if(iStat %in% c("count.neutral","neutral","count.uninf","uninf") == FALSE){
          expect_equal(coef(e.BT, statistic = iStat),
                       coef(e.BT, statistic = iStat, strata = FALSE, cumulative = TRUE), tol = 1e-5)
          test <- coef(e.BT, statistic = iStat, strata = TRUE, cumulative = TRUE)
          if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf")){
              expect_equal(as.double(coef(e.BT, statistic = iStat)), as.double(sum(test)), tol = 1e-5)
          }else if(iStat != "winRatio"){
              expect_equal(as.double(coef(e.BT, statistic = iStat)), as.double(sum(test * e.BT@weightStrata)), tol = 1e-5)
          }
      }
      expect_equal(coef(e.BT, statistic = iStat),
                   coef(e.BT, statistic = iStat, strata = FALSE, cumulative = FALSE), tol = 1e-5)
      test <- coef(e.BT, statistic = iStat, strata = TRUE, cumulative = FALSE)
      if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf")){
          expect_equal(as.double(coef(e.BT, statistic = iStat)), as.double(sum(test)), tol = 1e-5)
      }else if(iStat != "winRatio"){
          expect_equal(as.double(coef(e.BT, statistic = iStat)), as.double(sum(test * e.BT@weightStrata)), tol = 1e-5)
      }

      if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf") == FALSE){
          if(iStat %in% c("neutral","uninf") == FALSE){
              expect_equal(coef(e.BT, statistic = iStat, resampling = TRUE),
                           coef(e.BT, statistic = iStat, strata = FALSE, cumulative = TRUE, resampling = TRUE), tol = 1e-5)
              test <- coef(e.BT, statistic = iStat, strata = TRUE, cumulative = TRUE, resampling = TRUE)
              if(iStat != "winRatio"){
                  expect_equal(as.double(coef(e.BT, statistic = iStat, resampling = TRUE)),
                               as.double(rowSums(test * e.BT@weightStrataResampling)), tol = 1e-5)
              }
          }
          expect_equal(coef(e.BT, statistic = iStat, resampling = TRUE),
                       coef(e.BT, statistic = iStat, strata = FALSE, cumulative = FALSE, resampling = TRUE), tol = 1e-5)
          test <- coef(e.BT, statistic = iStat, strata = TRUE, cumulative = FALSE, resampling = TRUE)
          if(iStat != "winRatio"){
              expect_equal(as.double(coef(e.BT, statistic = iStat, resampling = TRUE)),
                           as.double(rowSums(test * e.BT@weightStrataResampling)), tol = 1e-5)
          }
      }
  }

  ## ** two endpoints
  e.BT2 <- BuyseTest(treatment ~ bin(Y_1) + cont(Y_2) + stratum, data = dt.strata, trace = FALSE,
                     method.inference = "bootstrap", n.resampling = 10)

  for(iStat in all.stat){ ## iStat <- "favorable"
      if(iStat %in% c("count.neutral","neutral","count.uninf","uninf") == FALSE){
          expect_equal(coef(e.BT2, statistic = iStat),
                       coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = TRUE), tol = 1e-5)
          test <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = TRUE)
          if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf")){
              expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(colSums(test)), tol = 1e-5)
          }else if(iStat != "winRatio"){
              expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(colSums(BuyseTest:::.colMultiply_cpp(test, e.BT2@weightStrata))), tol = 1e-5)
          }
      }

      test <- coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = FALSE)
      if(iStat %in% c("count.neutral","neutral","count.uninf","uninf")){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(test), tol = 1e-5)
      }else if(iStat != "winRatio"){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(cumsum(test)), tol = 1e-5)
      }
      test <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = FALSE)
      if(iStat %in% c("count.neutral","count.uninf")){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(colSums(test)), tol = 1e-5)
      }else if(iStat %in% c("count.favorable","count.unfavorable")){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(cumsum(colSums(test))), tol = 1e-5)
      }else if(iStat %in% c("neutral","uninf")){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(colSums(BuyseTest:::.colMultiply_cpp(test, e.BT2@weightStrata))), tol = 1e-5)
      }else if(iStat != "winRatio"){
          expect_equal(as.double(coef(e.BT2, statistic = iStat)), as.double(cumsum(colSums(BuyseTest:::.colMultiply_cpp(test, e.BT2@weightStrata)))), tol = 1e-5)
      }

      if(iStat %in% c("count.favorable","count.unfavorable","count.neutral","count.uninf") == FALSE){
          if(iStat %in% c("neutral","uninf") == FALSE){
              expect_equal(coef(e.BT2, statistic = iStat, resampling = TRUE),
                           coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = TRUE, resampling = TRUE), tol = 1e-5)
              test <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = TRUE, resampling = TRUE)
              if(iStat != "winRatio"){
                  expect_equal(as.double(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                               c(as.double(rowSums(test[,,1] * e.BT2@weightStrataResampling)), as.double(rowSums(test[,,2] * e.BT2@weightStrataResampling))), tol = 1e-5)
              }
          }
          test <- coef(e.BT2, statistic = iStat, strata = FALSE, cumulative = FALSE, resampling = TRUE)
          if(iStat %in% c("neutral","uninf")){
              expect_equal(unname(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                           unname(test), tol = 1e-5)
          }else if(iStat != "winRatio"){
              expect_equal(unname(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                           unname(BuyseTest:::.rowCumSum_cpp(test)), tol = 1e-5)
          }
          test <- coef(e.BT2, statistic = iStat, strata = TRUE, cumulative = FALSE, resampling = TRUE)
          if(iStat %in% c("neutral","uninf")){
              expect_equal(as.double(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                           c(as.double(rowSums(test[,,1] * e.BT2@weightStrataResampling)), as.double(rowSums((test[,,2]) * e.BT2@weightStrataResampling))), tol = 1e-5)
          }else if(iStat != "winRatio"){
              expect_equal(as.double(coef(e.BT2, statistic = iStat, resampling = TRUE)),
                           c(as.double(rowSums(test[,,1] * e.BT2@weightStrataResampling)), as.double(rowSums((test[,,1]+test[,,2]) * e.BT2@weightStrataResampling))), tol = 1e-5)
          }
      }
  }

})

## * Strata (CMH weight)
test_that("strata (CMH weights)",{
    ## GS.WINS <- WINS::win.stat(data = dt.strata, summary.print = FALSE,
    ##                           ep_type = c("bin"),
    ##                           arm.name = c("T","C"), tau = 0, priority = 1,
    ##                           alpha = 0.05, digit = 3, censoring_adjust = "No",
    ##                           weight = "MH-type", pvalue = "two-sided")
    ## c(GS.WINS$Win_statistic$Net_Benefit[1], GS.WINS$Win_statistic$Win_Ratio[1])
    GS.WINS <- c("NB" = -0.00916137, "WR" = 0.96538417)
    GS.stats <- mantelhaen.test(xtabs(data = dt.strata, ~ treatment + Y_1 + stratum))

    ## ** check equivalence WINS
    e.BT <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "CMH", trace = FALSE)
    ls.eBT <- by(dt.strata, INDICES = dt.strata$stratum, FUN = function(iData){
        BuyseTest(treatment ~ bin(Y_1), data = iData, trace = FALSE)
    })
    ecount.favorable <- coef(e.BT, statistic = "count.favorable", strata = TRUE, simplify = FALSE)[,1]
    ecount.unfavorable <- coef(e.BT, statistic = "count.unfavorable", strata = TRUE, simplify = FALSE)[,1]
    ecount.strata <- table(dt.strata$stratum)
    ecount.pairs <- e.BT@n.pairs
    weight <- ecount.pairs/ecount.strata

    ## *** point estimate
    test <- list(FA = weighted.mean( ecount.favorable/ecount.pairs, weight),
                 UN = weighted.mean( ecount.unfavorable/ecount.pairs, weight),
                 NB = weighted.mean( ecount.favorable/ecount.pairs, weight) - weighted.mean( ecount.unfavorable/ecount.pairs, weight),
                 WR = weighted.mean( ecount.favorable/ecount.pairs, weight) / weighted.mean( ecount.unfavorable/ecount.pairs, weight)
                 )
    expect_equal(unname(GS.WINS["NB"]), test$NB, tol = 1e-5)
    expect_equal(unname(GS.WINS["WR"]), test$WR, tol = 1e-5)
    expect_equal(unname(GS.stats$estimate), test$WR, tol = 1e-5)
    expect_equal(as.double(e.BT@Delta[,c("favorable","unfavorable","netBenefit","winRatio")]), as.double(test), tol = 1e-5)

})

## * Strata (variance weights)
test_that("strata (variance weights)",{

    e.BTopt1 <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "var-favorable",
                          method.inference = "none", trace = FALSE)
    e.BTopt2 <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "var-unfavorable",
                          method.inference = "none", trace = FALSE)
    e.BTopt3 <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "var-netBenefit",
                          method.inference = "none", trace = FALSE)
    e.BTopt4 <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "var-winRatio",
                          method.inference = "none", trace = FALSE)

    test.optimal <- c(FA = as.double(coef(e.BTopt1, statistic = "favorable")),
                      UN = as.double(coef(e.BTopt2, statistic = "unfavorable")),
                      NB = as.double(coef(e.BTopt3, statistic = "netBenefit")),
                      WR = as.double(coef(e.BTopt4, statistic = "winRatio")))
  
    ls.eBT <- by(dt.strata, INDICES = dt.strata$stratum, FUN = function(iData){
        BuyseTest(treatment ~ bin(Y_1), data = iData, trace = FALSE, method.inference = "u-statistic")
    })
    ecount.favorable <- sapply(ls.eBT, function(iBT){coef(iBT, statistic = "count.favorable", strata = TRUE, simplify = FALSE)[,1]})
    ecount.unfavorable <- sapply(ls.eBT, function(iBT){coef(iBT, statistic = "count.unfavorable", strata = TRUE, simplify = FALSE)[,1]})
    ecount.strata <- table(dt.strata$stratum)
    ecount.pairs <- sapply(ls.eBT, function(iBT){iBT@n.pairs})

    GS.optimal <- c(FA = weighted.mean( ecount.favorable/ecount.pairs,
                                       w = do.call(rbind,lapply(ls.eBT, confint, statistic = "favorable", order.Hprojection = 1))[,"se"]^(-2)),
                    UN = weighted.mean( ecount.unfavorable/ecount.pairs,
                                       w = do.call(rbind,lapply(ls.eBT, confint, statistic = "unfavorable", order.Hprojection = 1))[,"se"]^(-2)),
                    NB = weighted.mean( (ecount.favorable-ecount.unfavorable)/ecount.pairs,
                                       w = do.call(rbind,lapply(ls.eBT, confint, statistic = "netBenefit", order.Hprojection = 1))[,"se"]^(-2)),
                    WR = weighted.mean( ecount.favorable/ecount.unfavorable,
                                       w = do.call(rbind,lapply(ls.eBT, confint, statistic = "winRatio", order.Hprojection = 1))[,"se"]^(-2))
                    )

    expect_equal(test.optimal, GS.optimal, tol = 1e-5)

    expect_equal(coef(e.BTopt4, cumulative = FALSE, strata = FALSE, statistic = "winRatio"),
                 coef(e.BTopt4, cumulative = TRUE, strata = FALSE, statistic = "winRatio"),
                 tol = 1e-6)

    ##
    ## e.BTopt1 <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, pool.strata = "var-favorable",
    ##                       method.inference = "bootstrap", trace = FALSE)
    ## summary(e.BTopt1)
})
##----------------------------------------------------------------------
### test-BuyseTest-strata.R ends here
