### test-BuyseTest-strata.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan  5 2023 (11:45) 
## Version: 
## Last-Updated: Jan  5 2023 (12:40) 
##           By: Brice Ozenne
##     Update #: 9
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

    library(testthat)
}

context("Check BuyseTest without strata")


## * simulate data
n <- 100

set.seed(10)
dt.strata <- cbind(id = 1:n,
                   simBuyseTest(n,
                                n.strata = 3,
                                argsBin = list(name = "Y_1"),
                                argsCont = NULL,
                                argsTTE = NULL,
                                names.strata = "stratum"))
dt.strata$Y_1 <- as.numeric(dt.strata$Y_1)-1
dt.strata$stratum <- as.numeric(dt.strata$stratum)
dt.strata <- dt.strata[order(dt.strata$stratum),]

## * No strata
test_that("no strata (equivalence glm)",{
    ## GS <- WINS::win.stat(data = dt.strata, summary.print = FALSE,
    ##                      ep_type = c("bin"),
    ##                      arm.name = c("T","C"), tau = 0, priority = 1,
    ##                      alpha = 0.05, digit = 3, censoring_adjust = "No",
    ##                      weight = "unstratified", pvalue = "two-sided")
    ## c(GS$Win_statistic$Win_Ratio[1],GS$Win_statistic$Net_Benefit[1])

    e.glmLOGIT <- glm(Y_1 ~ treatment, data = dt.strata, family = binomial(link="logit"))
    e.glmID <- glm(Y_1 ~ treatment, data = dt.strata, family = binomial(link="identity"))

    GS <- c(WR = as.double(exp(coef(e.glmLOGIT)["treatmentT"])),
            NB = as.double(coef(e.glmID)["treatmentT"])
            )

    e.BT <- BuyseTest(treatment ~ bin(Y_1), data = dt.strata, trace = FALSE)
    test <- c(WR = as.double(coef(e.BT, statistic = "winRatio")),
              NB = as.double(coef(e.BT, statistic = "netBenefit"))
              )

    expect_equal(test, GS, tol = 1e-5)
})

## * Strata (historical weight)
test_that("strata (historical weight)",{

    ## GS <- WINS::win.stat(data = dt.strata, summary.print = FALSE,
    ##                      ep_type = c("bin"),
    ##                      arm.name = c("T","C"), tau = 0, priority = 1,
    ##                      alpha = 0.05, digit = 3, censoring_adjust = "No",
    ##                      weight = "equal", pvalue = "two-sided")
    ## c(GS$Win_statistic$Win_Ratio[1],GS$Win_statistic$Net_Benefit[1])
    GS <- c("WR" = 1.33187773, "NB" = 0.06859206)

    e.BT <- BuyseTest(treatment ~ bin(Y_1) + stratum, data = dt.strata, trace = FALSE)
    ls.eBT <- by(dt.strata, INDICES = dt.strata$stratum, FUN = function(iData){
        BuyseTest(treatment ~ bin(Y_1), data = iData, trace = FALSE)
    })

    ## point estimate
    test <- c(WR = as.double(coef(e.BT, statistic = "winRatio")),
              NB = as.double(coef(e.BT, statistic = "netBenefit"))
              )
    expect_equal(test, GS, tol = 1e-5)

    count.favorable <- sapply(ls.eBT, function(iBT){iBT@count.favorable})
    count.unfavorable <- sapply(ls.eBT, function(iBT){iBT@count.unfavorable})
    n.pairs <- sapply(ls.eBT, function(iBT){iBT@n.pairs})
    weight <- n.pairs/sum(n.pairs)

    expect_equal(unname(test["NB"]),
                 weighted.mean(sapply(ls.eBT, coef, statistic = "netBenefit"), weight))
    expect_equal(unname(test["WR"]),
                 sum(count.favorable)/sum(count.unfavorable))

    ## iid
    test.iid <- do.call(rbind,lapply(1:3, function(iS){(ls.eBT[[iS]]@iidAverage$favorable - ls.eBT[[iS]]@iidAverage$unfavorable)*weight[[iS]]}))

    GS.iid <- e.BT@iidAverage$favorable - e.BT@iidAverage$unfavorable
    expect_equal(test.iid, GS.iid, tol = 1e-5)


    e.BT@covariance

    GS.iid <- e.BT@iidAverage$favorable/sum(count.unfavorable) - e.BT@iidAverage$unfavorable*(sum(count.favorable)/sum(count.unfavorable)^2) 
    crossprod(GS.iid)

    0.00147339/sum(count.unfavorable)^2 + 0.001064895*sum(count.favorable)^2/sum(count.unfavorable)^4 + 2*0.001135432*sum(count.favorable)/sum(count.unfavorable)^3
    0.00147339 0.001064895 -0.001135432 0.004809149 0.1495202
})


## * Strata (historical weight)
test_that("strata (MH weights)",{
ecount.favorable <- coef(e.BT, statistic = "count.favorable", stratified = TRUE)[,1]
ecount.unfavorable <- coef(e.BT, statistic = "count.unfavorable", stratified = TRUE)[,1]
count.strata <- table(ddd$stratum)
n.pairs <- e.BT@n.pairs

res_mix_equalU$Win_statistic$Win_Ratio
coef(e.BT, statistic = "winRatio")
sum(ecount.favorable)/sum(ecount.unfavorable)
## weighted.mean(ecount.favorable/ecount.unfavorable, n.pairs/sum(n.pairs))

res_mix_equalU$Win_statistic$Net_Benefit
coef(e.BT, statistic = "netBenefit")
sum(ecount.favorable)/sum(n.pairs) - sum(ecount.unfavorable)/sum(n.pairs)
weighted.mean(ecount.favorable/n.pairs - ecount.unfavorable/n.pairs, n.pairs/sum(n.pairs))
})


##----------------------------------------------------------------------
### test-BuyseTest-strata.R ends here
