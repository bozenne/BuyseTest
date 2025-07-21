### test-BuyseTest-sensitivity.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 21 2025 (09:41) 
## Version: 
## Last-Updated: Jul 21 2025 (09:50) 
##           By: Brice Ozenne
##     Update #: 5
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

context("Check sensitivity for BuyseTest object")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = FALSE,
                  method.inference = "u-statistic",
                  trace = 0)

## * simulate data (taken from the GPC book)

argsSurv <- list(name = c("OS","PFS"),
                 name.censoring = c("statusOS","statusPFS"),
                 scale.C = c(8.995655, 4.265128),
                 scale.T = c(13.76543, 7.884477),
                 shape.C = c(1.28993, 1.391015),
                 shape.T = c(1.275269, 1.327461),
                 scale.censoring.C = c(34.30562, 20.748712),
                 scale.censoring.T = c(27.88519, 17.484281),
                 shape.censoring.C = c(1.369449, 1.463876),
                 shape.censoring.T = c(1.490881, 1.835526))

argsTox <- list(name = "toxicity",
                p.C =  c(1.17, 2.92, 36.26, 39.18, 19.88, 0.59)/100,
                p.T = c(3.51, 4.09, 23.39, 47.37, 21.05, 0.59)/100,
                rho.T = 1, rho.C = 1)

set.seed(1)
dt.data <- simBuyseTest(n.T = 200, n.C = 200,
                        argsBin = argsTox,
                        argsCont = NULL,
                        argsTTE = argsSurv,
                        level.strata = c("M","F"), names.strata = "gender")
dt.data$toxicity.num <- as.numeric(dt.data$toxicity)


test_that("BuyseTest - sensitivity",{
  
    eRBB.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                         data=dt.data, scoring.rule = "Gehan", trace = FALSE)

    eRBB.Se <- sensitivity(eRBB.BT, threshold = list(1:2,c(0,5)),
                           band = TRUE, adj.p.value = TRUE, seed = 10, trace = FALSE)

    ## Gehan (to save time)
    GS <- data.frame("toxicity.num" = c(1, 2, 1, 2), 
                     "OS" = c(0, 0, 5, 5), 
                     "estimate" = c(-0.0395, 0.0569, -0.058575, 0.02695), 
                     "se" = c(0.05698582, 0.04022624, 0.05632815, 0.04002103), 
                     "lower.ci" = c(-0.1502393, -0.02213292, -0.16782072, -0.05149468), 
                     "upper.ci" = c(0.07221819, 0.13522617, 0.05209192, 0.10506416), 
                     "null" = c(0, 0, 0, 0), 
                     "p.value" = c(0.48866503, 0.15811355, 0.29949904, 0.5009029), 
                     "lower.band" = c(-0.16306712, -0.0314249, -0.18044898, -0.06069361), 
                     "upper.band" = c(0.08528717, 0.14434312, 0.06507041, 0.11418123), 
                     "adj.p.value" = c(0.667633, 0.24483912, 0.43493874, 0.68133707))
    ## Peron
    ## GS <- data.frame("toxicity.num" = c(1, 2, 1, 2), 
    ##                  "OS" = c(0, 0, 5, 5), 
    ##                  "estimate" = c(-0.01972149, 0.10896273, -0.03347608, 0.09511226), 
    ##                  "se" = c(0.05440853, 0.03824761, 0.05411638, 0.03880383), 
    ##                  "lower.ci" = c(-0.12573534, 0.03351983, -0.13877248, 0.01865006), 
    ##                  "upper.ci" = c(0.08673769, 0.18317093, 0.07256876, 0.17046839), 
    ##                  "null" = c(0, 0, 0, 0), 
    ##                  "p.value" = c(0.71707061, 0.0047093, 0.53648837, 0.01483855), 
    ##                  "lower.band" = c(-0.14136068, 0.02221489, -0.15426777, 0.00720646), 
    ##                  "upper.band" = c(0.10250434, 0.1940821, 0.08830157, 0.1815593), 
    ##                  "adj.p.value" = c(0.93032879, 0.01075334, 0.77963633, 0.03122924))

    expect_equivalent(eRBB.Se, GS, tol = 1e-3)

})

##----------------------------------------------------------------------
### test-BuyseTest-sensitivity.R ends here
