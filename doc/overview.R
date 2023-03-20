## chunk 2
library(BuyseTest)
library(data.table)

## chunk 3
library(survival)
head(veteran)

## chunk 4
utils::packageVersion("BuyseTest")

## chunk 5
sessionInfo()

## * Performing generalized pairwise comparisons (GPC) using the =BuyseTest= function

## chunk 6
BT <- BuyseTest(data = veteran, 
                endpoint = "time", 
                type = "timeToEvent", 
                treatment = "trt", 
                status = "status", 
                threshold = 20)

## chunk 7
BT.f <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status"),
                  data = veteran)

## chunk 8
BT.f@call <- list(); BT@call <- list();
testthat::expect_equal(BT.f,BT)

## ** Displaying the results

## chunk 9
summary(BT)

## chunk 10
summary(BT, percentage = FALSE)

## chunk 11
summary(BT, statistic = "winRatio")

## chunk 12
confint(BT, statistic = "winRatio")

## ** Stratified GPC

## chunk 13
ffstrata <- trt ~ tte(time, threshold = 20, status = "status") + celltype
BTstrata <- BuyseTest(ffstrata, data = veteran, trace = 0)

## chunk 14
summary(BTstrata, type.display = c("endpoint","threshold","strata",
                              "total","favorable","unfavorable","delta","Delta"))

## chunk 15
BTstrataCMH <- BuyseTest(ffstrata, data = veteran, trace = 0, pool.strata = "CMH")
summary(BTstrataCMH, type.display = c("endpoint","threshold","strata",
                                      "total","favorable","unfavorable","delta","Delta"))

## chunk 16
nStrata <- table(veteran$celltype, veteran$trt)
nStrata

## chunk 17
nPairStrata <- nStrata[,1]*nStrata[,2]
wStrata <- rbind(default = nPairStrata,
                 CMH = nPairStrata/rowSums(nStrata)) 
round(100*wStrata/rowSums(wStrata),2)

## chunk 18
confint(BTstrata, stratified = TRUE)

## ** Using multiple endpoints

## chunk 19
ff2 <- trt ~ tte(time, threshold = 20, status = "status") + cont(karno, threshold = 0)
BT.H <- BuyseTest(ff2, data = veteran, trace = 0)
summary(BT.H)

## chunk 20
BT.nH <- BuyseTest(ff2, hierarchical = FALSE, data = veteran, trace = 0)
summary(BT.nH)

## chunk 21
ff2w <- trt ~ tte(time, threshold = 20, status = "status", weight = 0.8)
ff2w <- update.formula(ff2w, . ~ . + cont(karno, threshold = 0, weight = 0.2))
BT.nHw <- BuyseTest(ff2w, hierarchical = FALSE, data = veteran, trace = 0)
summary(BT.nHw, print = FALSE)$table.print[,-13]

## chunk 22
confint(BuyseTest(trt ~ cont(karno, threshold = 0), data = veteran, trace = 0))

## chunk 23
confint(BT.nHw, cumulative = FALSE)

## chunk 24
BuyseMultComp(BT.nHw, cumulative = FALSE, endpoint = 1:2)

## ** What if smaller is better?

## chunk 25
ffop <- trt ~ tte(time, status = "status", threshold = 20, operator = "<0")
BTinv <- BuyseTest(ffop, data = veteran, trace = 0)
summary(BTinv)

## ** Stopping comparison for neutral pairs

## chunk 26
dt.sim <- data.table(Id = 1:2,
                     treatment = c("Yes","No"),
                     tumor = c("Yes","Yes"),
                     size = c(15,20))
dt.sim

## chunk 27
BT.pair <- BuyseTest(treatment ~ bin(tumor) + cont(size, operator = "<0"), data = dt.sim,
                     trace = 0, method.inference = "none")
summary(BT.pair)

## chunk 28
BT.pair2 <- BuyseTest(treatment ~ bin(tumor) + cont(size, operator = "<0"), data = dt.sim,
                     trace = 0, method.inference = "none", neutral.as.uninf = FALSE)
summary(BT.pair2)

## ** What about p-value and confidence intervals?

## chunk 29
BT.perm <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status"),
                     data = veteran, trace = 0, method.inference = "permutation",
                     seed = 10) 
summary(BT.perm)

## chunk 30
BT.boot <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status"),
                     data = veteran, trace = 0, method.inference = "bootstrap",
                     seed = 10) 
summary(BT.boot)

## chunk 31
BT.ustat <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status"),
                      data = veteran, trace = 0, method.inference = "u-statistic") 
summary(BT.ustat)

## chunk 32
set.seed(10)
sapply(1:10, function(i){mean(rbinom(1e4, size = 1, prob = 0.05))})

## ** Sensitivity analysis

## chunk 33
BTse.ustat <- sensitivity(BT.ustat, threshold = seq(0,500, length.out=10),
                          band = TRUE, trace = FALSE)
BTse.ustat[,c("time","estimate","se","lower.ci","upper.ci","null","lower.band","upper.band")]

## chunk 34
library(ggplot2)
autoplot(BTse.ustat)

## chunk 36
BTse.H <- sensitivity(BT.H, trace = FALSE,
                      threshold = list(time = seq(0,500,length = 10), karno = c(0,40,80)))
head(BTse.H)

## chunk 37
grid <- expand.grid(list("time_t20" = seq(0,500,length = 10), "karno" = c(0,40,80)))
cbind(head(grid)," " = "  ...   ",tail(grid))
BTse.H2 <-sensitivity(BT.H, threshold = grid, trace = FALSE)
range(BTse.H-BTse.H2)

## chunk 38
autoplot(BTse.H, col = NA)
##  alternative display:
## autoplot(BTse.H, position  = position_dodge(width = 15))

## * Getting additional inside: looking at the pair level
## ** Extracting the contribution of each pair to the statistic

## chunk 40
form <- trt ~ tte(time, threshold = 20, status = "status") + cont(karno)
BT.keep <- BuyseTest(form,
                     data = veteran, keep.pairScore = TRUE, 
                     trace = 0, method.inference = "none")

## chunk 41
getPairScore(BT.keep, endpoint = 1)

## chunk 42
veteran[c(70,1),]

## chunk 43
getPairScore(BT.keep, endpoint = 1)[, mean(favorable) - mean(unfavorable)]

## chunk 44
BT.keep

## ** Extracting the survival probabilities

## chunk 45
BuyseTest.options(keep.survival = TRUE, precompute = FALSE)
BT.keep2 <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status") + cont(karno),
                      data = veteran, keep.pairScore = TRUE, scoring.rule = "Peron",
                      trace = 0, method.inference = "none")

## chunk 46
outSurv <- getSurvival(BT.keep2, endpoint = 1, strata = 1)
str(outSurv)

## *** Computation of the score with only one censored event

## chunk 47
getPairScore(BT.keep2, endpoint = 1, rm.withinStrata = FALSE)[91]

## chunk 48
veteran[c(22,71),]

## chunk 49
iSurv <- outSurv$survTimeC[22,] 
iSurv

## chunk 50
Sc97 <- iSurv["survivalC_0"] 
Sc97

## chunk 51
iSurv <- outSurv$survTimeT[2,] ## survival at time 112+20
iSurv

## chunk 52
Sc132 <- iSurv["survivalC+threshold"] 
Sc132

## chunk 53
Sc132/Sc97

## *** Computation of the score with two censored events

## chunk 54
head(outSurv$survJumpT)

## chunk 55
getPairScore(BT.keep2, endpoint = 1, rm.withinStrata = FALSE)[148]

## chunk 56
veteran[c(10,72),]

## chunk 57
calcInt <- function(...){ ## no need for the functionnal derivative of the score 
    BuyseTest:::.calcIntegralSurv_cpp(..., 
                                      returnDeriv = FALSE, 
                                      derivSurv = matrix(0), 
                                      derivSurvD = matrix(0))
}

## chunk 58
denom <- as.double(outSurv$survTimeT[3,"survivalT_0"] * outSurv$survTimeC[10,"survivalC_0"])
M <- cbind("favorable" = -calcInt(outSurv$survJumpC, start = 100, 
                                  lastSurv = outSurv$lastSurv[2],
                                  lastdSurv = outSurv$lastSurv[1])/denom,
           "unfavorable" = -calcInt(outSurv$survJumpT, start = 87, 
                                    lastSurv = outSurv$lastSurv[1],
                                    lastdSurv = outSurv$lastSurv[2])/denom)
rownames(M) <- c("lowerBound","upperBound")
M

## chunk 59
outSurv$lastSurv

## * Dealing with missing values or/and right censoring 

## chunk 60
set.seed(10)
dt <- simBuyseTest(1e2, latent = TRUE, argsCont = NULL,
                   argsTTE = list(scale.T = 1/2, scale.C = 1,
                                  scale.censoring.C = 1, scale.censoring.T = 1))
dt[, eventtimeCensoring := NULL]
dt[, status1 := 1]
head(dt)

## chunk 61
100*dt[,mean(status==0)]

## chunk 62
BuyseTest(treatment ~ tte(eventtimeUncensored, status1, threshold = 0.5),
          data = dt,
          scoring.rule = "Gehan", method.inference = "none", trace = 0)

## ** Gehan's scoring rule

## chunk 63
e.G <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
          data = dt, scoring.rule = "Gehan", trace = 0)
summary(e.G, print=FALSE)$table.print

## ** Peron's scoring rule

## chunk 64
e.P <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
          data = dt, scoring.rule = "Peron", trace = 0)
summary(e.P, print=FALSE)$table.print

## chunk 65
dt[,.SD[which.max(eventtime)],by="treatment"]

## chunk 66
library(prodlim)
e.prodlim <- prodlim(Hist(eventtime, status) ~ treatment, data = dt)

## chunk 67
e.P1 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                  model.tte = e.prodlim,
                  data = dt, scoring.rule = "Peron", trace = 0)
summary(e.P1, print=FALSE)$table.print

## chunk 68
attr(e.prodlim, "iidNuisance") <- TRUE
e.P2 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                  model.tte = e.prodlim,
                  data = dt, scoring.rule = "Peron", trace = 0)
summary(e.P2, print=FALSE)$table.print

## chunk 69
library(survival)
e.survreg <- survreg(Surv(eventtime, status) ~ treatment, data = dt, 
                     dist = "weibull")
attr(e.survreg, "iidNuisance") <- TRUE

## chunk 70
e.P3 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                  model.tte = e.survreg,
                  data = dt, scoring.rule = "Peron", trace = 0)
summary(e.P3, print=FALSE)$table.print

## chunk 71
e.TTEM <- BuyseTTEM(e.survreg, treatment = "treatment", iid = TRUE, n.grid = 2500)
attr(e.TTEM, "iidNuisance") <- TRUE
str(e.TTEM$peron$jumpSurvHaz[[1]][[1]])

## chunk 72
e.P4 <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                  model.tte = e.TTEM,
                  data = dt, scoring.rule = "Peron", trace = 0)
summary(e.P4, print=FALSE)$table.print

## ** Correction via inverse probability-of-censoring weights (IPCW)

## chunk 73
BT <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                data = dt, keep.pairScore = TRUE, trace = 0,
                scoring.rule = "Gehan", method.inference = "none", correction.uninf = 2)
summary(BT)

## chunk 74
iScore <- getPairScore(BT, endpoint = 1)

## chunk 75
iScore[,.SD[1], 
       .SDcols = c("favorableC","unfavorableC","neutralC","uninfC"),
       by = c("favorable","unfavorable","neutral","uninf")]

## chunk 76
iScore[,1/mean(favorable + unfavorable + neutral)]

## ** Correction at the pair level

## chunk 77
BT <- BuyseTest(treatment ~ tte(eventtime, status, threshold = 0.5),
                data = dt, keep.pairScore = TRUE, trace = 0,
                scoring.rule = "Gehan", method.inference = "none", correction.uninf = TRUE)
summary(BT)

## chunk 78
iScore <- getPairScore(BT, endpoint = 1)
iScore[,.SD[1], 
       .SDcols = c("favorableC","unfavorableC","neutralC","uninfC"),
       by = c("favorable","unfavorable","neutral","uninf")]

## chunk 79
iScore[, .(favorable = weighted.mean(favorable, w = 1-uninf), 
           unfavorable = weighted.mean(unfavorable, w = 1-uninf), 
           neutral = weighted.mean(neutral, w = 1-uninf))]

## ** Note on the use of the corrections

## chunk 80
set.seed(10);n <- 250; 
df <- rbind(data.frame(group = "T1", time = rweibull(n, shape = 1, scale = 2), status = 1),
            data.frame(group = "T2", time = rweibull(n, shape = 2, scale = 1.8), status = 1))
df$censoring <- runif(NROW(df),0,2)
df$timeC <- pmin(df$time,df$censoring)
df$statusC <- as.numeric(df$time<=df$censoring)
plot(prodlim(Hist(time,status)~group, data = df)); title("complete data");
plot(prodlim(Hist(timeC,statusC)~group, data = df)); title("right-censored data");

## chunk 82
BuyseTest.options(method.inference = "none")
e.ref <- BuyseTest(group ~ tte(time,status), data = df, trace = FALSE)
s.ref <- summary(e.ref, print = FALSE)$table[1,c("favorable","unfavorable","neutral","uninf","Delta")]
s.ref

## chunk 83
e.correction <- BuyseTest(group ~ tte(timeC,statusC)+cont(time), data = df, trace = FALSE, correction.uninf = TRUE)
s.correction <- summary(e.correction, print = FALSE)$table[1,c("favorable","unfavorable","neutral","uninf","Delta")]

## chunk 84
e.Peron <- BuyseTest(group ~ tte(timeC,statusC), data = df, trace = FALSE)
s.Peron <- summary(e.Peron,print = FALSE)$table[1,c("favorable","unfavorable","neutral","uninf","Delta")]
rbind("reference" = s.ref,
      "no correction" = s.Peron,
      "correction" = s.correction)

## * Simulating data using =simBuyseTest=

## chunk 85
set.seed(10)
simBuyseTest(n.T = 5, n.C = 5)

## chunk 86
set.seed(10)
argsCont <- list(mu.T = c(5,5), mu.C = c(10,10), 
                 sigma.T = c(1,1), sigma.C = c(1,1),
                 name = c("tumorSize","score"))
dt <- simBuyseTest(n.T = 5, n.C = 5,
                   argsCont = argsCont)
dt

## * Power calculation using =powerBuyseTest=

## chunk 87
simFCT <- function(n.C, n.T){
     out <- rbind(cbind(Y=stats::rt(n.C, df = 5), group=0),
                  cbind(Y=stats::rt(n.T, df = 5) + 1/2, group=1))
     return(data.table::as.data.table(out))
}
simFCT(101,101)

## chunk 88
null <- c("netBenefit" = 0)

## chunk 89
powerW <- powerBuyseTest(sim = simFCT, method.inference = "u-statistic", null = null,
                         sample.size = c(5,10,20,30,50,100),                         
                         formula = group ~ cont(Y), 
                         n.rep = 1000, seed = 10, cpus = 6, trace = 0)

## chunk 90
summary(powerW)

## chunk 91
nW <- powerBuyseTest(sim = simFCT, method.inference = "u-statistic", null = null,
                     power = 0.8, max.sample.size = 10000,                     
                     formula = group ~ cont(Y), 
                     n.rep = 1000, seed = 10, cpus = 6, trace = 0)
summary(nW)

## * Modifying default options

## chunk 92
BuyseTest.options("trace")

## chunk 93
BuyseTest.options(trace = 0)

## chunk 94
BuyseTest.options(summary.display = list(c("endpoint","threshold","delta","Delta","information(%)")))
summary(BT)

## chunk 95
BuyseTest.options(reinitialise = TRUE)

## * References
