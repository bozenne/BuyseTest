### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: jan  7 2021 (10:08) 
##           By: Brice Ozenne
##     Update #: 180
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

context("Check that bugs that have been reported are fixed \n")


## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * Joris: jeudi 5 avril 2018 à 14:57
dt.sim <- data.table(
    ttt=c(rep(0,3),rep(1,3)),
    timeOS = c(10,20,30,15,20,35),
    eventOS = c(1,1,0,0,1,1),
    Mgrade.tox = -c(1,2,3,2,4,2)
)

test_that("number of pairs - argument neutral.as.uninf", {

    for(iCorrection in c(FALSE,TRUE)){ ## iCorrection <- TRUE ; iCorrection <- FALSE
        BT.T <- BuyseTest(ttt~TTE(timeOS,threshold=0,status=eventOS) + cont(Mgrade.tox,threshold=0),
                          data = dt.sim,
                          neutral.as.uninf = TRUE, scoring.rule = "Gehan", correction.uninf = iCorrection)
        BTS.T <- as.data.table(summary(BT.T, print = FALSE, percentage = FALSE)$table)
        BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,status=eventOS) + cont(Mgrade.tox,threshold=0),
                          data = dt.sim,
                          neutral.as.uninf = FALSE, scoring.rule = "Gehan", correction.uninf = iCorrection)
        BTS.F <- as.data.table(summary(BT.F, print = FALSE, percentage = FALSE)$table)

        ## neutral.as.uninf does not impact the results for first endpoint
        expect_equal(BTS.T[1,c("favorable","unfavorable","neutral","uninf","delta","Delta")],
                     BTS.F[1,c("favorable","unfavorable","neutral","uninf","delta","Delta")])

        ## check consistency of the number of pairs
        ## neutral.as.uninf = TRUE
        ## summary(BT.T)
        expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", favorable+unfavorable+neutral+uninf],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", total])
        expect_equal(BTS.T[endpoint == "timeOS" & strata == "global", neutral+uninf],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", total])
        expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", total],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", favorable+unfavorable+neutral+uninf])

        ## neutral.as.uninf = FALSE
        expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", favorable+unfavorable+neutral+uninf],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", total])
        expect_equal(BTS.F[endpoint == "timeOS" & strata == "global", uninf],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", total])
        expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", total],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", favorable+unfavorable+neutral+uninf])

        ## compared to known value
        if(iCorrection == FALSE){
            keep.col <- c("endpoint","threshold","strata","weight","total","favorable","unfavorable","neutral","uninf","delta","Delta")
            test <- as.data.table(summary(BT.T, print = FALSE)$table[,keep.col])
            GS <- data.table("endpoint" = c("timeOS", "timeOS", "Mgrade.tox", "Mgrade.tox"), 
                             "threshold" = c(1e-12, 1e-12, 1e-12, 1e-12), 
                             "strata" = c("global", "1", "global", "1"),
                             "weight" = c(1, 1, 1, 1), 
                             "total" = c(100.00000, 100.00000,  44.44444,  44.44444), 
                             "favorable" = c(44.44444, 44.44444, 22.22222, 22.22222), 
                             "unfavorable" = c(11.11111, 11.11111, 11.11111, 11.11111), 
                             "neutral" = c(11.11111, 11.11111, 11.11111, 11.11111), 
                             "uninf" = c(33.33333, 33.33333,  0.00000,  0.00000), 
                             "delta" = c(0.3333333, 0.3333333, 0.1111111, 0.1111111), 
                             "Delta" = c(0.3333333, NA, 0.4444444, NA))
            ##    butils::object2script(test)

            attr(test,"index") <- NULL
            expect_equal(test, GS, tol = 1e-6)
            ## class(BTS.T[["n.resampling"]])
            ## class(GS[["n.resampling"]])
        }
    }
})

## * Emeline T: samedi 26 mai 2018 à 14:39 (Version 1.0)
## ERROR: Error in xy.coords(x, y, setLab = FALSE) : 'x' and 'y' lengths differ
## butils:::object2script(data[175:325,], digits = 8)
data <- data.frame("X" = c(175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325), 
                   "trt" = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0), 
                   "time" = c(0.34972271, 0.15919528, 0.27021802, 0.95994001, 0.46312769, 0.31997409, 0.23637537, 0.89707932, 0.01708799, 0.39366592, 0.50017673, 0.446804, 0.50040844, 0.32638758, 0.6262995, 0.14755089, 0.12050248, 0.07458989, 0.04339593, 0.59982912, 1.41101136, 0.2675838, 0.52521586, 1.14631933, 0.74191272, 1.2196955, 0.02732667, 1.3869635, 0.75430971, 0.3780356, 0.42434206, 1.28254783, 0.65964535, 0.80568326, 1.09058069, 0.14099648, 1.30095204, 0.69223441, 1.33892841, 0.73062582, 0.28980283, 1.74724314, 0.85952631, 0.40828457, 1.26493484, 0.96396552, 0.75849828, 0.70308743, 1.71091642, 1.01266995, 0.29350899, 0.79999462, 0.90685983, 0.2697463, 0.92647206, 0.00936012, 0.69425291, 0.82894713, 0.28051478, 1.40047767, 0.83924557, 0.61605441, 0.56216195, 0.68796769, 1.83362936, 0.45955409, 1.381266, 1.34455702, 0.30326241, 2.42955884, 0.53467431, 1.00931952, 1.11490004, 0.72048666, 0.07125682, 0.34582823, 0.33357166, 0.47453535, 0.27259304, 0.60673207, 0.95520791, 0.05198433, 0.82662585, 1.15532297, 0.87506277, 1.37889663, 0.12846039, 0.68540728, 0.77377909, 0.81177511, 0.29095231, 2.02666276, 0.21531326, 0.45024274, 1.43151175, 0.46492612, 0.14985886, 0.22205914, 1.59582145, 0.76701798, 1.23825982, 0.33712561, 1.07747869, 0.06973708, 1.27342747, 0.42610371, 1.0686674, 2.03964558, 0.5787245, 1.05125486, 0.24393524, 1.02678662, 0.2725943, 0.59435986, 0.32627314, 0.39337226, 0.71167895, 0.58597973, 0.3605633, 1.24886565, 0.43183396, 0.75826836, 0.22063575, 0.28832416, 0.16407274, 0.91388552, 0.62053192, 2.46164696, 0.28193246, 0.33575549, 0.51327929, 0.90610562, 0.43071919, 1.392834, 0.69855789, 0.81717857, 0.46312768, 0.11466708, 0.42909682, 0.29334352, 0.76480274, 0.80197241, 0.40497033, 0.68113025, 0.98833506, 0.58629864, 0.00627822, 0.35254414, 0.52416901, 0.67108879, 0.49179438), 
                   "event" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                   stringsAsFactors = FALSE)

BT_tau0 <- BuyseTest(data=data,
                     treatment="trt",
                     endpoint="time",
                     type="timeToEvent",
                     threshold=as.numeric(0),
                     status="event",
                     scoring.rule="Peron",
                     method.inference = "none",
                     cpus=1,
                     trace = 0)

## * Brice: 09/06/18 6:51 (Tied event with tte endpoint)
## when computing the integral for peron with double censoring
## the ordering of the data modified the ouput
## this has been correct with version 1.4
data(veteran,package="survival")

test_that("ordering of tied event does not affect BuyseTest", {
    ## veteran2[veteran2$time==100,]
    BT.all <- BuyseTest(trt ~ tte(time, threshold = 0, status = "status"),
                        data = veteran, scoring.rule = "Peron", method.inference = "none", correction.uninf = FALSE)

    veteran1 <- veteran[order(veteran$time,veteran$status),c("time","status","trt")]
    ## veteran1[veteran2$time==100,]
    BT1.all <- BuyseTest(trt ~ tte(time, threshold = 0, status = "status"),
                         data = veteran1, scoring.rule = "Peron", method.inference = "none", correction.uninf = FALSE)

    veteran2 <- veteran[order(veteran$time,-veteran$status),c("time","status","trt")]
    ## ## veteran2[veteran2$time==100,]
    BT2.all <- BuyseTest(trt ~ tte(time, threshold = 0, status = "status"),
                         data = veteran2, scoring.rule = "Peron", method.inference = "none", correction.uninf = FALSE)

    ## effect of the ordering
    expect_equal(coef(BT.all, statistic = "winRatio"), coef(BT1.all, statistic = "winRatio"))
    expect_equal(coef(BT.all, statistic = "winRatio"), coef(BT2.all, statistic = "winRatio"))

    ## number of pairs
    expect_equal(as.double(BT.all@n.pairs), prod(table(veteran$trt)), tol = 1e-5)
    expect_equal(as.double(BT1.all@n.pairs), prod(table(veteran$trt)), tol = 1e-5)
    expect_equal(as.double(BT2.all@n.pairs), prod(table(veteran$trt)), tol = 1e-5)

    ## values of the pairs
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["favorable"]]>=0))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["favorable"]]<=1))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["unfavorable"]]>=0))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["unfavorable"]]<=1))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["neutral"]]>=0))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["neutral"]]<=1))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["uninformative"]]>=0))
    expect_true(all(getPairScore(BT.all, endpoint = 1)[["uninformative"]]<=1))

    expect_true(all(getPairScore(BT.all, endpoint = 1)[,favorable + unfavorable]<=1+1e-12)) ## tolerance

    ## survival
    ## getSurvival(BT.all, endpoint = 1, strata = 1)$lastSurv: only 0 so no uninformative paris
    expect_equal(as.double(coef(BT.all, statistic = "count.uninf", cumulative = FALSE)), 0.0, tol = 1e-12)
    
    ## result
    expect_equal(as.double(coef(BT.all, statistic = "winRatio")), 0.8384569, tol = 1e-5)
    

})

## * Brice: 26/09/18 x:xx (Multiple thresholds in Julien's simulations)

HR1 <- 0.65
TpsFin <- 60 #values for Taux.Censure 
HazC <- 0.1

set.seed(10)
HazT <- 0.1*(HR1)
n.Treatment <- 100
n.Control <- 100
n <- n.Treatment+n.Control
group <- c(rep(1, n.Treatment),rep(0, n.Control))

TimeEvent.Ctr <- rexp(n.Control,HazC)
TimeEvent.Tr <- rexp(n.Control,HazT)

TimeEvent<-c(TimeEvent.Tr,TimeEvent.Ctr)
Time.Cens<-runif(n,0,TpsFin)
Time<-pmin(Time.Cens,TimeEvent)
Event<-Time==TimeEvent
Event<-as.numeric(Event)

tab <- data.frame(group,Time,Event, stringsAsFactors = FALSE)

test_that("Multiple thresholds",{
    BuyseresPer <- BuyseTest(data=tab,
                             endpoint=c("Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time"),
                             treatment="group",
                             type=c("TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE"),
                             status=c("Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event"),
                             threshold=c(42,39,36,33,30,27,24,21,18,15,12,9,6,3,0),
                             n.resampling=500,
                             trace=0,
                             scoring.rule="Peron",
                             correction.uninf=F,
                             method.inference="none")

    resS <- as.data.table(summary(BuyseresPer, print = FALSE)$table)

    ## pairs are correctly transfered from one endpoint to another
    expect_equal(resS[strata == "global" & threshold > tail(threshold,1), neutral + uninf],
                 resS[strata == "global" & threshold < threshold[1], total], tol = 1e-2)

    ## butils::object2script(as.double(BuyseresPer@count.favorable), digit = 2)
    GS <- c(260.64, 35.93, 37.33, 147.32, 272.14, 263.6, 235.7, 213.21, 390.29, 408.73, 514.7, 514.34, 744.78, 865.21, 1095.26)
    expect_equal(as.double(coef(BuyseresPer, statistic = "count.favorable", cumulative = FALSE)), GS, tol = 1e-5)
    ## butils::object2script(as.double(BuyseresPer@count.unfavorable), digit = 2)
    GS <- c(0, 0, 6.97, 25.66, 43.89, 34.8, 46.38, 105.42, 199.85, 338.55, 407.72, 521.83, 548.02, 782.94, 938.8)
    expect_equal(as.double(coef(BuyseresPer, statistic = "count.unfavorable", cumulative = FALSE)), GS, tol = 1e-5)
    ## butils::object2script(as.double(BuyseresPer@count.neutral), digit = 2)
    GS <- c(9580.24173298, 9580.24173298, 9573.26856493, 9433.45578036, 9140.84153639, 8860.00546388, 8577.9282735, 8259.29654759, 7675.01664374, 6933.58435144, 6011.16786378, 4975.00117239, 3682.20846492, 2034.06256111, 0)
    expect_equal(as.double(coef(BuyseresPer, statistic = "count.neutral", cumulative = FALSE)), GS, tol = 1e-5)
    ## butils::object2script(as.double(BuyseresPer@count.uninf), digit = 2)
    GS <- c(159.12095011, 123.19041299, 85.85998481, 52.68680886, 29.27044937, 11.70817975, 11.70817975, 11.70817975, 5.85408987, 0, 0, 0, 0, 0, 0)
    expect_equal(as.double(coef(BuyseresPer, statistic = "count.uninf", cumulative = FALSE)), GS, tol = 1e-1)
    
    ## butils::object2script(as.double(BuyseresPer@delta.netBenefit), digit = 5)
    GS <- c(0.02606, 0.00359, 0.00304, 0.01217, 0.02282, 0.02288, 0.01893, 0.01078, 0.01904, 0.00702, 0.0107, -0.00075, 0.01968, 0.00823, 0.01565)
    expect_equal(as.double(coef(BuyseresPer, statistic = "netBenefit", cumulative = FALSE)), GS, tol = 1e-3)
    ## butils::object2script(as.double(BuyseresPer@delta.winRatio), digit = 5)
    GS <- c(Inf, Inf, 5.35344, 5.74093, 6.19986, 7.57457, 5.08161, 2.02241, 1.95291, 1.2073, 1.2624, 0.98564, 1.35904, 1.10508, 1.16666)
    expect_equal(as.double(coef(BuyseresPer, statistic = "winRatio", cumulative = FALSE)), GS, tol = 1e-3)
})

## * Brice: 30/10/18 4:36 Neutral pairs with 0 threshold
df <- data.frame("survie" = c(2.1, 4.1, 6.1, 8.1, 4, 6, 8, 10),
                 "event" = c(1, 1, 1, 0, 1, 0, 0, 1),
                 "group" = c(0, 0, 0, 0, 1, 1, 1, 1),
                 "score" = 1,
                 stringsAsFactors = FALSE)

test_that("1 TTE endpoint - Gehan (no correction)", {
    Peron <- BuyseTest(group ~ tte(survie, status = event, threshold = 0),
                       data = df, 
                       scoring.rule = "Peron", correction.uninf = FALSE)

    expect_equal(as.double(coef(Peron, statistic = "count.neutral", cumulative = FALSE)),0) ## should not be any neutral pair with a threshold of 0
})

## * Hickey, Graeme: 8 mars 2019 14:54 p-value permutation
## I have one question, which I hope you can help with.
## If using method.inference = “permutation”, the P-values are slightly different for the net benefit and win ratio summary methods.
## However, if you use using method.inference = “bootstrap”, the P-values are identical, as I would expect.
## Can you explain why they differ with the permutation test?

set.seed(1)
dt <- simBuyseTest(50)

test_that("same p.value (permutation test) for winRatio and net Benefit", {
    e.perm <- BuyseTest(treatment ~ bin(toxicity), data = dt,
                        method.inference = "permutation", n.resampling = 100, trace = 0)
    netBenefit.perm <- confint(e.perm, statistic = "netBenefit")
    winRatio.perm <- confint(e.perm, statistic = "winRatio")

    Delta.netBenefit <- coef(e.perm, statistic = "netBenefit")
    Delta.winRatio <- coef(e.perm, statistic = "winRatio")
    DeltaResampling.netBenefit <- e.perm@DeltaResampling[,1,"netBenefit"]
    DeltaResampling.winRatio <- e.perm@DeltaResampling[,1,"winRatio"]
    
    manual <- c(netBenefit = mean(abs(DeltaResampling.netBenefit) >= abs(Delta.netBenefit)),
                netBenefit.atanh = mean(abs(atanh(DeltaResampling.netBenefit)) >= abs(atanh(Delta.netBenefit))),
                winRatio = mean(abs(DeltaResampling.winRatio-1) >= abs(Delta.winRatio-1)),
                winRatio.log = mean(abs(log(DeltaResampling.winRatio)) >= abs(log(Delta.winRatio)))
                )

    expect_equal(netBenefit.perm[,"p.value"], winRatio.perm[,"p.value"])
    expect_equal(unname(manual["netBenefit"]), netBenefit.perm[,"p.value"])
    expect_equal(unname(manual["winRatio.log"]), winRatio.perm[,"p.value"])

    ## note CI are not agreeing with p-values
    confint(e.perm, statistic = "netBenefit", conf.level = 1-0.48)
    confint(e.perm, statistic = "winRatio", conf.level = 1-0.48)

})


## * Alice, Brouquet-Laglair: 3 avril 2019 p-value bootstrap
df <- rbind(data.frame(score = rep(1,5),
                       tox = 0,
                       group = 1,
                       stringsAsFactors = FALSE),
            data.frame(score = rep(0,5),
                       tox = 0,
                       group = 0,
                       stringsAsFactors = FALSE)
            )

test_that("BuyseTest without variability", {
    e.BT_ustat <- BuyseTest(group ~ bin(tox) + cont(score), data = df,
                            method.inference = "u-statistic", trace = 0)
    e.BT_boot <- BuyseTest(group ~ bin(tox) + cont(score), data = df,
                           method.inference = "studentized bootstrap",
                           n.resampling = 10, trace = 0)

    confintTempo <- confint(e.BT_ustat)
    expect_equal(unname(confintTempo[,"p.value"]),1:0)
    confintTempo <- suppressMessages(confint(e.BT_boot, transformation = FALSE))
    expect_equal(unname(confintTempo[,"p.value"]),1:0)
})

## * graemeleehickey (issue #2 on Github): 8 september 2019 p-value bootstrap
test_that("Boostrap - issue in the summary", {
    data(veteran,package="survival")
    BT.keep <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status") + cont(karno),
                         data = veteran, keep.pairScore = TRUE, scoring.rule = "Gehan", 
                         trace = 0, method.inference = "bootstrap", n.resampling = 20, seed = 10)
    expect_error(capture.output(summary(BT.keep, statistic = "winRatio")), regexp = NA) ## no error
})

## * graemeleehickey (issue #3 on Github): 22 september 2019 BuysePower
test_that("BuysePower - error in print", {
    simFCT <- function(n.C, n.T){
        out <- data.table(Y=rnorm(n.C+n.T),
                          T=c(rep(1,n.C),rep(0,n.T))
                          )
        return(out)
    }

    ## the error was when setting trace to 4
    tempo <- capture.output({
        xx <- powerBuyseTest(sim = simFCT, sample.sizeC = c(100), sample.sizeT = c(100), n.rep = 2,
                             formula = T ~ cont(Y), method.inference = "u-statistic", trace = 4,
                             seed = 10)
    })

    yy <- powerBuyseTest(sim = function(n.C, n.T){
        out <- data.table(Y=rnorm(n.C+n.T),
                          T=c(rep(1,n.C),rep(0,n.T))
                          )
        return(out)
    }, sample.sizeC = c(100), sample.sizeT = c(100), n.rep = 2,
    formula = T ~ cont(Y), method.inference = "u-statistic", trace = 0,
    seed = 10)

    expect_equal(xx,yy)

    ## xx <- powerBuyseTest(sim = simFCT,
    ##                      sample.sizeC = c(100),
    ##                      sample.sizeT = c(100),
    ##                      n.rep = 10,
    ##                      cpus = 3,
    ##                      formula = T ~ cont(Y),
    ##                      method.inference = "u-statistic",
    ##                      trace = 4)

})

## * brice ozenne: 11/13/19 4:11 hierachical in BuyseTest
test_that("BuyseTest - hierarchical", {
    data(veteran, package = "survival")
    BT.nH <- BuyseTest(trt ~ tte(time, threshold = 20, status = "status") + cont(karno, threshold = 0),
                       hierarchical = FALSE, data = veteran, 
                       method.inference = "none", trace = 0)
    expect_equal(coef(BT.nH),
                 c("time_20" = -0.08765836, "karno_1e-12" = -0.11898828),
                 tol = 1e-6)

})

## * graemeleehickey (issue #4 on Github): 6 october 2019 simBuyseTest
test_that("simBuyseTest - rate vs. scale", {
    scale <- 2

    args <- list(scale.T = scale, scale.Censoring.T = scale+1,
                 scale.C = scale, scale.Censoring.C = scale+1,
                 scale.CR =  scale)
    set.seed(10)
    test <- simBuyseTest(1e4, argsBin = NULL, argsCont = NULL, argsTTE = args,
                       latent = TRUE)
    
    expect_equal(scale,mean(test[treatment == "C", mean(eventtimeUncensored)]), tol = 1e-2)
    expect_equal(scale,mean(test[treatment == "T", mean(eventtimeUncensored)]), tol = 1e-2)

    expect_equal(scale+1,mean(test[treatment == "C", mean(eventtimeCensoring)]), tol = 1e-1)
    expect_equal(scale+1,mean(test[treatment == "T", mean(eventtimeCensoring)]), tol = 1e-1)
})

## * graemeleehickey (issue #6 on Github): 15 march 2020 powerBuyseTest

args <- list(scale.T = c((3:5) / 10), scale.Censoring.T = rep(1, 3))
simFCT <- function(n.C, n.T) {
  simBuyseTest(100, argsBin = NULL, argsCont = NULL, argsTTE = args)
}

test_that("powerBuyseTest - status vs. censoring", {
    valid <- powerBuyseTest(sim = simFCT, sample.size = c(100), n.rep = 2,
                            formula = treatment ~ tte(eventtime1, status = status1),
                            method.inference = "u-statistic",
                            scoring.rule = "Gehan", trace = 0)

    expect_error(powerBuyseTest(sim = simFCT, sample.size = c(100), n.rep = 2,
                                formula = treatment ~ tte(eventtime1, censoring = status1),
                                method.inference = "u-statistic",
                                scoring.rule = "Gehan"))

    valid <- capture.output(powerBuyseTest(sim = simFCT, sample.size = c(100), n.rep = 2,
                                           formula = treatment ~ tte(eventtime1, status = status1),
                                           method.inference = "u-statistic",
                                           scoring.rule = "Gehan", trace = 4))
})

## * brice ozenne : 04/26/20 2:36 uninformative pairs Peron
dt.prodlim <- rbind(data.table(treat=0,
                               time = c(1:8,rep(9,12)),
                               status = c(rep(1,8),rep(0,12))
                               ),
                    data.table(treat=1,
                               time = c(1:8,rep(9,12)),
                               status = c(0,rep(1,7),rep(0,12))
                               ))

e.prodlim <- prodlim(Hist(time, status) ~ treat, data = dt.prodlim)
## plot(e.prodlim)

dt.sim <- data.table(treat = c(0:1), time = 8, status = 0)

test_that("uniformative pair after last observation",{
    ## warning because only uninformative
    expect_warning(e.BP <- BuyseTest(treat ~ tte(time, status, threshold=2),
                                     model.tte = e.prodlim, data = dt.sim, method.inference = "none"))
    expect_equal(as.double(e.BP@count.neutral), 0)
    expect_equal(as.double(e.BP@count.uninf), 1)
})



## * brice ozenne : 10/08/20 3:26  last time tie (event/censor)
## butils::object2script(mydata, digit = 3)
dt <- data.table("bras" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                 "OS" = c(0.427, 1.708, 2.004, 2.792, 3.088, 3.384, 3.417, 3.647, 3.778, 3.844, 5.092, 5.355, 5.453, 6.012, 6.209, 6.209, 6.307, 6.702, 7.786, 8.049, 8.739, 9.461, 11.367, 11.728, 11.925, 11.991, 12.648, 12.746, 13.042, 13.338, 13.436, 13.666, 13.798, 16.097, 16.097, 0.854, 1.84, 3.055, 3.515, 4.172, 5.059, 5.158, 5.223, 5.519, 5.585, 6.307, 6.34, 6.373, 6.767, 6.899, 6.965, 7.129, 7.589, 7.589, 7.589, 7.753, 8.18, 9.133, 9.198, 9.855, 10.315, 11.498, 13.141, 13.239, 13.305, 13.568, 13.929, 15.21, 16.459, 19.087, 20.237, 20.532, 21.846, 22.273, 26.445, 27.989), 
                 "etat" = c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1))

test_that("last time is a tie with both event and censor",{
    test <- BuyseTest(bras ~ tte(OS, status = etat),
                      data = dt, method.inference = "u-statistic", scoring.rule = "Peron",
                      trace = 0)
    expect_equal(as.double(c(coef(test,"count.favorable"),coef(test,"count.unfavorable"),coef(test,"count.neutral"))),
                 c(892.6111, 520.3092,   0.0000 ), tol = 1e-3)

    ## dt[c(1,36)]
                 
})

## * brice ozenne : 10/12/20 9:46 only censored event in one group
test_that("one group with only censoring, one group with no censoring",{
    dt <- data.table("treatment" = c(rep("C",10),rep("T",10)),
                     "time" = c(1:10,1:10),
                     "status" = c(rep(1,10),rep(0,10)))

    e.Peron <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0),
                         data = dt, scoring.rule = "Peron")
    expect_equal(as.double(coef(e.Peron,"netBenefit")),1)

    dt2 <- data.table("treatment" = c(rep("C",10),rep("T",10)),
                      "time" = c(1:10,1:10),
                      "status" = c(c(rep(1,9),0),rep(0,10)))

    e2.Peron <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0),
                          data = dt2, scoring.rule = "Peron")
    expect_equal(as.double(coef(e2.Peron,"netBenefit")),0.9)

    dt3 <- data.table("treatment" = c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T"), 
                     "time" = c(0.302, 0.307, 0.336, 0.347, 0.348, 0.459, 0.494, 0.525, 0.587, 0.588, 0.098, 0.116, 0.180, 0.229, 0.306, 0.318, 0.452, 0.485, 1.025, 1.339), 
                     "status" = c(0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    e3.Peron <- BuyseTest(treatment ~ tte(time, status = status, threshold = 0),
                          data = dt3, scoring.rule = "Peron")
    expect_equal(as.double(coef(e3.Peron,"netBenefit")),0.733333333)
})
