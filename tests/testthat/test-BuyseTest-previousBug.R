### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: okt 30 2018 (16:37) 
##           By: Brice Ozenne
##     Update #: 84
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

    for(iCorrection in c(FALSE,TRUE)){ ## iCorrection <- TRUE
        BT.T <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                          data = dt.sim,
                          neutral.as.uninf = TRUE, method.tte = "Gehan", correction.uninf = iCorrection)
        BTS.T <- as.data.table(summary(BT.T, print = FALSE, percentage = FALSE)$table)

        BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                          data = dt.sim,
                          neutral.as.uninf = FALSE, method.tte = "Gehan", correction.uninf = iCorrection)
        BTS.F <- as.data.table(summary(BT.F, print = FALSE, percentage = FALSE)$table)

        ## neutral.as.uninf does not impact the results for first endpoint
        expect_equal(BTS.T[1,],BTS.F[1,])

        ## check consistency of the number of pairs
        ## neutral.as.uninf = TRUE
        ## summary(BT.T)
        expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total])
        expect_equal(BTS.T[endpoint == "timeOS" & strata == "global", n.neutral+n.uninf],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total])
        expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total],
                     BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf])

        ## neutral.as.uninf = FALSE
        expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total])
        expect_equal(BTS.F[endpoint == "timeOS" & strata == "global", n.uninf],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total])
        expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total],
                     BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf])

        ## compared to known value
        if(iCorrection == FALSE){
            test <- as.data.table(summary(BT.T, print = FALSE)$table)
            GS <- data.table("endpoint" = c("timeOS", "timeOS", "Mgrade.tox", "Mgrade.tox"), 
                             "threshold" = c(1e-12, 1e-12, 1e-12, 1e-12), 
                             "strata" = c("global", "1", "global", "1"), 
                             "pc.total" = c(100.00000, 100.00000,  44.44444,  44.44444), 
                             "pc.favorable" = c(44.44444, 44.44444, 22.22222, 22.22222), 
                             "pc.unfavorable" = c(11.11111, 11.11111, 11.11111, 11.11111), 
                             "pc.neutral" = c(11.11111, 11.11111, 11.11111, 11.11111), 
                             "pc.uninf" = c(33.33333, 33.33333,  0.00000,  0.00000), 
                             "delta" = c(0.3333333, 0.3333333, 0.1111111, 0.1111111), 
                             "Delta" = c(0.3333333, NA, 0.4444444, NA), 
                             "CIinf.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                             "CIsup.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                             "p.value" = as.numeric(c(NA, NA, NA, NA)), 
                             "n.resampling" = as.numeric(c(NA, NA, NA, NA)))
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
                   "event" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

BT_tau0 <- BuyseTest(data=data,
                     treatment="trt",
                     endpoint="time",
                     type="timeToEvent",
                     threshold=as.numeric(0),
                     censoring="event",
                     method.tte="Peron",
                     method.inference = "none",
                     cpus=1)

## * Brice: 09/06/18 6:51 (Tied event with tte endpoint)
## when computing the integral for peron with double censoring
## the ordering of the data modified the ouput
## this has been correct with version 1.4
data(veteran,package="survival")

test_that("ordering of tied event does not affect BuyseTest", {
    ## veteran2[veteran2$time==100,]
    BT.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                        data = veteran, method.tte = "Peron", method.inference = "none", correction.uninf = FALSE)

    veteran1 <- veteran[order(veteran$time,veteran$status),c("time","status","trt")]
    ## veteran1[veteran2$time==100,]
    BT1.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                         data = veteran1, method.tte = "Peron", method.inference = "none", correction.uninf = FALSE)

    veteran2 <- veteran[order(veteran$time,-veteran$status),c("time","status","trt")]
    ## ## veteran2[veteran2$time==100,]
    BT2.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                         data = veteran2, method.tte = "Peron", method.inference = "none", correction.uninf = FALSE)

    ## effect of the ordering
    expect_equal(BT.all@Delta.winRatio, BT1.all@Delta.winRatio)
    expect_equal(BT.all@Delta.winRatio, BT2.all@Delta.winRatio)

    ## number of pairs
    expect_equal(BT.all@n.pairs, prod(table(veteran$trt)), tol = 1e-5)
    expect_equal(BT1.all@n.pairs, prod(table(veteran$trt)), tol = 1e-5)
    expect_equal(BT2.all@n.pairs, prod(table(veteran$trt)), tol = 1e-5)

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
    expect_equal(as.double(BT.all@count.uninf), 0.0, tol = 1e-12)
    
    ## result
    expect_equal(as.double(BT.all@Delta.winRatio), 0.8384569, tol = 1e-5)
    

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

tab<-data.frame(group,Time,Event)

test_that("Multiple thresholds",{
    BuyseresPer <- BuyseTest(data=tab,
                             endpoint=c("Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time","Time"),
                             treatment="group",
                             type=c("TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE","TTE"),
                             censoring=c("Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event","Event"),
                             threshold=c(42,39,36,33,30,27,24,21,18,15,12,9,6,3,0),
                             n.resampling=500,
                             trace=0,
                             method.tte="Peron",
                             correction.uninf=F,
                             method.inference="none")

    resS <- as.data.table(summary(BuyseresPer, print = FALSE)$table)

    ## pairs are correctly transfered from one endpoint to another
    expect_equal(resS[strata == "global" & threshold > tail(threshold,1), pc.neutral +  pc.uninf],
                 resS[strata == "global" & threshold < threshold[1], pc.total], tol = 1e-2)

    ## butils::object2script(as.double(BuyseresPer@count.favorable), digit = 2)
    GS <- c(260.64, 35.93, 37.33, 147.32, 272.14, 263.6, 235.7, 213.21, 390.29, 408.73, 514.7, 514.34, 744.78, 865.21, 1095.26)
    expect_equal(as.double(BuyseresPer@count.favorable), GS, tol = 1e-1)
    ## butils::object2script(as.double(BuyseresPer@count.unfavorable), digit = 2)
    GS <- c(0, 0, 6.97, 25.66, 43.89, 34.8, 46.38, 105.42, 199.85, 338.55, 407.72, 521.83, 548.02, 782.94, 938.8)
    expect_equal(as.double(BuyseresPer@count.unfavorable), GS, tol = 1e-1)
    ## butils::object2script(as.double(BuyseresPer@count.neutral), digit = 2)
    GS <- c(9617.63, 9611.09, 9596.16, 9448.36, 9149.12, 8863.32, 8581.24, 8262.61, 7676.67, 6933.58, 6011.17, 4975, 3682.21, 2034.06, 0)
    expect_equal(as.double(BuyseresPer@count.neutral), GS, tol = 1e-1)
    ## butils::object2script(as.double(BuyseresPer@count.uninf), digit = 2)
    GS <- c(121.73, 92.35, 62.96, 37.78, 20.99, 8.4, 8.4, 8.4, 4.2, 0, 0, 0, 0, 0, 0)
    expect_equal(as.double(BuyseresPer@count.uninf), GS, tol = 1e-1)
    ## butils::object2script(as.double(BuyseresPer@delta.netBenefit), digit = 5)
    GS <- c(0.02606, 0.00359, 0.00304, 0.01217, 0.02282, 0.02288, 0.01893, 0.01078, 0.01904, 0.00702, 0.0107, -0.00075, 0.01968, 0.00823, 0.01565)
    expect_equal(as.double(BuyseresPer@delta.netBenefit), GS, tol = 1e-3)
    ## butils::object2script(as.double(BuyseresPer@delta.winRatio), digit = 5)
    GS <- c(Inf, Inf, 5.35344, 5.74093, 6.19986, 7.57457, 5.08161, 2.02241, 1.95291, 1.2073, 1.2624, 0.98564, 1.35904, 1.10508, 1.16666)
    expect_equal(as.double(BuyseresPer@delta.winRatio), GS, tol = 1e-3)
})


## * Brice: 12/10/18 3:02 (Wscheme)

BuyseTest_buildWscheme <- BuyseTest:::buildWscheme
## BuyseTest_buildWscheme <- buildWscheme

endpoint <- c("time","time","time")
threshold <- c(3:1)
D <- length(endpoint)
type <- rep(3, D)
D.TTE <- sum(type==3)

test_that("Wscheme: 3 times the same endpoint",{
    Wtest <- BuyseTest_buildWscheme(method.tte = 1,
                                    endpoint = endpoint,
                                    D.TTE = D.TTE,
                                    D = D,
                                    type = type,
                                    n.strata = 1,
                                    threshold = threshold)

    ## butils::object2script(Wtest$Wscheme)
    GS <- list(Wscheme = matrix(c(0, NA, NA, 0, 0, NA, 0, 0, 0),
                                nrow = 3, ncol = 3,
                                dimnames = list(c("weigth of time(3)", "weigth of time(2)", "weigth of time(1)"),
                                                c("for time(3)", "for time(2)", "for time(1)")) ),
               endpoint.UTTE = "time",
               index.UTTE = c(0,0,0),
               D.UTTE = 1,
               reanalyzed = c(TRUE, TRUE, FALSE) )

    expect_equal(Wtest[c("Wscheme","endpoint.UTTE","index.UTTE","D.UTTE","reanalyzed")], GS)
})

endpoint <- c("time","time1","time","time","time2","time1")
threshold <- c(6:1)
D <- length(endpoint)
type <- rep(3, D)
D.TTE <- sum(type==3)

test_that("Wscheme: 6 tte endpoint",{
    Wtest <- BuyseTest_buildWscheme(method.tte = 1,
                                    endpoint = endpoint,
                                    D.TTE = D.TTE,
                                    D = D,
                                    type = type,
                                    n.strata = 1,
                                    threshold = threshold)

    ## butils::object2script(Wtest)
    GS <- list(Wscheme = matrix(c(0, NA, NA, NA, NA, NA, 1, 0, NA, NA, NA, NA, 0, 1, 0, NA, NA, NA, 0, 1, 0, 0, NA, NA, 1, 1, 1, 1, 0, NA, 1, 0, 1, 1, 1, 0), 
                                nrow = 6, 
                                ncol = 6, 
                                dimnames = list(c("weigth of time(6)", "weigth of time1(5)", "weigth of time(4)", "weigth of time(3)", "weigth of time2(2)", "weigth of time1(1)"),
                                                c("for time(6)", "for time1(5)", "for time(4)", "for time(3)", "for time2(2)", "for time1(1)")) 
                                ),
               endpoint.UTTE = c("time","time1","time2"),
               index.UTTE = c(0,1,0,0,2,1),
               D.UTTE = 3,
               reanalyzed = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) )

    expect_equal(Wtest[c("Wscheme","endpoint.UTTE","index.UTTE","D.UTTE","reanalyzed")], GS)
})

endpoint <- c("time","bin","bin","time","bin","time")
threshold <- c(6:1)
D <- length(endpoint)
type <- 1+(endpoint=="time")*2
D.TTE <- sum(type==3)

test_that("Wscheme: 6 mixed endpoint",{
    Wtest <- BuyseTest_buildWscheme(method.tte = 1,
                                    endpoint = endpoint,
                                    D.TTE = D.TTE,
                                    D = D,
                                    type = type,
                                    n.strata = 1,
                                    threshold = threshold)

    ## butils::object2script(Wtest)
    GS <- list(Wscheme = matrix(c(0, NA, NA, NA, NA, NA, 1, 0, NA, NA, NA, NA, 1, 1, 0, NA, NA, NA, 0, 1, 1, 0, NA, NA, 1, 1, 1, 1, 0, NA, 0, 1, 1, 0, 1, 0), 
                                nrow = 6, 
                                ncol = 6, 
                                dimnames = list(c("weigth of time(6)", "weigth of bin(5)", "weigth of bin(4)", "weigth of time(3)", "weigth of bin(2)", "weigth of time(1)"),
                                                c("for time(6)", "for bin(5)", "for bin(4)", "for time(3)", "for bin(2)", "for time(1)")) 
                                ),
               endpoint.UTTE = "time",
               index.UTTE = c(0,-1,-1,0,-1,0),
               D.UTTE = 1,
               reanalyzed = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) )

    expect_equal(Wtest[c("Wscheme","endpoint.UTTE","index.UTTE","D.UTTE","reanalyzed")], GS)
})


## * Brice: 30/10/18 4:36 Neutral pairs with 0 threshold
df <- data.frame("survie" = c(2.1, 4.1, 6.1, 8.1, 4, 6, 8, 10),
                 "event" = c(1, 1, 1, 0, 1, 0, 0, 1),
                 "group" = c(0, 0, 0, 0, 1, 1, 1, 1),
                 "score" = 1)

test_that("1 TTE endpoint - Gehan (no correction)", {
    Peron <- BuyseTest(group ~ tte(survie, censoring = event, threshold = 0),
                       data = df, 
                       method.tte = "Peron", correction.uninf = FALSE)

    expect_equal(as.double(Peron@count.neutral),0) ## should not be any neutral pair with a threshold of 0
})
