### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: sep 10 2018 (09:36) 
##           By: Brice Ozenne
##     Update #: 45
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
BuyseTest.options(check = FALSE,
                  keep.individualScore = TRUE,
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
    BT.T <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = TRUE, method.tte = "Gehan")
    BTS.T <- as.data.table(summary(BT.T, print = FALSE, percentage = FALSE)$table)
    
    BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = FALSE, method.tte = "Gehan")
    BTS.F <- as.data.table(summary(BT.F, print = FALSE, percentage = FALSE)$table)


    
    expect_equal(BTS.T[1,],BTS.F[1,])

    expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf],
                 BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total])
    expect_equal(BTS.T[endpoint == "timeOS" & strata == "global", n.neutral+n.uninf],
                 BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total])
    expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.total],
                 BTS.T[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf])

    expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf],
                 BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total])

    expect_equal(BTS.F[endpoint == "timeOS" & strata == "global", n.uninf],
                 BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total])

    expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.total],
                 BTS.F[endpoint == "Mgrade.tox" & strata == "global", n.favorable+n.unfavorable+n.neutral+n.uninf])

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

## * Joris: jeudi 31 mai 2018
## Beware of a bug in this new CRAN version (1.3.2)
## pvalues computed through bootstrap and stratified bootstrap seem to be reversed, while associated confidence intervals are correct
if(FALSE){
    dt <- data.table(
        ttt = rep(c(0,0,0,0,0,1,1,1,1,1),10),
        y1 = rep(c(1,2,3,4,5,0,3,4,5,6),10),
        strat = rep(1:5,each=20)
    )

    BT.boot <- BuyseTest(data = dt, formula = ttt ~ cont(y1, threshold=1),
                         method.inference = "bootstrap",
                         n.resampling = 500)

    summary(BT.boot)

    BT.perm <- BuyseTest(data = dt, formula = ttt ~ cont(y1, threshold=1),
                         method.inference = "permutation",
                         n.resampling = 500)
    summary(BT.perm)
}


## * Brice: 09/06/18 6:51 (Tied event with tte endpoint)
## when computing the integral for peron with double censoring
## the ordering of the data modified the ouput
## this has been correct with version 1.4
data(veteran,package="survival")

test_that("ordering of tied event does not affect BuyseTest", {
    ## veteran2[veteran2$time==100,]
    BT.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                        data = veteran, method.inference = "none")

    veteran1 <- veteran[order(veteran$time,veteran$status),c("time","status","trt")]
    ## veteran1[veteran2$time==100,]
    BT1.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                         data = veteran1, method.inference = "none")

    veteran2 <- veteran[order(veteran$time,-veteran$status),c("time","status","trt")]
    ## ## veteran2[veteran2$time==100,]
    BT2.all <- BuyseTest(trt ~ tte(time, threshold = 0, censoring = "status"),
                         data = veteran2, method.inference = "none")

    expect_equal(as.double(BT.all@Delta.winRatio), 0.8384072, tol = 1e-5)
    expect_equal(BT.all@Delta.winRatio, BT1.all@Delta.winRatio)
    expect_equal(BT.all@Delta.winRatio, BT2.all@Delta.winRatio)
})
