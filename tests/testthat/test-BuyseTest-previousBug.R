### test-BuyseTest-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2018 (16:46) 
## Version: 
## Last-Updated: maj 26 2018 (17:44) 
##           By: Brice Ozenne
##     Update #: 31
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
                  keep.comparison = TRUE,
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
                      neutral.as.uninf = TRUE)
    BTS.T <- as.data.table(summary(BT.T, print = FALSE)$table)
    
    BT.F <- BuyseTest(ttt~TTE(timeOS,threshold=0,censoring=eventOS) + cont(Mgrade.tox,threshold=0),
                      data = dt.sim,
                      neutral.as.uninf = FALSE)
    BTS.F <- as.data.table(summary(BT.F, print = FALSE)$table)

    expect_equal(BTS.T[1,],BTS.F[1,])
    expect_equal(BTS.T[endpoint == "Mgrade.tox" & strata == "global", pc.total],
                 BTS.T[endpoint == "timeOS" & strata == "global", pc.neutral+pc.uninf])
    expect_equal(BTS.F[endpoint == "Mgrade.tox" & strata == "global", pc.total],
                 BTS.F[endpoint == "timeOS" & strata == "global", pc.uninf])

    GS <- data.table("endpoint" = c("timeOS", "timeOS", "Mgrade.tox", "Mgrade.tox"), 
                     "threshold" = c(1e-12, 1e-12, 1e-12, 1e-12), 
                     "strata" = c("global", "1", "global", "1"), 
                     "pc.total" = c(100.00000, 100.00000,  33.33333,  33.33333), 
                     "pc.favorable" = c(50.00000, 50.00000, 16.66667, 16.66667), 
                     "pc.unfavorable" = c(16.66667, 16.66667, 11.11111, 11.11111), 
                     "pc.neutral" = c(11.111111, 11.111111,  5.555556,  5.555556), 
                     "pc.uninf" = c(22.22222, 22.22222,  0.00000,  0.00000), 
                     "delta" = c(0.33333333, 0.33333333, 0.05555556, 0.05555556), 
                     "Delta" = c(0.3333333, NA, 0.3888889, NA), 
                     "CIinf.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                     "CIsup.Delta" = as.numeric(c(NA, NA, NA, NA)), 
                     "p.value" = as.numeric(c(NA, NA, NA, NA)),
                     "n.resampling" = as.numeric(c(NA, NA, NA, NA)))
    ##    butils::object2script(BTS.T)

    attr(BTS.T,"index") <- NULL
    expect_equal(BTS.T, GS, tol = 1e-6)
    ## class(BTS.T[["n.resampling"]])
    ## class(GS[["n.resampling"]])
    
})

## * Joris: jeudi 5 avril 2018 à 14:57
if(FALSE){
    system.time(
        BT_Gehan <- BuyseTest(treatment = "trt",
                              strata = "celltype",
                              type = "TTE",
                              endpoint = "time",
                              threshold = 0,
                              censoring = "status", 
                              data=veteran, method.tte = "Gehan",
                              n.resampling = 1e4)
    )
    1
    confint(BT_Gehan)
    dt.sim <- data.table(
        ttt = c(0,0,0,0,0,1,1,1,1,1),
        y = c(1,2,3,4,5,2,3,4,5,6),
        end = c(1,1,1,1,1,1,1,1,1,1),
        time = c(10,11,12,13,14,11,12,13,14,15)
    )

    BT <- BuyseTest(ttt ~ cont(y,threshold=4) + TTE(time,censoring=end) + cont(y,threshold=3),
                    data = dt.sim, method.inference = "permutation", method.tte = "Peron", n.resampling = 1000)
    summary(BT, statistic = "netChance")
 ##        Generalized pairwise comparison with 3 prioritized endpoints

 ## > statistic       : net chance of a better outcome (delta: endpoint specific, Delta: global) 
 ## > null hypothesis : Delta == 0 
 ## > permutation test: 1000 samples, confidence level 0.95 
 ## > treatment groups: 0 (control) vs. 1 (treatment) 
 ## > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

 ## > results
 ## endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
 ##        y         4   100        12           0      88     0  0.12  0.12   [-0.04;0.281]   0.184 
 ##     time     1e-12    88        48          24      16     0  0.24  0.36    [-0.32;1.04]   0.297 
 ##        y         3    16         0           0      16     0  0.00  0.36    [-0.32;1.04]   0.297 

  xx <- summary(BT,statistic="winRatio")
## Generalized pairwise comparison with 3 prioritized endpoints

##  > statistic       : win ratio (delta: endpoint specific, Delta: global) 
##  > null hypothesis : Delta == 1 
##  > permutation test: 1000 samples, confidence level 0.95 
##  > treatment groups: 0 (control) vs. 1 (treatment) 
##  > censored pairs  : imputation using Kaplan Meier stratified by treatment group 

##  > results
##  endpoint threshold total favorable unfavorable neutral uninf delta Delta CI [2.5 ; 97.5] p.value 
##         y         4   100        12           0      88     0   Inf   Inf                         
##      time     1e-12    88        48          24      16     0     2   2.5    [1.65;8.167]   0.118 
##         y         3    16         0           0      16     0         2.5    [1.65;8.167]   0.118
}

## ** table Comparison
## $y_4
##    strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1       1       6       1                   1                   1         0           0       1             0
## 2       1       6       2                   1                   2         0           0       1             0
## 3       1       6       3                   1                   3         0           0       1             0
## 4       1       6       4                   1                   4         0           0       1             0
## 5       1       6       5                   1                   5         0           0       1             0
## 6       1       7       1                   2                   1         0           0       1             0
## 7       1       7       2                   2                   2         0           0       1             0
## 8       1       7       3                   2                   3         0           0       1             0
## 9       1       7       4                   2                   4         0           0       1             0
## 10      1       7       5                   2                   5         0           0       1             0
## 11      1       8       1                   3                   1         0           0       1             0
## 12      1       8       2                   3                   2         0           0       1             0
## 13      1       8       3                   3                   3         0           0       1             0
## 14      1       8       4                   3                   4         0           0       1             0
## 15      1       8       5                   3                   5         0           0       1             0
## 16      1       9       1                   4                   1         1           0       0             0
## 17      1       9       2                   4                   2         0           0       1             0
## 18      1       9       3                   4                   3         0           0       1             0
## 19      1       9       4                   4                   4         0           0       1             0
## 20      1       9       5                   4                   5         0           0       1             0
## 21      1      10       1                   5                   1         1           0       0             0
## 22      1      10       2                   5                   2         1           0       0             0
## 23      1      10       3                   5                   3         0           0       1             0
## 24      1      10       4                   5                   4         0           0       1             0
## 25      1      10       5                   5                   5         0           0       1             0

## $`time_1e-12`
##    strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1       1       6       1                   1                   1      1.00        0.00       0          0.00
## 2       1       6       2                   1                   2      1.00        0.00       0          0.00
## 3       1       6       3                   1                   3      0.75        0.00       0          0.25
## 4       1       6       4                   1                   4      0.50        0.25       0          0.25
## 5       1       6       5                   1                   5      0.25        0.50       0          0.25
## 6       1       7       1                   2                   1      1.00        0.00       0          0.00
## 7       1       7       2                   2                   2      1.00        0.00       0          0.00
## 8       1       7       3                   2                   3      0.00        0.00       1          0.00
## 9       1       7       4                   2                   4      0.00        1.00       0          0.00
## 10      1       7       5                   2                   5      0.00        1.00       0          0.00
## 11      1       8       1                   3                   1      1.00        0.00       0          0.00
## 12      1       8       2                   3                   2      1.00        0.00       0          0.00
## 13      1       8       3                   3                   3      1.00        0.00       0          0.00
## 14      1       8       4                   3                   4      0.00        0.00       1          0.00
## 15      1       8       5                   3                   5      0.00        1.00       0          0.00
## 16      1       9       2                   4                   2      1.00        0.00       0          0.00
## 17      1       9       3                   4                   3      1.00        0.00       0          0.00
## 18      1       9       4                   4                   4      1.00        0.00       0          0.00
## 19      1       9       5                   4                   5      0.00        0.00       1          0.00
## 20      1      10       3                   5                   3      1.00        0.00       0          0.00
## 21      1      10       4                   5                   4      1.00        0.00       0          0.00
## 22      1      10       5                   5                   5      1.00        0.00       0          0.00

## $y_3
##   strata index.1 index.0 indexWithinStrata.1 indexWithinStrata.0 favorable unfavorable neutral uninformative
## 1      1       7       3                   2                   3         0        0.00    1.00             0
## 2      1       8       4                   3                   4         0        0.00    1.00             0
## 3      1       9       5                   4                   5         0        0.00    1.00             0
## 4      1       6       3                   1                   3         0        0.00    0.25             0
## 5      1       6       4                   1                   4         0        0.00    0.25             0
## 6      1       6       5                   1                   5         0        0.25    0.00             0


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

## * running time (hidden)
## tps <- system.time(
##     BT.T <- .BuyseTest(alternative = "two.sided",
##                        censoring = NULL,
##                        correctionTTE = FALSE,
##                        cpus = 1,
##                        data = dt.sim,
##                        endpoint = c("timeOS","Mgrade.tox"),
##                        index.survivalM1 = numeric(0),
##                        keep.comparison = FALSE,
##                        method.tte = 0,
##                        method.inference = "none",
##                        neutral.as.uninf = TRUE,
##                        operator = c(">0",">0"),
##                        seed = 10,
##                        strata = NULL,
##                        threshold = c(1e-12,1e-12),
##                        threshold.TTEM1 = numeric(0),
##                        trace = 0,
##                        treatment = "ttt",
##                        type = c(1,1),
##                        Wscheme = matrix(nrow = 0, ncol = 0)
##                        )
## )    
## tps  

######################################################################
### test-BuyseTest-previousBug.R ends here


