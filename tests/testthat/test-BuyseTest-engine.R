### test-BuyseTest-engine.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2020 (09:46) 
## Version: 
## Last-Updated: okt  4 2021 (20:13) 
##           By: Brice Ozenne
##     Update #: 20
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

context("Check BuyseTest without strata")

## * Settings
n.patients <- c(60,65)
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = NULL,
                       argsCont = NULL,
                       argsTTE = list(scale.T = 1:2, scale.Censoring.T = rep(1,2)))

## * Compare
test_that("TTE with decreasing thresholds",{
    iFormula <- treatment ~ tte(eventtime1,status1,5) + tte(eventtime1,status1,1) + tte(eventtime1,status1,0.5) + tte(eventtime1,status1,0)

    BuyseTest.options(engine = "GPC_cpp")
    e.BT1 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")

    BuyseTest.options(engine = "GPC2_cpp")
    e.BT2 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")

    ## range(e.BT1@iidAverage$favorable - e.BT2@iidAverage$favorable)
    ## range(e.BT1@iidAverage$unfavorable - e.BT2@iidAverage$unfavorable)
    ## range(e.BT1@iidAverage$neutral - e.BT2@iidAverage$neutral)
    
    expect_equal(confint(e.BT1)[,"estimate"],confint(e.BT2)[,"estimate"], tol = 1e-6)
    expect_equal(confint(e.BT1)[,"se"],confint(e.BT2)[,"se"], tol = 1e-2)
    ## expected difference because GPC2 compute the influence function over all pairs,
    ## while GPC only over pairs with informative scores.
    ## the difference is expected to be small though in large samples

    GS <- matrix(c(0, 0.13766049, 0.07000947, 0.08471151, 0, 0.11214375, 0.12563388, 0.13305876, 0, -0.08529558, -0.17549094, -0.17591136, 0, 0.34748762, 0.30731249, 0.33423682, 0, 0, 0, 0, 1, 0.22552443, 0.5786024, 0.52634371), 
                 nrow = 4, 
                 ncol = 6, 
                 dimnames = list(c("eventtime1_5", "eventtime1_1", "eventtime1_0.5", "eventtime1_1e-12"),c("estimate", "se", "lower.ci", "upper.ci", "null", "p.value")) 
                 ) 

    test <- confint(e.BT2)
    attr(test,"n.resampling") <- NULL
    attr(test,"iid") <- NULL
    attr(test,"transform") <- NULL
    attr(test,"backtransform") <- NULL
    expect_equal(as.data.frame(GS), test, tol = 1e-6)
})

test_that("different TTE with decreasing thresholds",{
    iFormula <- treatment ~ tte(eventtime1,status1,1) + tte(eventtime2,status2,1) + tte(eventtime1,status1,0.25) + tte(eventtime2,status2,0)

    BuyseTest.options(engine = "GPC_cpp")
    e.BT1 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")

    BuyseTest.options(engine = "GPC2_cpp")
    e.BT2 <- BuyseTest(iFormula, data = dt.sim,
                       method.inference = "u-statistic", scoring.rule = "Peron")
    expect_equal(confint(e.BT1)[,"estimate"],confint(e.BT2)[,"estimate"], tol = 1e-6)
    expect_equal(confint(e.BT1)[,"se"],confint(e.BT2)[,"se"], tol = 1e-2)
    ## expected difference because GPC2 compute the influence function over all pairs,
    ## while GPC only over pairs with informative scores.
    ## the difference is expected to be small though in large samples

    GS <- matrix(c(0.13766049, 0.29168224, 0.26716412, 0.26881354, 0.11214375, 0.14029539, 0.14634712, 0.14672071, -0.08529558, -0.00013978, -0.03506022, -0.03436758, 0.34748762, 0.53772357, 0.52461749, 0.52668769, 0, 0, 0, 0, 0.22552443, 0.05010665, 0.08231598, 0.08140938), 
                 nrow = 4, 
                 ncol = 6, 
                 dimnames = list(c("eventtime1_1", "eventtime2_1", "eventtime1_0.25", "eventtime2_1e-12"),c("estimate", "se", "lower.ci", "upper.ci", "null", "p.value")) 
                 ) 
    test <- confint(e.BT2)
    attr(test,"n.resampling") <- NULL
    attr(test,"iid") <- NULL
    attr(test,"transform") <- NULL
    attr(test,"backtransform") <- NULL
    expect_equal(as.data.frame(GS), test, tol = 1e-3)
})


## * dataset [save]
## dt.sim <- data.table("treatment" = c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T"), 
##            "eventtime1" = c(0.79577437, 0.29261394, 0.29295204, 0.01077697, 0.32989749, 0.46464367, 0.95386903, 0.64348019, 0.16468683, 0.96397945, 0.15145199, 0.01248151, 0.19325083, 0.66066663, 0.37595327, 0.58180543, 0.56149305, 0.63643639, 1.34491304, 0.09811647, 0.15531358, 0.24577818, 0.91206455, 0.08622315, 0.00437068, 0.07262749, 1.89017379, 0.32409615, 0.99894561, 0.43884796, 0.41673022, 0.01915872, 0.04266630, 0.00129060, 0.08802724, 0.15854922, 0.19085512, 0.13225449, 0.17847595, 0.00011246, 0.14967014, 0.04105886, 0.60220956, 0.10690356, 0.44929205, 0.06879102, 0.29051506, 0.32497984, 0.60555591, 0.00485320, 0.92640116, 0.24089830, 1.01894511, 0.59803620, 0.43285432, 0.46751140, 0.19660899, 0.31900979, 0.00757993, 0.25049825, 0.06965880, 0.66302185, 0.22487687, 1.58865687, 0.69883175, 1.08300892, 1.09701618, 0.85784765, 0.57017628, 0.63918564, 0.46688273, 0.13235940, 0.38690290, 0.15446529, 0.05053404, 0.11874805, 2.62835120, 0.26261163, 0.65918371, 0.20573156, 0.80062460, 0.33799949, 0.14834002, 0.11642816, 0.33428564, 0.43329451, 0.62381704, 0.41169080, 0.08002308, 1.90258085, 0.19554009, 0.27228695, 0.00136648, 1.33447157, 0.10077969, 0.12835260, 0.39954960, 0.29344588, 0.88313950, 1.09824202, 0.26829154, 0.16479848, 0.63080594, 0.06606621, 2.08985421, 0.32385948, 0.37545326, 0.13866077, 0.26051812, 0.27175045, 1.35024580, 0.27024855, 0.01775182, 0.19560968, 0.52775497, 0.57430011, 0.27207542, 0.16156373, 0.34555520, 0.01327954, 0.79489001, 0.11052999, 0.64695253, 0.18732135, 0.30728657), 
##            "status1" = c(1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0), 
##            "eventtime2" = c(0.72615286, 0.30990753, 0.19734068, 1.45015313, 0.01678544, 0.05411686, 0.44599463, 1.31406643, 0.11086857, 0.35896983, 0.04642124, 0.01738564, 0.12670825, 0.35440325, 0.79098364, 0.08579868, 0.00819823, 0.10740522, 0.50854831, 0.01129737, 0.31479059, 0.14026615, 0.37035992, 0.01556742, 0.07427124, 0.40681971, 0.78088914, 0.12198793, 2.14035911, 0.16382409, 0.34727401, 0.29651620, 0.47276796, 0.02235366, 0.36557684, 0.14715453, 1.35100756, 0.26630413, 0.00338837, 0.98729447, 0.75771650, 0.23588531, 0.35851747, 0.74237587, 2.21618578, 2.94382345, 0.23658775, 0.22663702, 0.12846880, 0.49238925, 0.23567446, 0.01570052, 1.47169237, 1.09048580, 0.81635260, 0.69766535, 1.00493089, 0.42227316, 0.68362323, 1.02917180, 0.22932350, 1.09688596, 0.75255164, 0.51641702, 0.59559315, 1.05595486, 1.20180822, 0.91293408, 0.53043209, 0.13241632, 0.34454510, 1.64199075, 0.22009462, 0.61963308, 0.20214388, 3.42796265, 0.26680520, 0.31651855, 0.91888102, 0.13361218, 0.40457320, 2.05480419, 1.62592766, 2.86376850, 0.36863418, 0.15910733, 0.53935052, 0.61094415, 0.00220048, 0.35833974, 1.35170549, 1.41819754, 0.96531610, 1.10756136, 0.03641272, 1.38989406, 0.00927533, 0.18481723, 1.58930474, 0.53308180, 0.01821104, 0.55704295, 0.29852972, 0.37228051, 0.07973612, 2.05818451, 0.00878831, 1.08116307, 0.14750010, 0.58478828, 0.05780919, 0.14660226, 0.16776685, 0.87629874, 0.21061788, 0.10621545, 0.02334416, 0.52714246, 0.01245405, 0.97092680, 0.73916365, 0.18460891, 0.21003983, 0.34469510, 0.31734420), 
##            "status2" = c(1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0))
##----------------------------------------------------------------------
### test-BuyseTest-engine.R ends here
