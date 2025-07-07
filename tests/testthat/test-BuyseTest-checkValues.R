if(FALSE){
    library(mvtnorm)
    library(testthat)
    library(BuyseTest)
    library(data.table)
}

context("Check BuyseTest without strata")

## * Settings
n.patients <- c(90,100)
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  pool.strata = "Buyse",
                  trace = 0)

## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = list(c(0.5,0.5),c(0.25,0.75))),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(scale.T = 1:3, scale.censoring.T = rep(1,3)))
## butils::object2script(dt.sim)
dt.sim$status1.noC <- 1

dtS.sim <- rbind(cbind(dt.sim, strata = 1),
                 cbind(dt.sim, strata = 2),
                 cbind(dt.sim, strata = 3))


## * Binary endpoint
## ** No strata
test_that("BuyseTest - binary (no strata)", {
    BT.bin <- BuyseTest(treatment ~ bin(toxicity1),
                        data = dt.sim)

    BT2 <- BuyseTest(data = dt.sim,
                     endpoint = "toxicity1",
                     treatment = "treatment",
                     type = "bin")
    
    ## *** test against fixed value
    test <- list(favorable = as.double(coef(BT.bin, statistic = "count.favorable", cumulative = FALSE)),
                 unfavorable = as.double(coef(BT.bin, statistic = "count.unfavorable", cumulative = FALSE)),
                 neutral = as.double(coef(BT.bin, statistic = "count.neutral", cumulative = FALSE)),
                 uninf = as.double(coef(BT.bin, statistic = "count.uninf", cumulative = FALSE)),
                 favorable = as.double(coef(BT.bin, statistic = "favorable", cumulative = TRUE)),
                 unfavorable = as.double(coef(BT.bin, statistic = "unfavorable", cumulative = TRUE)),
                 netChange = as.double(coef(BT.bin, statistic = "netBenefit", cumulative = TRUE)),
                 winRatio = as.double(coef(BT.bin, statistic = "winRatio", cumulative = TRUE))
                 )

    GS <- list(favorable = c(1968) ,
               unfavorable = c(2478) ,
               neutral = c(4554) ,
               uninf = c(0) ,
               favorable = c(0.21866667) ,
               unfavorable = c(0.27533333) ,
               netChange = c(-0.05666667) ,
               winRatio = c(0.79418886) )
    ## butils::object2script(test, digit = 6)

    expect_equal(test, GS, tol = 1e-6, scale = 1)
    BT.bin@call <- list()
    BT2@call <- list()
    expect_equal(BT.bin,BT2)
    
    ## *** count pairs
    tableS <- model.tables(BT.bin, percentage = FALSE)
    expect_equal(tableS$total,
                 tableS$favorable + tableS$unfavorable + tableS$neutral + tableS$uninf)

})

## ** Strata
test_that("BuyseTest - binary (strata)", {

    BT.bin <- BuyseTest(treatment ~ bin(toxicity1) + strata,
                        data = dtS.sim)

    tableS <- model.tables(BT.bin, percentage = FALSE)
    dt.tableS <- as.data.table(tableS)
    
    ## *** count pairs
    expect_equal(dt.tableS[,total],
                 unname(dt.tableS[,favorable + unfavorable + neutral + uninf]
                 ))
    expect_equal(dt.tableS[,total], c(27000,9000,9000,9000))
    expect_equal(dt.tableS[,favorable], c(5904, 1968, 1968, 1968))
    expect_equal(dt.tableS[,unfavorable], c(7434, 2478, 2478, 2478))
    expect_equal(dt.tableS[,neutral], c(13662, 4554, 4554, 4554))
    expect_equal(dt.tableS[,uninf], c(0, 0, 0, 0))

    ## *** test summary statistic
    expect_equal(dt.tableS[,delta], c(-0.05666667, -0.05666667, -0.05666667, -0.05666667), tol = 1e-6)
    expect_equal(dt.tableS[,Delta], c(-0.05666667, NA, NA, NA), tol = 1e-6)
})

## * Continuous endpoint
## ** No strata
test_that("BuyseTest - continuous (no strata)", {
    BT.cont <- BuyseTest(treatment ~ cont(score1, 1) + cont(score2, 0),
                         data = dt.sim)
    
    BT2 <- BuyseTest(data = dt.sim,
                     endpoint = c("score1","score2"),
                     treatment = "treatment",
                     type = c("cont","cont"),
                     threshold = c(1,0)
                     )
    
    ## *** test against fixed value    
    test <- list(favorable = as.double(coef(BT.cont, statistic = "count.favorable", cumulative = FALSE)),
                 unfavorable = as.double(coef(BT.cont, statistic = "count.unfavorable", cumulative = FALSE)),
                 neutral = as.double(coef(BT.cont, statistic = "count.neutral", cumulative = FALSE)),
                 uninf = as.double(coef(BT.cont, statistic = "count.uninf", cumulative = FALSE)),
                 favorable = as.double(coef(BT.cont, statistic = "favorable", cumulative = TRUE)),
                 unfavorable = as.double(coef(BT.cont, statistic = "unfavorable", cumulative = TRUE)),
                 netChange = as.double(coef(BT.cont, statistic = "netBenefit", cumulative = TRUE)),
                 winRatio = as.double(coef(BT.cont, statistic = "winRatio", cumulative = TRUE))
                 )

    GS <- list(favorable = c(2196, 2142) ,
               unfavorable = c(2501, 2161) ,
               neutral = c(4303, 0) ,
               uninf = c(0, 0) ,
               favorable = c(0.244, 0.482) ,
               unfavorable = c(0.27788889, 0.518) ,
               netChange = c(-0.03388889, -0.036) ,
               winRatio = c(0.87804878, 0.93050193) )

    ## butils::object2script(test, digit = 6)
    BT.cont@call <- list()
    BT2@call <- list()
    expect_equal(test, GS, tol = 1e-6, scale = 1)
    expect_equal(BT.cont,BT2)

    ## *** count pairs
    tableS <- model.tables(BT.cont, percentage = FALSE)
    expect_equal(tableS$total,
                 tableS$favorable + tableS$unfavorable + tableS$neutral + tableS$uninf)
})

## ** Strata
test_that("BuyseTest - continuous (strata)", {

    BT.cont <- BuyseTest(treatment ~ cont(score1, 1) + cont(score2, 0) + strata,
                         data = dtS.sim)

    tableS <- model.tables(BT.cont, percentage = FALSE)
    dt.tableS <- as.data.table(tableS)

        ## *** count pairs
    expect_equal(dt.tableS[,total],
                 unname(dt.tableS[,favorable + unfavorable + neutral + uninf]
                 ))
    expect_equal(dt.tableS[,total], c(27000, 9000, 9000, 9000, 12909, 4303, 4303, 4303))
    expect_equal(dt.tableS[,favorable], c(6588, 2196, 2196, 2196, 6426, 2142, 2142, 2142))
    expect_equal(dt.tableS[,unfavorable], c(7503, 2501, 2501, 2501, 6483, 2161, 2161, 2161))
    expect_equal(dt.tableS[,neutral], c(12909, 4303, 4303, 4303, 0, 0, 0, 0))
    expect_equal(dt.tableS[,uninf], c(0, 0, 0, 0, 0, 0, 0, 0))

    ## *** test summary statistic
    expect_equal(dt.tableS[,delta], c(-0.03388889, -0.03388889, -0.03388889, -0.03388889, -0.00211111, -0.00211111, -0.00211111, -0.00211111), tol = 1e-6)
    expect_equal(dt.tableS[,Delta], c(-0.03388889, NA, NA, NA, -0.036, NA, NA, NA), tol = 1e-6)
})


## * Time to event endpoint
## ** No strata - same endpoint
for(method in c("Gehan","Peron")){ ## method <- "Gehan" ## method <- "Peron"
    test_that(paste0("BuyseTest - tte (same, ",method,", no strata)"),{ 

        BT.tte <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + tte(eventtime1, status1, threshold = 0.5) + tte(eventtime1, status1, threshold = 0.25),
                            data = dt.sim,
                            scoring.rule = method,
                            correction.uninf = FALSE
                            )
        
        BT.1tte <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0.25),
                            data = dt.sim,
                            scoring.rule = method,
                            correction.uninf = FALSE
                            )
        
        BT2 <- BuyseTest(data = dt.sim,
                         endpoint = c("eventtime1","eventtime1","eventtime1"),
                         status = c("status1","status1","status1"),
                         treatment = "treatment",
                         type = c("tte","tte","tte"),
                         threshold = c(1,0.5,0.25),
                         scoring.rule = method,
                         correction.uninf = FALSE
                         )

        ## *** compatibility between BuyseTests
        BT.tte@call <- list()
        BT2@call <- list()
        expect_equal(BT.tte, BT2)
        expect_equivalent(sum(coef(BT.tte, statistic = "count.favorable", cumulative = FALSE)),
                          coef(BT.1tte, statistic = "count.favorable", cumulative = FALSE))
        expect_equivalent(sum(coef(BT.tte, statistic = "count.unfavorable", cumulative = FALSE)),
                          coef(BT.1tte, statistic = "count.unfavorable", cumulative = FALSE))
        expect_equivalent(unname(coef(BT.tte, statistic = "count.neutral", cumulative = FALSE)[3]),
                          coef(BT.1tte, statistic = "count.neutral", cumulative = FALSE))
        expect_equivalent(unname(coef(BT.tte, statistic = "count.uninf", cumulative = FALSE)[3]),
                          coef(BT.1tte, statistic = "count.uninf", cumulative = FALSE))
        expect_equivalent(unname(coef(BT.tte, statistic = "netBenefit", cumulative = TRUE)[3]),
                          coef(BT.1tte, statistic = "netBenefit", cumulative = TRUE))
        expect_equivalent(unname(coef(BT.tte, statistic = "winRatio", cumulative = TRUE)[3]),
                          coef(BT.1tte, statistic = "winRatio", cumulative = TRUE))

        ## *** test against fixed value
        test <- list(favorable = as.double(coef(BT.tte, statistic = "count.favorable", cumulative = FALSE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "count.unfavorable", cumulative = FALSE)),
                     neutral = as.double(coef(BT.tte, statistic = "count.neutral", cumulative = FALSE)),
                     uninf = as.double(coef(BT.tte, statistic = "count.uninf", cumulative = FALSE)),
                     favorable = as.double(coef(BT.tte, statistic = "favorable", cumulative = TRUE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "unfavorable", cumulative = TRUE)),
                     netChange = as.double(coef(BT.tte, statistic = "netBenefit", cumulative = TRUE)),
                     winRatio = as.double(coef(BT.tte, statistic = "winRatio", cumulative = TRUE))
                     )

        if(method == "Gehan"){
            GS <- list(favorable = c(438, 719, 543) ,
                       unfavorable = c(325, 582, 500) ,
                       neutral = c(2284, 1569, 1084) ,
                       uninf = c(5953, 5367, 4809) ,
                       favorable = c(0.04866667, 0.12855556, 0.18888889) ,
                       unfavorable = c(0.03611111, 0.10077778, 0.15633333) ,
                       netChange = c(0.01255556, 0.02777778, 0.03255556) ,
                       winRatio = c(1.34769231, 1.27563396, 1.20824449) )
            ## butils::object2script(test, digit = 8)
            
        }else if(method == "Peron"){
            GS <- list(favorable = c(1289.0425448, 1452.9970531, 682.33602169) ,
                       unfavorable = c(2044.84933459, 908.62963327, 578.82862552) ,
                       neutral = c(5666.10812061, 3304.48143424, 2043.31678703) ,
                       uninf = c(0, 0, 0) ,
                       favorable = c(0.14322695, 0.30467107, 0.38048618) ,
                       unfavorable = c(0.22720548, 0.32816433, 0.39247862) ,
                       netChange = c(-0.08397853, -0.02349326, -0.01199244) ,
                       winRatio = c(0.6303851, 0.92841006, 0.96944434) )
        }

        expect_equal(test, GS, tolerance = 1e-6, scale = 1)

        ## *** count pairs
        tableS <- model.tables(BT.tte, percentage = FALSE)
        expect_equal(tableS$total,
                     tableS$favorable + tableS$unfavorable + tableS$neutral + tableS$uninf,
                     tolerance = 1e-1, scale = 1) ## inexact for Peron
                     
    })
}

## ** No strata - different endpoints
for(method in c("Gehan","Peron")){ ## method <- "Gehan" ## method <- "Peron"
    test_that(paste0("BuyseTest - tte (different, ",method,", no strata)"),{ 
    
        BT.tte <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + tte(eventtime2, status2, threshold = 0.5) + tte(eventtime3, status3, threshold = 0.25),
                            data = dt.sim, scoring.rule = method,
                            correction.uninf = FALSE)

        BT2 <- BuyseTest(data = dt.sim,
                         endpoint = c("eventtime1","eventtime2","eventtime3"),
                         status = c("status1","status2","status3"),
                         treatment = "treatment",
                         type = c("tte","tte","tte"),
                         threshold = c(1,0.5,0.25),
                         scoring.rule = method,
                         correction.uninf = FALSE
                         )

        test <- list(favorable = as.double(coef(BT.tte, statistic = "count.favorable", cumulative = FALSE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "count.unfavorable", cumulative = FALSE)),
                     neutral = as.double(coef(BT.tte, statistic = "count.neutral", cumulative = FALSE)),
                     uninf = as.double(coef(BT.tte, statistic = "count.uninf", cumulative = FALSE)),
                     favorable = as.double(coef(BT.tte, statistic = "favorable", cumulative = TRUE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "unfavorable", cumulative = TRUE)),
                     netChange = as.double(coef(BT.tte, statistic = "netBenefit", cumulative = TRUE)),
                     winRatio = as.double(coef(BT.tte, statistic = "winRatio", cumulative = TRUE))
                     )

        ## *** compatibility between BuyseTests
        BT.tte@call <- list()
        BT2@call <- list()
        expect_equal(BT.tte, BT2)

        ## *** test against fixed value
        if(method == "Gehan"){
            GS <- list(favorable = c(438, 620, 794) ,
                       unfavorable = c(325, 561, 361) ,
                       neutral = c(2284, 339, 73) ,
                       uninf = c(5953, 6717, 5828) ,
                       favorable = c(0.04866667, 0.11755556, 0.20577778) ,
                       unfavorable = c(0.03611111, 0.09844444, 0.13855556) ,
                       netChange = c(0.01255556, 0.01911111, 0.06722222) ,
                       winRatio = c(1.34769231, 1.19413093, 1.48516439) )
            ## butils::object2script(test, digit = 8)
            
        }else if(method == "Peron"){
            GS <- list(favorable = c(1289.0425448, 2318.38791489, 1231.91554493) ,
                       unfavorable = c(2044.84933459, 1529.8258322, 491.18260522) ,
                       neutral = c(5666.10812061, 867.93018367, 94.79622337) ,
                       uninf = c(0, 949.96418985, 0) ,
                       favorable = c(0.14322695, 0.40082561, 0.53770511) ,
                       unfavorable = c(0.22720548, 0.39718613, 0.45176197) ,
                       netChange = c(-0.08397853, 0.00363948, 0.08594314) ,
                       winRatio = c(0.6303851, 1.00916315, 1.19023986) )
        }

        ## *** count pairs
        tableS <- model.tables(BT.tte, percentage = FALSE)
        expect_equal(tableS$total,
                     tableS$favorable + tableS$unfavorable + tableS$neutral + tableS$uninf,
                     tolerance = 1e-1, scale = 1) ## inexact for Peron
    })
}

## ** Strata - same endpoint
for(method in c("Gehan","Peron")){ ## method <- "Peron"  ## method <- "Gehan"
    test_that(paste0("BuyseTest - tte (same, ",method,", strata)"),{ 
    
        BT.tte <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 1) + tte(eventtime1, status1, threshold = 0.5) + tte(eventtime1, status1, threshold = 0.25) + strata,
                            data = dtS.sim, scoring.rule = method)

        ## *** test against fixed value
        test <- list(favorable = as.double(coef(BT.tte, statistic = "count.favorable", strata = TRUE, cumulative = FALSE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "count.unfavorable", strata = TRUE, cumulative = FALSE)),
                     neutral = as.double(coef(BT.tte, statistic = "count.neutral", strata = TRUE, cumulative = FALSE)),
                     uninf = as.double(coef(BT.tte, statistic = "count.uninf", strata = TRUE, cumulative = FALSE)),
                     favorable = as.double(coef(BT.tte, statistic = "favorable", strata = FALSE, cumulative = TRUE)),
                     unfavorable = as.double(coef(BT.tte, statistic = "unfavorable", strata = FALSE, cumulative = TRUE)),
                     netChange = as.double(coef(BT.tte, statistic = "netBenefit", strata = FALSE, cumulative = TRUE)),
                     winRatio = as.double(coef(BT.tte, statistic = "winRatio", strata = FALSE, cumulative = TRUE))
                     )

        
        if(method == "Gehan"){
            GS <- list(favorable = c(438, 438, 438, 719, 719, 719, 543, 543, 543) ,
                       unfavorable = c(325, 325, 325, 582, 582, 582, 500, 500, 500) ,
                       neutral = c(2284, 2284, 2284, 1569, 1569, 1569, 1084, 1084, 1084) ,
                       uninf = c(5953, 5953, 5953, 5367, 5367, 5367, 4809, 4809, 4809) ,
                       favorable = c(0.04867, 0.12856, 0.18889) ,
                       unfavorable = c(0.03611, 0.10078, 0.15633) ,
                       netChange = c(0.01256, 0.02778, 0.03256) ,
                       winRatio = c(1.34769, 1.27563, 1.20824) )

        } else if(method == "Peron"){
            GS <- list(favorable = c(1289.04254, 1289.04254, 1289.04254, 1452.99705, 1452.99705, 1452.99705, 682.33602, 682.33602, 682.33602) ,
                       unfavorable = c(2044.84933, 2044.84933, 2044.84933, 908.62963, 908.62963, 908.62963, 578.82863, 578.82863, 578.82863) ,
                       neutral = c(5666.10812, 5666.10812, 5666.10812, 3304.48143, 3304.48143, 3304.48143, 2043.31679, 2043.31679, 2043.31679) ,
                       uninf = c(0, 0, 0, 0, 0, 0, 0, 0, 0) ,
                       favorable = c(0.14323, 0.30467, 0.38049) ,
                       unfavorable = c(0.22721, 0.32816, 0.39248) ,
                       netChange = c(-0.08398, -0.02349, -0.01199) ,
                       winRatio = c(0.63039, 0.92841, 0.96944) )
            ## butils::object2script(test, digit = 5)
        }
        expect_equal(GS, test, tol = 1e-4, scale = 1)
        
        ## *** same result for each pair
        tableS <- model.tables(BT.tte, percentage = FALSE)
        expect_equal(tableS[tableS$strata=="1","Delta"],tableS[tableS$strata=="2","Delta"])
        expect_equal(tableS[tableS$strata=="1","Delta"],tableS[tableS$strata=="3","Delta"])
        expect_equal(tableS[tableS$strata=="1","Delta"],tableS[tableS$strata=="3","Delta"])
        
        ## *** count pairs
        dt.tableS <- as.data.table(tableS)[strata == "global"]
        expect_equal(dt.tableS[,total],
                     unname(dt.tableS[,favorable + unfavorable + neutral + uninf]),
                     tolerance = 1e-1, scale = 1) ## inexact for Peron
})
}

## * Mixed endpoints 
for(method in c("Gehan","Peron")){ ## method <- "Peron" ## method <- "Gehan"
    test_that(paste0("BuyseTest - mixed (",method,", no strata)"),{ 
    
        BT.mixed <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0.5) + cont(score1, 1) + bin(toxicity1) + tte(eventtime1, status1, threshold = 0.25) + cont(score1, 0.5),
                              data = dt.sim, scoring.rule = method)

        BT2 <- BuyseTest(data=dt.sim,
                         endpoint=c("eventtime1","score1","toxicity1","eventtime1","score1"),
                         status=c("status1","..NA..","..NA..","status1","..NA.."),
                         treatment="treatment",
                         type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                         threshold=c(0.5,1,NA,0.25,0.5),
                         scoring.rule=method)
  
        ## *** compatibility between BuyseTests
        BT.mixed@call <- list()
        BT2@call <- list()
        expect_equal(BT.mixed, BT2)

        ## *** test against fixed value
        test <- list(favorable = as.double(coef(BT.mixed, statistic = "count.favorable", cumulative = FALSE)),
                     unfavorable = as.double(coef(BT.mixed, statistic = "count.unfavorable", cumulative = FALSE)),
                     neutral = as.double(coef(BT.mixed, statistic = "count.neutral", cumulative = FALSE)),
                     uninf = as.double(coef(BT.mixed, statistic = "count.uninf", cumulative = FALSE)),
                     favorable = as.double(coef(BT.mixed, statistic = "favorable", cumulative = TRUE)),
                     unfavorable = as.double(coef(BT.mixed, statistic = "unfavorable", cumulative = TRUE)),
                     netChange = as.double(coef(BT.mixed, statistic = "netBenefit", cumulative = TRUE)),
                     winRatio = as.double(coef(BT.mixed, statistic = "winRatio", cumulative = TRUE))
                     )

        if(method == "Gehan"){
            GS <- list(favorable = c(1157, 1753, 751, 134, 373) ,
                       unfavorable = c(907, 1806, 949, 129, 323) ,
                       neutral = c(1569, 3377, 1677, 277, 718) ,
                       uninf = c(5367, 0, 0, 1137, 0) ,
                       favorable = c(0.12855556, 0.32333333, 0.40677778, 0.42166667, 0.46311111) ,
                       unfavorable = c(0.10077778, 0.30144444, 0.40688889, 0.42122222, 0.45711111) ,
                       netChange = c(0.02777778, 0.02188889, -0.00011111, 0.00044444, 0.006) ,
                       winRatio = c(1.27563396, 1.07261334, 0.99972693, 1.00105513, 1.01312591) )
            ## butils::object2script(test, digit = 8)
            
        }else if(method == "Peron"){
            GS <- list(favorable = c(2742.0395979, 792.80301972, 403.03891763, 160.70305305, 134.38721963) ,
                       unfavorable = c(2953.47896786, 896.93725328, 407.50415506, 142.85049401, 122.54879121) ,
                       neutral = c(3304.48143424, 1614.74116124, 804.19808854, 500.64454148, 243.70853064) ,
                       uninf = c(0, 0, 0, 0, 0) ,
                       favorable = c(0.30467107, 0.39276029, 0.43754239, 0.45539829, 0.4703302) ,
                       unfavorable = c(0.32816433, 0.42782402, 0.47310226, 0.48897454, 0.50259107) ,
                       netChange = c(-0.02349326, -0.03506373, -0.03555987, -0.03357625, -0.03226087) ,
                       winRatio = c(0.92841006, 0.91804169, 0.92483682, 0.93133333, 0.93581089) )
        }
        
        expect_equal(test, GS, tolerance = 1e-6, scale = 1)

        ## *** count pairs
        tableS <- model.tables(BT.mixed, percentage = FALSE)
        expect_equal(tableS$total,
                     tableS$favorable + tableS$unfavorable + tableS$neutral + tableS$uninf)

    })
}


test_that("ordering does not matter", {
    BT.mixed1 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0.25) + cont(score1, 1),
                           data = dt.sim, scoring.rule = method)
    BT.mixed2 <- BuyseTest(treatment ~ tte(eventtime1, status1, threshold = 0.5) +  tte(eventtime1, status1, threshold = 0.25) + cont(score1, 1),
                           data = dt.sim, scoring.rule = method)
    expect_equal(coef(BT.mixed2, statistic = "netBenefit")[2:3], coef(BT.mixed1, statistic = "netBenefit"), tol = 1e-6)
    expect_equal(coef(BT.mixed2, statistic = "winRatio")[2:3], coef(BT.mixed1, statistic = "winRatio"), tol = 1e-6)
})

test_that(paste0("BuyseTest - Peron scoring rule with 2 TTE, one without censoring"),{ 
    ## 1 continuous
    ## 2 Gehan left-censoring
    ## 3 Gehan right-censoring
    ## 4 Peron right-censoring survival
    ## 5 Peron right-censoring competing risks
    
    BT.mixed <- BuyseTest(treatment ~ tte(eventtime2, status2, threshold = 0.5) + tte(eventtime1, status1.noC, threshold = 0),
                          data = dt.sim, scoring.rule = "Peron")
    expect_equal(unname(attr(BT.mixed@scoring.rule,"method.score")), c("SurvPeron","continuous"))
    ## summary(BT.mixed)
    BT.mixed <- BuyseTest(treatment ~ tte(eventtime1, status1.noC, threshold = 0) + tte(eventtime2, status2, threshold = 0.5),
                          data = dt.sim, scoring.rule = "Peron")
    ## summary(BT.mixed)
    expect_equal(unname(attr(BT.mixed@scoring.rule,"method.score")), c("continuous","SurvPeron"))
    
})

## * Left censoring
test_that("BuyseTest - left vs. right censoring", {

    BT.left <- BuyseTest(treatment ~ tte(eventtime1, status = status1, censoring = "left"),
                         data = dt.sim,
                         scoring.rule = "Gehan")

    expect_equal(as.double(coef(BT.left)), 0.09488889, tol = 1e-6)

    BT.left <- BuyseTest(treatment ~ tte(eventtime1, status = status1, censoring = "left"),
                         data = dt.sim,
                         scoring.rule = "Gehan",
                         correction.uninf = TRUE)

    expect_equal(as.double(coef(BT.left)), 0.1768116, tol = 1e-6)
})

## * Gaussian endpoint

## ** uncorrelated
df.1 <- data.frame(mean = 0:1, sd = 1, treatment = c("C","T"))
df.2 <- data.frame(mean = 0:1, sd = c(2,0.5), treatment = c("C","T"))
df.3 <- rbind(df.1,df.2)

test_that("BuyseTest - uncorrelated gaussians", {

    GS.1 <- 1 - pnorm(0, mean = 1, sd = sqrt(2))
    ## GS.1 - mean( rnorm(n.GS, mean = 1)> rnorm(n.GS, mean = 0))
    BTG.1 <- BuyseTest(treatment ~ gaus(mean, sd),
                       data = df.1, method.inference = "none")
    expect_equal(GS.1,as.double(coef(BTG.1, statistic = "favorable")),tol=1e-6)

    GS.2 <- 1 - pnorm(0, mean = 1, sd = sqrt(4.25))
    ## GS.2 - mean( rnorm(n.GS, mean = 1, sd = 0.5)> rnorm(n.GS, mean = 0, sd = 2))
    BTG.2 <- BuyseTest(treatment ~ gaus(mean = mean, std = sd),
                       data = df.2, method.inference = "none")
    expect_equal(GS.2,as.double(coef(BTG.2, statistic = "favorable")),tol=1e-6)

    GS.3 <- mean(c(GS.1, (1 - pnorm(0, mean = 1, sd = sqrt(5))), 1 - pnorm(0, mean = 1, sd = sqrt(1.25)), GS.2))
    ## GS.3 - mean(c(GS.1,mean( rnorm(n.GS, mean = 1)> rnorm(n.GS, mean = 0, sd = 2)), mean( rnorm(n.GS, mean = 1, sd = 0.5)> rnorm(n.GS, mean = 0)),GS.2))
    BTG.3 <- BuyseTest(treatment ~ gaus(mean = mean, std = sd),
                       data = df.3, method.inference = "none")
    expect_equal(GS.3,as.double(coef(BTG.3, statistic = "favorable")),tol=1e-6)

})

## ** correlated


complement <- function(rho, n) {## generate a dataset with given correlation
    ## adapted from
    ## https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
    x <- rnorm(n) 
    y <- rnorm(n) 
    y.perp <- residuals(lm(x ~ y))
    z <- rho * sd(y.perp) * y + sqrt(1 - rho^2) * sd(y) * y.perp
    return(list(Y=as.double(y),X=as.double(z)))
}
## cor(complement(rho = 0.5, n = 10))

df.1$iid <- complement(rho = 0.5, n = 10)
df.2$iid <- complement(rho = 0.5, n = 10)
df.3 <- rbind(df.1,df.2)


test_that("BuyseTest - correlated gaussians", {

    GS.1 <- 1 - pnorm(0, mean = 1, sd = sqrt(1))
    ## GS.1 - mean(apply(mvtnorm::rmvnorm(n.GS, mean = 0:1, sigma = matrix(c(1,0.5,0.5,1),2,2)),1, FUN = function(x){x[2]>x[1]}))
    BTG.1 <- BuyseTest(treatment ~ gaus(mean, sd, iid),
                       data = df.1, method.inference = "none")
    expect_equal(GS.1,as.double(coef(BTG.1, statistic = "favorable")),tol=1e-6)
    
    GS.2 <- 1 - pnorm(0, mean = 1, sd = sqrt(3.25)) ## 2^2+0.5^2-2*0.5*0.5*2
    ## GS.2 - mean(apply(mvtnorm::rmvnorm(10*n.GS, mean = 0:1, sigma = matrix(c(0.5^2,0.5,0.5,2^2),2,2)),1, FUN = function(x){x[2]>x[1]}))
    BTG.2 <- BuyseTest(treatment ~ gaus(mean = mean, std = sd, iid),
                       data = df.2, method.inference = "none")
    expect_equal(GS.2,as.double(coef(BTG.2, statistic = "favorable")),tol=1e-6)

    GS.3 <- mean(c(GS.1,
              1 - pnorm(0, mean = 1, sd = sqrt(1.25-cor(df.1$iid[[1]],df.2$iid[[2]]))),
              1 - pnorm(0, mean = 1, sd = sqrt(5-4*cor(df.1$iid[[2]],df.2$iid[[1]]))),
              GS.2))
    ## GS.3 - c(GS.1,mean( rnorm(n.GS, mean = 1)> rnorm(n.GS, mean = 0, sd = 2)), mean( rnorm(n.GS, mean = 1, sd = 0.5)> rnorm(n.GS, mean = 0)),GS.2)
    BTG.3 <- BuyseTest(treatment ~ gaus(mean = mean, std = sd, iid = iid),
                       data = df.3, method.inference = "none")
    expect_equal(GS.3,as.double(coef(BTG.3, statistic = "favorable")),tol=1e-6)

})


## * dataset [save]
## dt.sim <- data.table("treatment" = c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T"), 
           ## "eventtime1" = c(0.29972933, 0.17098301, 0.03202016, 0.44202235, 0.10311930, 1.15511106, 0.56124439, 1.21925282, 0.17187895, 0.29113518, 0.09773182, 0.67653288, 0.03861652, 0.07063795, 0.01458732, 0.71823701, 1.08732315, 0.15199073, 0.60176965, 0.41266664, 0.70669922, 0.43573838, 0.07363507, 0.96531658, 0.05755952, 0.94544071, 3.92245778, 0.26717898, 0.16922697, 0.73171108, 1.56385321, 0.00937406, 0.16954569, 0.09197029, 0.09036225, 0.18974451, 0.04073623, 0.02378406, 0.56634940, 1.61507125, 0.49139404, 0.19115956, 0.13740882, 0.36936734, 0.36919469, 1.41367871, 0.94269045, 0.01191205, 0.37297697, 0.19502322, 0.01296422, 0.31343847, 0.00213360, 0.28171661, 0.17320335, 0.06992269, 0.03277328, 0.21255628, 1.43421433, 0.50712777, 0.24571909, 1.00813698, 0.51160661, 0.09173387, 0.28063656, 1.14177394, 0.12593593, 1.58859472, 0.07077964, 0.11041468, 0.09741926, 0.56342077, 0.23108781, 0.76598116, 0.02193362, 0.14312356, 1.36059222, 0.85553186, 0.38761972, 0.05592164, 0.24080708, 2.23146741, 0.65659820, 0.12662146, 0.33644115, 0.93738422, 0.93216642, 0.80139621, 0.65390255, 0.60241389, 0.34299720, 0.66186296, 1.10529116, 0.14865979, 0.12501623, 0.04451988, 0.48423927, 0.92904199, 0.32060823, 0.20941169, 0.29373301, 0.99816128, 1.33980963, 0.16543365, 1.22099704, 0.03737215, 0.16298912, 0.32335369, 0.39027702, 0.10348081, 1.03796020, 0.47513692, 0.24106903, 0.45926525, 0.49608224, 1.44827529, 0.52186516, 0.68353467, 0.01981440, 0.18416592, 0.97426659, 1.77382739, 0.33398520, 0.19615994, 0.03780470, 0.17649501, 0.22901002, 0.02323039, 0.20845366, 0.52986278, 0.74053528, 0.27117162, 0.19489030, 0.66019467, 0.88323068, 0.32137872, 0.17473734, 0.10029520, 0.08490402, 0.34625873, 1.92253508, 1.24734174, 0.20656446, 1.47517308, 0.00019922, 0.33644227, 0.26031585, 0.24064544, 0.87270237, 0.50428180, 0.55966656, 1.09623527, 0.00653237, 0.51073806, 0.36994855, 0.74641533, 0.44120419, 0.98441248, 0.27874909, 0.29785716, 0.19272977, 0.03585685, 0.07027667, 0.00518237, 0.13138687, 0.03262746, 0.26673138, 0.22325116, 0.71796943, 0.29428682, 0.74450076, 0.29965883, 0.17469397, 1.73676014, 1.38930578, 1.61992553, 0.73321636, 0.79600567, 0.04142438, 0.94565307, 0.00825042, 0.65877918, 0.76745178, 1.11121647, 1.58163545, 0.10784914, 0.94274529, 0.05602611, 0.59380396, 1.25969953), 
           ## "status1" = c(1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1), 
           ## "eventtime2" = c(0.13556029, 0.53306031, 0.75447141, 0.25613517, 2.07382035, 0.03757998, 0.14045294, 1.43718090, 0.25834648, 0.00990976, 0.09679214, 0.10656134, 0.12146487, 1.65973920, 0.36003623, 0.27730941, 0.38657147, 0.24456645, 1.45490109, 1.36873057, 0.59120666, 0.13334589, 0.01156500, 0.53612156, 0.50263578, 0.83841401, 1.13769949, 0.32597204, 2.28409336, 0.24155223, 0.14974690, 0.44926621, 0.04444326, 1.74047606, 0.14653971, 0.05550446, 0.91919142, 0.19709401, 0.13809616, 0.19533969, 2.16819532, 1.63080095, 0.63950588, 0.39308933, 0.93342527, 1.44372187, 0.07228017, 1.65850544, 1.65081603, 1.59301263, 0.56652696, 1.18547005, 0.17771856, 0.88895104, 0.55326678, 0.53893584, 0.06138524, 0.41325058, 0.50743982, 0.56196957, 0.05072848, 0.78399042, 0.14126094, 0.37339708, 1.71804695, 0.61959578, 0.37048513, 0.19876601, 1.13166471, 0.16526419, 1.00895604, 0.27660263, 0.15692162, 0.56680821, 1.02953170, 0.15395316, 0.18412961, 0.35121113, 1.71637364, 0.37027203, 0.05331582, 0.41455140, 0.40164440, 0.40714141, 1.60638089, 0.42633103, 0.21886920, 0.12911882, 0.21075684, 0.41380614, 0.13020199, 0.83162531, 0.33213999, 0.25378188, 0.03565690, 1.79972143, 0.49513339, 0.85519650, 0.95797393, 1.18930068, 1.52944416, 0.21211345, 0.36342043, 1.12946317, 0.11842668, 1.50611081, 0.47826400, 0.58815796, 0.20995225, 0.25050953, 0.38504902, 2.57865824, 2.37486593, 0.37757152, 0.11404643, 0.05407206, 0.42755586, 0.06360704, 0.04317937, 0.45965630, 0.40623887, 0.21847145, 0.39437507, 0.88480211, 1.40718306, 0.64707974, 0.08332118, 0.36962127, 0.60152779, 0.39706135, 0.55125693, 0.36913746, 1.42278678, 0.69311190, 1.01065256, 1.08925374, 1.34066288, 0.59957988, 0.04203430, 2.77233260, 3.28708257, 1.73709539, 0.45768357, 0.32263242, 0.29657430, 0.02366551, 0.20247683, 1.35654772, 0.00694441, 1.38201424, 0.89090216, 0.88823543, 1.41377148, 0.37135459, 0.36557318, 1.90512208, 0.31316393, 1.10058790, 0.36843826, 1.04621615, 0.99875000, 0.12788404, 0.36530394, 0.05811976, 2.05009814, 0.51824171, 0.87219406, 0.13617999, 1.00594646, 0.74437044, 0.00258926, 0.57609633, 0.39368111, 0.39772202, 0.31094959, 0.37548816, 2.17934168, 0.99261368, 0.25028018, 0.04431970, 0.77118728, 1.56589807, 2.07293061, 0.90534207, 1.07834985, 0.16480664, 0.14750491, 0.30542754, 0.19788267, 0.07055950), 
           ## "status2" = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1), 
           ## "eventtime3" = c(0.25956955, 0.50870505, 0.38243703, 0.48507356, 0.58090695, 0.09143226, 0.11595455, 2.45142917, 0.31529917, 0.01063624, 0.23269836, 0.00314356, 0.35192745, 0.05702643, 0.77426085, 0.04011495, 0.47305223, 1.32361908, 0.49734166, 0.08057535, 0.31222287, 0.49705636, 0.78467962, 0.85062628, 0.01585449, 0.84971402, 0.14097533, 0.64007436, 0.47504948, 0.47065190, 0.91619564, 3.14863908, 0.72449497, 0.34146462, 0.06298503, 1.26862569, 0.07311503, 0.26950937, 0.24296576, 1.46229570, 0.44175144, 0.08437995, 0.11765742, 1.26484624, 0.13972311, 0.88368353, 0.10077329, 0.30071004, 0.18178031, 0.40319616, 0.19262871, 0.84156278, 1.01319195, 1.09334295, 0.06925393, 0.08496785, 0.38289653, 0.77851078, 2.90459458, 0.70906019, 1.44433835, 0.31710947, 0.83804625, 0.52672195, 1.39708324, 0.05738464, 0.00424703, 0.02745054, 0.15640178, 0.60252617, 0.45624187, 0.03877660, 0.26583575, 1.93489936, 0.16157491, 0.16150214, 3.12000133, 1.15754730, 0.20733374, 1.36244884, 0.85908195, 1.17088649, 1.04190785, 0.58512798, 0.17684563, 0.39304759, 0.50360868, 0.25826671, 1.36782193, 0.79286184, 0.54019913, 0.54899883, 0.01927732, 0.83191354, 2.95844611, 0.66324356, 0.37850024, 1.01325887, 0.68367717, 0.16975714, 1.07644784, 2.05425366, 0.76593812, 0.93194348, 0.46623093, 2.96814573, 0.12555074, 1.85279179, 0.91838000, 2.96795061, 0.20853482, 0.55747755, 1.16290689, 0.25204838, 1.45458273, 0.47887218, 2.14245439, 0.73046914, 1.21973505, 0.24528169, 0.60239017, 1.79625436, 0.09840920, 0.42368372, 0.07741995, 1.28366723, 0.51361326, 1.44102172, 0.23235485, 0.00745763, 1.46408645, 0.28432717, 0.50773758, 0.12780739, 1.62522052, 0.60240232, 0.53551248, 0.37865307, 1.43088588, 0.49141086, 0.27257514, 1.28147177, 1.11686803, 0.37442988, 1.00084367, 1.78079525, 0.36024791, 1.50952573, 0.36300718, 0.73043847, 0.25946183, 0.86342025, 0.86724760, 0.65525025, 0.34944216, 0.78352676, 0.76614068, 0.03508025, 1.10827027, 0.13490347, 2.82395488, 0.42936653, 0.15014156, 0.82605928, 0.38453132, 1.19652345, 0.54175957, 0.40951641, 0.25130183, 2.10913985, 3.90959749, 0.83906640, 0.35827788, 0.82174584, 0.43750343, 0.72346693, 2.07799650, 0.03194980, 0.02397542, 0.84753338, 0.39459503, 1.40010494, 1.05098332, 2.16823693, 1.17526902, 1.21647314, 0.20328870, 0.08513324, 0.20774038, 0.14752052), 
           ## "status3" = c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0), 
           ## "toxicity1" = c("2", "1", "1", "2", "1", "2", "2", "1", "1", "1", "2", "2", "2", "1", "2", "1", "2", "2", "2", "2", "2", "2", "2", "2", "1", "2", "2", "2", "2", "1", "1", "2", "1", "1", "2", "2", "2", "1", "2", "2", "2", "2", "1", "2", "2", "2", "2", "1", "1", "1", "2", "1", "2", "1", "2", "1", "1", "2", "2", "2", "2", "2", "1", "2", "1", "1", "1", "1", "1", "1", "2", "1", "2", "2", "1", "1", "2", "1", "2", "2", "2", "2", "1", "1", "1", "2", "2", "2", "1", "1", "2", "2", "2", "2", "1", "2", "1", "2", "1", "2", "1", "1", "2", "2", "2", "1", "1", "2", "2", "2", "1", "2", "1", "2", "1", "2", "2", "1", "2", "2", "1", "2", "1", "1", "2", "2", "2", "2", "2", "2", "1", "2", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "1", "1", "1", "1", "2", "2", "2", "2", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "1", "1", "2", "1", "1", "1", "2", "2", "1", "1", "2", "1", "2", "1", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "1", "2", "2", "1", "2", "2"), 
           ## "toxicity2" = c("2", "2", "2", "2", "2", "1", "2", "2", "1", "2", "1", "2", "2", "2", "2", "2", "1", "2", "1", "1", "2", "2", "2", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "1", "2", "2", "2", "2", "2", "2", "2", "2", "1", "1", "1", "2", "2", "2", "2", "1", "2", "2", "2", "2", "1", "2", "2", "2", "2", "2", "1", "2", "1", "2", "2", "2", "1", "2", "2", "1", "2", "2", "2", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "2", "1", "1", "1", "2", "2", "2", "1", "2", "2", "1", "2", "2", "2", "2", "2", "2", "2", "2", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "1", "1", "1", "2", "1", "1", "1", "1", "2", "2", "2", "2", "1", "2", "2", "1", "2", "1", "2", "2", "2", "1", "2", "1", "2", "2", "2", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "1", "2", "2", "2", "2", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "1", "1", "1", "2", "1", "1", "1", "1", "2", "2", "1", "2", "2", "2", "2", "1", "2", "2"), 
           ## "score1" = c( 1.37248589,  2.19299163,  1.73483525, -0.00752315,  1.88739264,  1.41874690,  1.29806497,  1.10816840,  1.04006544,  2.47152734, -0.75581445,  0.90738612,  2.23780927,  1.42039507,  1.16518174,  0.17216119,  1.20579932,  0.00104613,  1.02850830,  1.72188935,  1.20350882,  0.75061290,  3.07140367,  0.01723101,  1.24266243,  1.00337130, -0.55113598,  2.22269088,  0.62567229,  1.03200894,  0.76908950,  0.72202409, -0.33953529,  0.16539723,  1.33126684,  0.13868705,  1.44159579,  1.73077039,  2.37832330,  1.48123808,  0.43160998,  0.60855246,  1.64288621,  1.27620613,  0.74186695,  2.46582983, -0.97327590,  1.92010988,  0.85092865,  2.23738748,  1.52738054,  2.49049279,  0.60756862,  1.18796458,  0.98872659, -0.24592590,  1.94396826,  1.45941533, -0.65322897,  0.22658973,  0.31024850,  2.67183286,  0.90395013,  0.20399559,  0.77525623, -0.26017269, -0.42640232, -1.65363714, -0.03865936,  1.57982435,  3.24926361,  0.65269903,  1.35892133,  0.43765887,  0.03994969,  1.34253659, -0.04648803,  1.55333508,  0.71948130,  0.74599877,  1.99124934,  0.33964643,  1.69329156,  2.52728405,  1.23695659,  2.17387684, -1.48614865,  2.49742470,  1.97382087,  0.51588451,  1.27685067,  0.78477715,  2.70827066,  2.53843123,  0.76378630, -1.07643183,  0.39430703,  1.95235963, -0.18769287,  4.81258034,  2.16396750, -0.52338460, -1.51833505,  0.29247077,  0.71256712,  0.56469169,  0.65692123,  0.96068912,  1.88696599,  0.64005160,  0.27104573,  2.75174562,  0.91396141,  2.10636302,  0.98082216, -0.49346018,  3.70063662,  0.25630558,  2.06519498,  0.96791827,  0.46004031, -0.92564357,  2.00783138,  0.72076520, -0.24956586,  2.24849113,  0.80778662,  1.91197623, -1.15826968,  2.28961053,  0.57189003,  0.74499711,  2.32715824,  1.83869012, -0.46500068,  0.38374730,  1.10610584, -0.84401676, -0.61642285,  1.22767670,  1.92121670,  0.66879656,  2.28723403,  1.05726080, -1.20393313,  0.47708965,  2.08078531,  0.75685866,  0.09055344,  0.13170049,  1.86947504,  0.31999040,  1.17321454,  0.84056196,  1.79349943,  2.69435049,  2.23996869,  1.02943674,  2.65721450,  2.13122314, -0.40241063,  1.15677046,  1.86750860,  0.96676484,  1.95406456,  2.13009673,  1.39369731,  1.37680833,  0.74825709,  1.50409489, -0.32980347,  0.88170055,  0.36675153,  1.67769134,  0.70211824,  1.62792225,  2.67843073, -1.40000254, -0.34984289,  2.20200790,  2.29384272,  0.94040565,  0.31493765, -0.45343930,  2.29412656,  2.17616909,  1.86404508, -1.29673252, -0.39181001, -0.38587675), 
           ## "score2" = c( 2.65317519,  2.19804984,  3.39958837,  2.36418421,  3.64077145,  3.13527144,  0.78355959,  2.90098701,  1.38803990,  3.13141377,  0.84382164,  0.89827149,  2.57784739,  2.18087366,  0.01207711,  3.36059519,  0.82579966,  1.06209054,  3.26831405,  3.36400522,  0.15778253,  1.43518251,  0.34003090,  3.25256983,  0.36760006,  3.05366327, -0.76870569,  0.21894658,  2.25633914,  2.00545663,  2.74176995, -0.68560240,  3.10162133,  3.21802045,  3.50102058,  0.81409682,  0.57813789,  1.51272481,  2.58998056,  1.54575954,  1.78942788,  2.93706409,  1.41297510,  2.49658599,  1.73979822,  3.38222640,  2.56692197,  2.07484513,  1.08979248,  1.04274682,  0.72159509,  2.00701403,  1.81918381,  0.86052608,  2.74381724,  2.62444167,  0.58176694,  2.93314119,  3.15064738,  2.47002921,  2.87963124,  1.07675302,  2.25137934,  0.74625923,  1.84460737,  2.51732094,  2.42306429,  3.15885851,  3.90794404,  2.12231655,  2.95211161,  3.38124150,  1.31147508,  1.99527444,  2.70690665,  2.63814472,  0.21252668,  1.79456144,  2.15728030,  2.77809395,  2.44579281,  1.93244981,  0.51950279,  1.79247109,  1.68026738,  3.28908680,  2.07278306,  1.99467467,  2.81879144,  2.95246510,  0.06090360,  2.87884707,  1.59755053,  1.11761304,  2.61240767,  3.13087171,  0.36334251,  3.08359626,  1.91447252,  2.56312532,  3.09462984,  2.72640935,  0.60034286,  0.71220484,  1.03158452,  0.26331424,  2.83169736,  1.24104689,  4.15507700,  1.57674112,  1.99435006,  2.40036173,  1.63432220,  1.57085616,  2.24469713,  4.21977625,  2.31426981,  1.43621493,  0.51430384,  2.41529859,  1.47193632,  1.64952938,  2.73062194,  1.44337583,  0.55023946,  0.94800357,  1.31716138,  1.38705882,  3.89078004,  1.85458947,  1.36057628,  1.58882395,  2.34279188,  0.86675325,  2.06607030,  2.03799977,  2.92106475,  4.13496564,  1.40017307,  3.24444646,  1.92054298,  3.18175155,  4.18614406,  2.40617493,  1.26164090,  0.04351330,  0.04995430,  1.05900218,  3.19778677,  1.37576061,  1.86713034,  1.98277908,  1.54037072,  3.47293679,  4.16931961, -1.00143131,  0.22801442,  1.63516078,  2.37453378,  0.76592207,  2.47489342,  3.23305746,  1.40687134,  2.75226244,  2.61204551,  1.77134263,  1.29516065,  3.22951605,  3.51851513,  3.34509418,  1.87630847,  2.93390245,  2.01949685,  2.83085211,  1.69550533,  1.94071872,  3.00739670,  0.62220491,  0.99579727,  1.97671970,  3.54674727,  1.24063675,  3.00679266,  2.27421492,  1.91961764,  1.09644750,  2.47314510,  3.38848898,  2.03474609,  0.70741497), 
           ## "score3" = c(1.1422912, 3.7273816, 3.3474971, 3.8424655, 2.9674466, 1.8834731, 1.9977386, 1.9200766, 2.9129226, 1.6749344, 2.7190415, 1.9429009, 4.8992433, 4.5210321, 3.2763423, 1.7971510, 2.1940079, 2.4692713, 3.7387694, 5.6624053, 2.4008662, 2.5411284, 1.4538581, 3.9670348, 0.6024646, 4.3625804, 3.2195859, 2.3472134, 3.2848838, 3.1957958, 2.0305307, 1.8458276, 2.9608437, 2.6592635, 4.0677380, 2.8616322, 3.8583077, 4.1206752, 2.1229097, 2.2782372, 3.5029077, 1.6258354, 3.0368000, 3.5862009, 3.6839435, 2.4312439, 3.2703556, 2.4580630, 3.4026097, 2.8134101, 3.6427508, 3.6714750, 2.5989603, 5.1650258, 1.6047065, 2.8292999, 2.8305656, 2.9670166, 5.3885198, 2.0099881, 2.8989458, 1.0858509, 2.4511256, 3.7138949, 4.8744887, 3.7946575, 2.9052863, 1.6529439, 3.2264245, 3.0241093, 1.1952988, 3.7405965, 2.6342226, 3.9380437, 2.7353587, 2.8336375, 2.9961175, 3.4351701, 2.3799796, 2.7542375, 2.8883925, 3.5208775, 2.7418243, 2.9288784, 4.1494674, 0.5699271, 1.2895520, 3.2585224, 3.4980213, 2.9506039, 4.4946761, 0.4182285, 3.0610335, 3.3277493, 3.1286587, 2.5031202, 3.4321744, 4.4126444, 0.9265798, 3.6894032, 2.7362238, 2.6441475, 2.5188452, 3.2227355, 5.4299103, 4.4962981, 2.2827918, 2.5329461, 3.6629235, 5.3001769, 3.3275101, 3.0638635, 1.8604391, 4.1804102, 3.0413885, 1.7864064, 3.0731958, 2.7426754, 3.2668064, 4.3877243, 3.1930796, 3.5923166, 2.1700255, 3.3925733, 3.3848676, 4.0510447, 4.1557975, 1.9655620, 2.7455319, 4.2736843, 4.5025446, 3.5904095, 2.3693145, 3.7923495, 3.1253846, 3.3227550, 2.5544168, 3.7668439, 1.5964970, 1.8239532, 3.5115965, 4.3167653, 5.3929130, 2.9322306, 2.7962149, 4.9964344, 4.5432178, 1.7716624, 5.6444880, 3.7570137, 4.0553452, 3.9579449, 3.8338791, 2.8307244, 2.5955664, 3.2956016, 2.7681255, 2.5519218, 3.0223108, 3.0444673, 3.4807212, 3.6356199, 0.9992576, 2.3093445, 2.8693548, 2.6557482, 2.9475872, 3.0961239, 3.2664071, 3.5547935, 4.2354459, 3.2851334, 2.5412071, 3.6234780, 2.2760049, 4.6194187, 2.3833266, 3.3845411, 5.1145213, 1.9714326, 2.1132120, 4.2711460, 1.3949146, 4.1222734, 5.1584386, 3.4282466, 4.2011787, 4.0316901, 3.6538742, 5.0120818), 
           ## "status1.noC" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
