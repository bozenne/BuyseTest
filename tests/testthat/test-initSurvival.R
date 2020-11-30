if(FALSE){
    library(testthat)
    library(BuyseTest)
    library(data.table)
    library(prodlim)
}

context("Check computation")
BuyseTest.options(keep.survival = TRUE,
                  trace = 0,
                  precompute = FALSE,
                  method.inference = "none")

## does the same as prodlim:::predict.prodlim
## except continue the survival at 0 after last event if it is 0 at last event
predictPL2 <- function(object, times, ...){
    pred <- predict(object, times = times, ...)

    index.plus <- which(times>max(object$time))
    if( (length(index.plus)>0) && (utils::tail(object$surv,1) == 0) ){
        pred[index.plus] <- 0
    }

    return(pred)
}

## * No strata

## ** settings
dataT <- data.table(time = 1:5,
                    treatment = "T",
                    status1 = c(1,0,1,1,1),
                    status2 = c(1,0,1,1,1),
                    status3 = c(1,1,1,1,1))
dataC <- data.table(time = c(1:5-0.1,5,5),
                    treatment = "C",
                    status1 = c(1,1,0,1,0,0,0),
                    status2 = c(1,1,0,1,0,1,1),
                    status3 = c(1,1,1,1,1,1,1))
data <- rbind(dataC, dataT)
seqThreshold <- c(0,0.5,1.1)

## ** tests
for(iData in 1:3){ ## iData <- 1

    data[, status := .SD[[paste0("status",iData)]]]

    e.survC <- prodlim(Hist(time, status) ~ 1, data = data[treatment=="C"])
    e.survT <- prodlim(Hist(time, status) ~ 1, data = data[treatment=="T"])
    ## plot(e.survT)
    ## plot(e.survC)

    for(iThreshold in seqThreshold){ ## iThreshold <- 0
        
        ## *** Compute survival
        form <- as.formula(paste0("treatment ~ tte(time, status = status, threshold = ",iThreshold,")"))
        outBT <- BuyseTest(form, data = data) ## outBT <- BuyseTest(form, data = data, method.inference = "u-statistic")
        outSurv <- getSurvival(outBT, endpoint = 1, strata = 1, unlist = TRUE)

        if(iData==3){next} ## survival not used because no censoring
        
        ## *** Check at jump times
        test_that("initSurvival (jump times, no strata)",{
            ## correct jump times
            expect_equal(e.survC$time[e.survC$hazard>0],
                         outSurv$survJumpC[,"time"])
            expect_equal(e.survT$time[e.survT$hazard>0],
                         outSurv$survJumpT[,"time"])

            ## correct survival
            expect_equal(predictPL2(e.survC, times = outSurv$survJumpT[,"time"] + iThreshold),
                         outSurv$survJumpT[,"survival"])
            expect_equal(predictPL2(e.survT, times = outSurv$survJumpC[,"time"] + iThreshold),
                         outSurv$survJumpC[,"survival"])

            ## correct jumps in survival
            expect_equal(diff(predictPL2(e.survC, times = c(0,outSurv$survJumpC[,"time"]))),
                         outSurv$survJumpC[,"dSurvival"])
            expect_equal(diff(predictPL2(e.survT, times = c(0,outSurv$survJumpT[,"time"]))),
                         outSurv$survJumpT[,"dSurvival"])

            
        })
        
        ## *** Check at observation time
    
        test_that("initSurvival (observation times, no strata)",{
            ## correct time
            expect_equal(data[treatment=="C",time],
                         outSurv$survTimeC[,"time"])

            expect_equal(data[treatment=="T",time],
                         outSurv$survTimeT[,"time"])

            ## correct survival
            expect_equal(predictPL2(e.survC, times = data[treatment=="C",time] - outBT@threshold),
                         outSurv$survTimeC[,"survivalC-threshold"])
            expect_equal(predictPL2(e.survC, times = data[treatment=="C",time]),
                         outSurv$survTimeC[,"survivalC_0"])
            expect_equal(predictPL2(e.survC, times = data[treatment=="C",time] + iThreshold),
                         outSurv$survTimeC[,"survivalC+threshold"])

            expect_equal(predictPL2(e.survT, times = data[treatment=="C",time] - outBT@threshold),
                         outSurv$survTimeC[,"survivalT-threshold"])
            expect_equal(predictPL2(e.survT, times = data[treatment=="C",time]),
                         outSurv$survTimeC[,"survivalT_0"])
            expect_equal(predictPL2(e.survT, times = data[treatment=="C",time] + iThreshold),
                         outSurv$survTimeC[,"survivalT+threshold"])

            ## correct survival (treatment observation times)
            expect_equal(predictPL2(e.survC, times = data[treatment=="T",time] - outBT@threshold),
                         outSurv$survTimeT[,"survivalC-threshold"])
            expect_equal(predictPL2(e.survC, times = data[treatment=="T",time]),
                         outSurv$survTimeT[,"survivalC_0"])
            expect_equal(predictPL2(e.survC, times = data[treatment=="T",time] + iThreshold),
                         outSurv$survTimeT[,"survivalC+threshold"])

            expect_equal(predictPL2(e.survT, times = data[treatment=="T",time] - outBT@threshold),
                         outSurv$survTimeT[,"survivalT-threshold"])
            expect_equal(predictPL2(e.survT, times = data[treatment=="T",time]),
                         outSurv$survTimeT[,"survivalT_0"])
            expect_equal(predictPL2(e.survT, times = data[treatment=="T",time] + iThreshold),
                         outSurv$survTimeT[,"survivalT+threshold"])
        })
    
    }
}
    
## * With strata

## ** settings
threshold <- 0.001

dataT <- data.table(time = 1:5,
                    treatment = "T",
                    status = c(1,0,1,1,1)
                    )
dataC <- data.table(time = c(1:5-0.1,5,5),
                    treatment = "C",
                    status = c(1,1,0,1,0,0,0)
                    )
data <- rbind(dataC, dataT)
dataStrata <- rbind(cbind(data, strata = 1),
                    cbind(data, strata = 2)
                    )

## ** tests
outBT <- BuyseTest(treatment ~ tte(time, status = status) + strata,
                   data = dataStrata)
outSurv <- getSurvival(outBT, endpoint = 1, unlist = TRUE)

for(iStrata in 1:2){ ## iStrata <- 1
    iE.survC <- prodlim(Hist(time, status) ~ 1, data = dataStrata[treatment=="C" & strata == iStrata])
    iE.survT <- prodlim(Hist(time, status) ~ 1, data = dataStrata[treatment=="T" & strata == iStrata])

    ## *** Control at jump times
    test_that("initSurvival (jump times, strata)",{    
        ## correct jump times
        expect_equal(iE.survC$time[iE.survC$hazard>0],
                     outSurv$survJumpC[[iStrata]][,"time"])
        expect_equal(iE.survT$time[iE.survT$hazard>0],
                     outSurv$survJumpT[[iStrata]][,"time"])

        ## correct jumps in survival
        expect_equal(diff(c(1,iE.survC$surv[iE.survC$hazard>0])),
                     outSurv$survJumpC[[iStrata]][,"dSurvival"])
        expect_equal(diff(c(1,iE.survT$surv[iE.survT$hazard>0])),
                     outSurv$survJumpT[[iStrata]][,"dSurvival"])

        ## correct survival
        expect_equal(predictPL2(iE.survT, times = iE.survC$time[iE.survC$hazard>0] - outBT@threshold),
                     outSurv$survJumpC[[iStrata]][,"survival"])
        expect_equal(predictPL2(iE.survC, times = iE.survT$time[iE.survT$hazard>0]),
                     outSurv$survJumpT[[iStrata]][,"survival"])
    })
    
    ## *** Control at observation time
    test_that("initSurvival (observation times, strata)",{        
        ## correct time
        expect_equal(data[treatment=="C",time],
                     outSurv$survTimeC[[iStrata]][,"time"])

        expect_equal(data[treatment=="T",time],
                     outSurv$survTimeT[[iStrata]][,"time"])

        ## correct survival (control observation times)
        expect_equal(predictPL2(iE.survC, times = data[treatment=="C",time] - outBT@threshold),
                     outSurv$survTimeC[[iStrata]][,"survivalC-threshold"])
        expect_equal(predictPL2(iE.survC, times = data[treatment=="C",time]),
                     outSurv$survTimeC[[iStrata]][,"survivalC_0"])
        expect_equal(predictPL2(iE.survC, times = data[treatment=="C",time]),
                     outSurv$survTimeC[[iStrata]][,"survivalC+threshold"])

        expect_equal(predictPL2(iE.survT, times = data[treatment=="C",time] - outBT@threshold),
                     outSurv$survTimeC[[iStrata]][,"survivalT-threshold"])
        expect_equal(predictPL2(iE.survT, times = data[treatment=="C",time]),
                     outSurv$survTimeC[[iStrata]][,"survivalT_0"])
        expect_equal(predictPL2(iE.survT, times = data[treatment=="C",time]),
                     outSurv$survTimeC[[iStrata]][,"survivalT+threshold"])

        ## correct survival (treatment observation times)
        expect_equal(predictPL2(iE.survC, times = data[treatment=="T",time] - outBT@threshold),
                     outSurv$survTimeT[[iStrata]][,"survivalC-threshold"])
        expect_equal(predictPL2(iE.survC, times = data[treatment=="T",time]),
                     outSurv$survTimeT[[iStrata]][,"survivalC_0"])
        expect_equal(predictPL2(iE.survC, times = data[treatment=="T",time]),
                     outSurv$survTimeT[[iStrata]][,"survivalC+threshold"])

        expect_equal(predictPL2(iE.survT, times = data[treatment=="T",time] - outBT@threshold),
                     outSurv$survTimeT[[iStrata]][,"survivalT-threshold"])
        expect_equal(predictPL2(iE.survT, times = data[treatment=="T",time]),
                     outSurv$survTimeT[[iStrata]][,"survivalT_0"])
        expect_equal(predictPL2(iE.survT, times = data[treatment=="T",time]),
                     outSurv$survTimeT[[iStrata]][,"survivalT+threshold"])
    })
}



