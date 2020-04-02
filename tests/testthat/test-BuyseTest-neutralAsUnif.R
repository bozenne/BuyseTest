### test-BT_neutralAsUninf.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:50) 
## Version: 
## last-updated: apr  2 2020 (16:57) 
##           By: Brice Ozenne
##     Update #: 58
#----------------------------------------------------------------------
## 
### Commentary: Check whether the option neutral.as.uninf is working
## this option allows to stop the analysis of the neutral pairs
## instead of looking at endpoints with lower priority.
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(testthat)
    library(BuyseTest)
    library(data.table)
}

context("Check that the option neutral.as.uninf in BuyseTest is working correctly \n")

## * settings
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * generate data
## two survival endpoints and no censoring
## for the first endpoint all observations are neutral
## for the second endpoint all observations are unfavorable
dt.data <- data.table(ID = 1:4,
                      memory = c(10,10,10,10),
                      time = c(10,20,12,32),
                      status = c(1,1,1,1),
                      treat = c("Y","N","Y","N"))
n.data <- NROW(dt.data)

## same with some NA
dt.dataNA  <- copy(dt.data)
dt.dataNA[1,memory := NA]

## * neutral.as.uninf = TRUE (default)
## the neutral observations are analysed using the following endpoints
test_that("continue after NA (no NA)", {
    BT.TRUE <- BuyseTest(treat ~ cont(memory) + TTE(time, status, threshold = 0),
                         data = dt.data,
                         neutral.as.uninf = TRUE)
    
    expect_equal(as.double(coef(BT.TRUE, statistic = "count.favorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT.TRUE, statistic = "count.unfavorable", cumulative = FALSE)),c(0,4))
    expect_equal(as.double(coef(BT.TRUE, statistic = "count.neutral", cumulative = FALSE)),c(4,0))
    expect_equal(as.double(coef(BT.TRUE, statistic = "count.uninf", cumulative = FALSE)),c(0,0))
    
})

test_that("continue after NA (NA)", {
    BT.TRUE_NA <- BuyseTest(treat ~ Cont(memory) + TTE(time, status, threshold = 0),
                            data = dt.dataNA, scoring.rule = "Peron",
                            neutral.as.uninf = TRUE)
    ## summary(BT.TRUE_NA, percentage = FALSE)
    expect_equal(as.double(coef(BT.TRUE_NA, statistic = "count.favorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT.TRUE_NA, statistic = "count.unfavorable", cumulative = FALSE)),c(0,4))
    expect_equal(as.double(coef(BT.TRUE_NA, statistic = "count.neutral", cumulative = FALSE)),c(2,0))
    expect_equal(as.double(coef(BT.TRUE_NA, statistic = "count.uninf", cumulative = FALSE)),c(2,0))
    
})

## * neutral.as.uninf = FALSE
## the neutral observations are not analysed using the following endpoints
test_that("stop after NA (no NA)", {
    BT.FALSE <- BuyseTest(treat ~ Cont(memory) + TTE(time, status, threshold = 0),
                          data = dt.data,
                          neutral.as.uninf = FALSE)
    
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.favorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.unfavorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.neutral", cumulative = FALSE)),c(4,0))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.uninf", cumulative = FALSE)),c(0,0))
    
})

test_that("stop after NA (NA)", {
    BT.FALSE <- BuyseTest(treat ~ Cont(memory) + TTE(time, status, threshold = 0),
                          data = dt.dataNA,
                          neutral.as.uninf = FALSE)
    
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.favorable", cumulative = FALSE)),c(0,0))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.unfavorable", cumulative = FALSE)),c(0,2))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.neutral", cumulative = FALSE)),c(2,0))
    expect_equal(as.double(coef(BT.FALSE, statistic = "count.uninf", cumulative = FALSE)),c(2,0))
    
})

#----------------------------------------------------------------------
### test-BT_neutralAsUninf.R ends here
