### test-BT_neutralAsUninf.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:50) 
## Version: 
## last-updated: apr 15 2018 (14:22) 
##           By: Brice Ozenne
##     Update #: 29
#----------------------------------------------------------------------
## 
### Commentary: Check whether the option neutralAsUninf is working
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
}

context("Check that the option neutralAsUninf in BuyseTest is working correctly \n")

## * settings
BuyseTest.option(n.permutation = 0, trace = 0)

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


## * neutralAsUninf = TRUE (default)
## the neutral observations are analysed using the following endpoints
test_that("continue after NA (no NA)", {
    BT.TRUE <- BuyseTest(treat ~ cont(memory) + TTE(time, 0, status),
                         data = dt.data,
                         neutralAsUninf = TRUE)
    
    expect_equal(as.double(BT.TRUE@count_favorable),c(0,0))
    expect_equal(as.double(BT.TRUE@count_unfavorable),c(0,4))
    expect_equal(as.double(BT.TRUE@count_neutral),c(4,0))
    expect_equal(as.double(BT.TRUE@count_uninf),c(0,0))
    
})

test_that("continue after NA (NA)", {
    BT.TRUE_NA <- BuyseTest(treat ~ Cont(memory) + TTE(time, 0, status),
                            data = dt.dataNA,
                            neutralAsUninf = TRUE)

    expect_equal(as.double(BT.TRUE_NA@count_favorable),c(0,0))
    expect_equal(as.double(BT.TRUE_NA@count_unfavorable),c(0,4))
    expect_equal(as.double(BT.TRUE_NA@count_neutral),c(2,0))
    expect_equal(as.double(BT.TRUE_NA@count_uninf),c(2,0))
    
})

## * neutralAsUninf = FALSE
## the neutral observations are not analysed using the following endpoints
test_that("stop after NA (no NA)", {
    BT.FALSE <- BuyseTest(treat ~ Cont(memory) + TTE(time, 0, status),
                          data = dt.data,
                          neutralAsUninf = FALSE)
    
    expect_equal(as.double(BT.TRUE@count_favorable),c(0,0))
    expect_equal(as.double(BT.TRUE@count_unfavorable),c(0,4))
    expect_equal(as.double(BT.TRUE@count_neutral),c(4,0))
    expect_equal(as.double(BT.TRUE@count_uninf),c(0,0))
    
})

test_that("stop after NA (NA)", {
    BT.FALSE <- BuyseTest(treat ~ Cont(memory) + TTE(time, 0, status),
                          data = dt.dataNA,
                          neutralAsUninf = FALSE)
    
    expect_equal(as.double(BT.FALSE@count_favorable),c(0,0))
    expect_equal(as.double(BT.FALSE@count_unfavorable),c(0,2))
    expect_equal(as.double(BT.FALSE@count_neutral),c(2,0))
    expect_equal(as.double(BT.FALSE@count_uninf),c(2,0))
    
})

#----------------------------------------------------------------------
### test-BT_neutralAsUninf.R ends here
