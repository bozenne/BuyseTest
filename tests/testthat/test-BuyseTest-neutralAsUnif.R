### test-BT_neutralAsUninf.R --- 
#----------------------------------------------------------------------
## author: Brice
## created: maj 12 2017 (14:50) 
## Version: 
## last-updated: apr  1 2018 (11:43) 
##           By: Brice Ozenne
##     Update #: 17
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

verboseContext("Check that the option neutralAsUninf in BuyseTest is working correctly \n")

## * settings
BuyseTest.option(n.permutation = 0, trace = 0)

## * generate data
## two survival endpoints and no censoring
## for the first endpoint all observations are neutral
## for the second endpoint all observations are unfavorable
dt.data <- data.table(ID = 1:4,
                      time1 = c(10,10,10,10),
                      time2 = c(10,20,12,32),
                      status = c(1,1,1,1),
                      treat = c("Y","N","Y","N"))
n.data <- NROW(dt.data)

## same with some NA
dt.dataNA  <- copy(dt.data)
dt.dataNA[1,time1 := NA]


## * neutralAsUninf = TRUE (default)
## the neutral observations are analysed using the following endpoints
test_that("stop after neutral", {
    BT.TRUE <- BuyseTest(treat ~ TTE(time1, 0, status) + TTE(time2, 0, status),
                         data = dt.data,
                         neutralAsUninf = TRUE)
    
    expect_equal(as.double(BT.TRUE@count_favorable)[1],0)
    expect_equal(as.double(BT.TRUE@count_unfavorable)[1],0)
    expect_equal(as.double(BT.TRUE@count_neutral)[1],4)
    expect_equal(as.double(BT.TRUE@count_uninf)[1],0)
    
    expect_equal(as.double(BT.TRUE@count_favorable)[2],0)
    expect_equal(as.double(BT.TRUE@count_unfavorable)[2],4)
    expect_equal(as.double(BT.TRUE@count_neutral)[2],0)
    expect_equal(as.double(BT.TRUE@count_uninf)[2],0)
})

test_that("stop after neutral", {
    suppressWarnings(
        BT.TRUE_NA <- BuyseTest(treat ~ Continous(time1) + TTE(time2, 0, status),
                                data=dt.dataNA,
                                neutralAsUninf = TRUE)
    )
    expect_equal(as.double(BT.TRUE_NA@count_favorable)[1],0)
    expect_equal(as.double(BT.TRUE_NA@count_unfavorable)[1],0)
    expect_equal(as.double(BT.TRUE_NA@count_neutral)[1],2)
    expect_equal(as.double(BT.TRUE_NA@count_uninf)[1],2)
    
    expect_equal(as.double(BT.TRUE_NA@count_favorable)[2],0)
    expect_equal(as.double(BT.TRUE_NA@count_unfavorable)[2],2)
    expect_equal(as.double(BT.TRUE_NA@count_neutral)[2],0)
    expect_equal(as.double(BT.TRUE_NA@count_uninf)[2],0)
})

# }}}


# {{{ neutralAsUninf = FALSE
# the neutral observations are not analysed using the following endpoints
BT.FALSE <- BuyseTest(treat ~ TTE(time1, 0, status) + TTE(time2, 0, status),
                      data=dt.data,
                      neutralAsUninf = FALSE)
test_that("continue after neutral", {
    expect_equal(as.double(BT.FALSE@count_favorable)[1],0)
    expect_equal(as.double(BT.FALSE@count_unfavorable)[1],0)
    expect_equal(as.double(BT.FALSE@count_neutral)[1],4)
    expect_equal(as.double(BT.FALSE@count_uninf)[1],0)
    
    expect_equal(as.double(BT.FALSE@count_favorable)[2],0)
    expect_equal(as.double(BT.FALSE@count_unfavorable)[2],0)
    expect_equal(as.double(BT.FALSE@count_neutral)[2],0)
    expect_equal(as.double(BT.FALSE@count_uninf)[2],0)
})

suppressWarnings(
BT.FALSE_NA <- BuyseTest(treat ~ Continous(time1) + TTE(time2, 0, status),
                         data=dt.dataNA,
                         neutralAsUninf = FALSE)
)
test_that("stop after neutral", {
    expect_equal(as.double(BT.FALSE_NA@count_favorable)[1],0)
    expect_equal(as.double(BT.FALSE_NA@count_unfavorable)[1],0)
    expect_equal(as.double(BT.FALSE_NA@count_neutral)[1],2)
    expect_equal(as.double(BT.FALSE_NA@count_uninf)[1],2)
    
    expect_equal(as.double(BT.FALSE_NA@count_favorable)[2],0)
    expect_equal(as.double(BT.FALSE_NA@count_unfavorable)[2],0)
    expect_equal(as.double(BT.FALSE_NA@count_neutral)[2],0)
    expect_equal(as.double(BT.FALSE_NA@count_uninf)[2],0)
})

# }}}


#----------------------------------------------------------------------
### test-BT_neutralAsUninf.R ends here
