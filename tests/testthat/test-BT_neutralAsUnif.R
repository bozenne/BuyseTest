verboseContext("Check neutralAsUnif")


Fpairs <- data.table(ID = 1:4,
                     time1 = c(10,10,10,10),
                     time2 = c(10,20,12,32),
                     status = c(1,1,1,1),
                     treat = c("Y","N","Y","N")) 


test_that("continue after neutral", {
  BT <- BuyseTest(treat ~ TTE(time1, 0, status) + TTE(time2, 0, status), data=Fpairs, n.bootstrap = 0)
  expect_equal(as.double(BT@count_neutral)[2],0)
})

test_that("stop at neutral", {
  BT <- BuyseTest(treat ~ TTE(time1, 0, status) + TTE(time2, 0, status), data=Fpairs, n.bootstrap = 0, neutralAsUninf = FALSE)
  expect_equal(as.double(BT@count_favorable+BT@count_unfavorable+BT@count_neutral+BT@count_uninf)[2],0)
})