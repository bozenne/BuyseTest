verboseContext("Check specific examples")

#### 1- test 2 pairs ####
cat("* check results with 2 pairs in each group \n")

## a)
gehan4 <- data.table(ID = 1:4,
                     pair = 1,
                     time = c(10,20,12,32),
                     cens = c(0,1,1,1),
                     treat = c("control","6-MP","control","6-MP")) 



test_that("2 pairs - Gehan",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,method="Gehan")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),2)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),2)
 
  expect_equal(as.double(BT@delta$netChance),-1/2)
  expect_equal(as.double(BT@delta$winRatio),0)
})

test_that("2 pairs - Peto",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,method="Peto")
  
  expect_equal(as.double(BT@count_favorable),1/3)
  expect_equal(as.double(BT@count_unfavorable),3)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),2/3)

  expect_equal(as.double(BT@delta$netChance),-2/3)
  expect_equal(as.double(BT@delta$winRatio),0.1)
})

test_that("2 pairs - Efron",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,method="Efron")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),4)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),0)
  
  expect_equal(as.double(BT@delta$netChance),-1)
  expect_equal(as.double(BT@delta$winRatio),0)
})

test_that("2 pairs - Peron",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,method="Peron")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),4)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),0)
  
  expect_equal(as.double(BT@delta$netChance),-1)
  expect_equal(as.double(BT@delta$winRatio),0)
  
})

summary(survival::survfit(survival::Surv(time=time, event=cens) ~ 1, data=gehan4))

#### 2- more patients? ####
cat("* check results with more patients \n")

## a) lambda.C = 0.5
set.seed(10)
TpsFin <- 1 # 0.75
lambda.T <- 0.5
n.Treatment <- 10
n.Control <- 10
n <- n.Treatment + n.Control
group <- c(rep(1, n.Treatment),rep(0, n.Control))

lambda.C <- 0.5
TimeEvent <- c(rexp(n.Treatment,rate=lambda.T),
             rexp(n.Control,rate=lambda.C))
Time.Cens <- runif(n,0,TpsFin)
Time <-pmin(Time.Cens,TimeEvent)
Event <- Time == TimeEvent
Event <- as.numeric(Event)
tab <- data.frame(group,Time,Event)


test_that("lambdaC = 0.5 - Gehan",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Gehan",seed=11)
  
  expect_equal(as.double(BT@count_favorable),6)
  expect_equal(as.double(BT@count_unfavorable),8)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),86)
  
  expect_equal(as.double(BT@delta$netChance),-0.02)
  expect_equal(as.double(BT@delta$winRatio),6/14)
})

test_that("lambdaC = 0.5 - Peto",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peto",seed=11)
  
  expect_equal(as.double(BT@count_favorable),40.95, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),41.67, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),17.38, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),-0.00712, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0.4956924, tolerance = 1e-3)
})

test_that("lambdaC = 0.5 - Efron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Efron",seed=11)
  
  expect_equal(as.double(BT@count_favorable),11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),1, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),76.78, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),0, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0.5, tolerance = 1e-3)
  
})

test_that("lambdaC = 0.5 - Peron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peron",seed=11)
 
  expect_equal(as.double(BT@count_favorable), 11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable), 11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral), 0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf), 77.78, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),0, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0.5, tolerance = 1e-3)
})

## b) lambda.C = 1
lambda.C <- 1
TimeEvent<-c(rexp(n.Treatment,rate=lambda.T),
             rexp(n.Control,rate=lambda.C))
Time.Cens<-runif(n,0,TpsFin)
Time<-pmin(Time.Cens,TimeEvent)
Event<-Time==TimeEvent
Event<-as.numeric(Event)
tab<-data.frame(group,Time,Event)

test_that("lambdaC = 1 - Gehan",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Gehan",seed=11)
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),9)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),91)
  
  expect_equal(as.double(BT@delta$netChance),-0.09)
  expect_equal(as.double(BT@delta$winRatio),0)
})

test_that("lambdaC = 1 - Peto",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peto",seed=11)
  
  expect_equal(as.double(BT@count_favorable),36, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),42, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),22, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),-0.06, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0.4615385, tolerance = 1e-3)
})

test_that("lambdaC = 1 - Efron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Efron",seed=11)
  expect_equal(as.double(BT@count_favorable),0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),55, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),1, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),44, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),-0.55, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0, tolerance = 1e-3)
  
})

test_that("lambdaC = 1 - Peron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,method="Peron",seed=11)
  
  expect_equal(as.double(BT@count_favorable), 0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable), 10, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral), 0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf), 90, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta$netChance),-0.1, tolerance = 1e-3)
  expect_equal(as.double(BT@delta$winRatio),0, tolerance = 1e-3)
})

#### 3- Threshold ####
test_that("details - threshold",{
  expect_equal(1, BuyseTest:::initThreshold(threshold=1, type=3, D=1, method="Peto", endpoint="time"))
  expect_equal(1e-12, BuyseTest:::initThreshold(threshold=0, type=3, D=1, method="Peto", endpoint="time"))
})
