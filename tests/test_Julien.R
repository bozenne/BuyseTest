#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Validation Julien
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(BuyseTest)
library(data.table)
library(testthat)

#### 1- test 2 pairs ####

## a)
gehan4 <- data.table(ID = 1:4,
                     pair = 1,
                     time = c(10,20,12,32),
                     cens = c(0,1,1,1),
                     treat = c("control","6-MP","control","6-MP")) 



test_that("2 pairs - Gehan",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,trace=0,method="Gehan")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),2)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),2)
 
  expect_equal(as.double(BT@delta),-1/2)
  
})

test_that("2 pairs - Peto",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,trace=0,method="Peto")
  
  expect_equal(as.double(BT@count_favorable),1/3)
  expect_equal(as.double(BT@count_unfavorable),3)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),2/3)

  expect_equal(as.double(BT@delta),-2/3)
  
})

test_that("2 pairs - Efron",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,trace=0,method="Efron")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),4)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),0)
  
  expect_equal(as.double(BT@delta),-1)
  
})

test_that("2 pairs - Peron",{
  BT <- BuyseTest(data=gehan4,endpoint="time",treatment="treat",
                  type="TTE",censoring="cens",threshold=0,n.bootstrap=0,trace=0,method="Peron")
  
  expect_equal(as.double(BT@count_favorable),0)
  expect_equal(as.double(BT@count_unfavorable),4)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),0)
  
  expect_equal(as.double(BT@delta),-1)
  
})

summary(survival::survfit(survival::Surv(time=time, event=cens) ~ 1, data=gehan4))


## b)
set.seed(10)
TpsFin<- 1 # 0.75
lambda.T <-0.5
n.Treatment <- 10
n.Control <- 10
n<-n.Treatment+n.Control
group<-c(rep(1, n.Treatment),rep(0, n.Control))

lambda.C <- 0.5
TimeEvent<-c(rexp(n.Treatment,rate=lambda.T),
             rexp(n.Control,rate=lambda.C))
Time.Cens<-runif(n,0,TpsFin)
Time<-pmin(Time.Cens,TimeEvent)
Event<-Time==TimeEvent
Event<-as.numeric(Event)
tab<-data.frame(group,Time,Event)


test_that("2 pairs - Gehan",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,trace=0,method="Gehan",seed=11)
  
  expect_equal(as.double(BT@count_favorable),6)
  expect_equal(as.double(BT@count_unfavorable),8)
  expect_equal(as.double(BT@count_neutral),0)
  expect_equal(as.double(BT@count_uninf),86)
  
  expect_equal(as.double(BT@delta),-0.02)
})

test_that("2 pairs - Peto",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,trace=0,method="Peto",seed=11)
  
  expect_equal(as.double(BT@count_favorable),40.95, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),41.67, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),17.38, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta),-0.00712, tolerance = 1e-3)
})

test_that("2 pairs - Efron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,trace=0,method="Efron",seed=11)
  
  expect_equal(as.double(BT@count_favorable),11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable),11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral),1, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf),76.78, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta),0, tolerance = 1e-3)
  
})

test_that("2 pairs - Peron",{
  BT <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                  type="TTE",censoring="Event",threshold=0.1,n.bootstrap=1,trace=0,method="Peron",seed=11)
 
  expect_equal(as.double(BT@count_favorable), 11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_unfavorable), 11.11, tolerance = 1e-3)
  expect_equal(as.double(BT@count_neutral), 0, tolerance = 1e-3)
  expect_equal(as.double(BT@count_uninf), 77.78, tolerance = 1e-3)
  
  expect_equal(as.double(BT@delta), 0, tolerance = 1e-3)
})

### lambda.C = 1
lambda.C <- 1
TimeEvent<-c(rexp(n.Treatment,rate=lambda.T),
             rexp(n.Control,rate=lambda.C))
Time.Cens<-runif(n,0,TpsFin)
Time<-pmin(Time.Cens,TimeEvent)
Event<-Time==TimeEvent
Event<-as.numeric(Event)
tab<-data.frame(group,Time,Event)


Buyseres.Gehan <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                            type="TTE",censoring="Event",threshold=0.1,n.bootstrap=0,trace=2,method="Gehan",seed=11)
Buyseres.Peto <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                           type="TTE",censoring="Event",threshold=0.1,n.bootstrap=0,trace=2,method="Peto",seed=11)
Buyseres.Efron <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                            type="TTE",censoring="Event",threshold=0.1,n.bootstrap=0,trace=2,method="Efron",seed=11)
Buyseres.Peron <- BuyseTest(data=tab,endpoint="Time",treatment="group",
                            type="TTE",censoring="Event",threshold=0.1,n.bootstrap=0,trace=2,method="Peron",seed=11)

summary(Buyseres.Gehan)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100            0              9          0       91 -0.09 -0.09           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100           19             16          4       61  0.03  0.03           0
summary(Buyseres.Peto)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100           36             42          0       22 -0.06 -0.06           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta  Delta n.bootstrap
# 1 global      100        41.43          38.44          4    16.13 0.0298 0.0298           0
summary(Buyseres.Efron)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100            0             55          1       44 -0.55 -0.55           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta  Delta n.bootstrap
# 1 global      100        27.51           67.3          4     1.19 -0.398 -0.398           0
summary(Buyseres.Peron)
# trata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100            0             10          0       90  -0.1  -0.1           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta  Delta n.bootstrap
# 1 global      100        27.51          23.74          4    44.75 0.0377 0.0377           0


#### check in detail ####
test_that("details - threshold",{
  expect_equal(1, BuyseTest:::initThreshold(threshold=1, type=3, D=1, method="Peto", endpoint="time"))
  expect_equal(1e-12, BuyseTest:::initThreshold(threshold=0, type=3, D=1, method="Peto", endpoint="time"))
})

test <- FALSE
if(test){ #### to be finished

SurvivalPeto_th1 <- BuyseTest:::initSurvival(M.Treatment=gehan4[gehan4$treat=="6-MP","time",drop=FALSE], 
                                         M.Control=gehan4[gehan4$treat=="control","time",drop=FALSE], 
                                         M.delta_Treatment=gehan4[gehan4$treat=="6-MP","cens",drop=FALSE], 
                                         M.delta_Control=gehan4[gehan4$treat=="control","cens",drop=FALSE], 
                                         endpoint="time",
                                         D.TTE=1, type=3, threshold=1, method="Peto")

SurvivalPeto_th0 <- BuyseTest:::initSurvival(M.Treatment=gehan4[gehan4$treat=="6-MP","time",drop=FALSE], 
                                        M.Control=gehan4[gehan4$treat=="control","time",drop=FALSE], 
                                        M.delta_Treatment=gehan4[gehan4$treat=="6-MP","cens",drop=FALSE], 
                                        M.delta_Control=gehan4[gehan4$treat=="control","cens",drop=FALSE], 
                                        endpoint="time",
                                        D.TTE=1, type=3, threshold=1, method="Peto")

resAllPairsPeto <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=gehan4$time[gehan4$treat=="6-MP"],
                                                    Control=gehan4$time[gehan4$treat=="control"],
                                                    deltaT=gehan4$cens[gehan4$treat=="6-MP"],
                                                    deltaC=gehan4$cens[gehan4$treat=="control"],
                                                    matKMT=SurvivalPeto_th0$list_survivalT[[1]], matKMC=SurvivalPeto_th0$list_survivalC[[1]],
                                                    type=1,threshold=0)
# type = 0 Gehan
# type = 1 Peto
# type = 2 Efron

SurvivalEfron_th1 <- BuyseTest:::initSurvival(M.Treatment=gehan4[gehan4$treat=="6-MP","time",drop=FALSE], 
                                             M.Control=gehan4[gehan4$treat=="control","time",drop=FALSE], 
                                             M.delta_Treatment=gehan4[gehan4$treat=="6-MP","cens",drop=FALSE], 
                                             M.delta_Control=gehan4[gehan4$treat=="control","cens",drop=FALSE], 
                                             endpoint="time",
                                             D.TTE=1, type=3, threshold=1, method="Efron", callMethod = "initSurvival")

SurvivalEfron_th0 <- BuyseTest:::initSurvival(M.Treatment=gehan4[gehan4$treat=="6-MP","time",drop=FALSE], 
                                             M.Control=gehan4[gehan4$treat=="control","time",drop=FALSE], 
                                             M.delta_Treatment=gehan4[gehan4$treat=="6-MP","cens",drop=FALSE], 
                                             M.delta_Control=gehan4[gehan4$treat=="control","cens",drop=FALSE], 
                                             endpoint="time",
                                             D.TTE=1, type=3, threshold=1, method="Efron", callMethod = "initSurvival")

resAllPairsEfron <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=gehan4$time[gehan4$treat=="6-MP"],
                                                     Control=gehan4$time[gehan4$treat=="control"],
                                                     deltaT=gehan4$cens[gehan4$treat=="6-MP"],
                                                     deltaC=gehan4$cens[gehan4$treat=="control"],
                                                     matKMT=SurvivalEfron_th0$list_survivalT[[1]], matKMC=SurvivalEfron_th0$list_survivalC[[1]],
                                                     type=2,threshold=0)


#### on one pair
Rcpp::sourceCpp(file.path(path,"test/test_calcOnePair.cpp"),rebuild=TRUE)

index_T <- 1
index_C <- 1

resOnePairPeto <-  Test_calcOnePair_TTEOutcome_cpp(endpoint_T=gehan4$time[gehan4$treat=="6-MP"][index_T],
                                                   endpoint_C=gehan4$time[gehan4$treat=="control"][index_C], 
                                                   delta_T=gehan4$cens[gehan4$treat=="6-MP"][index_T],
                                                   delta_C=gehan4$cens[gehan4$treat=="control"][index_C],
                                                   survival_T=rbind(SurvivalPeto_th0$list_survivalT[[1]][index_T,]),
                                                   survival_C=rbind(SurvivalPeto_th0$list_survivalC[[1]][index_C,]),
                                                   threshold=0, Wpair=1, iter_pair=c(0,0),type=1)

index_T <- 1
index_C <- 1

resOnePairPeto <-  Test_calcOnePair_TTEOutcome_cpp(endpoint_T=gehan4$time[gehan4$treat=="6-MP"][index_T],
                                                   endpoint_C=gehan4$time[gehan4$treat=="control"][index_C], 
                                                   delta_T=gehan4$cens[gehan4$treat=="6-MP"][index_T],
                                                   delta_C=gehan4$cens[gehan4$treat=="control"][index_C],
                                                   survival_T=rbind(SurvivalEfron_th0$list_survivalT[[1]]),
                                                   survival_C=rbind(SurvivalEfron_th0$list_survivalC[[1]]),
                                                   threshold=0, Wpair=1, iter_pair=c(index_T-1,index_C-1),type=2)




#### debug Peron
threshold <- 1
SurvivalPeto_th0 <- BuyseTest:::initSurvival(M.Treatment=tab[tab$group==1,"Time",drop=FALSE], 
                                              M.Control=tab[tab$group==0,"Time",drop=FALSE], 
                                              M.delta_Treatment=tab[tab$group==1,"Event",drop=FALSE], 
                                              M.delta_Control=tab[tab$group==0,"Event",drop=FALSE], 
                                              endpoint="Time",
                                              D.TTE=1, type=3, threshold=threshold, method="Peto", callMethod = "initSurvival")

plot(tab[tab$group==1,"Time"],tab[tab$group==1,"Event"],col="red",ylim=c(0,1.1))
points(tab[tab$group==0,"Time"],tab[tab$group==0,"Event"]+0.1,col="blue")
abline(v=0.11)



#### 3 another ####

x<-matrix(c(0.08,0.33,0.46,0.465,0.54,0.544,0.73,0.81,0.91,1.33 ,1,0,1,1,0,0,1,0,0,0 ,rep(1,10)),nrow=10,ncol=3,dimnames=list(1:10,c("Time","Event","group")))
y<-matrix(c(0.02,0.035,0.22,0.27,0.352,0.369,0.473,0.523,0.946,1.27 ,1,1,1,0,0,1,0,1,1,0 ,rep(0,10)),nrow=10,ncol=3,dimnames=list(11:20,c("Time","Event","group")))
tab<-rbind(x,y)
tab<-as.data.frame(tab)

seedGPC<-20150115
BuyseresE<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=1000,trace=0,method="Efron",seed=seedGPC)
summary(BuyseresE)
plot(survfit(Surv(time=Time, event=Event) ~ group, data=tab))
summary(survfit(Surv(time=Time, event=Event) ~ group, data=tab))


seedGPC<-20150115
BuyseresE<- BuyseTest(data=tab,endpoint="Time",treatment="group",
                      type="TTE",censoring="Event",threshold=0,n.bootstrap=1,trace=0,method="Efron",seed=seedGPC)
summary(BuyseresE)
}