#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Simulation et test sur donneees par pairs : donnees binaires
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path <- "E:/Creation_package/Package_BuyseTest/BuyseTest" # path to the uncompressed tar.gz file

Rcpp::sourceCpp(file.path(path,"test/test_calcAllPairs.cpp"),rebuild=TRUE)

#### 1- Binary oucome ####

#### base ####
data_testBin0 <- data.frame(treatment=c(0,1),control=c(0,1))

resC0_testBin0 <- Test_calcAllPairs_BinaryOutcome_cpp(Treatment=data_testBin0$treatment,
                                                     Control=data_testBin0$control)
# count  (favorable / unfavorable / neutral / uninf) : 1 1 2 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 2 2 0 0
# # index_neutralT : 0 1 
# # index_neutralC : 0 1 
# index size (wNeutral/ wUninf) : 0 0

resC1_testBin0 <- Test_calcAllPairs_BinaryOutcome_cpp(Treatment=c(NA,data_testBin0$treatment[-1]),
                                                      Control=data_testBin0$control)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 1 2
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 2 2
# # index_neutralT : 1 
# # index_neutralC : 1 
# # index_uninfT : 0 0 
# # index_uninfC : 0 1 
# index size (wNeutral/ wUninf) : 0 0

resC2_testBin0 <- Test_calcAllPairs_BinaryOutcome_cpp(Treatment=c(data_testBin0$treatment[-2],NA),
                                                      Control=data_testBin0$control)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 1 2
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 2 2
# # index_neutralT : 0 
# # index_neutralC : 0 
# # index_uninfT : 1 1 
# # index_uninfC : 0 1 
# index size (wNeutral/ wUninf) : 0 0

#### large ####
n.Treatment_testBin1 <- 50
n.Control_testBin1 <- 50
prob.Treatment_testBin1 <- c(0.5,0.75)
prob.Control_testBin1 <- c(0.5,0.25)

set.seed(10)
data_testBin1 <- data.frame(treatment=c(rep(1,n.Treatment_testBin1),rep(0,n.Treatment_testBin1)))
data_testBin1$endpoint1 <- c(rbinom(n.Treatment_testBin1,size=1,prob=prob.Treatment_testBin1[1]),
                            rbinom(n.Control_testBin1,size=1,prob=prob.Control_testBin1[1]))
data_testBin1$endpoint2 <- c(rbinom(n.Control_testBin1,size=1,prob=prob.Treatment_testBin1[2]),
                            rbinom(n.Control_testBin1,size=1,prob=prob.Control_testBin1[2]))


resC1_testBin1 <- Test_calcAllPairs_BinaryOutcome_cpp(Treatment=data_testBin1$endpoint1[data_testBin1$treatment==1],
                                                     Control=data_testBin1$endpoint1[data_testBin1$treatment==0],
                                                     trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 660 560 1280 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1280 1280 0 0
# index size (wNeutral/ wUninf) : 0 0

#### subset
resC1_testBin2 <- Test_calcSubsetPairs_BinaryOutcome_cpp(Treatment=data_testBin1$endpoint1[data_testBin1$treatment==1],
                                                         Control=data_testBin1$endpoint1[data_testBin1$treatment==0],
                                                         index_neutralT=resC1_testBin1$index_neutralT, index_neutralC=resC1_testBin1$index_neutralC, nNeutral_pairs=length(resC1_testBin1$index_neutralC),
                                                         index_uninfT=resC1_testBin1$index_uninfT, index_uninfC=resC1_testBin1$index_uninfC, nUninf_pairs=length(resC1_testBin1$index_uninfC),
                                                         trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1280 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1280 1280 0 0
# index size (wNeutral/ wUninf) : 0 0
quantile(resC1_testBin2$index_neutralT==resC1_testBin1$index_neutralT)
quantile(resC1_testBin2$index_neutralC==resC1_testBin1$index_neutralC)


resC2_testBin2 <- Test_calcSubsetPairs_BinaryOutcome_cpp(Treatment=data_testBin1$endpoint1[data_testBin1$treatment==1],
                                                         Control=data_testBin1$endpoint1[data_testBin1$treatment==0],
                                                         index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                         index_uninfT=expand.grid(0:(n.Treatment_testBin1-1),0:(n.Control_testBin1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testBin1-1),0:(n.Control_testBin1-1))[,2], nUninf_pairs=n.Treatment_testBin1*n.Control_testBin1,
                                                         trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 660 560 1280 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1280 1280 0 0
# index size (wNeutral/ wUninf) : 0 0
quantile(sort(resC2_testBin2$index_neutralT)==sort(resC1_testBin1$index_neutralT))
quantile(sort(resC2_testBin2$index_neutralC)==sort(resC1_testBin1$index_neutralC))

resC3_testBin2 <- Test_calcSubsetPairs_BinaryOutcome_cpp(Treatment=data_testBin1$endpoint2[data_testBin1$treatment==1],
                                                         Control=data_testBin1$endpoint2[data_testBin1$treatment==0],
                                                         index_neutralT=resC1_testBin1$index_neutralT, index_neutralC=resC1_testBin1$index_neutralC, nNeutral_pairs=length(resC1_testBin1$index_neutralC),
                                                         index_uninfT=resC1_testBin1$index_uninfT, index_uninfC=resC1_testBin1$index_uninfC, nUninf_pairs=length(resC1_testBin1$index_uninfC),
                                                         trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 725 75 480 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 480 480 0 0
# index size (wNeutral/ wUninf) : 0 0

#### subset weighted
resC1_testBin3 <- Test_calcSubsetPairs_WeightedBinaryOutcome_cpp(Treatment=data_testBin1$endpoint1[data_testBin1$treatment==1],
                                                         Control=data_testBin1$endpoint1[data_testBin1$treatment==0],
                                                         index_neutralT=resC1_testBin1$index_neutralT, index_neutralC=resC1_testBin1$index_neutralC, nNeutral_pairs=length(resC1_testBin1$index_neutralC),
                                                         index_uninfT=resC1_testBin1$index_uninfT, index_uninfC=resC1_testBin1$index_uninfC, nUninf_pairs=length(resC1_testBin1$index_uninfC),
                                                         Wpairs=rep(0.3,length(resC1_testBin1$index_neutralT+length(resC1_testBin1$index_uninfT))),trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 384 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1280 1280 0 0
# index size (wNeutral/ wUninf) : 0 1280
quantile(resC1_testBin3$index_neutralT==resC1_testBin1$index_neutralT)
quantile(resC1_testBin3$index_neutralC==resC1_testBin1$index_neutralC)

resC2_testBin3 <- Test_calcSubsetPairs_WeightedBinaryOutcome_cpp(Treatment=data_testBin1$endpoint1[data_testBin1$treatment==1],
                                                         Control=data_testBin1$endpoint1[data_testBin1$treatment==0],
                                                         index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                         index_uninfT=expand.grid(0:(n.Treatment_testBin1-1),0:(n.Control_testBin1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testBin1-1),0:(n.Control_testBin1-1))[,2], nUninf_pairs=n.Treatment_testBin1*n.Control_testBin1,
                                                         Wpairs=rep(0.3,n.Treatment_testBin1*n.Control_testBin1),trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 198=(660*0.3) 168=(560*0.3) 384 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1280 1280 0 0
# index size (wNeutral/ wUninf) : 0 1280
quantile(sort(resC2_testBin3$index_neutralT)==sort(resC1_testBin1$index_neutralT))
quantile(sort(resC2_testBin3$index_neutralC)==sort(resC1_testBin1$index_neutralC))

resC3_testBin3 <- Test_calcSubsetPairs_WeightedBinaryOutcome_cpp(Treatment=data_testBin1$endpoint2[data_testBin1$treatment==1],
                                                                 Control=data_testBin1$endpoint2[data_testBin1$treatment==0],
                                                                 index_neutralT=resC1_testBin1$index_neutralT, index_neutralC=resC1_testBin1$index_neutralC, nNeutral_pairs=length(resC1_testBin1$index_neutralC),
                                                                 index_uninfT=resC1_testBin1$index_uninfT, index_uninfC=resC1_testBin1$index_uninfC, nUninf_pairs=length(resC1_testBin1$index_uninfC),
                                                                 Wpairs=rep(0.3,length(resC1_testBin1$index_neutralT+length(resC1_testBin1$index_uninfT))),trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 217.5=(725*0.3) 22.5=(75*0.3) 144 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 480 480 0 0
# index size (wNeutral/ wUninf) : 0 480

#### 2- Continuous oucome ####

#### base ####
data_testCont0 <- data.frame(treatment=c(-0.5,0,0.5),control=c(-0.5,0,0.5))

resC0_testCont0 <- Test_calcAllPairs_ContinuousOutcome_cpp(Treatment=data_testCont0$treatment,
                                                      Control=data_testCont0$control,
                                                      threshold=1)
# count  (favorable / unfavorable / neutral / uninf) : 1 1 7 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 7 7 0 0
# # index_neutralT : 0 0 1 1 1 2 2 
# # index_neutralC : 0 1 0 1 2 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC1_testCont0 <- Test_calcAllPairs_ContinuousOutcome_cpp(Treatment=data_testCont0$treatment,
                                                           Control=data_testCont0$control,
                                                           threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 3 3 3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 3 3 0 0
# # index_neutralT : 0 1 2 
# # index_neutralC : 0 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC1_testCont0 <- Test_calcAllPairs_ContinuousOutcome_cpp(Treatment=c(NA,data_testCont0$treatment[-1]),
                                                      Control=data_testCont0$control,
                                                      threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 3 1 2 3
# index size (neutralT / neutralC / uninfT/ uninfC) : 2 2 3 3
# # index_neutralT : 1 2 
# # index_neutralC : 1 2 
# # index_uninfT : 0 0 0 
# # index_uninfC : 0 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC2_testCont0 <- Test_calcAllPairs_ContinuousOutcome_cpp(Treatment=c(data_testCont0$treatment[-2],NA),
                                                      Control=data_testCont0$control,
                                                      threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 2 2 2 3
# index size (neutralT / neutralC / uninfT/ uninfC) : 2 2 3 3
# # index_neutralT : 0 1 
# # index_neutralC : 0 2 
# # index_uninfT : 2 2 2 
# # index_uninfC : 0 1 2 
# index size (wNeutral/ wUninf) : 0 0


#### large ####
n.Treatment_testCont1 <- 50
n.Control_testCont1 <- 50
mu.Treatment_testCont1 <- c(2,4)
mu.Control_testCont1 <- c(2,0)

set.seed(10)
data_testCont1 <- data.frame(treatment=c(rep(1,n.Treatment_testCont1),rep(0,n.Treatment_testCont1)))
data_testCont1$endpoint1 <- c(rnorm(n.Treatment_testCont1,mean=mu.Treatment_testCont1[1]),
                              rnorm(n.Control_testCont1,mean=mu.Control_testCont1[1]))
data_testCont1$endpoint2 <- c(rnorm(n.Control_testCont1,mean=mu.Treatment_testCont1[2]),
                              rnorm(n.Control_testCont1,mean=mu.Control_testCont1[2]))


resC1_testCont1 <- Test_calcAllPairs_ContinuousOutcome_cpp(Treatment=data_testCont1$endpoint1[data_testCont1$treatment==1],
                                                      Control=data_testCont1$endpoint1[data_testCont1$treatment==0],
                                                      threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 367 816 1317 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1317 1317 0 0
# index size (wNeutral/ wUninf) : 0 0

#### subset
resC1_testCont2 <- Test_calcSubsetPairs_ContinuousOutcome_cpp(Treatment=data_testCont1$endpoint1[data_testCont1$treatment==1],
                                                         Control=data_testCont1$endpoint1[data_testCont1$treatment==0],
                                                         index_neutralT=resC1_testCont1$index_neutralT, index_neutralC=resC1_testCont1$index_neutralC, nNeutral_pairs=length(resC1_testCont1$index_neutralC),
                                                         index_uninfT=resC1_testCont1$index_uninfT, index_uninfC=resC1_testCont1$index_uninfC, nUninf_pairs=length(resC1_testCont1$index_uninfC),
                                                         threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1317 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1317 1317
# index size (wNeutral/ wUninf) : 0 0
quantile(resC1_testCont2$index_neutralT==resC1_testCont1$index_neutralT)
quantile(resC1_testCont2$index_neutralC==resC1_testCont1$index_neutralC)

resC2_testCont2 <- Test_calcSubsetPairs_ContinuousOutcome_cpp(Treatment=data_testCont1$endpoint1[data_testCont1$treatment==1],
                                                         Control=data_testCont1$endpoint1[data_testCont1$treatment==0],
                                                         index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                         index_uninfT=expand.grid(0:(n.Treatment_testCont1-1),0:(n.Control_testCont1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testCont1-1),0:(n.Control_testCont1-1))[,2], nUninf_pairs=n.Treatment_testCont1*n.Control_testCont1,
                                                         threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 367 816 1317 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1317 1317 0 0
# index size (wNeutral/ wUninf) : 0 0
quantile(sort(resC2_testCont2$index_neutralT)==sort(resC1_testCont1$index_neutralT))
quantile(sort(resC2_testCont2$index_neutralC)==sort(resC1_testCont1$index_neutralC))

resC3_testCont2 <- Test_calcSubsetPairs_ContinuousOutcome_cpp(Treatment=data_testCont1$endpoint2[data_testCont1$treatment==1],
                                                         Control=data_testCont1$endpoint2[data_testCont1$treatment==0],
                                                         index_neutralT=resC1_testCont1$index_neutralT, index_neutralC=resC1_testCont1$index_neutralC, nNeutral_pairs=length(resC1_testCont1$index_neutralC),
                                                         index_uninfT=resC1_testCont1$index_uninfT, index_uninfC=resC1_testCont1$index_uninfC, nUninf_pairs=length(resC1_testCont1$index_uninfC),
                                                         threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 1306 0 11 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 11 11 0 0
# index size (wNeutral/ wUninf) : 0 0

#### subset weighted
resC1_testCont3 <- Test_calcSubsetPairs_WeightedContinuousOutcome_cpp(Treatment=data_testCont1$endpoint1[data_testCont1$treatment==1],
                                                                 Control=data_testCont1$endpoint1[data_testCont1$treatment==0],
                                                                 index_neutralT=resC1_testCont1$index_neutralT, index_neutralC=resC1_testCont1$index_neutralC, nNeutral_pairs=length(resC1_testCont1$index_neutralC),
                                                                 index_uninfT=resC1_testCont1$index_uninfT, index_uninfC=resC1_testCont1$index_uninfC, nUninf_pairs=length(resC1_testCont1$index_uninfC),
                                                                 Wpairs=rep(0.3,length(resC1_testCont1$index_neutralT+length(resC1_testCont1$index_uninfT))),threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 395.1=(1317*0.3) 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1317 1317 0 0
# index size (wNeutral/ wUninf) : 0 1317
quantile(resC1_testCont3$index_neutralT==resC1_testCont1$index_neutralT)
quantile(resC1_testCont3$index_neutralC==resC1_testCont1$index_neutralC)

resC2_testCont3 <- Test_calcSubsetPairs_WeightedContinuousOutcome_cpp(Treatment=data_testCont1$endpoint1[data_testCont1$treatment==1],
                                                                 Control=data_testCont1$endpoint1[data_testCont1$treatment==0],
                                                                 index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                                 index_uninfT=expand.grid(0:(n.Treatment_testCont1-1),0:(n.Control_testCont1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testCont1-1),0:(n.Control_testCont1-1))[,2], nUninf_pairs=n.Treatment_testCont1*n.Control_testCont1,
                                                                 Wpairs=rep(0.3,n.Treatment_testCont1*n.Control_testCont1),threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 110.1=(367*0.3) 244.8=(816*0.3) 395.1=(1317*0.3) 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1317 1317 0 0
# index size (wNeutral/ wUninf) : 0 1317
quantile(sort(resC2_testCont3$index_neutralT)==sort(resC1_testCont1$index_neutralT))
quantile(sort(resC2_testCont3$index_neutralC)==sort(resC1_testCont1$index_neutralC))

resC3_testCont3 <- Test_calcSubsetPairs_WeightedContinuousOutcome_cpp(Treatment=data_testCont1$endpoint2[data_testCont1$treatment==1],
                                                                 Control=data_testCont1$endpoint2[data_testCont1$treatment==0],
                                                                 index_neutralT=resC1_testCont1$index_neutralT, index_neutralC=resC1_testCont1$index_neutralC, nNeutral_pairs=length(resC1_testCont1$index_neutralC),
                                                                 index_uninfT=resC1_testCont1$index_uninfT, index_uninfC=resC1_testCont1$index_uninfC, nUninf_pairs=length(resC1_testCont1$index_uninfC),
                                                                 Wpairs=rep(0.3,length(resC1_testCont1$index_neutralT+length(resC1_testCont1$index_uninfT))),threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 391.8=(1306*0.3) 0 3.3=(11*0.3) 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 11 11 0 0
# index size (wNeutral/ wUninf) : 0 11

#### 3- TTE oucome (Gehan) ####

#### base ####
data_testTTE_Gehan0 <- data.frame(treatment=c(1,2,2.5),control=c(1,2,2.5))

resC0_testTTE_Gehan0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan0$treatment,
                                                    Control=data_testTTE_Gehan0$control,
                                                    deltaT=c(1,1,1,1),deltaC=c(1,1,1,1), matKMT=matrix(), matKMC=matrix(),type=0,
                                                    threshold=1)
# count  (favorable / unfavorable / neutral / uninf) : 2 2 5 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 5 5 0 0
# # index_neutralT : 0 1 1 2 2 
# # index_neutralC : 0 1 2 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC1_testTTE_Gehan0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan0$treatment,
                                                    Control=data_testTTE_Gehan0$control,
                                                    deltaT=c(1,1,1,1),deltaC=c(1,1,1,1), matKMT=matrix(), matKMC=matrix(),type=0,
                                                    threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 3 3 3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 3 3 0 0
# # index_neutralT : 0 1 2 
# # index_neutralC : 0 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC2_testTTE_Gehan0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan0$treatment,
                                                    Control=data_testTTE_Gehan0$control,
                                                    deltaT=c(1,1,1,1),deltaC=c(0,0,0,0), matKMT=matrix(), matKMC=matrix(),type=0,
                                                    threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 3 0 6
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 6 6
# # index_uninfT : 0 1 1 2 2 2 
# # index_uninfC : 0 0 1 0 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC3_testTTE_Gehan0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan0$treatment,
                                                    Control=data_testTTE_Gehan0$control,
                                                    deltaT=c(0,0,0,0),deltaC=c(1,1,1,1), matKMT=matrix(), matKMC=matrix(),type=0,
                                                    threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 3 0 0 6
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 6 6
# # index_uninfT : 0 0 0 1 1 2 
# # index_uninfC : 0 1 2 1 2 2 
# index size (wNeutral/ wUninf) : 0 0

resC4_testTTE_Gehan0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan0$treatment,
                                                    Control=data_testTTE_Gehan0$control,
                                                    deltaT=c(0,0,0,0),deltaC=c(0,0,0,0), matKMT=matrix(), matKMC=matrix(),type=0,
                                                    threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 9
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 9 9
# # index_uninfT : 0 0 0 1 1 1 2 2 2 
# # index_uninfC : 0 1 2 0 1 2 0 1 2 
# index size (wNeutral/ wUninf) : 0 0


#### large ####
n.Treatment_testTTE1 <- 50
n.Control_testTTE1 <- 50
lambda.Treatment_testTTE1 <- c(0.75,0.5)
lambda.Control_testTTE1 <- c(0.75,5)
lambda.Censoring_testTTE1 <- c(0.5,0.5)

set.seed(10)
data_testTTE_Gehan1 <- data.frame(treatment=c(rep(1,n.Treatment_testTTE1),rep(0,n.Treatment_testTTE1)))
data_testTTE_Gehan1$EventTime1 <- c(rexp(n.Treatment_testTTE1,rate=lambda.Treatment_testTTE1[1]),
                              rexp(n.Control_testTTE1,rate=lambda.Control_testTTE1[1]))
data_testTTE_Gehan1$EventTime2 <- c(rexp(n.Control_testTTE1,rate=lambda.Treatment_testTTE1[2]),
                              rexp(n.Control_testTTE1,rate=lambda.Control_testTTE1[2]))
data_testTTE_Gehan1$CensoringTime1 <- c(rexp(n.Treatment_testTTE1,rate=lambda.Censoring_testTTE1[1]),
                                  rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[1]))
data_testTTE_Gehan1$CensoringTime2 <- c(rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[2]),
                                  rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[2]))

data_testTTE_Gehan1$endpoint1 <- apply(data_testTTE_Gehan1[,c("EventTime1","CensoringTime1")],1,min)
data_testTTE_Gehan1$event1 <- as.numeric(apply(data_testTTE_Gehan1[,c("EventTime1","CensoringTime1")],1,which.min)==1)
data_testTTE_Gehan1$endpoint2 <- apply(data_testTTE_Gehan1[,c("EventTime2","CensoringTime2")],1,min)
data_testTTE_Gehan1$event2 <- as.numeric(apply(data_testTTE_Gehan1[,c("EventTime2","CensoringTime2")],1,which.min)==1)

resC1_testTTE_Gehan1 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==1],
                                                         Control=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==0],
                                                         deltaT=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==1],
                                                         deltaC=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==0],
                                                         matKMT=matrix(), matKMC=matrix(),type=0,
                                                         threshold=0.25,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 487 624 208 1181
# index size (neutralT / neutralC / uninfT/ uninfC) : 208 208 1181 1181
# index size (wNeutral/ wUninf) : 0 0

#### subset
resC1_testTTE_Gehan2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==1],
                                                         Control=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==0],threshold=0.25,
                                                         deltaT=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==1],
                                                         deltaC=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==0],
                                                         matKMT=matrix(), matKMC=matrix(),
                                                         index_neutralT=resC1_testTTE_Gehan1$index_neutralT, index_neutralC=resC1_testTTE_Gehan1$index_neutralC, nNeutral_pairs=length(resC1_testTTE_Gehan1$index_neutralC),
                                                         index_uninfT=resC1_testTTE_Gehan1$index_uninfT, index_uninfC=resC1_testTTE_Gehan1$index_uninfC, nUninf_pairs=length(resC1_testTTE_Gehan1$index_uninfC),
                                                         Wpairs=numeric(0),threshold_M1=0,matKMT_M1=matrix(), matKMC_M1=matrix(),type=0,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 208 1181
# index size (neutralT / neutralC / uninfT/ uninfC) : 208 208 1181 1181
# index size (wNeutral/ wUninf) : 0 0

resC2_testTTE_Gehan2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==1],
                                                            Control=data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==0],threshold=0.25,
                                                            deltaT=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==1],
                                                            deltaC=data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==0],
                                                            matKMT=matrix(), matKMC=matrix(),
                                                            index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                            index_uninfT=expand.grid(0:(n.Treatment_testTTE1-1),0:(n.Control_testTTE1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testTTE1-1),0:(n.Control_testTTE1-1))[,2], nUninf_pairs=n.Treatment_testTTE1*n.Control_testTTE1,
                                                            Wpairs=numeric(0),threshold_M1=0,matKMT_M1=matrix(), matKMC_M1=matrix(),type=0,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 487 624 208 1181
# index size (neutralT / neutralC / uninfT/ uninfC) : 208 208 1181 1181
# index size (wNeutral/ wUninf) : 0 0

resC3_testTTE_Gehan2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Gehan1$endpoint2[data_testTTE_Gehan1$treatment==1],
                                                            Control=data_testTTE_Gehan1$endpoint2[data_testTTE_Gehan1$treatment==0],threshold=0.25,
                                                            deltaT=data_testTTE_Gehan1$event2[data_testTTE_Gehan1$treatment==1],
                                                            deltaC=data_testTTE_Gehan1$event2[data_testTTE_Gehan1$treatment==0],
                                                            matKMT=matrix(), matKMC=matrix(),
                                                            index_neutralT=resC1_testTTE_Gehan1$index_neutralT, index_neutralC=resC1_testTTE_Gehan1$index_neutralC, nNeutral_pairs=length(resC1_testTTE_Gehan1$index_neutralC),
                                                            index_uninfT=resC1_testTTE_Gehan1$index_uninfT, index_uninfC=resC1_testTTE_Gehan1$index_uninfC, nUninf_pairs=length(resC1_testTTE_Gehan1$index_uninfC),
                                                            Wpairs=numeric(0),threshold_M1=0,matKMT_M1=matrix(), matKMC_M1=matrix(),type=0,trace=1)

# count  (favorable / unfavorable / neutral / uninf) : 790 24 160 415
# index size (neutralT / neutralC / uninfT/ uninfC) : 160 160 415 415
# index size (wNeutral/ wUninf) : 0 0


#### 4- TTE oucome (Peto) ####
data("veteran",package="survival")

#### base ####
data_testTTE_Peto0 <- data.frame(treatment=c(1,2,2.5),control=c(1,2,2.5))
resKMT_tempo <- survival::survfit(survival::Surv(veteran$time/50,veteran$trt)~1) # compute the survival over the controls with common Kaplan Meier estimator.  
survestKMT_tempo <- stats::stepfun(c(resKMT_tempo$time,max(veteran$time/50)), c(1,resKMT_tempo$surv,NA)) # step interpolation of the survival function

listKMT_Peto0 <- list(cbind(survestKMT_tempo(c(1,2,2.5)-1),survestKMT_tempo(c(1,2,2.5)),survestKMT_tempo(c(1,2,2.5)+1)))
listKMC_Peto0 <- list(cbind(survestKMT_tempo(c(1,2,2.5)-1),survestKMT_tempo(c(1,2,2.5)),survestKMT_tempo(c(1,2,2.5)+1)))

resC0_testTTE_Peto0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto0$treatment,
                                                         Control=data_testTTE_Peto0$control,
                                                         deltaT=c(1,1,1,1),deltaC=c(1,1,1,1), matKMT=listKMT_Peto0[[1]], matKMC=listKMC_Peto0[[1]],type=1,
                                                         threshold=1)

# count  (favorable / unfavorable / neutral / uninf) : 2 2 5 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 5 5 0 0
# # index_neutralT : 0 1 1 2 2 
# # index_neutralC : 0 1 2 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC1_testTTE_Peto0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto0$treatment,
                                                         Control=data_testTTE_Peto0$control,
                                                         deltaT=c(1,1,1,1),deltaC=c(1,1,1,1), matKMT=listKMT_Peto0[[1]], matKMC=listKMC_Peto0[[1]],type=1,
                                                         threshold=0)
# count  (favorable / unfavorable / neutral / uninf) : 3 3 3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 3 3 0 0
# # index_neutralT : 0 1 2 
# # index_neutralC : 0 1 2 
# index size (wNeutral/ wUninf) : 0 0

resC2_testTTE_Peto0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto0$treatment,
                                                         Control=data_testTTE_Peto0$control,
                                                         deltaT=c(1,1,1,1),deltaC=c(0,0,0,0), matKMT=listKMT_Peto0[[1]], matKMC=listKMC_Peto0[[1]],type=1,
                                                         threshold=1)

# count  (favorable / unfavorable / neutral / uninf) : 0.0862373 7.69315 0 1.22061
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 7 7
# # index_uninfT : 0 1 1 1 2 2 2 
# # index_uninfC : 0 0 1 2 0 1 2 
# size (weights/ index_weights) : 7 0
# # weights : 0.221903 0.320573 0.12681 0.0519751 0.257764 0.15692 0.0846656


resC4_testTTE_Peto0 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto0$treatment,
                                                        Control=data_testTTE_Peto0$control,
                                                        deltaT=c(0,0,0,0),deltaC=c(0,0,0,0), matKMT=listKMT_Peto0[[1]], matKMC=listKMC_Peto0[[1]],type=1,
                                                        threshold=1)
# count  (favorable / unfavorable / neutral / uninf) : 3.88969 3.88969 0 1.22061
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 9 9
# # index_uninfT : 0 0 0 1 1 1 2 2 2 
# # index_uninfC : 0 1 2 0 1 2 0 1 2 
# size (weights/ index_weights) : 9 0
# # weights : 0.221903 0.160287 0.128882 0.160287 0.12681 0.104447 0.128882 0.104447 0.0846656 

#### large ####
n.Treatment_testTTE1 <- 50
n.Control_testTTE1 <- 50
lambda.Treatment_testTTE1 <- c(0.75,0.5)
lambda.Control_testTTE1 <- c(0.75,5)
lambda.Censoring_testTTE1 <- c(0.5,0.5)

set.seed(10)
data_testTTE_Peto1 <- data.frame(treatment=c(rep(1,n.Treatment_testTTE1),rep(0,n.Treatment_testTTE1)))
data_testTTE_Peto1$EventTime1 <- c(rexp(n.Treatment_testTTE1,rate=lambda.Treatment_testTTE1[1]),
                                    rexp(n.Control_testTTE1,rate=lambda.Control_testTTE1[1]))
data_testTTE_Peto1$EventTime2 <- c(rexp(n.Control_testTTE1,rate=lambda.Treatment_testTTE1[2]),
                                    rexp(n.Control_testTTE1,rate=lambda.Control_testTTE1[2]))
data_testTTE_Peto1$CensoringTime1 <- c(rexp(n.Treatment_testTTE1,rate=lambda.Censoring_testTTE1[1]),
                                        rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[1]))
data_testTTE_Peto1$CensoringTime2 <- c(rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[2]),
                                        rexp(n.Control_testTTE1,rate=lambda.Censoring_testTTE1[2]))

data_testTTE_Peto1$endpoint1 <- apply(data_testTTE_Peto1[,c("EventTime1","CensoringTime1")],1,min)
data_testTTE_Peto1$event1 <- as.numeric(apply(data_testTTE_Peto1[,c("EventTime1","CensoringTime1")],1,which.min)==1)
data_testTTE_Peto1$endpoint2 <- apply(data_testTTE_Peto1[,c("EventTime2","CensoringTime2")],1,min)
data_testTTE_Peto1$event2 <- as.numeric(apply(data_testTTE_Peto1[,c("EventTime2","CensoringTime2")],1,which.min)==1)

listKMT_Peto1 <- list()
listKMC_Peto1 <- list()
threshold_Peto1 <- c(1,1)

for(iter_endpoint in 1:2){
  resKM_tempo <- survival::survfit(survival::Surv(data_testTTE_Peto1[,c("endpoint1","endpoint2")[iter_endpoint]],
                              data_testTTE_Peto1[,c("event1","event2")[iter_endpoint]])~1)
  survestKM_tempo <- stats::stepfun(resKM_tempo$time, c(1, resKM_tempo$surv))
  
  time_treatment <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("endpoint1","endpoint2")[iter_endpoint]]
  listKMT_Peto1[[iter_endpoint]]  <- cbind(survestKM_tempo(time_treatment-threshold_Peto1[iter_endpoint]),
                                     survestKM_tempo(time_treatment),
                                     survestKM_tempo(time_treatment+threshold_Peto1[iter_endpoint])
  )
  rownames(listKMT_Peto1[[iter_endpoint]]) <- time_treatment
  
  time_control <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("endpoint1","endpoint2")[iter_endpoint]]
  listKMC_Peto1[[iter_endpoint]]  <- cbind(survestKM_tempo(time_control-threshold_Peto1[iter_endpoint]),
                                     survestKM_tempo(time_control),
                                     survestKM_tempo(time_control+threshold_Peto1[iter_endpoint])
  )
  rownames(listKMC_Peto1[[iter_endpoint]]) <- time_control
  
}


resC1_testTTE_Peto1 <- Test_calcAllPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto1$endpoint1[data_testTTE_Peto1$treatment==1],
                                                         Control=data_testTTE_Peto1$endpoint1[data_testTTE_Peto1$treatment==0],
                                                         deltaT=data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==1],
                                                         deltaC=data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==0],
                                                         matKMT=listKMT_Peto1[[1]], matKMC=listKMC_Peto1[[1]],type=1,
                                                         threshold=1,trace=1)
# count  (favorable / unfavorable / neutral / uninf) : 616.493 711.047 536 636.46
# index size (neutralT / neutralC / uninfT/ uninfC) : 536 536 1480 1480
# size (weights/ index_weights) : 1480 0

#### subset
resC1_testTTE_Peto2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto1$endpoint1[data_testTTE_Gehan1$treatment==1],
                                                            Control=data_testTTE_Peto1$endpoint1[data_testTTE_Gehan1$treatment==0],threshold=1,
                                                            deltaT=data_testTTE_Peto1$event1[data_testTTE_Gehan1$treatment==1],
                                                            deltaC=data_testTTE_Peto1$event1[data_testTTE_Gehan1$treatment==0],
                                                            matKMT=listKMT_Peto1[[1]], matKMC=listKMC_Peto1[[1]],
                                                            index_neutralT=resC1_testTTE_Peto1$index_neutralT, index_neutralC=resC1_testTTE_Peto1$index_neutralC, nNeutral_pairs=length(resC1_testTTE_Peto1$index_neutralC),
                                                            index_uninfT=resC1_testTTE_Peto1$index_uninfT, index_uninfC=resC1_testTTE_Peto1$index_uninfC, nUninf_pairs=length(resC1_testTTE_Peto1$index_uninfC),
                                                            Wpairs=rep(1,length(resC1_testTTE_Peto1$index_neutralC)+length(resC1_testTTE_Peto1$index_uninfC)),threshold_M1=1,matKMT_M1=listKMT_Peto1[[1]], matKMC_M1=listKMC_Peto1[[1]],type=1,trace=1)


# count  (favorable / unfavorable / neutral / uninf) : 0 0 536 1480
# index size (neutralT / neutralC / uninfT/ uninfC) : 536 536 1480 1480
# size (weights/ index_weights) : 2016 2016

resC1_testTTE_Peto2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto1$endpoint1[data_testTTE_Gehan1$treatment==1],
                                                           Control=data_testTTE_Peto1$endpoint1[data_testTTE_Gehan1$treatment==0],threshold=1,
                                                           deltaT=data_testTTE_Peto1$event1[data_testTTE_Gehan1$treatment==1],
                                                           deltaC=data_testTTE_Peto1$event1[data_testTTE_Gehan1$treatment==0],
                                                           matKMT=listKMT_Peto1[[1]], matKMC=listKMC_Peto1[[1]],
                                                           index_neutralT=numeric(0), index_neutralC=numeric(0), nNeutral_pairs=0,
                                                           index_uninfT=expand.grid(0:(n.Treatment_testTTE1-1),0:(n.Control_testTTE1-1))[,1], index_uninfC=expand.grid(0:(n.Treatment_testTTE1-1),0:(n.Control_testTTE1-1))[,2], nUninf_pairs=n.Treatment_testTTE1*n.Control_testTTE1,
                                                           Wpairs=rep(1,n.Treatment_testTTE1*n.Control_testTTE1),threshold_M1=1,matKMT_M1=matrix(), matKMC_M1=matrix(),type=1,trace=1)

# count  (favorable / unfavorable / neutral / uninf) : 616.493 711.047 536 636.46
# index size (neutralT / neutralC / uninfT/ uninfC) : 536 536 1480 1480
# size (weights/ index_weights) : 2016 2016

resC3_testTTE_Peto2 <- Test_calcSubsetPairs_TTEOutcome_cpp(Treatment=data_testTTE_Peto1$endpoint2[data_testTTE_Peto1$treatment==1],
                                                            Control=data_testTTE_Peto1$endpoint2[data_testTTE_Peto1$treatment==0],threshold=1,
                                                            deltaT=data_testTTE_Peto1$event2[data_testTTE_Peto1$treatment==1],
                                                            deltaC=data_testTTE_Peto1$event2[data_testTTE_Peto1$treatment==0],
                                                            matKMT=listKMT_Peto1[[2]], matKMC=listKMC_Peto1[[2]],
                                                            index_neutralT=resC1_testTTE_Peto1$index_neutralT, index_neutralC=resC1_testTTE_Peto1$index_neutralC, nNeutral_pairs=length(resC1_testTTE_Peto1$index_neutralC),
                                                            index_uninfT=resC1_testTTE_Peto1$index_uninfT, index_uninfC=resC1_testTTE_Peto1$index_uninfC, nUninf_pairs=length(resC1_testTTE_Peto1$index_uninfC),
                                                            Wpairs=c(rep(1,length(resC1_testTTE_Peto1$index_neutralT)),resC1_testTTE_Peto1$w),threshold_M1=1,matKMT_M1=matrix(), matKMC_M1=matrix(),type=1,trace=1)

# count  (favorable / unfavorable / neutral / uninf) : 509.616 36.2431 324.721 301.88
# index size (neutralT / neutralC / uninfT/ uninfC) : 569 569 1173 1173
# size (weights/ index_weights) : 1742 1742
