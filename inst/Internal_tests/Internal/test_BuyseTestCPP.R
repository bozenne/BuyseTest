#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Simulation et test sur donneees par pairs : donnees binaires
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path <- "E:/Creation_package/Package_BuyseTest/BuyseTest" # path to the uncompressed tar.gz file

Rcpp::sourceCpp(file.path(path,"src/FCT_BuyseTest.cpp"),rebuild=TRUE)

#### 1- Binary endpoint ####

#### base ####
data_testBin0 <- data.frame(treatment=c(0,1),control=c(0,1),strataT=1,strataC=1)


res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testBin0$treatment), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testBin0$treatment), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])

# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 1                 1                 2                 0 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(c(NA,data_testBin0$treatment[-1])), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(c(NA,data_testBin0$treatment[-1])), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.25              1.00              0.00              1.00              2.00 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(c(data_testBin0$treatment[-2],NA)), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(c(data_testBin0$treatment[-2],NA)), Control=cbind(data_testBin0$control), threshold=NA, type=1, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testBin0$strataT==1))-1), strataC=list(which((data_testBin0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.25              0.00              1.00              1.00              2.00 

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
data_testBin1$strata <- 1


res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testBin1$endpoint1[data_testBin1$treatment==1]), 
                           Control=cbind(data_testBin1$endpoint1[data_testBin1$treatment==0]), 
                           threshold=NA, type=1, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testBin1$endpoint1[data_testBin1$treatment==1]), 
                           Control=cbind(data_testBin1$endpoint1[data_testBin1$treatment==0]), 
                           threshold=NA, type=1, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.04            660.00            560.00           1280.00              0.00 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testBin1$endpoint2[data_testBin1$treatment==1]), 
                           Control=cbind(data_testBin1$endpoint2[data_testBin1$treatment==0]), 
                           threshold=NA, type=1, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testBin1$endpoint2[data_testBin1$treatment==1]), 
                           Control=cbind(data_testBin1$endpoint2[data_testBin1$treatment==0]), 
                           threshold=NA, type=1, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.5            1400.0             150.0             950.0               0.0

res <- BuyseTest_Gehan_cpp(Treatment=as.matrix(data_testBin1[data_testBin1$treatment==1,c("endpoint1","endpoint2")]), 
                           Control=as.matrix(data_testBin1[data_testBin1$treatment==0,c("endpoint1","endpoint2")]), 
                           threshold=c(NA,NA), type=c(1,1), delta_Treatment=matrix(), delta_Control=matrix(),
                           D=2, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=as.matrix(data_testBin1[data_testBin1$treatment==1,c("endpoint1","endpoint2")]), 
                           Control=as.matrix(data_testBin1[data_testBin1$treatment==0,c("endpoint1","endpoint2")]), 
                           threshold=c(NA,NA), type=c(1,1), delta_Treatment=matrix(), delta_Control=matrix(),
                           D=2, returnIndex=TRUE, 
                           strataT=list(which(data_testBin1[data_testBin1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testBin1[data_testBin1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta1             delta2   count_favorable1   count_favorable2 count_unfavorable1 count_unfavorable2     count_neutral1 
# 0.04               0.26             660.00             725.00             560.00              75.00            1280.00 
# count_neutral2       count_uninf1       count_uninf2 
# 480.00               0.00               0.00 

#### 2- Continuous endpoint ####

#### base ####
data_testCont0 <- data.frame(treatment=c(-0.5,0,0.5),control=c(-0.5,0,0.5),strataT=1,strataC=1)

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testCont0$treatment), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testCont0$treatment), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 1                 1                 7                 0 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testCont0$treatment), Control=cbind(data_testCont0$control), threshold=0, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testCont0$treatment), Control=cbind(data_testCont0$control), threshold=0, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 3                 3                 3                 0 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(c(NA,data_testCont0$treatment[-1])), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(c(NA,data_testCont0$treatment[-1])), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.1111111         1.0000000         0.0000000         5.0000000         3.0000000 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(c(data_testCont0$treatment[-3],NA)), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(c(data_testCont0$treatment[-3],NA)), Control=cbind(data_testCont0$control), threshold=1, type=2, 
                           delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.1111111         0.0000000         1.0000000         5.0000000         3.0000000 

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
data_testCont1$strata <- 1


res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testCont1$endpoint1[data_testCont1$treatment==1]), 
                           Control=cbind(data_testCont1$endpoint1[data_testCont1$treatment==0]), 
                           threshold=1, type=2, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testCont1$endpoint1[data_testCont1$treatment==1]), 
                           Control=cbind(data_testCont1$endpoint1[data_testCont1$treatment==0]), 
                           threshold=1, type=2, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.1796          367.0000          816.0000         1317.0000            0.0000

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testCont1$endpoint2[data_testCont1$treatment==1]), 
                           Control=cbind(data_testCont1$endpoint2[data_testCont1$treatment==0]), 
                           threshold=1, type=2, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testCont1$endpoint2[data_testCont1$treatment==1]), 
                           Control=cbind(data_testCont1$endpoint2[data_testCont1$treatment==0]), 
                           threshold=1, type=2, delta_Treatment=matrix(), delta_Control=matrix(),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.9924         2481.0000            0.0000           19.0000            0.0000 

res <- BuyseTest_Gehan_cpp(Treatment=as.matrix(data_testCont1[data_testCont1$treatment==1,c("endpoint1","endpoint2")]), 
                           Control=as.matrix(data_testCont1[data_testCont1$treatment==0,c("endpoint1","endpoint2")]), 
                           threshold=c(1,1), type=c(2,2), delta_Treatment=matrix(), delta_Control=matrix(),
                           D=2, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0)
unlist(res[1:5])
res <- BuyseTest_PetoEfron_cpp(Treatment=as.matrix(data_testCont1[data_testCont1$treatment==1,c("endpoint1","endpoint2")]), 
                           Control=as.matrix(data_testCont1[data_testCont1$treatment==0,c("endpoint1","endpoint2")]), 
                           threshold=c(1,1), type=c(2,2), delta_Treatment=matrix(), delta_Control=matrix(),
                           D=2, returnIndex=TRUE, 
                           strataT=list(which(data_testCont1[data_testCont1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testCont1[data_testCont1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=0,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=numeric(0),list_survivalT=list(),list_survivalC=list(),Efron=FALSE)
unlist(res[1:5])
# delta1             delta2   count_favorable1   count_favorable2 count_unfavorable1 count_unfavorable2     count_neutral1 
# -0.1796             0.5224           367.0000          1306.0000           816.0000             0.0000          1317.0000 
# count_neutral2       count_uninf1       count_uninf2 
# 11.0000             0.0000             0.0000 


#### 3- Time To event endpoint ####

#### base (Gehan) ####
data_testTTE_Gehan0 <- data.frame(treatment=c(1,2,2.5),control=c(1,2,2.5),strataT=1,strataC=1)

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan0$treatment), Control=cbind(data_testTTE_Gehan0$control), threshold=1, type=3, 
                           delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(1,1,1,1)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 2                 2                 5                 0 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan0$treatment), Control=cbind(data_testTTE_Gehan0$control), threshold=0, type=3, 
                           delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(1,1,1,1)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 3                 3                 3                 0 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan0$treatment), Control=cbind(data_testTTE_Gehan0$control), threshold=0, type=3, 
                           delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(0,0,0,0)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.3333333         0.0000000         3.0000000         0.0000000         6.0000000 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan0$treatment), Control=cbind(data_testTTE_Gehan0$control), threshold=0, type=3, 
                           delta_Treatment=cbind(c(0,0,0,0)),delta_Control=cbind(c(1,1,1,1)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.3333333         3.0000000         0.0000000         0.0000000         6.0000000 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan0$treatment), Control=cbind(data_testTTE_Gehan0$control), threshold=0, type=3, 
                           delta_Treatment=cbind(c(0,0,0,0)),delta_Control=cbind(c(0,0,0,0)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 0                 0                 0                 9 

#### large (Gehan) ####
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
data_testTTE_Gehan1$strata <- 1


res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==1]), 
                           Control=cbind(data_testTTE_Gehan1$endpoint1[data_testTTE_Gehan1$treatment==0]), 
                           threshold=0.25, type=3, 
                           delta_Treatment=cbind(data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==1]), 
                           delta_Control=cbind(data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==0]),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.0548          487.0000          624.0000          208.0000         1181.0000 

res <- BuyseTest_Gehan_cpp(Treatment=cbind(data_testTTE_Gehan1$endpoint2[data_testTTE_Gehan1$treatment==1]), 
                           Control=cbind(data_testTTE_Gehan1$endpoint2[data_testTTE_Gehan1$treatment==0]), 
                           threshold=0.25, type=3, 
                           delta_Treatment=cbind(data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==1]), 
                           delta_Control=cbind(data_testTTE_Gehan1$event1[data_testTTE_Gehan1$treatment==0]),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=1)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.3328          916.0000           84.0000          272.0000         1228.0000 

res <- BuyseTest_Gehan_cpp(Treatment=as.matrix(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==1,c("endpoint1","endpoint2")]), 
                           Control=as.matrix(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==0,c("endpoint1","endpoint2")]), 
                           threshold=c(0.25,0.25), type=c(3,3), 
                           delta_Treatment=as.matrix(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==1,c("event1","event2")]), 
                           delta_Control=as.matrix(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==0,c("event1","event2")]),
                           D=2, returnIndex=TRUE, 
                           strataT=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testTTE_Gehan1[data_testTTE_Gehan1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=2)
unlist(res[1:5])
# delta1             delta2   count_favorable1   count_favorable2 count_unfavorable1 count_unfavorable2     count_neutral1 
# -0.0548             0.3064           487.0000           790.0000           624.0000            24.0000           208.0000 
# count_neutral2       count_uninf1       count_uninf2 
# 160.0000          1181.0000           415.0000 


#### base (Peto) ####
data("veteran",package="survival")
data_testTTE_Peto0 <- data.frame(treatment=c(1,2,2.5),control=c(1,2,2.5),strataT=1,strataC=1)
resKMT_tempo <- survival::survfit(survival::Surv(veteran$time/50,veteran$trt)~1) # compute the survival over the controls with common Kaplan Meier estimator.  
survestKMT_tempo <- stats::stepfun(c(resKMT_tempo$time,max(veteran$time/50)), c(1,resKMT_tempo$surv,NA)) # step interpolation of the survival function

list_survivalT_Peto0 <- list(cbind(survestKMT_tempo(c(1,2,2.5)-1),survestKMT_tempo(c(1,2,2.5)),survestKMT_tempo(c(1,2,2.5)+1)))
list_survivalC_Peto0 <- list(cbind(survestKMT_tempo(c(1,2,2.5)-1),survestKMT_tempo(c(1,2,2.5)),survestKMT_tempo(c(1,2,2.5)+1)))

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto0$treatment), Control=cbind(data_testTTE_Peto0$control), threshold=1, type=3, 
                           delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(1,1,1,1)),
                           D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1,
                           Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=c(0),list_survivalT=list_survivalT_Peto0,list_survivalC=list_survivalC_Peto0,Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 2                 2                 5                 0 

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto0$treatment), Control=cbind(data_testTTE_Peto0$control), threshold=0, type=3, 
                               delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(1,1,1,1)),
                               D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1,
                               Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=c(0),list_survivalT=list_survivalT_Peto0,list_survivalC=list_survivalC_Peto0,Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0                 3                 3                 3                 0 

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto0$treatment), Control=cbind(data_testTTE_Peto0$control), threshold=1, type=3, 
                               delta_Treatment=cbind(c(1,1,1,1)),delta_Control=cbind(c(0,0,0,0)),
                               D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1,
                               Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=c(0),list_survivalT=list_survivalT_Peto0,list_survivalC=list_survivalC_Peto0,Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.84521277        0.08623725        7.69315218        0.00000000        1.22061057 

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto0$treatment), Control=cbind(data_testTTE_Peto0$control), threshold=1, type=3, 
                               delta_Treatment=cbind(c(0,0,0,0)),delta_Control=cbind(c(1,1,1,1)),
                               D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1,
                               Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=c(0),list_survivalT=list_survivalT_Peto0,list_survivalC=list_survivalC_Peto0,Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.84521277        7.69315218        0.08623725        0.00000000        1.22061057 

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto0$treatment), Control=cbind(data_testTTE_Peto0$control), threshold=1, type=3, 
                               delta_Treatment=cbind(c(0,0,0,0)),delta_Control=cbind(c(0,0,0,0)),
                               D=1, returnIndex=TRUE, strataT=list(which((data_testCont0$strataT==1))-1), strataC=list(which((data_testCont0$strataC==1))-1), n_strata=1, n_TTE=1,
                               Wscheme=matrix(),threshold_TTEM1=-1,index_TTEM1=c(0),list_survivalT=list_survivalT_Peto0,list_survivalC=list_survivalC_Peto0,Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.000000          3.889695          3.889695          0.000000          1.220611 

#### large (Peto) ####
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
data_testTTE_Peto1$strata <- 1

list_survivalT_Peto1 <- list()
list_survivalC_Peto1 <- list()
threshold_Peto1 <- c(1,0.5)

for(iter_endpoint in 1:2){
  for(iter_threshold in 1:2){
 
    resKM_tempo <- survival::survfit(survival::Surv(data_testTTE_Peto1[,c("endpoint1","endpoint2")[iter_endpoint]],
                                                  data_testTTE_Peto1[,c("event1","event2")[iter_endpoint]])~1)
  survestKM_tempo <- stats::stepfun(resKM_tempo$time, c(1, resKM_tempo$surv))
  
  time_treatment <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("endpoint1","endpoint2")[iter_endpoint]]
  list_survivalT_Peto1[[2*(iter_endpoint-1)+iter_threshold]]  <- cbind(survestKM_tempo(time_treatment-threshold_Peto1[iter_threshold]),
                                           survestKM_tempo(time_treatment),
                                           survestKM_tempo(time_treatment+threshold_Peto1[iter_threshold])
  )
  rownames(list_survivalT_Peto1[[2*(iter_endpoint-1)+iter_threshold]]) <- time_treatment
  
  time_control <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("endpoint1","endpoint2")[iter_endpoint]]
  list_survivalC_Peto1[[2*(iter_endpoint-1)+iter_threshold]]  <- cbind(survestKM_tempo(time_control-threshold_Peto1[iter_threshold]),
                                           survestKM_tempo(time_control),
                                           survestKM_tempo(time_control+threshold_Peto1[iter_threshold])
  )
  rownames(list_survivalC_Peto1[[2*(iter_endpoint-1)+iter_threshold]]) <- time_control
  }
}
names(list_survivalT_Peto1) <- c(paste("endpoint1_th",1:2,sep=""),paste("endpoint2_th",1:2,sep=""))
names(list_survivalC_Peto1) <- c(paste("endpoint1_th",1:2,sep=""),paste("endpoint2_th",1:2,sep=""))

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto1$endpoint1[data_testTTE_Peto1$treatment==1]), 
                           Control=cbind(data_testTTE_Peto1$endpoint1[data_testTTE_Peto1$treatment==0]), 
                           threshold=1, type=3, 
                           delta_Treatment=cbind(data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==1]), 
                           delta_Control=cbind(data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==0]),
                           D=1, returnIndex=TRUE, 
                           strataT=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"strata"]==1)-1), 
                           strataC=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"strata"]==1)-1),
                           n_strata=1, n_TTE=1,
                           Wscheme=matrix(),index_TTEM1=0,threshold_TTEM1=-1,list_survivalT=list_survivalT_Peto1[1:2],list_survivalC=list_survivalC_Peto1[1:2],Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# -0.03782172      616.49263787      711.04692705      536.00000000      636.46043508

res <- BuyseTest_PetoEfron_cpp(Treatment=cbind(data_testTTE_Peto1$endpoint2[data_testTTE_Peto1$treatment==1]), 
                               Control=cbind(data_testTTE_Peto1$endpoint2[data_testTTE_Peto1$treatment==0]), 
                               threshold=1, type=3, 
                               delta_Treatment=cbind(data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==1]), 
                               delta_Control=cbind(data_testTTE_Peto1$event1[data_testTTE_Peto1$treatment==0]),
                               D=1, returnIndex=TRUE, 
                               strataT=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"strata"]==1)-1), 
                               strataC=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"strata"]==1)-1),
                               n_strata=1, n_TTE=1,
                               Wscheme=matrix(),index_TTEM1=0,threshold_TTEM1=-1,list_survivalT=list_survivalT_Peto1[c(1,4)],list_survivalC=list_survivalC_Peto1[c(1,4)],Efron=FALSE)
unlist(res[1:5])
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.206382        756.226015        240.270980        661.000000        842.503005
# delta   count_favorable count_unfavorable     count_neutral       count_uninf 
# 0.02555611      574.65584994      510.76557184      661.00000000      753.57857822 

res <- BuyseTest_PetoEfron_cpp(Treatment=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("endpoint1","endpoint2")]), 
                               Control=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("endpoint1","endpoint2")]), 
                               threshold=c(1,1), type=c(3,3), 
                               delta_Treatment=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("event1","event2")]), 
                               delta_Control=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("event1","event2")]),
                               D=2, returnIndex=TRUE, 
                               strataT=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"strata"]==1)-1), 
                               strataC=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"strata"]==1)-1),
                               n_strata=1, n_TTE=2,
                               Wscheme=cbind(c(1,0)),index_TTEM1=c(-1,-1),threshold_TTEM1=c(0,0),
                               list_survivalT=list_survivalT_Peto1[c(1,3)],
                               list_survivalC=list_survivalC_Peto1[c(1,3)],Efron=FALSE)
unlist(res[1:5])
# delta1             delta2   count_favorable1   count_favorable2 count_unfavorable1 count_unfavorable2     count_neutral1 
# -0.03782172         0.18934921       616.49263787       509.61610364       711.04692705        36.24307430       536.00000000 
# count_neutral2       count_uninf1       count_uninf2 
# 324.72114914       636.46043508       301.88010799

res <- BuyseTest_PetoEfron_cpp(Treatment=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("endpoint1","endpoint1")]), 
                               Control=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("endpoint1","endpoint1")]), 
                               threshold=c(1,0.5), type=c(3,3), 
                               delta_Treatment=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,c("event1","event1")]), 
                               delta_Control=as.matrix(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,c("event1","event1")]),
                               D=2, returnIndex=TRUE, 
                               strataT=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"strata"]==1)-1), 
                               strataC=list(which(data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"strata"]==1)-1),
                               n_strata=1, n_TTE=2,
                               Wscheme=cbind(c(0,0)),index_TTEM1=c(-1,0),threshold_TTEM1=c(0,1),
                               list_survivalT=list_survivalT_Peto1[c(1,2)],
                               list_survivalC=list_survivalC_Peto1[c(1,2)],Efron=FALSE)
unlist(res[1:5])
# delta1             delta2   count_favorable1   count_favorable2 count_unfavorable1 count_unfavorable2     count_neutral1 
# -3.782172e-02       3.460274e-04       6.164926e+02       7.623161e+01       7.110469e+02       7.536654e+01       5.360000e+02 
# count_neutral2       count_uninf1       count_uninf2 
# 3.560000e+02       6.364604e+02       4.608313e+02 



#### TO DO ####
# verifier la validite des resultats de la section (large) PETO