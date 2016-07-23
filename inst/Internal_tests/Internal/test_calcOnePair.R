#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%                                            %%%%%
#%%% check the validity of the calcCPPfunctions #%%%%
#%%%                                            %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path <- "E:/Creation_package/Package_BuyseTest/BuyseTest" # path to the uncompressed tar.gz file
Rcpp::sourceCpp(file.path(path,"test/test_calcOnePair.cpp"))

#### 1- binary outcome ####

## non weigthed
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=1,endpoint_C=0,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=0,endpoint_C=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=1,endpoint_C=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=0,endpoint_C=0,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=c(1,0,1,0),endpoint_C=c(0,1,1,0),Wpair=c(1,1,1,1),iter_pair=rep(-1,4))
# count  (favorable / unfavorable / neutral / uninf) : 1 1 2 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 2 2 0 0
# # index_neutralT : 2 3 
# # index_neutralC : 2 3 
# index size (wNeutral/ wUninf) : 0 0

## weighted
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=1,endpoint_C=0,Wpair=0.3,iter_pair=0)
# count  (favorable / unfavorable / neutral / uninf) : 0.3 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=0,endpoint_C=1,Wpair=0.3,iter_pair=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0.3 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=1,endpoint_C=1,Wpair=0.3,iter_pair=2)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0.3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 2 
res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=0,endpoint_C=0,Wpair=0.3,iter_pair=3)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0.3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 3 

res <- Test_calcOnePair_BinaryOutcome_cpp(endpoint_T=c(1,0,1,0),endpoint_C=c(0,1,1,0),Wpair=c(0.3,0.3,0.3,0.3),iter_pair=0:3)
# count  (favorable / unfavorable / neutral / uninf) : 0.3 0.3 0.6 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 2 2 0 0
# # index_neutralT : 2 3 
# # index_neutralC : 2 3 
# index size (wNeutral/ wUninf) : 2 0
# # index_wUNeutral : 2 3 

#### 2- continuous outcome ####

## non weighted
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,threshold=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=-5.5,endpoint_C=0.3,threshold=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,threshold=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,threshold=1,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,threshold=0,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,threshold=0,Wpair=1,iter_pair=-1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=c(5.5,-5.5,5.5,5.5),endpoint_C=c(0.3,0.3,5.3,4.5),
                                           threshold=1,Wpair=c(1,1,1,1,1),iter_pair=rep(-1,5))
# count  (favorable / unfavorable / neutral / uninf) : 2 1 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 2 
# # index_neutralC : 2 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=c(5.5,5.5),endpoint_C=c(5.3,5.5),
                                           threshold=0,Wpair=c(1,1),iter_pair=0:1)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 1 
# # index_neutralC : 1 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 1 

## weighted
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,threshold=1,Wpair=0.4,iter_pair=0)
# count  (favorable / unfavorable / neutral / uninf) : 0.4 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=-5.5,endpoint_C=0.3,threshold=1,Wpair=0.4,iter_pair=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0.4 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,threshold=1,Wpair=0.4,iter_pair=2)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0.4 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 2 
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,threshold=1,Wpair=0.4,iter_pair=3)
# count  (favorable / unfavorable / neutral / uninf) : 0.4 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,threshold=0,Wpair=0.4,iter_pair=0)
# count  (favorable / unfavorable / neutral / uninf) : 0.4 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,threshold=0,Wpair=0.4,iter_pair=1)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0.4 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 1 

res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=c(5.5,-5.5,5.5,5.5),endpoint_C=c(0.3,0.3,5.3,4.5),
                                           threshold=1,Wpair=c(0.3,0.3,0.3,0.3,0.3),iter_pair=0:3)
# count  (favorable / unfavorable / neutral / uninf) : 0.6 0.3 0.3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 2 
# # index_neutralC : 2 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 2 
res <- Test_calcOnePair_ContinuousOutcome_cpp(endpoint_T=c(5.5,5.5),endpoint_C=c(5.3,5.5),
                                           threshold=0,Wpair=c(0.3,0.3),iter_pair=0:1)
# count  (favorable / unfavorable / neutral / uninf) : 0.3 0 0.3 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 1 
# # index_neutralC : 1 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 1 

#### 3- TTE outcome ####

#### Gehan ####

## no censoring
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=0.5,endpoint_C=5.3,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=1,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 0 
# # index_neutralC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=c(5.5,0.5,5.5,5.5),endpoint_C=c(0.3,5.3,5.3,4.5),
                                           delta_T=rep(1,4),delta_C=rep(1,4),survival_T=matrix(),survival_C=matrix(),
                                           threshold=1,Wpair=c(1,1,1,1),iter_pair=0:3,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 2 1 1 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 1 1 0 0
# # index_neutralT : 2 
# # index_neutralC : 2 
# index size (wNeutral/ wUninf) : 1 0
# # index_wUNeutral : 2 

## treatment censored
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=0.5,endpoint_C=5.3,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 1 0 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=0,delta_C=1,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=c(5.5,0.5,5.5,5.5),endpoint_C=c(0.3,5.3,5.3,4.5),
                                    delta_T=rep(0,4),delta_C=rep(1,4),survival_T=matrix(),survival_C=matrix(),
                                    threshold=1,Wpair=c(1,1,1,1),iter_pair=0:3,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 2 0 0 2
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 2 2
# # index_uninfT : 1 2 
# # index_uninfC : 1 2 
# index size (wNeutral/ wUninf) : 0 2
# # index_wUninf : 1 2 

## control censored
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=0.5,endpoint_C=5.3,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 0 0
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 0 0
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=1,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=c(5.5,0.5,5.5,5.5),endpoint_C=c(0.3,5.3,5.3,4.5),
                                    delta_T=rep(1,4),delta_C=rep(0,4),survival_T=matrix(),survival_C=matrix(),
                                    threshold=1,Wpair=c(1,1,1,1),iter_pair=0:3,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 1 0 3
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 3 3
# # index_uninfT : 0 2 3 
# # index_uninfC : 0 2 3 
# index size (wNeutral/ wUninf) : 0 3
# # index_wUninf : 0 2 3 

## treatment and control censored
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=0.5,endpoint_C=5.3,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=1,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=0,delta_C=0,survival_T=matrix(),survival_C=matrix(),threshold=0,Wpair=1,iter_pair=-1,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 1
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 1 1
# # index_uninfT : 0 
# # index_uninfC : 0 
# index size (wNeutral/ wUninf) : 0 0

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=c(5.5,0.5,5.5,5.5),endpoint_C=c(0.3,5.3,5.3,4.5),
                                    delta_T=rep(0,4),delta_C=rep(0,4),survival_T=matrix(),survival_C=matrix(),
                                    threshold=1,Wpair=c(1,1,1,1),iter_pair=0:3,type=0)
# count  (favorable / unfavorable / neutral / uninf) : 0 0 0 4
# index size (neutralT / neutralC / uninfT/ uninfC) : 0 0 4 4
# # index_uninfT : 0 1 2 3 
# # index_uninfC : 0 1 2 3 
# index size (wNeutral/ wUninf) : 0 4
# # index_wUninf : 0 1 2 3 

#### Peto ####
resKMT_tempo <- survival::survfit(survival::Surv(c(4,3,11.5,1,6,2,3),c(1,1,1,0,1,1,0))~1) # compute the survival over the controls with common Kaplan Meier estimator.  
survestKMT_tempo <- stats::stepfun(c(resKMT_tempo$time,11.6), c(1,resKMT_tempo$surv,NA)) # step interpolation of the survival function

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(-0.7,0.3,1.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 1 0 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=10.3,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(9.3,10.3,11.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 1 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(3.5,4.5,5.5))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 1 0 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : -1 -1 1 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=0,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 1 0 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=1,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.5,5.5,6.5))),threshold=0,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : -1 -1 1 0

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(-0.7,0.3,1.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0.555556 0.222222 0 1
# (C=0,T=1) >= tau : # 1-survestKMT_tempo(4.5)/survestKMT_tempo(0.3) # survestKMT_tempo(6.5)/survestKMT_tempo(0.3)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=10.3,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(9.3,10.3,11.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 1 0 0
# (C=0,T=1) <= -tau
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(3.5,4.5,5.5))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 0.5 0 1
# (C=0,T=1) >= tau : # 1-survestKMT_tempo(4.5)/survestKMT_tempo(4.5) # survestKMT_tempo(6.5)/survestKMT_tempo(4.5)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 0.5 0 1
# (C=0,T=1) >= tau : # 1-survestKMT_tempo(4.5)/survestKMT_tempo(4.3) # survestKMT_tempo(6.5)/survestKMT_tempo(4.3)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=0,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 0.5 0 1
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=1,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.5,5.5,6.5))),threshold=0,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 0 0.5 0 1

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(-0.7,0.3,1.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
# output (probaF / probaUF / test.neutral / test.uninf) : 1 0 0 0
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=10.3,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(9.3,10.3,11.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(3.5,4.5,5.5))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=0,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=0,delta_C=1,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.5,5.5,6.5))),threshold=0,Wpair=1,iter_pair=-1,type=1)

res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=0.3,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(-0.7,0.3,1.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=10.3,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(9.3,10.3,11.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=4.5,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(3.5,4.5,5.5))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=1,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.3,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.3,5.3,6.3))),threshold=0,Wpair=1,iter_pair=-1,type=1)
res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=5.5,endpoint_C=5.5,delta_T=0,delta_C=0,survival_T=rbind(survestKMT_tempo(c(4.5,5.5,6.5))), survival_C=rbind(survestKMT_tempo(c(4.5,5.5,6.5))),threshold=0,Wpair=1,iter_pair=-1,type=1)






# index_T <- 2
# index_C <- 6
# endpoint_T <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"endpoint1"][index_T]
# endpoint_C <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"endpoint1"][index_C]
# delta_T <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==1,"event1"][index_T]
# delta_C <- data_testTTE_Peto1[data_testTTE_Peto1$treatment==0,"event1"][index_C]
# survival_T <- listKMT_Peto1[[2]][index_T,]
# survival_C <- listKMT_Peto1[[2]][index_C,]
# 
# res <- Test_calcOnePair_TTEOutcome_cpp(endpoint_T=endpoint_T,endpoint_C=endpoint_C,delta_T=delta_T,delta_C=delta_C,
#                                        survival_T=rbind(survival_T), survival_C=rbind(survival_C),
#                                        threshold=0.25,Wpair=1,iter_pair=-1,type=1)
# delta>tau
# (0,1)
# 
# survival_C[3]/survival_T[2]
# 1-survival_C[1]/survival_T[2]

#### Efron ####
load("E:/Creation_package/Package_BuyseTest/BuyseTest/test/test_Efron.RData")
iter_T <- 1
iter_C <- 1

res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                       test.Efron$list_survivalT[iter_T,,drop=FALSE], test.Efron$list_survivalC[iter_C,,drop=FALSE],
                                       0.25, 1, -1,0, 1)
res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                       test.Efron$list_survivalT, test.Efron$list_survivalC,
                                       0.25, 1, c(iter_T-1,iter_C-1),1, 1)
res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                       test.Efron$list_survivalT, test.Efron$list_survivalC,
                                       0.25, 1, c(iter_T-1,iter_C-1),2, 1)


iter_T <- 10
iter_C <- 1

res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                       test.Efron$list_survivalT[iter_T,,drop=FALSE], test.Efron$list_survivalC[iter_C,,drop=FALSE],
                                       0.25, 1, -1,0, 1)
res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                       test.Efron$list_survivalT, test.Efron$list_survivalC,
                                       0.25, 1, c(iter_T-1,iter_C-1),2, 1)
output (probaF / probaUF / test.neutral / test.uninf) : 0.293257 0.575947 0 1


for(iter_T in 1:500){
  for(iter_C in 1:500){
    res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                           test.Efron$list_survivalT, test.Efron$list_survivalC,
                                           0.25, 1, c(iter_T-1,iter_C-1),2, 0)
    if(any(res$res[1:2]>1) || any(res$res[1:2]<0) || sum(res$res[1:2]>1)){
      res <- Test_calcOnePair_TTEOutcome_cpp(test.Efron$M.Treatment[iter_T], test.Efron$M.Control[iter_C], test.Efron$M.delta_Treatment[iter_T], test.Efron$M.delta_Control[iter_C],
                                             test.Efron$list_survivalT, test.Efron$list_survivalC,
                                             0.25, 1, c(iter_T-1,iter_C-1),2, 1)
    }
      
  }
}

test.Efron$list_survivalT[iter_T,"SurvivalC_TimeT-threshold"]
test.Efron$list_survivalT[iter_T,"SurvivalC_TimeT-0"]
test.Efron$list_survivalT[iter_T,"SurvivalC_TimeT-threshold"]


#####
integral = 0
Sintegral  <- rep(NA,nrow(test.Efron$list_survivalC))

if(test.Efron$list_survivalC[1,"time_control(ordered)"]>test.Efron$M.Control[iter_C]){
  Sintegral[1] <-  test.Efron$list_survivalC[1,"SurvivalT_TimeC-threshold(ordered)"]*(test.Efron$list_survivalC[1,"SurvivalC_TimeC_0(ordered)"]-1)
  integral = integral + Sintegral[1]
}


for(iter_time in 2:nrow(test.Efron$list_survivalC)){
  
  if(test.Efron$list_survivalC[iter_time,"time_control(ordered)"]>test.Efron$M.Control[iter_C]){
    Sintegral[iter_time] <-  test.Efron$list_survivalC[iter_time,"SurvivalT_TimeC-threshold(ordered)"]*(test.Efron$list_survivalC[iter_time,"SurvivalC_TimeC_0(ordered)"]-test.Efron$list_survivalC[iter_time-1,"SurvivalC_TimeC_0(ordered)"])
    integral = integral + Sintegral[iter_time]
  }
  
} 

integral/(0.816645 * 0.633109)
Sintegral


#### Julien ####