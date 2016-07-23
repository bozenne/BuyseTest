#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Simulation et test sur donneees par pairs 
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path <- "E:/Creation_package/Package_BuyseTest/BuyseTest" # path to the uncompressed tar.gz file

Rcpp::sourceCpp(file.path(path,"src/FCT_BuyseTest.cpp"),rebuild=TRUE)
source(file.path(path,"R/FCT_buyseTest.R"))
source(file.path(path,"R/OBJET_buyseTest.R"))
source(file.path(path,"R/FCT_buyseInit.R"))
# require(BuyseTest)

#### 1- binary endpoint ####
n.Treatment_testBin <- 500
n.Control_testBin <- 500
prob.Treatment_testBin <- c(0.5,0.75)
prob.Control_testBin <- c(0.5,0.25)

set.seed(10)
data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)

system.time(
BuyseTest_testBin <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                               treatment="treatment", type=c("bin","bin"),
                               n.bootstrap=1000)
)
# utilisateur     système      écoulé 
# 12.59        1.08       13.75 
summary_BuyseTest_testBin <- summary(BuyseTest_testBin)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta  Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global   100.00        24.30          25.70      49.99        0 -0.014 -0.014     -0.0742       0.050        1000   0.651   
# 3 global    49.99        28.81           2.91      18.28        0  0.259  0.245      0.1750       0.313        1000   0.001 **

#### 2- survival outcome ####

#### data ####

n.Treatment_testTTE <- 500
n.Control_testTTE <- 500
lambda.Treatment_testTTE <- c(0.75,0.5)
lambda.Control_testTTE <- c(0.75,5)
lambda.Censoring_testTTE <- c(0.5,0.5)

set.seed(10)
data_testTTE <- data.frame(treatment=c(rep(1,n.Treatment_testTTE),rep(0,n.Treatment_testTTE)))
data_testTTE$EventTime1 <- c(rexp(n.Treatment_testTTE,rate=lambda.Treatment_testTTE[1]),
                             rexp(n.Control_testTTE,rate=lambda.Control_testTTE[1]))
data_testTTE$EventTime2 <- c(rexp(n.Control_testTTE,rate=lambda.Treatment_testTTE[2]),
                             rexp(n.Control_testTTE,rate=lambda.Control_testTTE[2]))
data_testTTE$CensoringTime1 <- c(rexp(n.Treatment_testTTE,rate=lambda.Censoring_testTTE[1]),
                                 rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[1]))
data_testTTE$CensoringTime2 <- c(rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[2]),
                                 rexp(n.Control_testTTE,rate=lambda.Censoring_testTTE[2]))
data_testTTE$strata <- rbinom(n.Treatment_testTTE+n.Control_testTTE,size=4,prob=0.5)

data_testTTE$endpoint1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,min)
data_testTTE$event1 <- apply(data_testTTE[,c("EventTime1","CensoringTime1")],1,which.min)==1
data_testTTE$event1 <- as.numeric(data_testTTE$event1)
data_testTTE$endpoint2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,min)
data_testTTE$event2 <- apply(data_testTTE[,c("EventTime2","CensoringTime2")],1,which.min)==1
data_testTTE$event2 <- as.numeric(data_testTTE$event2)

#### Gehan ####

#### outcome 1 
BuyseTest_Gehan1_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Gehan",
                               treatment="treatment",censoring=c("event1","event1","event1"),
                               type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Gehan1_3TTE <- summary(BuyseTest_Gehan1_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.unif    delta   Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global   100.00        10.90          13.11      20.59   55.41 -0.02210 -0.0221          NA          NA           0      NA NA
# 3 global    75.99         3.72           4.53      15.66   52.09 -0.00811 -0.0302          NA          NA           0      NA NA
# 5 global    67.74         5.32           6.12       8.95   47.36 -0.00803 -0.0383          NA          NA           0      NA NA

BuyseTest_Gehan1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Gehan",
                               treatment="treatment",censoring="event1",
                               type="TTE",threshold=0.25)
BuyseTestS_Gehan1_1TTE <- summary(BuyseTest_Gehan1_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.unif   delta   Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global      100        19.93          23.76       8.95   47.36 -0.0383 -0.0383          NA          NA           0      NA NA

BuyseTestS_Gehan1_3TTE$nb$Delta[BuyseTestS_Gehan1_3TTE$nb$strata=="global"][3]
BuyseTestS_Gehan1_1TTE$nb$Delta[BuyseTestS_Gehan1_1TTE$nb$strata=="global"][1]


#### outcome 2
BuyseTest_Gehan2_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint2","endpoint2","endpoint2"),method="Gehan",
                                   treatment="treatment",censoring=c("event2","event2","event2"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Gehan2_3TTE <- summary(BuyseTest_Gehan2_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.unif delta Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global   100.00        34.35           0.18      28.70   36.77 0.342 0.342          NA          NA           0      NA NA
# 3 global    65.47        11.07           0.43      22.25   31.71 0.106 0.448          NA          NA           0      NA NA
# 5 global    53.97        14.33           1.32      13.97   24.34 0.130 0.578          NA          NA           0      NA NA

BuyseTest_Gehan2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Gehan",
                                   treatment="treatment",censoring="event2",
                                   type="TTE",threshold=0.25)
BuyseTestS_Gehan2_1TTE <- summary(BuyseTest_Gehan2_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.unif delta Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global      100        59.75           1.94      13.97   24.34 0.578 0.578          NA          NA           0      NA NA

# check number of pairs
rowSums(BuyseTestS_Gehan2_3TTE$nb[BuyseTestS_Gehan2_3TTE$nb$strata=="global",c("n.neutral","n.uninf")])
BuyseTestS_Gehan2_3TTE$nb[BuyseTestS_Gehan2_3TTE$nb$strata=="global",c("n.total")]

# check equal survival
BuyseTestS_Gehan2_3TTE$nb$Delta[BuyseTestS_Gehan2_3TTE$nb$strata=="global"][3]
BuyseTestS_Gehan2_1TTE$nb$Delta[BuyseTestS_Gehan2_1TTE$nb$strata=="global"][1]

# bootstrap
system.time(
BuyseTest_Gehan2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Gehan",
                                   treatment="treatment",censoring="event2",
                                   type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)
)
# utilisateur     système      écoulé 
# 0.28        0.01        8.44 
BuyseTestS_Gehan2_1TTE <- summary(BuyseTest_Gehan2_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.unif delta Delta CIinf.Delta CIsup.Delta n.bootstrap p.value    
# 1 global      100        59.75           1.94      13.97   24.34 0.578 0.578       0.524       0.631        1000       0 ***

#### Peto ####
BuyseTest_Peto1_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Peto",
                                   treatment="treatment",censoring=c("event1","event1","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peto1_3TTE <- summary(BuyseTest_Peto1_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf    delta   Delta n.bootstrap
# 1 global   100.00        27.36          31.12      20.59    20.94 -0.03760 -0.0376           0
# 3 global    41.52         5.39           6.06      15.66    14.41 -0.00671 -0.0443           0
# 5 global    30.07         6.65           7.08       8.95     7.39 -0.00434 -0.0486           0

BuyseTest_Peto1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Peto",
                                   treatment="treatment",censoring="event1",
                                   type="TTE",threshold=0.25)
BuyseTestS_Peto1_1TTE <- summary(BuyseTest_Peto1_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100         39.4          44.26       8.95     7.39 -0.0486 -0.0486           0

# check number of pairs
rowSums(BuyseTestS_Peto1_3TTE$nb[BuyseTestS_Peto1_3TTE$nb$strata=="global",c("n.neutral","n.uninf")])
BuyseTestS_Peto1_3TTE$nb[BuyseTestS_Peto1_3TTE$nb$strata=="global",c("n.total")]

# check equal survival
BuyseTestS_Peto1_3TTE$nb$Delta[BuyseTestS_Peto1_3TTE$nb$strata=="global"][3]
BuyseTestS_Peto1_1TTE$nb$Delta[BuyseTestS_Peto1_1TTE$nb$strata=="global"][1]

# bootstrap
system.time(
BuyseTest_Peto1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Peto",
                                  treatment="treatment",censoring="event1",
                                  type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)
)
# utilisateur     système      écoulé 
# 0.17        0.03       48.26 
BuyseTestS_Peto1_1TTE <- summary(BuyseTest_Peto1_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta CIinf.Delta CIsup.Delta n.bootstrap p.value 
# 1 global      100         39.4          44.26       8.95     7.39 -0.0486 -0.0486      -0.112      0.0136        1000   0.134 

#### outcome 2
BuyseTest_Peto2_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint2","endpoint2","endpoint2"),method="Peto",
                                   treatment="treatment",censoring=c("event2","event2","event2"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peto2_3TTE <- summary(BuyseTest_Peto2_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta Delta n.bootstrap
# 1 global   100.00        54.73           2.16      28.70    14.41 0.5260 0.526           0
# 3 global    43.10         9.79           0.89      22.25    10.17 0.0889 0.615           0
# 5 global    32.42        10.81           2.02      13.97     5.62 0.0879 0.703           0

BuyseTest_Peto2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Peto",
                                   treatment="treatment",censoring="event2",
                                   type="TTE",threshold=0.25)
BuyseTestS_Peto2_1TTE <- summary(BuyseTest_Peto2_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100        75.33           5.08      13.97     5.62 0.703 0.703           0

# check equal survival
BuyseTestS_Peto2_3TTE$nb$Delta[BuyseTestS_Peto2_3TTE$nb$strata=="global"][3]
BuyseTestS_Peto2_1TTE$nb$Delta[BuyseTestS_Peto2_1TTE$nb$strata=="global"][1]


# bootstrap
BuyseTest_Peto2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Peto",
                                  treatment="treatment",censoring="event2",
                                  type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)
BuyseTestS_Peto2_1TTE <- summary(BuyseTest_Peto2_1TTE)

#### Efron ####
BuyseTest_Efron1_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Efron",
                                  treatment="treatment",censoring=c("event1","event1","event1"),
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Efron1_3TTE <- summary(BuyseTest_Efron1_3TTE)
#### new
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf    delta   Delta n.bootstrap
# 1 global   100.00        27.05          31.25      20.59    21.12 -0.04200 -0.0420           0
# 3 global    41.71         5.38           6.19      15.66    14.48 -0.00810 -0.0501           0
# 5 global    30.13         6.66           7.09       8.95     7.43 -0.00429 -0.0544           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf    delta   Delta n.bootstrap
# 1 global   100.00        27.04          31.29      20.59    21.08 -0.04250 -0.0425           0
# 3 global    41.67         5.38           6.19      15.66    14.45 -0.00805 -0.0505           0
# 5 global    30.10         6.66           7.05       8.95     7.45 -0.00393 -0.0544           0

BuyseTest_Efron1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Efron",
                                  treatment="treatment",censoring="event1",
                                  type="TTE",threshold=0.25)
BuyseTestS_Efron1_1TTE <- summary(BuyseTest_Efron1_1TTE)
#### new
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        39.09          44.53       8.95     7.43 -0.0544 -0.0544           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        39.08          44.52       8.95     7.45 -0.0544 -0.0544           0

# check number of pairs
rowSums(BuyseTestS_Efron1_3TTE$nb[BuyseTestS_Efron1_3TTE$nb$strata=="global",c("n.neutral","n.uninf")])
BuyseTestS_Efron1_3TTE$nb[BuyseTestS_Efron1_3TTE$nb$strata=="global",c("n.total")]

# check equal survival
BuyseTestS_Efron1_3TTE$nb$Delta[BuyseTestS_Efron1_3TTE$nb$strata=="global"][3] # hum correspondance pas exacte
BuyseTestS_Efron1_1TTE$nb$Delta[BuyseTestS_Efron1_1TTE$nb$strata=="global"][1] # hum correspondance pas exacte

# bootstrap
system.time(
BuyseTest_Efron1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Efron",
                                  treatment="treatment",censoring="event1",
                                  type="TTE",threshold=0.25,n.bootstrap=1000,cpus=1)
)
# utilisateur     système      écoulé 
# 0.45        0.06      185.81 
BuyseTestS_Efron1_1TTE <- summary(BuyseTest_Efron1_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta CIinf.Delta CIsup.Delta n.bootstrap p.value   
# 1 global      100        39.09          44.53       8.95     7.43 -0.0544 -0.0544     -0.0831     -0.0257        1000   0.001 **
  
#### outcome 2
BuyseTest_Efron2_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint2","endpoint2","endpoint2"),method="Efron",
                                  treatment="treatment",censoring=c("event2","event2","event2"),
                                  type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Efron2_3TTE <- summary(BuyseTest_Efron2_3TTE)
#### new
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta Delta n.bootstrap
# 1 global   100.00        60.00           0.27      28.70    11.04 0.5970 0.597           0
# 3 global    39.74         9.74           0.62      22.25     7.12 0.0912 0.689           0
# 5 global    29.37        10.20           1.72      13.97     3.48 0.0848 0.773           0
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta Delta n.bootstrap
# 1 global   100.00        59.98           0.27      28.70    11.06 0.5970 0.597           0
# 3 global    39.75         9.74           0.62      22.25     7.13 0.0912 0.688           0
# 5 global    29.39        10.20           1.72      13.97     3.49 0.0848 0.773           0

BuyseTest_Efron2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Efron",
                                  treatment="treatment",censoring="event2",
                                  type="TTE",threshold=0.25)
BuyseTestS_Efron2_1TTE <- summary(BuyseTest_Efron2_1TTE)
#### old
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100        79.92           2.61      13.97     3.49 0.773 0.773           0
#### new
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100        79.94           2.61      13.97     3.48 0.773 0.773           0

# check equal survival
BuyseTestS_Efron2_3TTE$nb$Delta[BuyseTestS_Efron2_3TTE$nb$strata=="global"][3]
BuyseTestS_Efron2_1TTE$nb$Delta[BuyseTestS_Efron2_1TTE$nb$strata=="global"][1]

# bootstrap
BuyseTest_Efron2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Efron",
                                  treatment="treatment",censoring="event1",
                                  type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)
# BuyseTestS_Efron2_1TTE <- summary(BuyseTest_Efron2_1TTE)
# test.Efron <- list(M.Treatment=M.Treatment,
#                    M.Control=M.Control, 
#                    M.delta_Treatment=M.delta_Treatment, 
#                    M.delta_Control=M.delta_Control,
#                    list_survivalT=list_survivalT[[1]],
#                    list_survivalC=list_survivalC[[1]],
#                    iter_T=iter_T,
#                    iter_C=iter_C)
# save(test.Efron,file="E:/Creation_package/Package_BuyseTest/BuyseTest/test/test_Efron.RData")

#### Peron ####
BuyseTest_Peron1_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Peron",
                                   treatment="treatment",censoring=c("event1","event1","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peron1_3TTE <- summary(BuyseTest_Peron1_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf    delta   Delta n.bootstrap
# 1 global   100.00        27.05          30.65      20.59    21.72 -0.03600 -0.0360           0
# 3 global    42.31         5.38           5.59      15.66    15.68 -0.00205 -0.0381           0
# 5 global    31.34         6.66           5.84       8.95     9.88  0.00817 -0.0299           0

BuyseTest_Peron1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Peron",
                                   treatment="treatment",censoring="event1",
                                   type="TTE",threshold=0.25)
BuyseTestS_Peron1_1TTE <- summary(BuyseTest_Peron1_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        39.09          42.08       8.95     9.88 -0.0299 -0.0299           0

# check number of pairs
rowSums(BuyseTestS_Peron1_3TTE$nb[BuyseTestS_Peron1_3TTE$nb$strata=="global",c("n.neutral","n.uninf")])
BuyseTestS_Peron1_3TTE$nb[BuyseTestS_Peron1_3TTE$nb$strata=="global",c("n.total")]

# check equal survival
BuyseTestS_Peron1_3TTE$nb$Delta[BuyseTestS_Peron1_3TTE$nb$strata=="global"][3] # hum correspondance pas exacte
BuyseTestS_Peron1_1TTE$nb$Delta[BuyseTestS_Peron1_1TTE$nb$strata=="global"][1] # hum correspondance pas exacte

# bootstrap
system.time(
  BuyseTest_Peron1_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Peron",
                                     treatment="treatment",censoring="event1",
                                     type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)
)
# utilisateur     système      écoulé 
# 0.45        0.06      185.81 
BuyseTestS_Peron1_1TTE <- summary(BuyseTest_Peron1_1TTE)

#### outcome 2
BuyseTest_Peron2_3TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint2","endpoint2","endpoint2"),method="Peron",
                                   treatment="treatment",censoring=c("event2","event2","event2"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peron2_3TTE <- summary(BuyseTest_Peron2_3TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf  delta Delta n.bootstrap
# 1 global   100.00        60.00           0.27      28.70    11.04 0.5970 0.597           0
# 3 global    39.74         9.74           0.62      22.25     7.12 0.0912 0.689           0
# 5 global    29.37        10.20           1.71      13.97     3.49 0.0849 0.773           0
BuyseTestS_Peron2_3TTE <- summary(BuyseTest_Peron2_3TTE,type="nb")

BuyseTest_Peron2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint2",method="Peron",
                                   treatment="treatment",censoring="event2",
                                   type="TTE",threshold=0.25)
BuyseTestS_Peron2_1TTE <- summary(BuyseTest_Peron2_1TTE)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf delta Delta n.bootstrap
# 1 global      100        79.94           2.59      13.97     3.49 0.773 0.773           

# check equal survival
BuyseTestS_Peron2_3TTE$nb$Delta[BuyseTestS_Peron2_3TTE$nb$strata=="global"][3]
BuyseTestS_Peron2_1TTE$nb$Delta[BuyseTestS_Peron2_1TTE$nb$strata=="global"][1]

# bootstrap
BuyseTest_Peron2_1TTE <- BuyseTest(data=data_testTTE,endpoint="endpoint1",method="Peron",
                                   treatment="treatment",censoring="event1",
                                   type="TTE",threshold=0.25,n.bootstrap=1000,cpus=2)



#### repeted or different outcome ####
BuyseTest_Peron1_4TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint2","endpoint1"),method="Peron",
                                   treatment="treatment",censoring=c("event1","event2","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peron1_4TTE <- summary(BuyseTest_Peron1_4TTE)

data_testTTE$endpoint3 <- data_testTTE$endpoint1
data_testTTE$event3 <- data_testTTE$event1

BuyseTest_Peron1_4TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint1","endpoint1"),method="Peron",
                                   treatment="treatment",censoring=c("event1","event1","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peron1_4TTE <- summary(BuyseTest_Peron1_4TTE)

BuyseTest_Peron1_4TTE <- BuyseTest(data=data_testTTE,endpoint=c("endpoint1","endpoint3","endpoint1"),method="Peron",
                                   treatment="treatment",censoring=c("event1","event3","event1"),
                                   type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25))
BuyseTestS_Peron1_4TTE <- summary(BuyseTest_Peron1_4TTE)

#### 3- mixed outcome ####

n.Treatment_testCont <- 500
n.Control_testCont <- 500
mu.Treatment_testCont <- c(2,4)
mu.Control_testCont <- c(2,0)

set.seed(10)
data_testCont <- data.frame(treatment=c(rep(1,n.Treatment_testCont),rep(0,n.Treatment_testCont)))
data_testCont$endpoint1 <- c(rnorm(n.Treatment_testCont,mean=mu.Treatment_testCont[1]),
                             rnorm(n.Control_testCont,mean=mu.Control_testCont[1]))
data_testCont$endpoint2 <- c(rnorm(n.Control_testCont,mean=mu.Treatment_testCont[2]),
                             rnorm(n.Control_testCont,mean=mu.Control_testCont[2]))
data_testCont$strata <- rbinom(n.Treatment_testCont+n.Control_testCont,size=4,prob=0.5)

data_testMixed <- data.frame(data_testTTE[,c("treatment","strata","event1")],
                             endpointBin1=data_testBin[,"endpoint1"],
                             endpointCont1=data_testCont[,"endpoint1"],
                             endpointTTE1=data_testTTE[,"endpoint1"])


method <- "Peto" # "Peron" "Gehan" "Peto"

#### without strata ####
BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint="endpointTTE1",
                                 treatment="treatment",censoring="event1",
                                 type="timeToEvent",method=method,
                                 threshold=0.5)

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)
#### Gehan
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        14.62          17.64      15.66    52.09 -0.0302 -0.0302           0
#### Peto
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        32.75          37.18      15.66    14.41 -0.0443 -0.0443           0
#### Efron
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf   delta   Delta n.bootstrap
# 1 global      100        32.43          37.44      15.66    14.48 -0.0501 -0.0501           0

BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointTTE1","endpointCont1","endpointBin1"),
                                 treatment="treatment",censoring=c("event1",NA,NA),
                                 type=c("timeToEvent","continuous","binary"),method=method,
                                 threshold=c(0.5,1,NA))

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)
#### Gehan
# 1 global   100.00        14.62          17.64      15.66    52.09 -0.0302 -0.0302           0
# 3 global    67.74        15.48          17.11      35.15     0.00 -0.0163 -0.0466           0
# 5 global    35.15         8.02           9.57      17.56     0.00 -0.0155 -0.0621           0
#### Peto
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf     delta   Delta n.bootstrap
# 1 global   100.00        32.75          37.18      15.66    14.41 -0.044300 -0.0443           0
# 3 global    30.07         7.27           7.27      15.52     0.00 -0.000046 -0.0443           0
# 5 global    15.52         3.96           3.79       7.77     0.00  0.001710 -0.0426           0

BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointTTE1","endpointBin1","endpointCont1"),
                                 treatment="treatment",censoring=c("event1",NA,NA),
                                 type=c("timeToEvent","binary","continuous"),method=method,
                                 threshold=c(0.5,NA,1))

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)
#### Gehan
# 1 global   100.00        14.62          17.64      15.66    52.09 -0.03020 -0.0302           0
# 3 global    67.74        15.58          18.31      33.85     0.00 -0.02730 -0.0575           0
# 5 global    33.85         7.75           8.54      17.56     0.00 -0.00786 -0.0654           0
#### Peto
# 1 global   100.00        32.75          37.18      15.66    14.41 -4.43e-02 -0.0443           0
# 3 global    30.07         7.73           7.32      15.01     0.00  4.10e-03 -0.0402           0
# 5 global    15.01         3.62           3.62       7.77     0.00  1.19e-05 -0.0401           0

BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointBin1","endpointCont1","endpointTTE1"),
                                 treatment="treatment",censoring=c(NA,NA,"event1"),strata="strata",
                                 type=c("binary","continuous","timeToEvent"),method=method,
                                 threshold=c(NA,1,0.5))

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)
# strata pc.total pc.favorable pc.unfavorable pc.neutral pc.uninf    delta   Delta n.bootstrap
# 1  global   100.00        23.87          26.30      49.83     0.00 -0.02430 -0.0243           0
# 7  global    49.83        11.48          12.51      25.84     0.00 -0.01030 -0.0346           0
# 13 global    26.95         8.13           9.08       4.19     5.55 -0.00947 -0.0441           0




BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointTTE1","endpointCont1","endpointBin1"),
                                 treatment="treatment",censoring=c("event1",NA,NA),n.bootstrap=100,
                                 type=c("timeToEvent","continuous","binary"),method=method,
                                 threshold=c(0.5,1,NA))

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)




BuyseTest_testMixte <- BuyseTest(data=data_testMixed,
                                 endpoint=c("endpointTTE1","endpointBin1","endpointTTE1","endpointCont1","endpointTTE1"),
                                 treatment="treatment",censoring=c("event1",NA,"event1",NA,"event1"),
                                 type=c("timeToEvent","binary","timeToEvent","continuous","timeToEvent"),
                                 threshold=c(1,NA,0.5,0.5,0.25),method="Peron")

summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte)
summary_BuyseTest_testMixte <- summary(BuyseTest_testMixte,type="nb")

