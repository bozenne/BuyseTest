#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Test BuyseTest with strata
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# check whether the number of pairs displayed in the summary is ok (e.g. total = favorable + unfavorable + neutral + uninf)
# check whether when using identical values across strata, the result is the same across strata
# check consistency with results from previous version (on common slots)

#### spec
library(BuyseTest) # butils:::package.source("BuyseTest", Rcode = TRUE, Ccode = TRUE)
library(testthat)
library(lava)
library(data.table)

precision <- 10^{-7}
n.patients <- c(90,100)
n.bootstrap <- 10

save <- NULL # TRUE to save results, FALSE to test, NULL to ignore
conv2df <- FALSE

#### 0- function - file FCT/FCT_check.R ####
validPairs <- function(BuyseRes, type = c("strata","sum")){
  
  BuyseSummary <- summary(BuyseRes, show = NULL)
  enpoint_threshold <- paste(BuyseSummary$nb$endpoint,BuyseSummary$nb$threshold, sep = "_")
  endpoints <- unique(enpoint_threshold)
  D <- length(endpoints)
  index_strata <- which(BuyseSummary$nb$strata!="global")
  levels.strata <- unique(BuyseSummary$nb$strata[index_strata])
  n.strata <- length(levels.strata) 
  
  diff <- NULL
  
  if("strata" %in% type){
    diff.strata <- matrix(NA,nrow=(1+n.strata)*D,ncol=5)
    colnames(diff.strata) <- c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")
    rownames(diff.strata) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    for(iter_endpoint in 1:D){
      index_endpoint <- which(enpoint_threshold==endpoints[iter_endpoint])  
      res_strata <- as.matrix(BuyseSummary$nb[intersect(index_strata,index_endpoint),c("n.total","n.favorable","n.unfavorable","n.neutral","n.uninf")])
      diff.strata[intersect(index_strata,index_endpoint),] <- t(apply(res_strata,1,function(x){x-res_strata[1,,drop=TRUE]}))
    }
    diff <- cbind(diff, diff.strata)
  }
  if("sum" %in% type){
    diff.sum <- cbind(sum = BuyseSummary$nb$n.total - (BuyseSummary$nb$n.favorable+BuyseSummary$nb$n.unfavorable+BuyseSummary$nb$n.neutral+BuyseSummary$nb$n.uninf))
    rownames(diff.sum) <- paste(BuyseSummary$nb$endpoint, " th=",BuyseSummary$nbthreshold ," strata=", BuyseSummary$nb$strata)
    diff <- cbind(diff, diff.sum)
  }
  return(diff)
}

expect_equalPairsBT <- function(BuyseRes1, BuyseRes2){
  count1 <- getCount(BuyseRes1) 
  count2 <- getCount(BuyseRes2)
  expect_equal(count1, count2)
}

expect_equalBT <- function(BuyseRes1, BuyseRes2, slots = NULL, trace = 1){
  if(is.null(slots)){slots <- setdiff(intersect(names(attributes(BuyseRes1)), names(attributes(BuyseRes2))), "class")}
  
  test.strata1 <- !("try-error" %in% class(try(BuyseRes1@strata, silent = TRUE)))
  test.strata2 <- !("try-error" %in% class(try(BuyseRes2@strata, silent = TRUE)))
  
  if (test.strata1 && length(BuyseRes1@strata) > 1) {
    if(test.strata2 == FALSE){
      
      if(!identical(BuyseRes1@delta, BuyseRes2@delta)){
        if(trace>0){cat("* reorder according to strata in BuyseRes 1 and 2 \n")}
        index.match <- match(BuyseRes1@delta[,1], BuyseRes2@delta[,1]) # could be wrong to consider only the first outcome
        BuyseRes2@delta <- BuyseRes2@delta[index.match,,drop=FALSE]
        BuyseRes2@count_favorable <- BuyseRes2@count_favorable[index.match,,drop=FALSE]
        BuyseRes2@count_unfavorable <- BuyseRes2@count_unfavorable[index.match,,drop=FALSE]
        BuyseRes2@count_neutral <- BuyseRes2@count_neutral[index.match,,drop=FALSE]
        BuyseRes2@count_uninf <- BuyseRes2@count_uninf[index.match,,drop=FALSE]
        BuyseRes2@delta_boot <- BuyseRes2@delta_boot[index.match,,,drop = FALSE]
        
        BuyseRes2@index_neutralT <- sort(BuyseRes2@index_neutralT)
        BuyseRes2@index_neutralC <- sort(BuyseRes2@index_neutralC)
        BuyseRes2@index_uninfT <- sort(BuyseRes2@index_uninfT)
        BuyseRes2@index_uninfC <- sort(BuyseRes2@index_uninfC)
        
        BuyseRes1@index_neutralT <- sort(BuyseRes1@index_neutralT)
        BuyseRes1@index_neutralC <- sort(BuyseRes1@index_neutralC)
        BuyseRes1@index_uninfT <- sort(BuyseRes1@index_uninfT)
        BuyseRes1@index_uninfC <- sort(BuyseRes1@index_uninfC)
      }
    }
    
  }
  
  test.error <- FALSE
  for(iterSlot in slots){
    res <- try(expect_equal(slot(BuyseRes1,  iterSlot), slot(BuyseRes2,  iterSlot)), silent = TRUE)
    if("try-error" %in% class(res) ){
      test.error <- TRUE
      cat("Differences in slot: ",iterSlot,"\n")
      if(trace == 2){cat(res[1])}
    }
  }
  if(test.error){stop("difference between the slots of the two BuyseRes objects \n")}
}


Vexpect_less_than <- function(x,y,...){
  sapply(x, function(X){expect_less_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_more_than <- function(x,y,...){
  sapply(x, function(X){expect_more_than(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_equal <- function(x,y,...){
  sapply(x, function(X){expect_equal(X,y,...)})
  return(invisible(TRUE))
}
Vexpect_NA <- function(x,...){
  sapply(x, function(X){expect_true(is.na(X),...)})
  return(invisible(TRUE))
}

#### 1- Check number of pairs ####
set.seed(10)
argsBin <- list(p.T = c(0.5,0.75))
argsCont <- list(mu.T = 1:3, sigma.T = rep(1,3))
argsTTE <- list(rates.T = 1:3, rates.Censor = rep(1,3))

#### binary #### 
data_Bin <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = NULL, argsTTE = NULL)
data_BinS <- rbind(cbind(data_Bin, strata = 1),
                   cbind(data_Bin, strata = 2),
                   cbind(data_Bin, strata = 3))
if(conv2df){data_BinS <- as.data.frame(data_BinS)}


BT_Bin1 <- BuyseTest(data=data_BinS,endpoint=c("Y_bin1","Y_bin2"),
                     treatment="Treatment", type=c("bin","bin"),strata="strata",
                     n.bootstrap=n.bootstrap,trace=0)
BT_Bin1

test_that("count pairs summary - Binary",{
  valTest <- as.double(validPairs(BT_Bin1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Binary",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Bin1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### continuous ####
set.seed(10)
data_Cont <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = argsCont, argsTTE = NULL)
data_ContS <- rbind(cbind(data_Cont, strata = 1),
                        cbind(data_Cont, strata = 2),
                        cbind(data_Cont, strata = 3))
if(conv2df){data_ContS <- as.data.frame(data_ContS)}

BT_Cont1 <- BuyseTest(data=data_ContS,endpoint=c("Y_cont1","Y_cont2"),
                                 treatment="Treatment", type=c("cont","cont"),threshold=c(1,1),strata="strata",
                                 n.bootstrap=n.bootstrap,trace=0)

summary(BT_Cont1)

test_that("count pairs summary - Continuous",{
  valTest <- as.double(validPairs(BT_Cont1, type = "sum"))
  expect_equal(valTest, rep(0, times = length(valTest)))
})
test_that("identical strata - Continuous",{
  valTest <- as.vector(na.omit(as.double(validPairs(BT_Cont1, type = "strata"))))
  expect_equal(valTest, rep(0, times = length(valTest)))
})

#### TTE ####
set.seed(10)
data_TTE <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = NULL, argsCont = NULL, argsTTE = argsTTE)
data_TTES <- rbind(cbind(data_TTE, strata = 1),
                   cbind(data_TTE, strata = 2),
                   cbind(data_TTE, strata = 3))
if(conv2df){data_TTES <- as.data.frame(data_TTES)}

## different endpoints
BT_TTE1 <- vector(length = 4, mode = "list")
names(BT_TTE1) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE1[[method]] <- BuyseTest(data=data_TTES,endpoint=c("Y_TTE1","Y_TTE2","Y_TTE3"),method=method,
                                 treatment="Treatment",censoring=c("event1","event2","event1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(0.75,0.5,0.25),
                                 n.bootstrap = n.bootstrap, trace=0)
  (BT_TTE1[[method]])
  
  test_that("count pairs summary - TTE",{
    valTest <- as.double(validPairs(BT_TTE1[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - TTE",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_TTE1[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

## same endpoint
BT_TTE2 <- vector(length = 4, mode = "list")
names(BT_TTE2) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_TTE2[[method]] <- BuyseTest(data=data_TTES,endpoint=c("Y_TTE1","Y_TTE1","Y_TTE1"),method=method,
                                 treatment="Treatment",censoring=c("event1","event1","event1"),strata="strata",
                                 type=c("TTE","TTE","TTE"),threshold=c(1,0.5,0.25),
                                 n.bootstrap = n.bootstrap, trace=0)
  (BT_TTE2[[method]])
  
  test_that("count pairs summary - TTE",{
    valTest <- as.double(validPairs(BT_TTE2[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - TTE",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_TTE2[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### mixed outcomes ####
set.seed(10)
data_Mix <- simulBT(n.T = n.patients[1], n.C = n.patients[2], argsBin = argsBin, argsCont = argsCont, argsTTE = argsTTE)
data_MixS <- rbind(cbind(data_Mix, strata = 1),
                   cbind(data_Mix, strata = 2),
                   cbind(data_Mix, strata = 3))
if(conv2df){data_MixS <- as.data.frame(data_MixS)}

BT_Mix <- vector(length = 4, mode = "list")
names(BT_Mix) <- c("Gehan","Peto","Efron","Peron")
for(method in c("Gehan","Peto","Efron","Peron")){
  BT_Mix[[method]] <- BuyseTest(data=data_MixS,
                                endpoint=c("Y_TTE1","Y_cont1","Y_bin1","Y_TTE1","Y_cont1"),method=method,
                                treatment="Treatment",censoring=c("event1",NA,NA,"event1",NA),strata="strata",
                                type=c("timeToEvent","continuous","binary","timeToEvent","continuous"),
                                threshold=c(0.5,1,NA,0.25,0.5),
                                n.bootstrap = n.bootstrap, trace=0)
  
  SBT_Mix <- summary(BT_Mix[[method]])
  
  test_that("count pairs summary - mixed",{
    valTest <- as.double(validPairs(BT_Mix[[method]], type = "sum"))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
  test_that("identical strata - mixed",{
    valTest <- as.vector(na.omit(as.double(validPairs(BT_Mix[[method]], type = "strata"))))
    expect_equal(valTest, rep(0, times = length(valTest)))
  })
}

#### Real data ####
data(veteran, package = "survival")

BT_veteran <- BuyseTest(data = veteran, endpoint = "time", treatment = "trt", strata = "celltype",
                        type = "timeToEvent", censoring = "status",threshold = 0, 
                        n.bootstrap = n.bootstrap, trace = 0)
BT_veteran

##### export
version <- packageVersion("BuyseTest")
dir <- paste0("tests/Results-version",version)

if(!is.null(save) && save == TRUE){
  results_strata <- list(OutcomeBin = list(data = data_BinS, BT = BT_Bin1),
                         OutcomeCont = list(data = data_ContS, BT = BT_Cont1),
                         OutcomeTTE = list(data = data_TTES, BT1 = BT_TTE1, BT2 = BT_TTE2), 
                         OutcomeMix = list(data = data_MixS, BT = BT_Mix),
                         veteran = list(data = veteran, BT = BT_veteran))
  if(dir.exists(dir) == FALSE){dir.create(dir)}
  save(results_strata, file = file.path(dir,"test_strata.RData"))
}else if(!is.null(save) && save == FALSE){
  load(file = file.path(dir,"test_strata.RData"))
  
  test_that("comparison with the previous version", {
    expect_equalBT(BT_Bin1, results_strata$OutcomeBin$BT)
    expect_equalBT(BT_Cont1, results_strata$OutcomeCont$BT)
    expect_equalBT(BT_TTE1$Gehan, results_strata$OutcomeTTE$BT1$Gehan)
    expect_equalBT(BT_TTE1$Peto, results_strata$OutcomeTTE$BT1$Peto)
    expect_equalBT(BT_TTE1$Efron, results_strata$OutcomeTTE$BT1$Efron)
    expect_equalBT(BT_TTE1$Peron, results_strata$OutcomeTTE$BT1$Peron)
    expect_equalBT(BT_TTE2$Gehan, results_strata$OutcomeTTE$BT2$Gehan)
    expect_equalBT(BT_TTE2$Peto, results_strata$OutcomeTTE$BT2$Peto)
    expect_equalBT(BT_TTE2$Efron, results_strata$OutcomeTTE$BT2$Efron)
    expect_equalBT(BT_TTE2$Peron, results_strata$OutcomeTTE$BT2$Peron)
    
    expect_equalBT(BT_Mix$Gehan, results_strata$OutcomeMix$BT$Gehan)
    expect_equalBT(BT_Mix$Peto, results_strata$OutcomeMix$BT$Peto)
    expect_equalBT(BT_Mix$Efron, results_strata$OutcomeMix$BT$Efron)
    expect_equalBT(BT_Mix$Peron, results_strata$OutcomeMix$BT$Peron)
  
    expect_equalBT(BT_veteran, results_strata$veteran$BT)
  })
}