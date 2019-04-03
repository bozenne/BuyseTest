#################################################################
#                                                               #
#         Extension of GPC to CR - Version 2 - 20181003         #      
#                                                               #
#################################################################

##### Goal : extension of the GPC method to take censoring into account (like in Peron 2016) but in the 
##### presence of competing risks (TTE data). This allows correcting the bias of Delta.
##### One outcome, no stratum

##### Version 2 : solves problem with cif non available to all times, that is transforms cif into a function

source("CalcAllPairsV2.R")
source("CalcStatisticV2.R")
library(survival)
library(cmprsk)

##### Main function

# ----------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------- 

GPCWithCR = function(endpoint, treatment, censoring, threshold, data, keep.PairScore = NULL, correction = NULL, 
                     print.res = NULL) {
  
  ## censoring : 0 (censoring), 1 (event of interest), 2 (competing event)
  
  if (levels(as.factor(data[,treatment]))[1] != 0 | levels(as.factor(data[,treatment]))[2] != 1) { 
    ## if treatment is not coded with 0 (control) and 1 (treatment)
    data[,treatment] = as.numeric(data[,treatment])-1
  }
  
  censoring.trt = data[data[,treatment]==1, censoring]
  censoring.ctrl = data[data[,treatment]==0, censoring]
  
  event.TimesT = data[data[,treatment] == 1, endpoint] # times of all events (including censoring) in treatment arm
  event.TimesC = data[data[,treatment] == 0, endpoint] # times of all events (including censoring) in control arm
  event1.TimesC = sort(unique(data[data[,treatment] == 0 & data[,censoring]==1, endpoint])) 
  # distinct times of event of interest in control arm
  cifTimeT = matrix(nrow = length(event.TimesT), ncol = 9)
  cifTimeC = matrix(nrow = length(event.TimesC), ncol = 9)
  cifJumpC = matrix(nrow = length(event1.TimesC), ncol = 4)
  cifTimeT[,1] = event.TimesT
  cifTimeC[,1] = event.TimesC
  cifJumpC[,1] = event1.TimesC
  
  cif = cuminc(ftime = data[,endpoint], fstatus = data[,censoring], group = data[,treatment], 
               cencode = 0)
  
  ###### create a function from cif to have cif values until infinity
  if (max(data[,censoring]) != 2) {
    cif3 = 0
    cif4 = 0
    sumC = cif[[1]]$est[length(cif[[1]]$est)] + cif3 
    sumT = cif[[2]]$est[length(cif[[2]]$est)] + cif4
    lastObsC = round(sumC,6) == 1
    lastObsT = round(sumT,6) == 1
    predCif1 = approxfun(x = unique(cif[[1]]$time), y = timepoints(cif,unique(cif[[1]]$time))$est[1,], yleft = 0,
                         yright = switch(as.character(lastObsC), "TRUE" = cif[[1]]$est[length(cif[[1]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
    predCif2 = approxfun(x = unique(cif[[2]]$time), y = timepoints(cif,unique(cif[[2]]$time))$est[2,], yleft = 0,
                         yright = switch(as.character(lastObsT), "TRUE" = cif[[2]]$est[length(cif[[2]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
    predCif3 = approxfun(x = unique(cif[[1]]$time), y = rep(0, length(unique(cif[[1]]$time))), yleft = 0,
                         yright = 0, f = 0, method = "constant")
    predCif4 = approxfun(x = unique(cif[[1]]$time), y = rep(0, length(unique(cif[[1]]$time))), yleft = 0,
                         yright = 0, f = 0, method = "constant")
  } else {
    sumC = cif[[1]]$est[length(cif[[1]]$est)] + cif[[3]]$est[length(cif[[3]]$est)] 
    sumT = cif[[2]]$est[length(cif[[2]]$est)] + cif[[4]]$est[length(cif[[4]]$est)]
    lastObsC = round(sumC,10) == 1
    lastObsT = round(sumT,10) == 1
    predCif1 = approxfun(x = unique(cif[[1]]$time), y = timepoints(cif,unique(cif[[1]]$time))$est[1,], yleft = 0,
                         yright = switch(as.character(lastObsC), "TRUE" = cif[[1]]$est[length(cif[[1]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
    predCif2 = approxfun(x = unique(cif[[2]]$time), y = timepoints(cif,unique(cif[[2]]$time))$est[2,], yleft = 0,
                         yright = switch(as.character(lastObsT), "TRUE" = cif[[2]]$est[length(cif[[2]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
    predCif3 = approxfun(x = unique(cif[[3]]$time), y = timepoints(cif,unique(cif[[3]]$time))$est[3,], yleft = 0,
                         yright = switch(as.character(lastObsC), "TRUE" = cif[[3]]$est[length(cif[[3]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
    predCif4 = approxfun(x = unique(cif[[4]]$time), y = timepoints(cif,unique(cif[[4]]$time))$est[4,], yleft = 0,
                         yright = switch(as.character(lastObsT), "TRUE" = cif[[4]]$est[length(cif[[4]]$est)], "FALSE" = NA),
                         f = 0, method = "constant")
  }
  
  if (threshold==0) {threshold = 10^(-12)}
  
  ###### cifTimeT
  cifTimeT[,2] = predCif1(event.TimesT - threshold)
  cifTimeT[,3] = predCif1(event.TimesT)
  cifTimeT[,4] = predCif1(event.TimesT + threshold)
  cifTimeT[,5] = predCif2(event.TimesT - threshold)
  cifTimeT[,6] = predCif2(event.TimesT)
  cifTimeT[,7] = predCif2(event.TimesT + threshold)
  cifTimeT[,8] = predCif3(event.TimesT)
  cifTimeT[,9] = predCif4(event.TimesT)
  
  ###### cifTimeC
  cifTimeC[,2] = predCif1(event.TimesC - threshold)
  cifTimeC[,3] = predCif1(event.TimesC)
  cifTimeC[,4] = predCif1(event.TimesC + threshold)
  cifTimeC[,5] = predCif2(event.TimesC - threshold)
  cifTimeC[,6] = predCif2(event.TimesC)
  cifTimeC[,7] = predCif2(event.TimesC + threshold)
  cifTimeC[,8] = predCif3(event.TimesC)
  cifTimeC[,9] = predCif4(event.TimesC)
  
  ###### cifJumpC
  cifJumpC[,2] = predCif2(event1.TimesC - threshold)
  cifJumpC[,3] = predCif2(event1.TimesC + threshold)
  cifJumpC[,4] = predCif1(event1.TimesC) - predCif1(event1.TimesC-10^(-12))
  
  ###### cif at tmax for each arm and each event
  lastCif1C = cif[[1]]$est[length(cif[[1]]$est)]
  lastCif1T = cif[[2]]$est[length(cif[[2]]$est)]
  lastCif2C = predCif3(max(event.TimesC))
  lastCif2T = predCif4(max(event.TimesT))
  # lastCif2C = cif[[3]]$est[length(cif[[3]]$est)]
  # lastCif2T = cif[[4]]$est[length(cif[[4]]$est)]
  
  # return(list(event.TimesT = event.TimesT, event.TimesC = event.TimesC, event1.TimesC = event1.TimesC,
  #             cifTimeT = cifTimeT, cifTimeC = cifTimeC, cifJumpC = cifJumpC))

  ###### Loop over the pairs
  if (is.null(correction)) {correction = F}
  if (is.null(keep.PairScore)) {keep.PairScore= F}
  res_list = CalcAllPairs_Peron(event.TimesT, event.TimesC, censoring.trt, censoring.ctrl,
                             threshold, cifTimeT, cifTimeC, cifJumpC, lastCif1C,
                             lastCif2C, lastCif1T, lastCif2T, keep.PairScore, correction)

  ###### Compute statistics
  stat = CalcStatistic(res_list$wins, res_list$losses, res_list$n.pairs)

  ###### Print results
  if(is.null(print.res)){print.res = T}
  if (print.res) {
    table = matrix(nrow = 1, ncol = 10)
    table = as.data.frame(table)
    colnames(table) = c("endpoint","threshold", "n.total","n.wins","n.losses","n.neutral.competing", 
                        "n.neutral.event","n.uninf", "weight", "delta")
    sum = res_list$wins + res_list$losses + res_list$neutral.competing + res_list$neutral.event + res_list$uninformative
    vec.to.print = c(endpoint, threshold, 100*res_list$n.pairs/res_list$n.pairs, round(100*res_list$wins/res_list$n.pairs,2), 
                     round(100*res_list$losses/res_list$n.pairs,2), round(100*res_list$neutral.competing/res_list$n.pairs,2),
                     round(100*res_list$neutral.event/res_list$n.pairs,2), 
                     round(100*res_list$uninformative/res_list$n.pairs,2), round(100*sum/res_list$n.pairs,2), 
                     round(stat$delta.netChance,4))
    table[1,] = vec.to.print 
    print(table, row.names = F)
  }
  
  return(c(res_list, stat))
}
