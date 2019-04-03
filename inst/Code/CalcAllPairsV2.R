#################################################################
#                                                               #
#         Extension of GPC to CR - Version 2 - 20181003         #      
#                                                               #
#################################################################

##### Goal : extension of the GPC method to take censoring into account (like in Peron 2016) but in the 
##### presence of competing risks (TTE data). This allows correcting the bias of Delta.
##### One outcome, no stratum

##### Version 2 : solves problem with cif non available to all times, that is transforms cif into a function

##### Pairwise correction 

# ----------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------- 

source("CalcOnePairV2.R")

correctionPairwise = function(countFav, countUnfav, countNeutralComp, countNeutralEvent, countUninf, score) {
  
  ## compute factors
  countInf = countFav + countUnfav + countNeutralComp + countNeutralEvent
  factorFav = countFav/countInf
  factorUnfav = countUnfav/countInf
  factorNeutralComp = countNeutralComp/countInf
  factorNeutralEvent = countNeutralEvent/countInf
  factors = c(factorFav = factorFav, factorUnfav = factorUnfav, factorNeutralComp = factorNeutralComp, 
              factorNeutralEvent = factorNeutralEvent)
  
  ## update global scores
  countFav = countFav + factorFav*countUninf
  countUnfav = countUnfav + factorUnfav*countUninf
  countNeutralComp = countNeutralComp + factorNeutralComp*countUninf
  countNeutralEvent = countNeutralEvent + factorNeutralEvent*countUninf
  countUninf = 0
  
  ## update individual weights of pairs partially or totally uninformative
  score[,9] = score[,3] + factorFav*score[,7]
  score[,10] = score[,4] + factorUnfav*score[,7]
  score[,11] = score[,5] + factorNeutralComp*score[,7]
  score[,12] = score[,6] + factorNeutralEvent*score[,7]
  
  return(list(countFav = countFav, countUnfav = countUnfav, countNeutralComp = countNeutralComp,
              countNeutralEvent = countNeutralEvent, countUninf = countUninf, factors = factors, score = score))
}

##### Computation of the 4 probabilities of ALL pairs, via Peron

# ----------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------- 

CalcAllPairs_Peron = function(endpoint.trt, endpoint.ctrl, delta_T, delta_C, tau, cifTimeT, cifTimeC, 
                              cifJumpC, lastCif1C, lastCif2C, lastCif1T, lastCif2T, keep.PairScore, correction) {
  
  n.trt = length(endpoint.trt)
  n.ctrl = length(endpoint.ctrl)
  n.pairs = n.trt*n.ctrl
  zeroPlus = 10^(-12)
  iter_pair = 0
  countFav = 0
  countUnfav = 0
  countNeutralComp = 0
  countNeutralEvent = 0
  countUninf = 0
  index_neutralCompT = c()
  index_neutralCompC = c()
  index_neutralEventT = c()
  index_neutralEventC = c()
  index_uninfT = c()
  index_uninfC = c()
  proba = rep(0,5) ## [1] favorable, [2] unfavorable, [3] neutral competing, 
                   ## [4] neutral event, [5] uninformative
  score = matrix(nrow = n.pairs, ncol = 12) 
  
  for (i in 1:n.trt) {
    for (j in 1:n.ctrl) {
      iter_pair = iter_pair + 1
      proba = CalcOnePair_Peron(endpoint.trt[i], endpoint.ctrl[j], delta_T[i], delta_C[j], 
                                tau, i, j, cifTimeT, cifTimeC, cifJumpC, 
                                lastCif1C, lastCif2C, lastCif1T, lastCif2T)
      countFav = countFav + proba[1]
      countUnfav = countUnfav + proba[2]
      
      if(proba[3] > zeroPlus) {
        countNeutralComp = countNeutralComp + proba[3]
        # index_neutralCompT = c(index_neutralCompT, i)
        # index_neutralCompC = c(index_neutralCompC, j)
      }
      if(proba[4] > zeroPlus) {
        countNeutralEvent = countNeutralEvent + proba[4]
        # index_neutralEventT = c(index_neutralEventT, i)
        # index_neutralEventC = c(index_neutralEventC, j)
      }
      if(proba[5] > zeroPlus) {
        countUninf = countUninf + proba[5]
        # index_uninfT = c(index_uninfT, i)
        # index_uninfC = c(index_uninfC, j)
      }
      if (keep.PairScore) {score[iter_pair,] = c(i, j, proba, sum(proba), proba[1:4])}
    }
  }
  
  if (correction) {
    correctedCounts = correctionPairwise(countFav, countUnfav, countNeutralComp, countNeutralEvent, countUninf, score)
    countFav = correctedCounts$countFav; countUnfav = correctedCounts$countUnfav; 
    countNeutralComp = correctedCounts$countNeutralComp; countNeutralEvent = correctedCounts$countNeutralEvent;
    countUninf = correctedCounts$countUninf; score = correctedCounts$score
  }
  
  if (keep.PairScore) {
    colnames(score) = c("indexT", "indexC", "Wins", "Losses", "Neutral_Competing", "Neutral_Event", "Uninformative", "Weight", 
                        "Wins.corrected", "Losses.corrected", "Neutral_Competing.corrected", "Neutral_Event.corrected")
    return(list(wins = countFav, losses = countUnfav, neutral.competing = countNeutralComp, neutral.event = countNeutralEvent, 
                uninformative = countUninf, n.pairs = iter_pair, PairScore = score))
  } else {
    return(list(wins = countFav, losses = countUnfav, neutral.competing = countNeutralComp, neutral.event = countNeutralEvent, 
                uninformative = countUninf, n.pairs = iter_pair))
  }
  
}
