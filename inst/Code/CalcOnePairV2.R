#################################################################
#                                                               #
#         Extension of GPC to CR - Version 2 - 20181003         #      
#                                                               #
#################################################################

##### Goal : extension of the GPC method to take censoring into account (like in Peron 2016) but in the presence
##### of competing risks (TTE data). This allows correcting the bias of Delta.
##### One outcome, no stratum

##### Version 2 : solves problem with cif non available to all times, that is transforms cif into a function

##### Computation of the integral and of the 4 probabilities of ONE pair, via Peron

# ----------------------------------------------------------------------------------------------------
# Calculate integrals for pairs (0,0) 
# ---------------------------------------------------------------------------------------------------- 

CalcIntegral = function(cif, start.val, stop.val, CIF.t, lastCIF, type) {
  ## cif[1] : jump times in control group
  ## cif[2-3] : cif of the treatment group at times-tau and times+tau
  ## cif[4] : jump in cif of control group at times
  ## start.val : value at which to start the integral
  ## stop.val : value at which to stop the integral
  ## CIF.t : cif of the treatment group at current treatment time value
  ## lastCIF : cif value at tmax 
  ## type : 1 to compute probabilities of being fav 
  ##        2 to compute probabilities of being defav
  ##        3 to compute probabilities of being neutral due to event of interest (integral with t+tau and xi)
  ##        4 to compute probabilities of being neutral due to event of interest (integral with t+tau and t-tau)
  
  integral = 0
  nJump = nrow(cif)
  
  if (nJump > 0) {
    if(type == 1) {
      for (i in 1:nJump) {
        if(is.na(cif[i,3])) {break}
        if(cif[i,1] > start.val) {
          integral = integral + (lastCIF - cif[i, 3])*cif[i, 4]
        }
      }
    } else if(type == 2) {
      for (i in 1:nJump) {
        lb = cif[i,2]
        if(is.na(cif[i,2])) {lb = lastCIF}
        if(cif[i,1] > start.val) {
          integral = integral + (lb - CIF.t)*cif[i, 4]
        }
      }
    } else if(type == 3) {
      for (i in 1:nJump) {
        lb = cif[i,3]
        if(is.na(cif[i,3])) {lb = lastCIF}
        if(cif[i,1] > start.val & cif[i,1] <= stop.val) {
          integral = integral + (lb - CIF.t)*cif[i, 4]
        }
      }
    } else if(type == 4) {
      for (i in 1:nJump) {
        lb = cif[i, 3]
        if(is.na(cif[i,3])) {lb = lastCIF}
        if(is.na(cif[i,2])) {break}
        if(cif[i,1] > start.val) {
          integral = integral + (lb - cif[i, 2])*cif[i, 4]
        }
      }
    }
  }
  
  return(integral)
  
}

# ---------------------------------------------------------------------------------------------
# Calculate probabilities to be favorable, unfavorable, neutral and uninformative for one pair
# ---------------------------------------------------------------------------------------------

CalcOnePair_Peron = function(endpoint_T, endpoint_C, delta_T, delta_C, tau, index_T, index_C, cifTimeT,
                             cifTimeC, cifJumpC, lastCif1C, lastCif2C, lastCif1T, lastCif2T) {
  
  ## cifTimeC and cifTimeT: cumulative incidence at control/treatment observation times
  ##        [1]    times 
  ##        [2-4]  cif of event of interest estimated in the control arm: time - tau, time, time + tau
  ##        [5-7]  cif of event of interest estimated in the treatment arm: time - tau, time, time + tau
  ##        [8]  cif of competing event estimated in the control arm: time
  ##        [9]  cif of competing event estimated in the treatment arm: time
  
  ## cifJumpC: cumulative incidence of event of interest in control group at jump times
  ##        [1]  jump times in control group (unique values in ascending order)
  ##        [2-3]  cif of the treatment group at times-tau and times+tau
  ##        [4]  d(cif) of control group estimated at time 

  diff = endpoint_T - endpoint_C
  Cif1T.t = cifTimeT[index_T,6]
  proba = rep(0, 5) ## [1] favorable, [2] unfavorable, [3] neutral competing, [4] neutral event, [5] uninformative
  denomC = 1 - cifTimeC[index_C, 3] - cifTimeC[index_C, 8]
  denomT = 1 - cifTimeT[index_T, 6] - cifTimeT[index_T, 9]

  if(delta_T == 2) {
    if(delta_C == 2) { # (2,2)
      proba[3] = 1 # systematically neutral competing
    } else if(delta_C == 1){ # (2,1)
      proba[1] = 1 # systematically favorable
    } else if(delta_C == 0) { # (2,0)
        proba[1] = (lastCif1C - cifTimeC[index_C, 3])/denomC
        # proba[2] = 0
        proba[3] = (lastCif2C - cifTimeC[index_C, 8])/denomC
        # proba[4] = 0 # since the treated patient had the competing event
        proba[5] = 1 - (proba[1] + proba[3])
    }
  } else if(delta_T == 1){
    if(delta_C == 2){ # (1,2)
      proba[2] = 1 # systematically defavorable
    } else if(delta_C == 1){ # (1,1)
      if(diff >= tau) {
        proba[1] = 1
      } else if(diff <= -tau) {
        proba[2] = 1
      } else { # |diff| < tau
        proba[4] = 1
      }
    } else if(delta_C == 0) { # (1,0)
      if(diff >= tau) {
        if(!is.na(cifTimeT[index_T,2])) {
          proba[1] = (cifTimeT[index_T,2] - cifTimeC[index_C,3])/denomC
        } else {
          proba[1] = (lastCif1C - cifTimeC[index_C,3])/denomC
        }
        if(!is.na(cifTimeT[index_T,4])) {
          proba[2] = (lastCif1C - cifTimeT[index_T,4] + lastCif2C - cifTimeC[index_C,8])/denomC
        } else {
          # proba[2] = 0
          proba[2] = (lastCif2C - cifTimeC[index_C,8])/denomC
        }
        if(!is.na(cifTimeT[index_T,4]) & !is.na(cifTimeT[index_T,2])) {
          proba[4] = (cifTimeT[index_T,4] - cifTimeT[index_T,2])/denomC
        } else if (is.na(cifTimeT[index_T,4]) & !is.na(cifTimeT[index_T,2])) {
          proba[4] = (lastCif1C - cifTimeT[index_T,2])/denomC
        } else {
          proba[4] = 0 
        }
      } else if(diff <= -tau) {
        proba[2] = 1
      } else { # |diff| < tau
        if(!is.na(cifTimeT[index_T,4])) {
          proba[2] = (lastCif1C - cifTimeT[index_T,4] + lastCif2C - cifTimeC[index_C,8])/denomC
          proba[4] = (cifTimeT[index_T,4] - cifTimeC[index_C,3])/denomC
        } else {
          # proba[2] = 0
          proba[2] = (lastCif2C - cifTimeC[index_C,8])/denomC
          proba[4] = (lastCif1C - cifTimeC[index_C,3])/denomC
        }
      }
      proba[5] = 1 - (proba[1] + proba[2] + proba[3] + proba[4])
    }
  } else { # delta_T == 0
    if(delta_C == 2) { # (0,2)
      # proba[1] = 0
      proba[2] = (lastCif1T - cifTimeT[index_T,6])/denomT
      proba[3] = (lastCif2T - cifTimeT[index_T,9])/denomT
      # proba[4] = 0
      # proba[5] = 0
      proba[5] = 1 - (proba[2] + proba[3])
    } else if(delta_C == 1) { # (0,1)
      if(diff >= tau) {
        proba[1] = 1
        # proba[2] = 0
        # proba[3] = 0
        # proba[4] = 0
        # proba[5] = 0
      } else if(diff <= -tau) {
        if(!is.na(cifTimeC[index_C,7])) {
          proba[1] = (lastCif1T - cifTimeC[index_C,7] + lastCif2T - cifTimeT[index_T,9])/denomT
        } else {
          # proba[1] = 0
          proba[1] = (lastCif2T - cifTimeT[index_T,9])/denomT
        }
        if(!is.na(cifTimeC[index_C,5])) {
          proba[2] = (cifTimeC[index_C,5] - cifTimeT[index_T,6])/denomT
        } else {
          proba[2] = (lastCif1T - cifTimeT[index_T,6])/denomT
        }
        if(!is.na(cifTimeC[index_C,7]) & !is.na(cifTimeC[index_C,5])) {
          proba[4] = (cifTimeC[index_C,7] - cifTimeC[index_C,5])/denomT
        } else if (is.na(cifTimeC[index_C,7]) & !is.na(cifTimeC[index_C,5])) {
          proba[4] = (lastCif1T - cifTimeC[index_C,5])/denomT
        } else {
          proba[4] = 0 
        }
      } else { # |diff| < tau
        if(!is.na(cifTimeC[index_C,7])) {
          proba[1] = (lastCif1T - cifTimeC[index_C,7] + lastCif2T - cifTimeT[index_T,9])/denomT
          proba[4] = (cifTimeC[index_C,7] - cifTimeT[index_T,6])/denomT
        } else {
          # proba[1] = 0
          proba[1] = (lastCif2T - cifTimeT[index_T,9])/denomT
          proba[4] = (lastCif1T - cifTimeT[index_T,6])/denomT
        }
      }
      proba[5] = 1 - (proba[1] + proba[2] + proba[3] + proba[4])
    } else if (delta_C == 0) { # (0,0)
      prob21 = (lastCif2T - cifTimeT[index_T,9])*(lastCif1C - cifTimeC[index_C,3])/(denomT*denomC)
      prob12 = (lastCif2C - cifTimeC[index_C,8])*(lastCif1T - cifTimeT[index_T,6])/(denomT*denomC)
      prob22 = (lastCif2C - cifTimeC[index_C,8])*(lastCif2T - cifTimeT[index_T,9])/(denomT*denomC)
      if(diff >= tau) {
        # lower bound of each integral
        intFav = CalcIntegral(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T.t, lastCif1T, 1)/(denomT*denomC)
        intDefav = CalcIntegral(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T.t, lastCif1T, 2)/(denomT*denomC)
        intNeutralEvent1 = CalcIntegral(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T.t, lastCif1T, 3)/(denomT*denomC) 
        intNeutralEvent2 = CalcIntegral(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T.t, lastCif1T, 4)/(denomT*denomC)
        if (!is.na(cifTimeT[index_T,2])) {
          proba[1] = ((cifTimeT[index_T,2] - cifTimeC[index_C,3])/denomC)*((lastCif1T - cifTimeT[index_T,6])/denomT) + 
            intFav + prob21          
        } else {
          proba[1] = ((lastCif1C - cifTimeC[index_C,3])/denomC)*((lastCif1T - cifTimeT[index_T,6])/denomT) + 
            intFav + prob21
        }
        proba[2] = intDefav + prob12
        proba[3] = prob22
        proba[4] = intNeutralEvent1 + intNeutralEvent2
      } else if(diff <= -tau) {
        intFav = CalcIntegral(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T.t, lastCif1T, 1)/(denomT*denomC)
        intDefav = CalcIntegral(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T.t, lastCif1T, 2)/(denomT*denomC)
        intNeutralEvent = CalcIntegral(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T.t, lastCif1T, 4)/(denomT*denomC) 
        proba[1] = intFav + prob21
        proba[2] = intDefav + prob12
        proba[3] = prob22
        proba[4] = intNeutralEvent
      } else { # |diff| < tau
        intFav = CalcIntegral(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T.t, lastCif1T, 1)/(denomT*denomC)
        intDefav = CalcIntegral(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T.t, lastCif1T, 2)/(denomT*denomC)
        intNeutralEvent1 = CalcIntegral(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T.t, lastCif1T, 3)/(denomT*denomC)
        intNeutralEvent2 = CalcIntegral(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T.t, lastCif1T, 4)/(denomT*denomC)
        proba[1] = intFav + prob21
        proba[2] = intDefav + prob12
        proba[3] = prob22
        proba[4] = intNeutralEvent1 + intNeutralEvent2
      }
      proba[5] = 1 - (proba[1] + proba[2] + proba[3] + proba[4])
    }
  }
  return(proba)
} # end of function

  

