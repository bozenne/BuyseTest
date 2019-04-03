#################################################################
#                                                               #
#         Extension of GPC to CR - Version 2 - 20181003         #      
#                                                               #
#################################################################

##### Goal : extension of the GPC method to take censoring into account (like in Peron 2016) but in the 
##### presence of competing risks (TTE data). This allows correcting the bias of Delta.
##### One outcome, no stratum

##### Version 2 : solves problem with cif non available to all times, that is transforms cif into a function

##### Computation of the Delta statistic - returns Delta and win ratio

# ----------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------- 

CalcStatistic = function(fav, defav, n.pairs) {
  delta_netChance = (fav - defav)/n.pairs
  delta_winRatio = fav/defav
  return (list(delta.netChance = delta_netChance, delta.winRatio = delta_winRatio))
}

