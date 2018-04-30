### BuyseRes-calcCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 29 2018 (17:42) 
## Version: 
## Last-Updated: apr 30 2018 (15:57) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcCIpermutation
calcCIpermutation <- function(Delta, Delta.permutation,
                              statistic, endpoint,
                              alternative, alpha){

    out <- list()
    
    ## ** number of permutations
    out$n.resampling_real <- rowSums(!is.na(Delta.permutation))
    
    ## ** null hypothesis
    null <- switch(statistic,
                   "netChance" = 0,
                   "winRatio" = 1)

    ## ** compute the quantiles for each endpoint of the cumulative proportions in favor of treatment (in the sample)

    out$Delta.CI <- sapply(endpoint,function(iE){ ## iE <- 1
        qDelta <- switch(alternative,
                         "two.sided" = stats::quantile(Delta.permutation[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                         "less" = c(stats::quantile(Delta.permutation[iE,], probs = alpha,na.rm = TRUE), Inf),
                         "greater" = c(-Inf,stats::quantile(Delta.permutation[iE,], probs = 1 - alpha,na.rm = TRUE))
                         )
        return(Delta[iE] + (qDelta-null))
    })
  
    ## ** p.value
    out$Delta.pvalue <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                           "two.sided" = sapply(endpoint, function(iE){
                               mean(abs(Delta[iE] - null) < abs(Delta.permutation[iE,] - null))
                           }),
                           "less" = sapply(endpoint, function(iE){
                               if(Delta[iE]<null){
                                   return(mean((Delta[iE] - null) > (Delta.permutation[iE,] - null)))
                               }else{
                                   return(1)
                               }
                           }),
                           "greater" = sapply(endpoint, function(iE){
                               if(Delta[iE]>null){
                                   return(mean((Delta[iE] - null) < (Delta.permutation[iE,] - null)))
                               }else{
                                   return(1)
                               }
                           })
                           )

    
    ## ** export  
    return(out)
}



##----------------------------------------------------------------------
### BuyseRes-calcCI.R ends here
