### BuyseRes-calcCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 29 2018 (17:42) 
## Version: 
## Last-Updated: maj 19 2018 (15:51) 
##           By: Brice Ozenne
##     Update #: 55
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
                                       return(mean((Delta[iE] - null) > (Delta.permutation[iE,] - null)))
                                   }),
                                   "greater" = sapply(endpoint, function(iE){
                                       return(mean((Delta[iE] - null) < (Delta.permutation[iE,] - null)))
                                   })
                               )
    
    ## ** export  
    return(out)
}

calcCIbootstrap <- function(Delta, Delta.permutation,                            
                            statistic, endpoint,
                            alternative, alpha, method, n.pairs){

    method <- match.arg(method, c("percentile","gaussian"))
    
    out <- list()
   
    ## ** number of permutations
    out$n.resampling_real <- rowSums(!is.na(Delta.permutation))

    ## ** null hypothesis
    null <- switch(statistic,
                   "netChance" = 0,
                   "winRatio" = 1)

    ## ** standard error
    vec.se <- sapply(endpoint, function(iE){
        stats::sd(Delta.permutation[iE,], na.rm = TRUE)            
    })
    
    ## ** Confidence intervals
    if(method == "percentile"){
        out$Delta.CI <- sapply(endpoint,function(iE){ ## iE <- 1
            qDelta <- switch(alternative,
                             "two.sided" = stats::quantile(Delta.permutation[iE,], probs = c(alpha/2,1 - alpha/2),na.rm = TRUE),
                             "less" = c(stats::quantile(Delta.permutation[iE,], probs = alpha,na.rm = TRUE), Inf),
                             "greater" = c(-Inf,stats::quantile(Delta.permutation[iE,], probs = 1 - alpha,na.rm = TRUE))
                             )
            return(qDelta)
        })
    }else if(method == "gaussian"){
        out$Delta.CI <- sapply(endpoint,function(iE){ ## iE <- 1
           
            qDelta <- switch(alternative,
                             "two.sided" = Delta[iE] + stats::qnorm(c(alpha/2,1 - alpha/2)) * vec.se[iE],
                             "less" = c(Delta[iE] + stats::qnorm(alpha) * vec.se[iE], Inf),
                             "greater" = c(-Inf,Delta[iE] + stats::qnorm(1-alpha) * vec.se[iE])
                             )
            return(qDelta)
        })
    }

    ## ** p.value
    if(method == "percentile"){
        out$Delta.pvalue <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                                   "two.sided" = sapply(endpoint, function(iE){
                                       if(Delta[iE]>null){
                                           return(mean(Delta.permutation[iE,] > null))
                                       }else if(Delta[iE]<null){
                                           return(mean(Delta.permutation[iE,] < null))
                                       }else{
                                           return(1)
                                       }
                                   }),
                                   "less" = sapply(endpoint, function(iE){
                                       return(mean(Delta.permutation[iE,] < null))
                                   }),
                                   "greater" = sapply(endpoint, function(iE){
                                       return(mean(Delta.permutation[iE,] > null))
                                   })
                                   )
    }else if(method == "gaussian"){
        out$Delta.p.value <- sapply(endpoint,function(iE){ ## iE <- 1
           
            pDelta <- switch(alternative,
                             "two.sided" = 2*(1-stats::pnorm(abs(Delta[iE]/vec.se[iE]))), ## 2*(1-pnorm(1.96))
                             "less" = stats::pnorm(Delta[iE]/vec.se[iE]), ## pnorm(1.96)
                             "greater" = 1-stats::pnorm(Delta[iE]/vec.se[iE])
                             )
            return(pDelta)
        })

    }
    return(out)
}



##----------------------------------------------------------------------
### BuyseRes-calcCI.R ends here
