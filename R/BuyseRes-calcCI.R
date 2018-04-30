### BuyseRes-calcCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 29 2018 (17:42) 
## Version: 
## Last-Updated: apr 29 2018 (17:42) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcCI
#' @rdname internal-computation
calcCI <- function(Delta,Delta.permutation,
                   endpoint,D,alternative,alpha,
                   n.resampling,cpus,trace){

    option <- BuyseTest.options()
        
    ## Confidence interval
    ## Delta.permutation <- apply(delta.permutation,c(2,3),sum) # sum of the proportion in favor of treatment over the strata for the permutation test.
    ## if (D > 1) {Delta.permutation <- apply(Delta.permutation,2,cumsum)} # cumulative sum over the endpoints to obtaine the cumulative proportion in favor of treatment for the permutation test
    n.resampling_real <- lapply(Delta.permutation, function(x){n.resampling - apply(is.na(x), 1, sum)})
  
    if (trace && any(unlist(n.resampling_real) < n.resampling)) {
        message <- "BuyseTest : for some samples, the permutation procedure failed and returned NA \n the sample will be ignore in the finale analysis \n"
        indexD <- lapply(n.resampling_real, function(x){which(x<n.resampling)})
    
        for (iter_endpoint in sort(unique(unlist(indexD)))) {
            message <- paste(message,
                             paste("endpoint ",iter_endpoint," (\"",endpoint[iter_endpoint],"\") : ",
                                   sum(is.na(Delta.permutation[[option$statistic]][iter_endpoint,])),
                                   " (",100*sum(is.na(Delta.permutation[[option$statistic]][iter_endpoint,]))/n.resampling,"%) failures  \n",sep = ""),
                             sep = " ")
        }
    
        if (cpus > 1 && ((n.resampling/cpus) %% 1 != 0)) {
            message <- paste(message,"probable cause : some iterations were not launched \n",
                             "because n.resampling (=",n.resampling,") is not a multiple of cpus ",cpus," \n",sep = "") 
        }else{
            message <- paste(message,"possible cause : resampling leaded to no control or no case \n",sep = "")
        }
    
        warning(message)
    
    }
    
    ## ** compute the quantiles for each endpoint of the cumulative proportions in favor of treatment (in the sample)
    Delta.quantile <- lapply(1:2, function(type){
        apply(Delta.permutation[[type]],1,function(x){stats::quantile(x,probs = c(alpha/2,1 - alpha/2),na.rm = TRUE)})
    })
    names(Delta.quantile) <- names(Delta)
  
    ## ** p. value (we don't really need to differientiate netChance and winRatio)
    testp.value <- switch(alternative, # test whether each sample is has a cumulative proportions in favor of treatment more extreme than the punctual estimate
                          "two.sided" = lapply(1:2, function(type){
                              apply(Delta.permutation[[type]], 2,function(x){abs(Delta[[type]] - c(0,0.5)[type]) > abs(x - c(0,0.5)[type])})
                          }),
                          "less" = lapply(1:2, function(type){
                              apply(Delta.permutation[[type]], 2,function(x){(Delta[[type]] - c(0,0.5)[type]) < (x - c(0,0.5)[type])})
                          }),
                          "more" =  lapply(1:2, function(type){
                              apply(Delta.permutation[[type]], 2,function(x){(Delta[[type]] - c(0,0.5)[type]) > (x - c(0,0.5)[type])})
                          })
                          )
  
    ## ** mean over the samples by endpoint
    ## (i.e. proportion of samples in which the statistic obtained after permutation is more extreme than the observed one)
    if (D == 1) {
        p.value <- lapply(testp.value, function(x){1 - mean(x, na.rm = TRUE)})
    }else{
        p.value <- lapply(testp.value, function(x){1 - rowMeans(x, na.rm = TRUE)})
    }
    names(p.value) <- names(Delta)
  
    ## ** export  
    res <- list()
    res$p.value <- p.value
    res$Delta.quantile <- Delta.quantile
    res$n.resampling <- n.resampling_real
    return(res)
}



##----------------------------------------------------------------------
### BuyseRes-calcCI.R ends here
