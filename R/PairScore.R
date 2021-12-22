### tablePairScore.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:54) 
## Version: 
## Last-Updated: Dec 21 2021 (17:01) 
##           By: Brice Ozenne
##     Update #: 137
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * pairScore2dt
## Convert output of .BuyseTest (list of vector) into a list of data.table
pairScore2dt <- function(pairScore,
                         level.treatment,
                         level.strata,
                         n.strata,
                         endpoint,
                         threshold,
                         restriction){
    
    ## Rcpp outputs vector: convert to matrix and rename
    name.tempo <- c("strata",
                    "index.C", "index.T", "index.pair",
                    "indexWithinStrata.C", "indexWithinStrata.T", 
                    "favorable","unfavorable","neutral","uninf",
                    "weight",
                    "favorableC","unfavorableC","neutralC","uninfC")
    p <- length(pairScore)
    
    pairScore2 <- vector(mode = "list", length = p)
    for(iL in 1:p){
        pairScore2[[iL]] <- data.table::as.data.table(matrix(pairScore[[iL]], ncol = 15, byrow = FALSE,
                                                             dimnames = list(NULL,name.tempo)))
        pairScore2[[iL]][, c("strata") := factor(.SD[["strata"]], levels = 0:(n.strata-1), labels = level.strata)] ## indexes start at 1 in R and not at 0 as in C++
        ## recall that indexes start at 1 in R and not at 0 as in C++
        pairScore2[[iL]][, c("index.C") := .SD$index.C + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        pairScore2[[iL]][, c("index.T") := .SD$index.T + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        pairScore2[[iL]][, c("index.pair") := .SD$index.pair + 1] 
        pairScore2[[iL]][, c("indexWithinStrata.T") := .SD$indexWithinStrata.T + 1]
        pairScore2[[iL]][, c("indexWithinStrata.C") := .SD$indexWithinStrata.C + 1]
        data.table::setkeyv(pairScore2[[iL]], c("index.T","index.C"))
    }
    names(pairScore2) <- paste0(endpoint,ifelse(!is.na(restriction),paste0("_r",restriction),""),ifelse(threshold>1e-12,paste0("_t",threshold),""))
    
    return(pairScore2)
}

