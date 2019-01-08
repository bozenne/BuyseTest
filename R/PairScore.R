### tablePairScore.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:54) 
## Version: 
## Last-Updated: jan  8 2019 (10:31) 
##           By: Brice Ozenne
##     Update #: 94
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
                         threshold){
    
    ## Rcpp outputs vector: convert to matrix and rename
    name.tempo <- c("strata",
                    "index.C", "index.T", 
                    "indexWithinStrata.C", "indexWithinStrata.T", 
                    "favorable","unfavorable","neutral","uninf",
                    "weight",
                    "favorableC","unfavorableC","neutralC","uninfC")
    pairScore2 <- lapply(pairScore, function(iC){ ## iC <- pairScore[[1]]
        iM <- data.table::as.data.table(matrix(iC, ncol = 14, byrow = FALSE,
                                               dimnames = list(NULL,name.tempo)))
        iM[, c("strata") := factor(.SD[["strata"]], levels = 0:(n.strata-1), labels = level.strata)] ## indexes start at 1 in R and not at 0 as in C++
        ## recall that indexes start at 1 in R and not at 0 as in C++
        iM[, c("index.C") := .SD$index.C + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("index.T") := .SD$index.T + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("indexWithinStrata.T") := .SD$indexWithinStrata.T + 1]
        iM[, c("indexWithinStrata.C") := .SD$indexWithinStrata.C + 1]
        return(iM[])
    })
    names(pairScore2) <- paste0(endpoint,"_",threshold)

    return(pairScore2)
}
