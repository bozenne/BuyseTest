### tableIndividualScore.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:54) 
## Version: 
## Last-Updated: sep  7 2018 (17:31) 
##           By: Brice Ozenne
##     Update #: 66
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * individualScore2dt
## Convert output of .BuyseTest (list of vector) into a list of data.table
individualScore2dt <- function(individualScore,
                               correction.tte,
                               level.treatment,
                               level.strata,
                               n.strata,
                               endpoint,
                               threshold,
                               indexT,
                               indexC){
    
    ## Rcpp outputs vector: convert to matrix and rename
    name.tempo <- c("strata",
                    "index.T", "index.C", 
                    "indexWithinStrata.T", "indexWithinStrata.C", 
                    "favorable","unfavorable","neutral","uninformative",
                    "weight",
                    "favorable.corrected","unfavorable.corrected","neutral.corrected")

    individualScore2 <- lapply(individualScore, function(iC){ ## iC <- individualScore[[1]]
        iM <- data.table::as.data.table(matrix(iC, ncol = 13, byrow = FALSE,
                                               dimnames = list(NULL,name.tempo)))
        iM[, c("strata") := factor(.SD[["strata"]], levels = 0:(n.strata-1), labels = level.strata)] ## indexes start at 1 in R and not at 0 as in C++
        ## recall that indexes start at 1 in R and not at 0 as in C++
        iM[, c("index.C") := indexC[.SD$index.C + 1]] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("index.T") := indexT[.SD$index.T + 1]] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("indexWithinStrata.T") := .SD$indexWithinStrata.C + 1]
        iM[, c("indexWithinStrata.C") := .SD$indexWithinStrata.T + 1]
        return(iM[])
    })
    names(individualScore2) <- paste0(endpoint,"_",threshold)


    return(individualScore2)
}
