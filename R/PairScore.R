### tablePairScore.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:54) 
## Version: 
## Last-Updated: mar 23 2020 (10:40) 
##           By: Brice Ozenne
##     Update #: 118
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
                    "index.C", "index.T", "index.pair",
                    "indexWithinStrata.C", "indexWithinStrata.T", 
                    "favorable","unfavorable","neutral","uninf",
                    "weight",
                    "favorableC","unfavorableC","neutralC","uninfC")
    pairScore2 <- lapply(pairScore, function(iC){ ## iC <- pairScore[[1]]
        iM <- data.table::as.data.table(matrix(iC, ncol = 15, byrow = FALSE,
                                               dimnames = list(NULL,name.tempo)))
        iM[, c("strata") := factor(.SD[["strata"]], levels = 0:(n.strata-1), labels = level.strata)] ## indexes start at 1 in R and not at 0 as in C++
        ## recall that indexes start at 1 in R and not at 0 as in C++
        iM[, c("index.C") := .SD$index.C + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("index.T") := .SD$index.T + 1] ## restaure position in the original dataset, not the datasets relative to T and C
        iM[, c("index.pair") := .SD$index.pair + 1] 
        iM[, c("indexWithinStrata.T") := .SD$indexWithinStrata.T + 1]
        iM[, c("indexWithinStrata.C") := .SD$indexWithinStrata.C + 1]
        data.table::setkeyv(iM, c("index.T","index.C"))
        return(iM)
    })
    names(pairScore2) <- paste0(endpoint,"_",threshold)

    return(pairScore2)
}

## * wsumPairScore
## cumulate over endpoint the scores
wsumPairScore <- function(pairScore, weight, n.endpoint){

    keep.col <- c("strata","index.C","index.T","index.pair","indexWithinStrata.C", "indexWithinStrata.T","favorableC","unfavorableC")
    old.col <- c("favorableC","unfavorableC")
    new.col <- c("favorable","unfavorable")

    out <- vector(mode = "list", length = n.endpoint)
    for(iE in 1:n.endpoint){ ## iE <- 2

        iTable <- data.table::copy(pairScore[[iE]][,.SD,.SDcols = keep.col])
        setnames(iTable, old = old.col, new = new.col)
        iTable[,c("favorable") := .SD$favorable * weight[iE]]
        iTable[,c("unfavorable") := .SD$unfavorable * weight[iE]]
        
        if(iE==1){
            out[[iE]] <- iTable
        }else{
            out[[iE]] <- data.table::copy(out[[iE-1]])
            out[[iE]][iTable$index.pair, c("favorable") := .SD$favorable + iTable$favorable]
            out[[iE]][iTable$index.pair, c("unfavorable") := .SD$unfavorable + iTable$unfavorable]
        }
    }

    return(out)
}
