### tableComparison.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 26 2018 (14:54) 
## Version: 
## Last-Updated: maj 26 2018 (17:09) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * tableComparison2Delta
## Compute the global statistics based on tableComparison - used to check the validity of tableComparison
tableComparison2Delta <- function(table, correct.tte){

     favorable <- unfavorable <- neutral <- uninformative <- NULL ## [:CRAN:] for CRAN check
    
    endpoint <- names(table)
    D <- length(endpoint)

    vec.favorable <- vector(mode = "numeric", length = D)
    vec.unfavorable <- vector(mode = "numeric", length = D)

    Delta.netChance <- vector(mode = "numeric", length = D)
    Delta.winRatio <- vector(mode = "numeric", length = D)

    for(iE in 1:D){ ## iE <- 1

        ## ** count number of pairs
        if(iE==1){
            n.pairs <- NROW(table[[iE]])
        }
        
        ## ** perform correction
        if(correct.tte){
            vec.tempo <- unlist(table[[iE]][,.(favorable = sum(favorable),
                                               unfavorable = sum(unfavorable),
                                               neutral = sum(neutral),
                                               uninformative = sum(uninformative))])
            factor <- sum(vec.tempo)/sum(vec.tempo[1:3])
            
            table[[iE]][, favorable := favorable * factor]
            table[[iE]][, unfavorable := unfavorable * factor]
            table[[iE]][, neutral := neutral * factor]
            table[[iE]][, uninformative := 0]
        }

        ## ** compute Delta
        if(iE==1){
            vec.favorable[iE] <- table[[iE]][,sum(favorable)]
            vec.unfavorable[iE] <- table[[iE]][,sum(unfavorable)]
        }else{
            vec.favorable[iE] <- table[[iE]][,sum(favorable)] + vec.favorable[iE-1]
            vec.unfavorable[iE] <- table[[iE]][,sum(unfavorable)] + vec.unfavorable[iE-1]
        }
        Delta.netChance[iE] <- (vec.favorable[iE]-vec.unfavorable[iE])/n.pairs
        Delta.winRatio[iE] <- vec.favorable[iE]/vec.unfavorable[iE]

    }



    ## ** export
    return(list(Delta.netChance = setNames(Delta.netChance, endpoint),
                Delta.winRatio = setNames(Delta.winRatio, endpoint)))
}

## * tableComparison2dt
## Convert output of .BuyseTest (list of vector) into a list of data.table
tableComparison2dt <- function(tableComparison,
                               level.treatment,
                               level.strata,
                               n.strata,
                               endpoint,
                               threshold,
                               indexT,
                               indexC){
    
    name.indexT <- paste0("index.",level.treatment[2])
    name.indexC <- paste0("index.",level.treatment[1])
    name.indexWT <- paste0("indexWithinStrata.",level.treatment[2])
    name.indexWC <- paste0("indexWithinStrata.",level.treatment[1])
    
    ## Rcpp outputs vector: convert to matrix and rename
    name.tempo <- c("strata",
                    name.indexT, name.indexC, 
                    name.indexWT, name.indexWC, 
                    "favorable","unfavorable","neutral","uninformative")

    tableComparison <- lapply(tableComparison, function(iC){
        iM <- data.table::as.data.table(matrix(iC, ncol = 9, byrow = FALSE,
                                               dimnames = list(NULL,name.tempo)))
        iM[, c("strata") := factor(.SD[["strata"]], levels = 0:(n.strata-1), labels = level.strata)] ## indexes start at 1 in R and not at 0 as in C++
        ## recall that indexes start at 1 in R and not at 0 as in C++
        if(!is.null(indexT)){
            iM[, c(name.indexT) := indexT[.SD[[1]]+1], .SDcols = name.indexT] ## restaure position in the original dataset, not the datasets relative to T and C
        }else{
            iM[, c(name.indexT) := .SD[[1]]+1, .SDcols = name.indexT]
        }
        if(!is.null(indexC)){
            iM[, c(name.indexC) := indexC[.SD[[1]]+1], .SDcols = name.indexC] ## restaure position in the original dataset, not the datasets relative to T and C
        }else{
            iM[, c(name.indexC) := .SD[[1]]+1, .SDcols = name.indexC]
        }
        iM[, c(name.indexWT) := .SD[[1]] + 1, .SDcols = name.indexWT]
        iM[, c(name.indexWC) := .SD[[1]] + 1, .SDcols = name.indexWC]
        ## data.table::setkeyv(iM, cols = c(name.indexT,name.indexC), verbose = FALSE)
        return(iM)
    })
    names(tableComparison) <- paste0(endpoint,"_",threshold)


    return(tableComparison)
}


##----------------------------------------------------------------------
### tableComparison.R ends here
