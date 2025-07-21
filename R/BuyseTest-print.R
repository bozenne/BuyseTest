## * Documentation - print function called by BuyseTest 
#' @name internal-print
#' @title internal functions for BuyseTest - display
#' @description Functions called by \code{\link{BuyseTest}} to display the settings.
#' @noRd
#' 
#' @author Brice Ozenne

## * Function printGeneral
printGeneral <- function(status,
                         D,
                         D.TTE,
                         data,
                         endpoint,
                         hierarchical,
                         level.strata,
                         level.treatment,
                         scoring.rule,
                         M.status,
                         method.score, 
                         neutral.as.uninf,
                         correction.uninf,
                         operator,
                         restriction,
                         strata,
                         threshold,
                         trace,
                         treatment,
                         type,
                         weightEndpoint,
                         Wscheme,
                         ...){

    if(!is.null(strata)){
        n.strata <- length(level.strata)
    }else{
        n.strata <- 1
    }

    ## ** Prepare
    ## endpoint
    name.col <- c("NA", "endpoint","type","operator","restriction","threshold","event")
    df.endpoint <- data.frame(matrix(NA, nrow = D, ncol = 7,
                                     dimnames = list(NULL, name.col)
                                     ), stringsAsFactors = FALSE)
    if(hierarchical){
        df.endpoint[,1] <- paste0("      ",1:D)
        names(df.endpoint)[1] <- "      priority"
    }else{
        df.endpoint[,1] <- paste0("      ",weightEndpoint)
        names(df.endpoint)[1] <- "      weight"
    }
    
    df.endpoint$endpoint <- endpoint
    if(any(type=="gaus")){
    df.endpoint$endpoint[type=="gaus"] <- paste0(df.endpoint$endpoint[type=="gaus"],",",status[type=="gaus"])
    }
    df.endpoint$type <- sapply(type,switch,
                               "bin"="binary",
                               "cont"="continuous",
                               "tte"="time to event",
                               "gaus"="gaussian")
    df.endpoint$operator <- ifelse(operator>0,"higher is favorable","lower is favorable")
    df.endpoint$threshold[type!="bin"] <- threshold[type!="bin"]
    df.endpoint$restriction <- restriction
    df.endpoint$event[type=="tte"] <- status[type=="tte"]
    
    
    ## add white space
    df.endpoint$endpoint <- paste0(df.endpoint$endpoint," ")
    df.endpoint$type <- paste0(df.endpoint$type," ")
    df.endpoint$operator <- paste0(df.endpoint$operator," ")
    if(all(threshold <= 1e-12)){
        df.endpoint$threshold <- NULL
    }else{
        df.endpoint$threshold <- ifelse(df.endpoint$threshold<=1e-12,NA,paste0(df.endpoint$threshold," "))
    }
    if(all(is.na(restriction))){
        df.endpoint$restriction <- NULL
    }else{
        df.endpoint$restriction <- ifelse(is.na(df.endpoint$restriction),NA,paste0(df.endpoint$restriction," "))
    }
    if(all(type!="tte")){
        df.endpoint$event <- NULL        
    }else{
        txt.eventType <- sapply(status[type=="tte"], function(iC){
            return(paste0(" (",paste(sort(unique(M.status[,iC])), collapse = " "),")"))
        })
        df.endpoint$event[type=="tte"] <- paste0(df.endpoint$event[type=="tte"],txt.eventType)
    }
    df.endpoint[is.na(df.endpoint)] <- ""
    
    ## ** Display
    cat("Settings \n")
    cat("   - 2 groups  ",if(D>1){" "},": Control = ",level.treatment[1]," and Treatment = ",level.treatment[2],"\n", sep = "")
    cat("   - ",D," endpoint",if(D>1){"s"},": \n", sep = "")
    print(df.endpoint, row.names = FALSE, quote = FALSE, right = FALSE)
    if(!is.null(strata) && attr(strata,"match")){
        txt.variable <- switch(as.character(length(strata)),
                               "1" = "variable",
                               "variables")        
        cat("   - ", n.strata, " clusters (",txt.variable,": ",paste(strata, collapse = " "),") \n", sep = "")
    }else if(n.strata>1){
        txt.variable <- switch(as.character(length(strata)),
                               "1" = "variable",
                               "variables")        
        cat("   - ", n.strata, " strata   : levels ",paste(level.strata, collapse = " ") , " (",txt.variable,": ",paste(strata, collapse = " "),") \n", sep = "")
    }
    if(D>1){
        cat("   - neutral pairs: ")
        if(all(neutral.as.uninf)){
            cat("re-analyzed using lower priority endpoints \n")
        }else if(all(!neutral.as.uninf)){
            cat("ignored at lower priority endpoints \n")
        }else{
            cat("re-analyzed using lower priority endpoints for endpoint ",
                paste(which(neutral.as.uninf), collapse = ", "),
                " \n                    otherwise ignored at lower priority endpoints \n",sep="")
        }
    }
    if(D.TTE>0){
        cat("   - right-censored pairs: ")
        n.CR <- sum(grep("2", txt.eventType))
        if(n.CR==D.TTE){
            txt.Peron <- "cif"
        }else if(n.CR==0){
            txt.Peron <- "survival"
        }else{
            txt.Peron <- "survival/cif"
        }

        switch(as.character(scoring.rule),
               "0" = cat("deterministic score or uninformative \n"),
               "1" = cat("probabilistic score based on the ",txt.Peron," curves \n",sep=""),
               "2" = cat("probabilistic score based on the ",txt.Peron," curves \n \t\t\t   (set to 0 beyond available follow-up) \n",sep="")
               )
    }
    ## if(trace>2){
    ##     if ( (scoring.rule == "3" || correction.uninf) && D > 1) {            
    ##         cat("   - Current contribution of a pair based on the weights computed at previous enpoints: \n")
    ##         print(Wscheme)
    ##     }
    ## }

    return(NULL)
}

## * Function printInference
printInference <- function(method.inference, strata, n.resampling, cpus, seed, ...){

    if(method.inference != "none"){

        ## method        
        if(attr(method.inference,"ustatistic")){
            if(!is.null(strata) && attr(strata,"match")){
                txt.type <- "variability of the estimate across strata"
            }else{
                txt.type <- "moments of the U-statistic"
            }
        }else if(attr(method.inference,"bootstrap")){
            txt.type <- paste0("non-parametric bootstrap with ",n.resampling," samples")
        }else if(method.inference == "varexact permutation"){
            txt.type <- paste0("permutation test with all possible permutations")
        }else if(attr(method.inference,"permutation")){
            txt.type <- paste0("permutation test with ",n.resampling," permutations")
        }
        if(any(!is.na(attr(method.inference,"resampling-strata")))){
            txt.type <- paste0(txt.type, " (stratified by \"",paste(attr(method.inference,"resampling-strata"),sep="\" \""),"\")")
        }

        ## display
        cat("Estimation of the estimator's distribution \n",
            "   - method: ",txt.type,"\n", sep = "")
        if(!attr(method.inference,"ustatistic")){
            cat("   - cpus  : ",cpus,"\n", sep = "")
            if (!is.null(seed)) {
                cat("   - seeds : ",seed, sep = "")       
            }
        }
    }

    return(NULL)
}
