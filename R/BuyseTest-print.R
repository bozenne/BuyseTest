## * Documentation - print function called by BuyseTest 
#' @name internal-print
#' @title internal functions for BuyseTest - display
#' @description Functions called by \code{\link{BuyseTest}} to display the settings.
#' 
#' @keywords internal
#' @author Brice Ozenne

## * Function printGeneral
#' @rdname internal-print
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
                         neutral.as.uninf,
                         correction.uninf,
                         operator,
                         strata,
                         threshold,
                         trace,
                         treatment,
                         type,
                         weight,
                         Wscheme,
                         ...){

    if(!is.null(strata)){
        n.strata <- length(level.strata)
    }else{
        n.strata <- 1
    }

    ## ** Prepare
    ## endpoint
    name.col <- c("NA", "endpoint","type","operator","threshold","event")
    df.endpoint <- data.frame(matrix(NA, nrow = D, ncol = 6,
                                     dimnames = list(NULL, name.col)
                                     ), stringsAsFactors = FALSE)
    if(hierarchical){
        df.endpoint[,1] <- paste0("      ",1:D)
        names(df.endpoint)[1] <- "      priority"
    }else{
        df.endpoint[,1] <- paste0("      ",weight)
        names(df.endpoint)[1] <- "      weight"
    }
    df.endpoint$endpoint <- endpoint
    df.endpoint$type <- c("binary","continuous","time to event")[type]
    df.endpoint$operator <- ifelse(operator>0,"higher is favorable","lower is favorable")
    df.endpoint$threshold[type!=1] <- threshold[type!=1]
    df.endpoint$event[type==3] <- status[type==3]
    
    
    ## add white space
    df.endpoint$endpoint <- paste0(df.endpoint$endpoint," ")
    df.endpoint$type <- paste0(df.endpoint$type," ")
    df.endpoint$operator <- paste0(df.endpoint$operator," ")
    df.endpoint$threshold <- ifelse(is.na(df.endpoint$threshold),NA,paste0(df.endpoint$threshold," "))

    if(all(type!=3)){
        df.endpoint$event <- NULL
        if(all(type!=2)){
            df.endpoint$threshold <- NULL
        }
    }else{
        txt.eventType <- sapply(status[type==3], function(iC){
            return(paste0(" (",paste(sort(unique(M.status[,iC])), collapse = " "),")"))
        })
        df.endpoint$event[type==3] <- paste0(df.endpoint$event[type==3],txt.eventType)
    }
    df.endpoint[is.na(df.endpoint)] <- ""
    
    ## ** Display
    cat("Settings \n")
    cat("   - 2 groups  ",if(D>1){" "},": Control = ",level.treatment[1]," and Treatment = ",level.treatment[2],"\n", sep = "")
    cat("   - ",D," endpoint",if(D>1){"s"},": \n", sep = "")
    print(df.endpoint, row.names = FALSE, quote = FALSE, right = FALSE)

    if(n.strata>1){
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
               "1" = cat("probabilistic score based on the ",txt.Peron," curves \n",sep="")
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
#' @rdname internal-print
printInference <- function(method.inference, n.resampling, cpus, seed, ...){

    if(method.inference != "none"){

        ## method        
        if(attr(method.inference,"ustatistic")){
            txt.type <- "moments of the U-statistic"
        }else if(attr(method.inference,"bootstrap")){
            txt.type <- paste0("non-parametric bootstrap with ",n.resampling," samples")
        }else if(attr(method.inference,"permutation")){
            txt.type <- paste0("permutation test with ",n.resampling," permutations")
        }
        if(!is.na(attr(method.inference,"resampling-strata"))){
            txt.type <- paste0(txt.type, " (stratified by \"",paste(attr(method.inference,"resampling-strata"),sep="\" \""),"\")")
        }

        ## display
        cat("Estimation of the estimator's distribution \n",
            "   - method: ",txt.type,"\n", sep = "")
        if(!attr(method.inference,"ustatistic")){
            cat("   - cpus  : ",cpus,"\n", sep = "")
            if (!is.null(seed)) {
                cat("   - seeds : ",paste(seq(seed,seed + cpus - 1), collapse = " "),"\n", sep = "")       
            }
        }
    }

    return(NULL)
}
