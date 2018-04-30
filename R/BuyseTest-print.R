## * Documentation - print function called by BuyseTest 
#' @name internal-print
#' @title internal functions for BuyseTest - display
#' @description Functions called by \code{\link{BuyseTest}} to display the settings.
#' 
#' @keywords internal

## * Function printGeneral
#' @rdname internal-print
printGeneral <- function(censoring,
                         D,
                         D.TTE,
                         data,
                         endpoint,
                         level.strata,
                         level.treatment,
                         method,
                         neutral.as.uninf,
                         operator,
                         strata,
                         threshold,
                         threshold.TTEM1,
                         treatment,
                         type,
                         Wscheme,
                         ...){

    D <- length(endpoint)
    D.TTE <- sum(type == 3)
    if(!is.null(strata)){
        n.strata <- length(level.strata)
    }else{
        n.strata <- 1
    }

    ## ** Prepare
    ## endpoint
    name.col <- c("priority", "endpoint","type","operator","threshold","censoring")
    df.endpoint <- data.frame(matrix(NA, nrow = D, ncol = 6,
                                     dimnames = list(NULL, name.col)
                                     ))
    df.endpoint[,1] <- paste0("     |",1:D)
    names(df.endpoint)[1] <- "     |priority"
    df.endpoint$endpoint <- endpoint
    df.endpoint$type <- c("binary","continuous","time to event")[type]
    df.endpoint$operator <- c("lower is favorable","higher is favorable")[1 + (operator == ">0")]
    df.endpoint$threshold[type!=1] <- threshold[type!=1]
    df.endpoint$censoring[type==3] <- censoring

    df.endpoint[is.na(df.endpoint)] <- ""
    if(all(type!=3)){
        df.endpoint$censoring <- NULL
        if(all(type!=2)){
            df.endpoint$threshold <- NULL
        }
    }
    df.endpoint <- cbind(df.endpoint, `|` = "|")
    
    ## threshold
    if(any(type==3)){
        if (length(threshold.TTEM1) > 0) {
            threshold.TTEM1.display <- threshold.TTEM1
            threshold.TTEM1.display[threshold.TTEM1.display < 0] <- +Inf
        }else{
            threshold.TTEM1.display <- +Inf
        }
        
        threshold.display <- rbind(sapply(1:D.TTE,
                                          function(x){paste(c("[",round(threshold.TTEM1.display[x],4),
                                                              " ; ",round(threshold[type == 3][x],4),
                                                              "] "), collapse = "")}))
        colnames(threshold.display) <- endpoint[type == 3]      
        rownames(threshold.display) <- "threshold interval"
    }
    
    ## ** Display
    cat("Settings (punctual estimation) \n")
    cat("   > reference: Control = ",level.treatment[1]," and Treatment = ",level.treatment[2],"\n", sep = "")
    cat("   > ",D," endpoint",if(D>1){"s"},": \n", sep = "")
    print(df.endpoint, row.names = FALSE, quote = FALSE, right = FALSE)
    if(n.strata>1){
        cat("   > ", n.strata, " strata with levels ", paste(level.strata,collapse = " "), "\n", sep = "")
    }
    cat("   > management of neutral pairs : ")
    if(neutral.as.uninf){
        cat("re-analyzed using endpoints of lower priority (if any) \n")
    }else{
        cat("ignore endpoints of lower priority \n")
    }
    if(any(type==3)){
        cat("   > management of censored survival pairs : ")
        switch(method,
               "Gehan" = cat("uninformative pairs \n"),
               "Peto" = cat("imputation using one survival curve estimated on all patients \n"),
               "Efron" = cat("imputation using different survival curve for control and treatment patients \n"),
               "Peron" = cat("imputation using different survival curve for control and treatment patients \n")
               ) 
        if (method %in% c("Peto","Efron","Peron")) {
            
            cat("   > weights of the pairs relatively to the enpoints : \n")
            print(Wscheme)
            
            cat("   > intervals thresholds for survival endpoints : \n")    
            print(threshold.display)
        }
    }
}

