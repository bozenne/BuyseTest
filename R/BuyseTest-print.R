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
                         method.tte,
                         neutral.as.uninf,
                         correction.uninf,
                         operator,
                         strata,
                         threshold,
                         trace,
                         treatment,
                         type,
                         Wscheme,
                         ...){


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
    df.endpoint[,1] <- paste0("      ",1:D)
    names(df.endpoint)[1] <- "      priority"
    df.endpoint$endpoint <- endpoint
    df.endpoint$type <- c("binary","continuous","time to event")[type]
    df.endpoint$operator <- c("lower is favorable","higher is favorable")[1 + (operator == ">0")]
    df.endpoint$threshold[type!=1] <- threshold[type!=1]
    df.endpoint$censoring[type==3] <- censoring[type==3]

    df.endpoint[is.na(df.endpoint)] <- ""
    ## add white space
    df.endpoint$endpoint <- paste0(df.endpoint$endpoint," ")
    df.endpoint$type <- paste0(df.endpoint$type," ")
    df.endpoint$operator <- paste0(df.endpoint$operator," ")
    df.endpoint$threshold <- paste0(df.endpoint$threshold," ")

    if(all(type!=3)){
        df.endpoint$censoring <- NULL
        if(all(type!=2)){
            df.endpoint$threshold <- NULL
        }
    }
    
    ## ** Display
    cat("Settings \n")
    cat("   - treatment groups: Control = ",level.treatment[1]," and Treatment = ",level.treatment[2],"\n", sep = "")
    cat("   - ",D," endpoint",if(D>1){"s"},": \n", sep = "")
    print(df.endpoint, row.names = FALSE, quote = FALSE, right = FALSE)
    if(n.strata>1){
        txt.variable <- switch(as.character(length(strata)),
                               "0" = " variable",
                               "variables")        
        cat("   - ", n.strata, " strata with levels: ",paste(level.strata, collapse = " ") , "\n", sep = "")
        cat("                ",txt.variable,": ",paste(strata, collapse = " ")," \n", sep = "")
    }
    cat("   - management of neutral pairs: ")
    if(neutral.as.uninf){
        cat("re-analyzed using endpoints of lower priority (if any) \n")
    }else{
        cat("ignore endpoints of lower priority \n")
    }
    if(D.TTE>0){
        cat("   - management of censored survival pairs: ")
        switch(as.character(method.tte),
               "0" = cat("uninformative pairs \n"),
               "1" = cat("use Kaplan Meier survival curves to compute the score \n")
               )
    }
    ## if(trace>2){
    ##     if ( (method.tte == "3" || correction.uninf) && D > 1) {            
    ##         cat("   - Current contribution of a pair based on the weights computed at previous enpoints: \n")
    ##         print(Wscheme)
    ##     }
    ## }

    return(NULL)
}

## * Function printInference
#' @rdname internal-print
printInference <- function(method.inference, n.resampling, cpus, seed, ...){

    txt.type <- switch(method.inference,
                       "asymptotic" = "moments of the U-statistic",
                       "bootstrap" = paste0("non parametric bootstrap with ",n.resampling," samples"),
                       "stratified bootstrap" = paste0("stratified non-parametric bootstrap with ",n.resampling," samples"),
                       "permutation" = paste0("permutation test with ", n.resampling, " permutations"),
                       "stratified permutation" = paste0("stratified permutation test with ", n.resampling, " permutations"))

    cat("Estimation of the asymptotic distribution \n",
        "   - method: ",txt.type,"\n", sep = "")
    if(method.inference != "asymptotic"){
        cat("   - cpus  : ",cpus,"\n", sep = "")
        if (!is.null(seed)) {
            cat("   - seeds : ",paste(seq(seed,seed + cpus - 1), collapse = " "),"\n", sep = "")       
        }
    }

    return(NULL)
}
