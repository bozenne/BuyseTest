### BuyseTest-check.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 27 2018 (23:32) 
## Version: 
## Last-Updated: feb 18 2021 (12:16) 
##           By: Brice Ozenne
##     Update #: 267
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * testArgs
##' @title Check Arguments Passed to BuyseTest
##' @description Check the validity of the argument passed the BuyseTest function by the user.
##'
##' @keywords internal
##' @author Brice Ozenne
testArgs <- function(name.call,
                     status,
                     correction.uninf,
                     cpus,
                     data,
                     endpoint,
                     formula,
                     iid,
                     iidNuisance,
                     keep.pairScore,
                     scoring.rule,
                     model.tte,
                     method.inference,
                     n.resampling,
                     strata.resampling,
                     hierarchical,
                     neutral.as.uninf,
                     operator,
                     censoring,
                     seed,
                     strata,
                     threshold,
                     trace,
                     treatment,
                     type,
                     weight,
                     ...){

    ## ** data
    if (!data.table::is.data.table(data)) {
        if(inherits(data,"function")){
            stop("Argument \'data\' is mispecified. \n",
                 "\'data\' cannot be a function. \n")
        }
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    if("..rowIndex.." %in% names(data)){
        stop("BuyseTest: Argument \'data\' must not contain a column \"..rowIndex..\". \n")
    }
    if("..NA.." %in% names(data)){
        stop("BuyseTest: Argument \'data\' must not contain a column \"..NA..\". \n")
    }
    if("..strata.." %in% names(data)){
        stop("BuyseTest: Argument \'data\' must not contain a column \"..strata..\". \n")
    }
    
    ## ** extract usefull quantities
    argnames <- c("treatment", "endpoint", "type", "threshold", "status", "strata")

    D <- length(endpoint) 
    D.TTE <- sum(type == 3) # number of time to event endpoints
    level.treatment <- levels(as.factor(data[[treatment]])) 
    if(is.null(strata)){
        n.strata <- 1
    }else{
        indexT <- which(data[[treatment]] == level.treatment[2])
        indexC <- which(data[[treatment]] == level.treatment[1])

        if(any(strata %in% names(data) == FALSE)){
            stop("Strata variable(s) \"",paste0(strata,collapse="\" \""),"\" not found in argument \'data\' \n")
        }
        
        strataT <- interaction(data[indexT,strata,with=FALSE], drop = TRUE, lex.order=FALSE,sep=".") 
        strataC <- interaction(data[indexC,strata,with=FALSE], drop = TRUE, lex.order=FALSE,sep=".") 
        level.strata <- levels(strataT)
        n.strata <- length(level.strata)
    }

    
    ## ** status
    if(length(status) != D){
        stop("BuyseTest: \'status\' does not match \'endpoint\' size. \n",
             "length(status): ",length(status),"\n",
             "length(endpoint) : ",D,"\n")
            
    }
    if(any(is.na(status))){
        stop("BuyseTest: \'status\' must not contain NA. \n")
    }
    index.pb <- which(status[type==3] == "..NA..")
    if(length(index.pb)>0){
        if(all(attr(censoring,"original")[index.pb] %in% names(data))){
            stop("BuyseTest: wrong specification of \'status\'. \n",
                 "\'status\' must indicate a variable in data for TTE endpoints. \n",
                 "\'censoring\' is used to indicate whether there is left or right censoring. \n",
                 "Consider changing \'censoring =\' into \'status =\' when in the argument \'formula\' \n")
        }else{        
            stop("BuyseTest: wrong specification of \'status\'. \n",
                 "\'status\' must indicate a variable in data for TTE endpoints. \n",
                 "TTE endoints: ",paste(endpoint[type==3],collapse=" "),"\n",
                 "proposed \'status\' for these endoints: ",paste(status[type==3],collapse=" "),"\n")
        }
    }
    if(any(status[type!=3] !="..NA..") ){
        stop("BuyseTest: wrong specification of \'status\'. \n",
             "\'status\' must be \"..NA..\" for binary or continuous endpoints. \n",
             "endoints : ",paste(endpoint[type!=3],collapse=" "),"\n",
             "proposed \'status\' for these endoints: ",paste(status[type!=3],collapse=" "),"\n")
    }
    Ustatus.TTE <- unique(status[type==3])

    if(is.null(strata)){
        if(any(sapply(Ustatus.TTE, function(iS){sum(data[[iS]]!=0)})==0)){
            warning("BuyseTest: time to event variables with only censored events \n")
        }   
    }else{
        strata.tempo <- data[[strata]]
        if(is.factor(strata.tempo)){strata.tempo <- droplevels(strata.tempo)} ## otherwise the next tapply statement generates NA when there are empty levels which leads to an error
        
        if(any(sapply(Ustatus.TTE, function(iS){tapply(data[[iS]], strata.tempo, function(iVec){sum(iVec!=0)})})==0)){
            warning("BuyseTest: time to event variables with only censored events in at least one strata \n")
        }
    }
    
    ## ** censoring
    if(any(is.na(censoring))){
        stop("BuyseTest: wrong specification of \'censoring\'. \n",
             "\'censoring\' must be \'as.character(NA)\', \"left\", or \"right\" \n",
             "incorrect \'censoring\' value(s): \"",paste(attr(censoring,"original")[is.na(censoring)], collapse = "\" \""),"\" \n")
    }
    if(any(censoring[type==3]==0)){
        stop("BuyseTest: wrong specification of \'censoring\'. \n",
             "\'censoring\' must be \"left\" or \"right\" for TTE endpoints \n")
    }

    ## ** cpus
    if(cpus>1){
        validInteger(cpus,
                     valid.length = 1,
                     min = 1,
                     max = parallel::detectCores(),
                     method = "BuyseTest")
    }

    ## ** scoring.rule
    ## must be before time to event endpoints
    if(is.na(scoring.rule)){
        stop("BuyseTest: wrong specification of \'scoring.rule\'. \n",
             "valid values: \"Gehan\" \"Gehan corrected\" \"Peron\" \"Peron corrected\". \n")
    }
    if(scoring.rule>0 && any(censoring>1)){
        warning("The Peron's scoring rule does not support left-censored endpoints \n",
                "For those endpoints, the Gehan's scoring rule will be used instead.")
    }

    ## ## ** model.tte
    if(!is.null(model.tte)){
        endpoint.UTTE <- unique(endpoint[type==3])
        D.UTTE <- length(endpoint.UTTE)

        if(!is.list(model.tte) || length(model.tte) != D.UTTE){
            stop("BuyseTest: argument \'model.tte\' must be a list containing ",D.UTTE," elements. \n",
                 "(one for each unique time to event endpoint). \n")
        }

        if(is.null(model.tte) || any(names(model.tte) != endpoint.UTTE)){
            stop("BuyseTest: argument \'model.tte\' must be a named list. \n",
                 "valid sequence of names: \"",paste0(endpoint.UTTE, collapse = "\" \""),"\" \n",
                 "proposed names: \"",paste0(names(model.tte), collapse = "\" \""),"\" \n")
        }

        valid.class <- setdiff(utils::methods(generic.function = "BuyseTTEM"), c("BuyseTTEM.formula"))
        vec.class  <- sapply(model.tte, function(iTTE){any(paste0("BuyseTTEM.",class(model.tte[[1]])) %in% valid.class)})
        if(any(vec.class == FALSE) ){
            stop("BuyseTest: argument \'model.tte\' must be a list of \"",paste0(gsub("BuyseTTEM\\.","",valid.class), collapse = "\", or \""),"\" objects. \n")
        }
        
        vec.predictors  <- sapply(model.tte, function(iTTE){identical(sort(all.vars(stats::update(stats::formula(model.tte[[1]]), "0~."))), sort(c(treatment,strata)))})
        if(any(vec.predictors == FALSE) ){
            stop("BuyseTest: argument \'model.tte\' must be a list of objects with \"",paste0(c(treatment,strata),collapse = "\" \""),"\" as predictors. \n")
        }        
    }
    
    ## ** data (endpoints)

    ## *** binary endpoints
    index.Bin <- which(type==1)
    if(length(index.Bin)>0){
        for(iBin in index.Bin){ ## iterY <- 1
            if(length(unique(na.omit(data[[endpoint[iBin]]])))>2){
                stop("Binary endpoint cannot have more than 2 levels. \n",
                     "endpoint: ",endpoint[iBin],"\n")
            }
            ## if(any(is.na(data[[endpoint[iBin]]]))){                
                ## warning("BuyseTest: endpoint ",endpoint[iBin]," contains NA \n")
            ## }
        }
    }

    ## *** continuous endpoints
    index.Cont <- which(type==2)
    if(length(index.Cont)>0){
        for(iCont in index.Cont){
            validNumeric(data[[endpoint[iCont]]],
                         name1 = endpoint[iCont],
                         valid.length = NULL,
                         refuse.NA =  FALSE,
                         method = "BuyseTest")

            ## if(any(is.na(data[[endpoint[iCont]]]))){                
                ## warning("BuyseTest: endpoint ",endpoint[iCont]," contains NA \n")
            ## }
        }
    }

    ## *** time to event endpoint
    index.TTE <- which(type==3)
    status.TTE <- status[type==3]
    if(length(index.TTE)>0){
        validNames(data,
                   name1 = "data",
                   required.values = status.TTE,
                   valid.length = NULL,
                   refuse.NULL = FALSE,
                   method = "BuyseTest")

        valid.values.status <- 0:2

        for(iTTE in index.TTE){
            validNumeric(data[[endpoint[iTTE]]],
                         name1 = endpoint[iTTE],
                         valid.length = NULL,
                         refuse.NA = TRUE,
                         method = "BuyseTest")
            validNumeric(unique(data[[status.TTE[which(index.TTE == iTTE)]]]),
                         name1 = status.TTE[which(index.TTE == iTTE)],
                         valid.values = valid.values.status,
                         valid.length = NULL,
                         method = "BuyseTest")
        }
    }

    ## ** endpoint
    validNames(data,
               name1 = "data",
               required.values = endpoint,
               valid.length = NULL,
               method = "BuyseTest")

    ## ** formula
    if(!is.null(formula) && any(name.call %in% argnames)){
        txt <- paste(name.call[name.call %in% argnames], collapse = "\' \'")
        warning("BuyseTest: argument",if(length(txt)>1){"s"}," \'",txt,"\' ha",if(length(txt)>1){"ve"}else{"s"}," been ignored. \n",
                "when specified, only argument \'formula\' is used. \n")
    }

    ## ** keep.pairScore
    validLogical(keep.pairScore,
                 valid.length = 1,
                 method = "BuyseTest")
 
    ## ** correction.uninf
    validInteger(correction.uninf,
                 valid.length = 1,
                 min = 0,
                 max = 3,
                 method = "BuyseTest")

    ## ** method.inference
    if(length(method.inference)!=1){
        stop("Argument \'method.inference\' must have length 1. \n")
    }
    if(method.inference != "u-statistic-bebu"){ ## asympototic bebu - hidden value only for debugging
        validCharacter(method.inference,
                       valid.length = 1,
                       valid.values = c("none","u-statistic","permutation", "studentized permutation", "bootstrap", "studentized bootstrap"),
                       method = "BuyseTest")
    }
    if(method.inference != "none" && any(table(data[[treatment]])<2) ){
        warning("P-value/confidence intervals will not be valid with only one observation. \n")
    }
    if(!is.na(attr(method.inference,"resampling-strata")) && any(attr(method.inference,"resampling-strata") %in% names(data) == FALSE)){
        stop("Incorrect value for argument \'strata.resampling\': must correspond to a column in argument \'data\'. \n")
    }
    if(!is.na(attr(method.inference,"resampling-strata")) && attr(method.inference,"permutation") && any(attr(method.inference,"resampling-strata") == treatment)){
        stop("Argument \'strata.resampling\' should not contain the variable used to form the treatment groups when using a permutation test. \n")
    }
    if(iid && correction.uninf > 0){
        warning("The current implementation of the asymptotic distribution is valid when using a correction. \n",
                "Standard errors / confidence intervals / p-values may not be correct. \n",
                "Consider using a resampling approach or checking the control of the type 1 error with powerBuyseTest. \n")
    }
    
    ## ** n.resampling
    if(method.inference %in% c("bootstrap","permutation","stratified bootstrap","stratified permutation")){
        validInteger(n.resampling,
                     valid.length = 1,
                     min = 1,
                     method = "BuyseTest")
    }
    
   ## ** hierarchical
    validLogical(hierarchical,
                 valid.length = 1,
                 method = "BuyseTest")

    ## ** neutral.as.uninf
    validLogical(neutral.as.uninf,
                 valid.length = D,
                 method = "BuyseTest")

    ## ** operator
    if(any(is.na(operator))){
        stop("BuyseTest: Wrong specification of \'operator\'. \n",
             "Should be either \"<0\" (lower is better) or \">0\" (higher is better)")
    }

    ## ** seed
    validInteger(seed,
                 valid.length = 1,
                 refuse.NULL = FALSE,
                 min = 1,
                 method = "BuyseTest")


     ## ** strata
    if (!is.null(strata)) {
        validNames(data,
                   name1 = "data",
                   required.values = strata,
                   valid.length = NULL,
                   method = "BuyseTest")

       
        if(length(level.strata) != length(levels(strataC)) || any(level.strata != levels(strataC))){
            stop("BuyseTest: wrong specification of \'strata\'. \n",
                 "different levels between Control and Treatment \n",
                 "levels(strataT) : ",paste(levels(strataT),collapse=" "),"\n",
                 "levels(strataC) : ",paste(levels(strataC),collapse=" "),"\n")
        }
    }

    ## ** threshold
    ## check numeric and no NA
    validNumeric(threshold,
                 valid.length = D,
                 min = 0,
                 refuse.NA = TRUE,
                 method = "BuyseTest")

    ## check threshold at 1/2 for binary endpoints
    if(any(threshold[type==1]!=1/2)){
        stop("BuyseTest: wrong specification of \'threshold\'. \n",
             "\'threshold\' must be 1/2 for binary endpoints (or equivalently NA) \n",
             "proposed \'threshold\' : ",paste(threshold[type==1],collapse=" "),"\n",
             "binary endpoint(s) : ",paste(endpoint[type==1],collapse=" "),"\n")
    }
    
    ## Check that the thresholds related to the same endoints are strictly decreasing
    ## is.unsorted(rev(2:1))
    ## is.unsorted(rev(1:2))

    vec.test <- tapply(threshold,endpoint, function(x){
        test.unsorted <- is.unsorted(rev(x))
        test.duplicated <- any(duplicated(x))
        return(test.unsorted+test.duplicated)
    })
    
    if(any(vec.test>0)){   
        stop("BuyseTest: wrong specification of \'endpoint\' or \'threshold\'. \n",
             "Endpoints must be used with strictly decreasing threshold when re-used with lower priority. \n",
             "Problematic endpoints: \"",paste0(names(vec.test)[vec.test>0], collapse = "\" \""),"\"\n")        
    }

    ## ** trace
    validInteger(trace,
                 valid.length = 1,
                 min = 0,
                 max = 2,
                 method = "BuyseTest")

    ## ** treatment
    validCharacter(treatment,
                   valid.length = 1,
                   method = "BuyseTest")

    validNames(data,
               name1 = "data",          
               required.values = treatment,
               valid.length = NULL,
               method = "BuyseTest")

    if (length(level.treatment) != 2) {
        stop("BuyseTest: wrong specification of \'treatment\'. \n",
             "The corresponding column in \'data\' must have exactly 2 levels. \n",
             "Proposed levels : ",paste(level.treatment,collapse = " "),"\n")
    }

    if(any(table(data[[treatment]])==0)){
        txt.stop <- names(which(table(data[[treatment]])==0))
        stop("BuyseTest: wrong specification of \'data\'. \n",
             "No observation taking level ",txt.stop," in the treatment variable. \n")
        
    }
    
    ## ** type
    if(any(type %in% 1:3 == FALSE)){
        txt <- type[type %in% 1:3 == FALSE]
        stop("BuyseTest: wrong specification of \'type\' \n",
             "valid values: \"binary\" \"continuous\" \"timetoevent\" \n",
             "incorrect values: \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    
    n.typePerEndpoint <- tapply(type, endpoint, function(x){length(unique(x))})
    if(any(n.typePerEndpoint>1)){
        message <- paste0("several types have been specified for endpoint(s) ",
                          paste0(unique(endpoint)[n.typePerEndpoint>1],collapse = ""),
                          "\n")        
        stop("BuyseTest: wrong specification of \'endpoint\' or \'type\' \n",message)
    }

    ## ** weight
    if(length(weight) != D){
            stop("BuyseTest: argument \'weight\' must have length the number of endpoints \n")
    }
    
    if(hierarchical){
        if(any(weight!=1) || any(is.na(weight))){
            stop("BuyseTest: all the weights must be 1 when using hierarchical GPC \n")
        }
    }
    
    ## ** export
    return(invisible(TRUE))
}



##----------------------------------------------------------------------
### BuyseTest-check.R ends here
