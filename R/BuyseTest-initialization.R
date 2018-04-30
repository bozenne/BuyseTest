## * Documentation initialization functions called by BuyseTest

#' @title internal functions for BuyseTest - initialization
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' @name internal-initialization
#' 
#' @details
#' 
#' \code{initializeArgs}: Normalize the argument 
#' \itemize{
#' \item method, neutral.as.uninf, keep.comparison, n.resampling, seed, cpus, trace: set to default value when not specified.
#' \item formula: call \code{initializeFormula} to extract arguments.
#' \item type: convert to numeric.
#' \item censoring: only keep censoring relative to TTE endpoint. Set to \code{NULL} if no TTE endpoint.
#' \item threshold: set default threshold to 1e-12 expect for binary variable where it is set to 1/2.
#' the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12.
#' When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' \item data: convert to data.table object.
#' \item method: convert to numeric.#' 
#' }
#' and create \code{Wscheme}. \cr \cr
#'
#' \code{initializeFormula}:  extract \code{treatment}, \code{type}, \code{endpoint}, \code{threshold}, \code{censoring}, \code{operator}, and \code{strata}
#' from the formula. \cr \cr
#'
#' \code{initializeData}: Divide the dataset into two, one relative to the treatment group and the other relative to the control group.
#' Merge the strata into one with the interaction variable.
#' Extract for each strata the index of the observations within each group.
#' Apply Efron correction (when method is set to Efron), i.e. set the last event to an observed event.
#' 
#' \code{initializeSurvival}: Compute the survival via KM.
#' 
#' @keywords function internal BuyseTest

## * initializeArgs
#' @rdname internal-initialization
initializeArgs <- function(alternative,
                           name.call,
                           censoring,
                           correctionTTE,
                           cpus,
                           data,
                           endpoint,
                           formula,
                           keep.comparison,
                           method,
                           method.inference,
                           n.resampling,
                           neutral.as.uninf,
                           operator,
                           option,
                           prob.alloc,
                           seed,
                           strata,
                           threshold,
                           trace,
                           treatment,
                           type){

    ## ** apply default options
    if(is.null(cpus)){ cpus <- option$cpus }
    if(is.null(keep.comparison)){ keep.comparison <- option$keep.comparison }
    if(is.null(method)){ method <- option$method }
    if(is.null(method.inference)){ method.inference <- option$method.inference }
    if(is.null(n.resampling)){ n.resampling <- option$n.resampling }
    if(is.null(neutral.as.uninf)){ neutral.as.uninf <- option$neutral.as.uninf }
    if(is.null(seed)){ seed <- option$seed }
    if(is.null(trace)){ trace <- option$trace }

    ## ** convert formula into separate arguments
    if(!missing(formula)){
        resFormula <- initializeFormula(formula)
        treatment <- resFormula$treatment
        type <- resFormula$type
        endpoint <- resFormula$endpoint
        threshold <- resFormula$threshold
        censoring <- resFormula$censoring
        operator <- resFormula$operator
        strata <- resFormula$strata
    }else{
        if(is.null(operator)){
            operator <- rep(">0",length(endpoint))
        }
        formula <- NULL
    }

    ## ** endpoint
    D <- length(endpoint) 
    
    ## ** type
    if(!is.numeric(type)){
        validType1 <- c("b","bin","binary")
        validType2 <- c("c","cont","continuous")
        validType3 <- c("t","tte","time","timetoevent")
        type <- tolower(type)

        type[grep(paste(validType1,collapse="|"), type)] <- "1" 
        type[grep(paste(validType2,collapse="|"), type)] <- "2" 
        type[grep(paste(validType3,collapse="|"), type)] <- "3" 
        type <- as.numeric(type) # type is an integer equal to 1 (binary endpoint), 2 (continuous endpoint) or 3 (time to event endpoint)
    }
    
    D.TTE <- sum(type == 3) # number of time to event endpoints

    ## ** method.inference
    method.inference <- tolower(method.inference)
    
    ## ** censoring
    if(D.TTE==0){
        if(!is.null(censoring) && all(is.na(censoring))){
            censoring <- NULL
        }
    }else if(length(censoring) == D){
        censoring <- censoring[type==3] # from now, censoring contains for each time to event endpoint the name of variable indicating censoring (0) or event (1)
    }

    ## ** threshold
    if(is.null(threshold)){
        threshold <- rep(10^{-12},D)  # if no treshold is proposed all threshold are by default set to 10^{-12}
        if(any(type==1)){threshold[type==1] <- 1/2} # except for threshold corresponding to binary endpoints that are set to NA.
    }else{
        if(any(is.na(threshold[type==1]))){
            index.tempo <- intersect(which(is.na(threshold)),which(type==1))
            threshold[index.tempo] <- 1/2
        }
        if(any(is.na(threshold[type!=1]))){
            index.tempo <- intersect(which(is.na(threshold)),which(type!=1))
            threshold[index.tempo] <- 10^{-12}
        }
        if(any(abs(stats::na.omit(threshold))<10^{-12})){
            threshold[which(abs(threshold)<10^{-12})] <- 10^{-12}
        }
    }

    ## ** data
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }

    ## ** method
    if(!is.numeric(method)){
        method <- switch(method,
                               "Gehan" = 0,
                               "Peto" = 1,
                               "Efron" = 2,
                               "Peron" = 3
                               )
    }

    if (D.TTE == 0) {
        method <- 0
        if ("method" %in% name.call && trace > 0) {
            message("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
        }
    }

    ## ** alternative
    alternative <- tolower(alternative)
    
    ## ** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }

    ## ** Wscheme
    if(D.TTE>0 && D>1){

        endpoint.TTE <- endpoint[type==3]
      
        Wscheme <- matrix(NA,nrow=D.TTE,ncol=D-1) # design matrix indicating to combine the weights obtained at differents TTE endpoints
        rownames(Wscheme) <- paste("weigth of ",endpoint.TTE,"(",threshold[type==3],")",sep="")
        colnames(Wscheme) <- paste("for ",endpoint[-1],"(",threshold[-1],")",sep="")
    
        index.survivalM1 <- rep(-1,D.TTE) # index of previous TTE endpoint (-1 if no previous TTE endpoint i.e. endpoint has not already been used)
        threshold.TTEM1 <- rep(-1,D.TTE) # previous threshold (-1 if no previous threshold i.e. endpoint has not already been used)
    
        iEndpoint.TTE <- if(type[1]==3){1}else{0} #  index.survivalM1 and threshold.TTEM1 are -1 even if the first endoint is a survival endpoint
    
        for(iEndpoint in 2:D){   
      
            if(type[iEndpoint]==3){ 
                iEndpoint.TTE <- iEndpoint.TTE + 1
            }
      
            ## select valid rows
            index.rowTTE <- which(paste(endpoint.TTE,threshold[type==3],sep="_") %in% paste(endpoint,threshold,sep="_")[1:(iEndpoint-1)])
            if(length(index.rowTTE)==0){next} # not yet TTE endpoints (no valid rows)
            Wscheme[index.rowTTE,iEndpoint-1] <- 0 # potential weights
      
            ## keep only the last repeated endpoints
            index.rowTTE <- sapply(unique(endpoint.TTE[index.rowTTE]),
                                   function(x){utils::tail(index.rowTTE[endpoint.TTE[index.rowTTE]==x],1)})
      
            ## if survival endpoint remove similar endpoints
            if(endpoint[iEndpoint] %in% endpoint.TTE[index.rowTTE]){
                index.survivalM1[iEndpoint.TTE] <- index.rowTTE[endpoint.TTE[index.rowTTE]==endpoint[iEndpoint]]-1
                threshold.TTEM1[iEndpoint.TTE] <- threshold[type==3][index.survivalM1[iEndpoint.TTE]+1]
                index.rowTTE <- setdiff(index.rowTTE,index.survivalM1[iEndpoint.TTE]+1)
            }
      
            ## update Wscheme
            if(length(index.rowTTE)>0){
                Wscheme[index.rowTTE,iEndpoint-1] <- 1
            }
      
        }

    }else{
        Wscheme <- matrix(nrow=0,ncol=0) #  factice matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
        index.survivalM1 <- numeric(0) # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
        threshold.TTEM1 <- numeric(0) # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }

    ## ** additional elements for print and object
    level.treatment <- levels(as.factor(data[[treatment]]))
    if(!is.null(strata)){
        level.strata <- levels(interaction(data[,strata,with=FALSE], drop = TRUE, lex.order = FALSE, sep = "."))
    }else{
        level.strata <- "1"
    }
    
    ## ** export
    return(list(
        alternative = alternative,
        name.call = name.call,
        censoring = censoring,
        correctionTTE = correctionTTE,
        cpus = cpus,
        data = data,
        endpoint = endpoint,
        formula = formula,
        index.survivalM1 = index.survivalM1,
        keep.comparison = keep.comparison,
        level.strata = level.strata,
        level.treatment = level.treatment,
        method = method,
        method.inference = method.inference,
        n.resampling = n.resampling,
        neutral.as.uninf = neutral.as.uninf,
        operator = operator,
        seed = seed,
        strata = strata,
        threshold = threshold,
        threshold.TTEM1 = threshold.TTEM1,
        trace = trace,
        treatment = treatment,
        type = type,
        Wscheme = Wscheme
    ))
}

## * initializeFormula
#' @rdname internal-initialization
initializeFormula <- function(x){

  validClass(x, valid.class = "formula")
    
  ## ** extract treatment
  treatment <- setdiff(all.vars(x), all.vars(stats::delete.response(stats::terms(x))))
  if(length(treatment)!=1){
    stop("initFormula: there must be exactly one response variable in formula\n",
         "number of response variables founded: ",length(treatment),"\n")
  }
  
  if(length(as.character(x))!=3){
    stop("initFormula: formula with unexpected length, as.character(x) should have length 3\n",
         "length founded: ",length(as.character(x)),"\n")
  }
  
    ## ** restrict to the right side of the formula
    x.rhs <- as.character(x)[3]
  
    ## remove all blancks
    x.rhs <- gsub("[[:blank:]]", "", x.rhs)
    vec.x.rhs <- unlist(strsplit(x.rhs, split = "+", fixed = TRUE))

    ## find all element in the vector corresponding to endpoints (i.e. ...(...) )
    ## \\w* any letter/number
    ## [[:print:]]* any letter/number/punctuation/space
    index.endpoint <- grep("^\\w*\\([[:print:]]*\\)$", vec.x.rhs)
    index.strata <- setdiff(1:length(vec.x.rhs), index.endpoint)

    ## ** strata variables
    if(length(index.strata)==0){
        strata <- NULL
    }else{
        strata <- vec.x.rhs[index.strata]
    }
    
    ## ** number of endpoint variables    
    vec.x.endpoint <- vec.x.rhs[index.endpoint]
    n.endpoint <- length(vec.x.endpoint)
    if(length(n.endpoint)==0){
        stop("initFormula: x must contains endpoints \n",
             "nothing of the form type(endpoint,threshold,censoring) found in the formula \n")
    }

    ## ** extract endpoints and additional arguments 
    threshold <- rep(NA, n.endpoint)
    censoring <- rep(NA, n.endpoint)
    endpoint <- rep(NA, n.endpoint)
    operator <- rep(">0", n.endpoint)
    validArgs <- c("endpoint","threshold","censoring","operator")

    ## split around parentheses
    ls.x.endpoint <- strsplit(vec.x.endpoint, split = "(", fixed = TRUE)

    type <- character(length = n.endpoint)
    for(iE in 1:n.endpoint){
        ## extract type
        type[iE] <- ls.x.endpoint[[iE]][1]
        
        ## get each argument
        iVec.args <- strsplit(gsub(")", replacement = "",ls.x.endpoint[[iE]][2]),
                              split = ",", fixed = TRUE)[[1]]
        n.args <- length(iVec.args)
        
        ## check size
        if(n.args==0){
            stop("initFormula: invalid formula \n",
                 vec.x.rhs[iE]," must contain the name of the endpoint between the parentheses \n"
                 )
        }        
        if(n.args>4){
            stop("initFormula: invalid formula \n",
                 x[iE]," has too many arguments (maximum 4: endpoint, threshold, censoring variable, operator) \n")
        }

        ## extract name of each argument
        iIndex.name <- grep("=",iVec.args)
        iArg <- gsub("^[[:print:]]*=", replacement = "", iVec.args)        
        iName <- rep(as.character(NA),n.args)
        
        ## use existing names
        if(length(iIndex.name)>0){
            iiName <- gsub("=[[:print:]]*$","",iVec.args[iIndex.name])
            iName[iIndex.name] <- iiName
            if(any(iiName %in% validArgs == FALSE)){
                stop("initFormula: invalid formula \n",
                     vec.x.rhs[iE]," contains arguments that are not endpoint, threshold, censoring \n")
            }
            if( any(duplicated(iiName)) ){
                stop("initFormula: invalid formula \n",
                     vec.x.rhs[iE]," contains arguments with the same name \n")
            }
        }else{
            iiName <- NULL
        }
        
        ## add missing names
        n.missingNames <- n.args - length(iiName) 
        if(n.missingNames>0){
            iName[setdiff(1:n.args,iIndex.name)] <- setdiff(validArgs,iiName)[1:n.missingNames]
        }

        ## extract arguments
        endpoint[iE] <- gsub("\"","",iArg[iName=="endpoint"])
        if("threshold" %in% iName){
            threshold[iE] <- as.numeric(iArg[iName=="threshold"])
        }
        if("censoring" %in% iName){
            censoring[iE] <- gsub("\"","",iArg[iName=="censoring"])
        }
        if("operator" %in% iName){
            operator[iE] <- gsub("\"","",iArg[iName=="operator"])
        }
    }

    ## ** export
    return(list(treatment = treatment,
                type = type,
                endpoint = endpoint,
                threshold = threshold,
                censoring = censoring,
                operator = operator,
                strata = strata))
}

## * initializeData
#' @rdname internal-initialization
initializeData <- function(data, treatment, endpoint, censoring, type, operator, method, strata){
    

    ## ** Treatment: extract the 2 levels
    level.treatment <- levels(as.factor(data[[treatment]])) # extraction of the levels of the treatment variable

    ## ** convert character/factor to numeric for binary endpoints
    name.bin <- endpoint[which(type %in% 1)]
    if(length(name.bin)>0){
        data.class <- sapply(data,class)
        
        for(iBin in name.bin){
            if(data.class[iBin] %in% c("numeric","integer") == FALSE){
                data[[iBin]] <- as.numeric(as.factor(data[[iBin]])) - 1
            }
        }
    }    

    ## ** operator
    operator.endpoint <- setNames(operator, endpoint)[!duplicated(endpoint)]
    name.negative <- names(operator.endpoint)[operator.endpoint=="<0"]
    if(length(name.negative)>0){
        name.negative.binary <- intersect(name.negative, endpoint[type==1])
        if(length(name.negative.binary)>0){
            data[, (name.negative.binary) := -.SD+1, .SDcols = name.negative.binary]
        }
        
        name.negative.other <- setdiff(name.negative, name.negative.binary)
        if(length(name.negative.other)){
            data[, (name.negative.other) := -.SD , .SDcols = name.negative.other]
        }
    }

    ## ** data: split the data according to the two levels
    indexT <- which(data[[treatment]] == level.treatment[2])
    indexC <- which(data[[treatment]] == level.treatment[1])
    dataT <- data[indexT, c(endpoint, strata, censoring), with = FALSE]
    dataC <- data[indexC, c(endpoint, strata, censoring), with = FALSE]

    n.Treatment <- NROW(dataT) # number of patient in the treatment arm
    n.Control <- NROW(dataC) # number of patient in the control arm
        
    ## ** strata
    if(!is.null(strata)){
        strataT <- interaction(dataT[,strata,with=FALSE], drop = TRUE, lex.order = FALSE, sep=".") # strata variable for the treatment arm transformed into factor
        strataC <- interaction(dataC[,strata,with=FALSE], drop = TRUE, lex.order = FALSE, sep=".") # strata variable for the control arm transformed into factor
        level.strata <- levels(strataT) # extraction of the levels of the strata variables

        strataT <- as.numeric(strataT) # strata variable for the treatment arm converted into numeric 
        strataC <- as.numeric(strataC) # strata variable for the control arm converted into numeric
    
    }else{ # if there is no strata variable the same strata is used for all patient
        strataT <- rep(1,n.Treatment) # all patient of the treatment arm are in the same strata : strata 1
        strataC <- rep(1,n.Control)  # all patient of the control arm are in the same strata : strata 1
        level.strata <- 1 # the only strata is strata 1
    }
    n.strata <- length(level.strata) # number of strata
    index.strataT <- lapply(1:n.strata,function(x){which(x==strataT)-1}) # for each strata, the index of the patients belonging to this strata in the treatment arm is stored in an element of the list. 
    index.strataC <- lapply(1:n.strata,function(x){which(x==strataC)-1}) # for each strata, the index of the patients belonging to this strata in the control arm is stored in an element of the list. 
                                        # Because of the minus one after the which, both indexes begin at 0. This is compulsory for C++ compatibility
    
    ## ** data: extract endpoint 
    M.Treatment <- as.matrix(dataT[,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
    M.Control <- as.matrix(dataC[,endpoint,with=FALSE]) # matrix of endpoints for the control arm

    ## ** data: censoring
    if(!is.null(censoring)){
        M.delta.Treatment <- as.matrix(dataT[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
        M.delta.Control <- as.matrix(dataC[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    }else{ # if the is no time to event variables
        M.delta.Treatment <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
        M.delta.Control <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }

    ## ** efron correction
    if(method==2){ # "Efron"
            for(iStrata in 1:n.strata){ ## iStrata <- 1

                index.typeTTE <- which(type==3)
                D.TTE <- length(index.typeTTE)
                Mstrata.Treatment <- M.Treatment[index.strataT[[iStrata]]+1,,drop=FALSE]
                Mstrata.Control <- M.Control[index.strataC[[iStrata]]+1,,drop=FALSE]
        
                ## set last observation for each TTE endpoint to non-censored
                for(iEndpoint.TTE in 1:D.TTE){ ## iEndpoint.TTE <- 1
                    iEndpoint <- index.typeTTE[iEndpoint.TTE]
                        
                    indexT_maxCensored <- which(Mstrata.Treatment[,iEndpoint] == max(Mstrata.Treatment[,iEndpoint])) ## cannot use which.max - not handlle correctly multiple times
                    M.delta.Treatment[index.strataT[[iStrata]][indexT_maxCensored]+1,iEndpoint.TTE] <- 1
                    indexC_maxCensored <- which(Mstrata.Control[,iEndpoint] == max(Mstrata.Control[,iEndpoint]))
                    M.delta.Control[index.strataC[[iStrata]][indexC_maxCensored]+1,iEndpoint.TTE] <- 1
                }
            }
    }
    
    ## ** export
    return(list(index.strataT = index.strataT,
                index.strataC = index.strataC,
                level.strata = level.strata,
                n.strata = n.strata,
                M.Treatment = M.Treatment,
                M.Control = M.Control,
                M.delta.Treatment = M.delta.Treatment,
                M.delta.Control = M.delta.Control))
}










## * initializeSurvival
## ** initializeSurvival_Peto
#' @rdname internal-initialization
initializeSurvival_Peto <- function(M.Treatment,
                                    M.Control,
                                    M.delta.Treatment,
                                    M.delta.Control,
                                    D.TTE,
                                    endpoint,
                                    type,
                                    threshold,
                                    index.strataT,
                                    index.strataC,
                                    n.strata){

    D.TTE <- sum(type==3)
    
    ## ** conversion to R index
    index.strataT <- lapply(index.strataT,function(x){x+1})
    index.strataC <- lapply(index.strataC,function(x){x+1})
    
    ## ** preparing 
    colnamesT <- paste("Survival_TimeT",c("-threshold","_0","+threshold"),sep="")
    colnamesC <- paste("Survival_TimeC",c("-threshold","_0","+threshold"),sep="")

    ## list of matrix of survival data for each endpoint in the treatment arm
    list.survivalT <- rep(list(matrix(NA,
                                      nrow=nrow(M.Treatment),
                                      ncol=3,
                                      dimnames=list(1:nrow(M.Treatment),colnamesT))),
                          times=D.TTE) 
    ## list of matrix of survival data for each endpoint in the control arm
    list.survivalC <- rep(list(matrix(NA,
                                      nrow=nrow(M.Control),
                                      ncol=3,
                                      dimnames=list(1:nrow(M.Control),colnamesC))),
                          times=D.TTE) 

    
    ## ** Compute survival in each strata 
    for(iStrata in 1:n.strata){

        Mstrata.Treatment <- M.Treatment[index.strataT[[iStrata]],,drop=FALSE]
        Mstrata.Control <- M.Control[index.strataC[[iStrata]],,drop=FALSE]
        Mstrata.delta.Treatment <- M.delta.Treatment[index.strataT[[iStrata]],,drop=FALSE]
        Mstrata.delta.Control <- M.delta.Control[index.strataC[[iStrata]],,drop=FALSE]
        
        for(iEndpoint.TTE in 1:D.TTE){
            
            time_treatment <- Mstrata.Treatment[,endpoint[type==3][iEndpoint.TTE]] # store the event times of the treatment arm for a given TTE endpoint
            time_control <- Mstrata.Control[,endpoint[type==3][iEndpoint.TTE]] # store the event times of the control arm for a given TTE endpoint
                
            ## prefer c(time_treatment,time_control) to time_all to avoid rounding error
            df.all <- data.frame(time = c(time_treatment,time_control),
                                 status = c(Mstrata.delta.Treatment[,iEndpoint.TTE], Mstrata.delta.Control[,iEndpoint.TTE]))

            ## *** estimate survival
            resKM_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.all)

            ## *** define interpolation function
            if(all(df.all$status[which(df.all$time==max(df.all$time))]==1)){ # step interpolation of the survival function
                survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                                    yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
                    
            }else{
                survestKM_tempo <- stats::approxfun(x=unique(sort(c(time_treatment,time_control))),y=resKM_tempo$surv,
                                                    yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is a death
            }

            ## *** apply interpolation function to event times (treatment group)
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],1] <- survestKM_tempo(time_treatment-threshold[type==3][iEndpoint.TTE]) # survival at t - tau, tau being the threshold
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],2] <- survestKM_tempo(time_treatment) # survival at t
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],3] <- survestKM_tempo(time_treatment+threshold[type==3][iEndpoint.TTE]) # survival at t + tau
                
            rownames(list.survivalT[[iEndpoint.TTE]])[index.strataT[[iStrata]]] <- time_treatment
                
            ## *** apply interpolation function to event times (control group)
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],1] <- survestKM_tempo(time_control-threshold[type==3][iEndpoint.TTE]) # survival at t - tau
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],2] <- survestKM_tempo(time_control) # survival at t
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],3] <- survestKM_tempo(time_control+threshold[type==3][iEndpoint.TTE]) # survival at t + tau
                
            rownames(list.survivalC[[iEndpoint.TTE]])[index.strataC[[iStrata]]] <- time_control    
        }        
        
    }

    ## ** export
    return(list(list.survivalT=list.survivalT,
                list.survivalC=list.survivalC))
    
}
## ** initializeSurvival_Peron
#' @rdname internal-initialization
initializeSurvival_Peron <- function(M.Treatment,
                                     M.Control,
                                     M.delta.Treatment,
                                     M.delta.Control,
                                     D.TTE,
                                     endpoint,
                                     type,
                                     threshold,
                                     index.strataT,
                                     index.strataC,
                                     n.strata){
    
    ## ** conversion to R index
    index.strataT <- lapply(index.strataT,function(x){x+1})
    index.strataC <- lapply(index.strataC,function(x){x+1})
    
    ## ** preparing 
    colnamesT <- c(paste("SurvivalT_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalC_TimeT",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalT_TimeT_0(ordered)","SurvivalC_TimeT-threshold(ordered)","SurvivalC_TimeT+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
        
    colnamesC <- c(paste("SurvivalC_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   paste("SurvivalT_TimeC",c("-threshold","_0","+threshold"),sep=""),
                   "SurvivalC_TimeC_0(ordered)","SurvivalT_TimeC-threshold(ordered)","SurvivalT_TimeC+threshold(ordered)",
                   "time_control(ordered)","events(ordered)") 
        
    list.survivalT <- rep(list(matrix(NA,
                                      nrow=nrow(M.Treatment),
                                      ncol=11,
                                      dimnames=list(1:nrow(M.Treatment),colnamesT))),
                          times=D.TTE) # list of matrix of survival data for each endpoint in the treatment arm
    list.survivalC <- rep(list(matrix(NA,
                                      nrow=nrow(M.Control),
                                      ncol=11,
                                      dimnames=list(1:nrow(M.Control),colnamesC))),
                          times=D.TTE) # list of matrix of survival data for each endpoint in the control arm
    
    ## ** Compute survival in each strata 
    for(iStrata in 1:n.strata){

        Mstrata.Treatment <- M.Treatment[index.strataT[[iStrata]],,drop=FALSE]
        Mstrata.Control <- M.Control[index.strataC[[iStrata]],,drop=FALSE]
        Mstrata.delta.Treatment <- M.delta.Treatment[index.strataT[[iStrata]],,drop=FALSE]
        Mstrata.delta.Control <- M.delta.Control[index.strataC[[iStrata]],,drop=FALSE]
        
        for(iEndpoint.TTE in 1:D.TTE){

            ## *** estimate survival
            ## treatment group
            df.treatment <- data.frame(time = Mstrata.Treatment[,endpoint[type==3][iEndpoint.TTE]], # store the event times of the treatment arm for a given TTE endpoint
                                       status = Mstrata.delta.Treatment[,iEndpoint.TTE])
                
            resKMT_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.treatment)
                                        # prefer sort(time_treatment) to resKMT_tempo$time to avoid rounding error

            ## control group
            df.control <- data.frame(time = Mstrata.Control[,endpoint[type==3][iEndpoint.TTE]], # store the event times of the treatment arm for a given TTE endpoint
                                     status = Mstrata.delta.Control[,iEndpoint.TTE])
                
            resKMC_tempo <- prodlim::prodlim(prodlim::Hist(time,status)~1, data = df.control) 

            ## *** define interpolation function (treatment group)
            ## treatment group
            if(all(Mstrata.delta.Treatment[which(df.treatment$time==max(df.treatment$time)),iEndpoint.TTE]==1)){ # step interpolation of the survival function (last = event)
                    
                    survestKMT_tempo <- stats::approxfun(x=unique(sort(df.treatment$time)),y=resKMT_tempo$surv,
                                                         yleft=1,yright=0,f=0, method = "constant") # 0 at infinity if last event is a death      
                    
            }else{ # step interpolation of the survival function (last = censoring)
                    
                survestKMT_tempo <- stats::approxfun(x=unique(sort(df.treatment$time)),y=resKMT_tempo$surv,
                                                     yleft=1,yright=NA,f=0, method = "constant") # NA at infinity if last event is censored
                    
            }

            ## control group
            if(all(Mstrata.delta.Control[which(df.control$time==max(df.control$time)),iEndpoint.TTE]==1)){  # step interpolation of the survival function (last = event)
                    
                survestKMC_tempo <- stats::approxfun(x=unique(sort(df.control$time)), y=resKMC_tempo$surv,
                                                     yleft=1,yright=0,f=0, method= "constant") # 0 at infinity if last event is a death
                    
            }else{ # step interpolation of the survival function (last = censoring)
                    
                survestKMC_tempo <- stats::approxfun(x=unique(sort(df.control$time)),y=resKMC_tempo$surv,
                                                     yleft=1,yright=NA,f=0, method= "constant") # NA at infinity if last event is censored      
                    
            }

            ## *** apply interpolation function to event times
            ## treatment group
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],1] <- survestKMT_tempo(df.treatment$time-threshold[type==3][iEndpoint.TTE]) # survival at t - tau, tau being the threshold
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],2] <- survestKMT_tempo(df.treatment$time) # survival at t
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],3] <- survestKMT_tempo(df.treatment$time+threshold[type==3][iEndpoint.TTE]) # survival at t + tau
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],4] <- survestKMC_tempo(df.treatment$time-threshold[type==3][iEndpoint.TTE]) # survival at t - tau
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],5] <- survestKMC_tempo(df.treatment$time) # survival at t
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],6] <- survestKMC_tempo(df.treatment$time+threshold[type==3][iEndpoint.TTE]) # survival at t + tau
                
            order_time_treatment <- order(df.treatment$time)
                
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],7:9] <- list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],c(2,4,6),drop=FALSE][order_time_treatment,,drop=FALSE] # St[t_treatment], Sc[t_treatment-threshold,t_treatment+threshold]
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],10] <-  df.treatment$time[order_time_treatment]
            list.survivalT[[iEndpoint.TTE]][index.strataT[[iStrata]],11] <-   Mstrata.delta.Treatment[order_time_treatment,iEndpoint.TTE] # to compute the integral (Efron)
                
            rownames(list.survivalT[[iEndpoint.TTE]])[index.strataT[[iStrata]]] <- df.treatment$time

            ## control group
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],1] <- survestKMC_tempo(df.control$time-threshold[type==3][iEndpoint.TTE]) # survival at t - tau
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],2] <- survestKMC_tempo(df.control$time) # survival at t
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],3] <- survestKMC_tempo(df.control$time+threshold[type==3][iEndpoint.TTE]) # survival at t + tau
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],4] <- survestKMT_tempo(df.control$time-threshold[type==3][iEndpoint.TTE]) # survival at t - tau, tau being the threshold
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],5] <- survestKMT_tempo(df.control$time) # survival at t
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],6] <- survestKMT_tempo(df.control$time+threshold[type==3][iEndpoint.TTE]) # survival at t + tau                                              
            order_time_control <- order(df.control$time)
                
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],7:9] <- list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],c(2,4,6),drop=FALSE][order_time_control,,drop=FALSE] # Sc[t_control], St[t_control-threshold,t_control+threshold]
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],10] <- df.control$time[order_time_control]
            list.survivalC[[iEndpoint.TTE]][index.strataC[[iStrata]],11] <- Mstrata.delta.Control[order_time_control,iEndpoint.TTE] # to compute the integral (Efron)
                
            rownames(list.survivalC[[iEndpoint.TTE]])[index.strataC[[iStrata]]] <- df.control$time    
                            
        }
    }
    
    return(list(list.survivalT=list.survivalT,
                list.survivalC=list.survivalC))
    
}






