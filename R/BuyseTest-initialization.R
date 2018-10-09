## * Documentation initialization functions called by BuyseTest

#' @title internal functions for BuyseTest - initialization
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' @name internal-initialization
#' 
#' @details
#' 
#' \code{initializeArgs}: Normalize the argument 
#' \itemize{
#' \item method.tte, neutral.as.uninf, keep.pairScore, n.resampling, seed, cpus, trace: set to default value when not specified.
#' \item formula: call \code{initializeFormula} to extract arguments.
#' \item type: convert to numeric.
#' \item censoring: only keep censoring relative to TTE endpoint. Set to \code{NULL} if no TTE endpoint.
#' \item threshold: set default threshold to 1e-12 expect for binary variable where it is set to 1/2.
#' the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12.
#' When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' \item data: convert to data.table object.
#' \item method.tte: convert to numeric.
#' }
#' and create \code{Wscheme}. \cr \cr
#'
#' \code{initializeFormula}:  extract \code{treatment}, \code{type}, \code{endpoint}, \code{threshold}, \code{censoring}, \code{operator}, and \code{strata}
#' from the formula. \cr \cr
#'
#' \code{initializeData}: Divide the dataset into two, one relative to the treatment group and the other relative to the control group.
#' Merge the strata into one with the interaction variable.
#' Extract for each strata the index of the observations within each group.
#' \code{initializeSurvival}: Compute the survival via KM.
#' 
#' @keywords function internal BuyseTest

## * initializeArgs
#' @rdname internal-initialization
initializeArgs <- function(alternative,
                           name.call,
                           censoring,
                           cpus = NULL,
                           data,
                           endpoint,
                           formula,
                           keep.pairScore = NULL,
                           method.tte = NULL,
                           correction.uninf = NULL,
                           model.tte,
                           method.inference = NULL,
                           n.resampling = NULL,
                           neutral.as.uninf = NULL,
                           operator,
                           option,
                           seed = NULL,
                           strata,
                           threshold,
                           trace = NULL,
                           treatment,
                           type){

    ## ** apply default options
    if(is.null(cpus)){ cpus <- option$cpus }
    if(is.null(keep.pairScore)){ keep.pairScore <- option$keep.pairScore }
    if(is.null(method.tte)){ method.tte <- option$method.tte }
    if(is.null(correction.uninf)){ correction.uninf <- option$correction.uninf }
    if(is.null(method.inference)){ method.inference <- option$method.inference }
    if(is.null(n.resampling)){ n.resampling <- option$n.resampling }
    if(is.null(neutral.as.uninf)){ neutral.as.uninf <- option$neutral.as.uninf }
    if(is.null(trace)){ trace <- option$trace }
    if(is.null(alternative)){ alternative <- option$alternative }

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
    if(is.null(strata) && length(grep("stratified ",method.inference))>0){ ## remove stratified if no strata variable
        method.inference <- gsub("stratified ","",method.inference)
    }

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
    
    ## ** method.tte
    ## WARNING: choices must be lower cases
    ##          remember to update check method.tte (in BuyseTest-check.R)
    method.tte <- tolower(method.tte)
    method.tte <- switch(method.tte,
                         "gehan" = 0,
                         "peron" = 1,
                         NA
                         )

    if (D.TTE == 0) {
        method.tte <- 0
        if ("method.tte" %in% name.call && trace > 0) {
            message("NOTE : there is no survival endpoint, \'method.tte\' argument is ignored \n")
        }
    }


    ## ** correction.uninf
    correction.uninf <- as.numeric(correction.uninf)

    ## ## ** model.tte
    if(method.tte > 0){
        if((!is.null(model.tte)) && (D.TTE == 1) && inherits(model.tte, "prodlim")){
            model.tte <- list(model.tte)
            names(model.tte) <- endpoint[type==3]
        }
    }else{
        model.tte <- NULL
    }

    ## ** alternative
    alternative <- tolower(alternative)
    
    ## ** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }

    ## ** export
    return(list(
        alternative = alternative,
        name.call = name.call,
        censoring = censoring,
        correction.uninf = correction.uninf,
        cpus = cpus,
        D = length(endpoint),
        D.TTE = sum(type == 3),
        data = data,
        endpoint = endpoint,
        formula = formula,
        keep.pairScore = keep.pairScore,
        keep.survival = option$keep.survival,
        method.tte = method.tte,
        model.tte = model.tte,
        method.inference = method.inference,
        n.resampling = n.resampling,
        neutral.as.uninf = neutral.as.uninf,
        operator = operator,
        seed = seed,
        strata = strata,
        threshold = threshold,
        trace = trace,
        treatment = treatment,
        type = type
    ))
}

## * initializeData
#' @rdname internal-initialization
initializeData <- function(data, type, endpoint, operator, strata, treatment){

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    n.strata <- length(strata)

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

    ## ** treatment
    level.treatment <- levels(as.factor(data[[treatment]]))
    data[ , c(treatment) := as.numeric(as.character(factor(.SD[[1]], levels = level.treatment, labels = 0:1))), .SDcols = treatment]

    ## ** strata
    if(is.null(strata)){
        level.strata <- "1"
        allstrata <- NULL
    }else {
        ## data[ , c(".allStrata") := paste0(strata[1],"=",.SD[[strata[1]]])]
        ## if(n.strata>1){
            ## for(iStrata in 2:n.strata){ ## add var=value
                ## data[ , c(".allStrata") := paste0(.SD$.allStrata,".",paste0(strata[iStrata],"=",.SD[[strata[iStrata]]]))]
            ## }
        ## }
        ## data[ , c(".allStrata") := as.factor(.SD$.allStrata)]

        data[ , c(".allStrata") := interaction(.SD, drop = TRUE, lex.order = FALSE, sep = "."), .SDcols = strata]
        level.strata <- levels(data[[".allStrata"]])        
        data[ , c(".allStrata") := as.numeric(.SD$.allStrata)] # convert to numeric
        allstrata <- ".allStrata"
    }

    ## ** n.obs
    n.obs <- data[,.N]
    if(!is.null(strata)){
        n.obsStrata <- data[,.N,by = allstrata][,setNames(.SD[[1]],.SD[[2]]),.SD = c("N",allstrata)]
    }else{
        n.obsStrata <- setNames(n.obs,level.strata)
    }
    
    ## ** export
    return(list(data = data,
                level.treatment = level.treatment,
                level.strata = level.strata,
                n.strata = length(level.strata),
                n.obs = n.obs,
                n.obsStrata = n.obsStrata,
                allstrata = allstrata))
}

## * buildWscheme
buildWscheme <- function(endpoint, D.TTE, D,
                         type, threshold){
    
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
    
    return(list(Wscheme = Wscheme,
                index.survivalM1 = index.survivalM1,
                threshold.TTEM1 = threshold.TTEM1)
           )
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
            thresholdTempo <- tryCatch(as.numeric(iArg[iName=="threshold"]),
                                       error = function(c){ "error" },
                                       warning = function(c){ "warning" }
                                       )
            if(thresholdTempo %in% c("error", "warning")){ ## maybe a variable was passed instead of a value
              thresholdTempo <- eval(expr = parse(text = iArg[iName=="threshold"]))
            }
            
            threshold[iE] <- thresholdTempo
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

## * initializeSurvival
#' @rdname internal-initialization
initializeSurvival_Peron <- function(data, dataT, dataC,
                                     model.tte,
                                     n.T, n.C,
                                     treatment,
                                     level.treatment,
                                     endpoint,
                                     censoring,
                                     D.TTE,
                                     type,
                                     strata,
                                     threshold,
                                     index.strataT,
                                     index.strataC,
                                     n.strata){

    threshold.TTE <- threshold[type==3]
    endpoint.TTE <- endpoint[type==3]

    ## ** conversion to R index
    index.strataT <- lapply(index.strataT,function(x){x+1})
    index.strataC <- lapply(index.strataC,function(x){x+1})
    
    ## ** estimate survival
    if(is.null(model.tte)){
        model.tte <- vector(length = D.TTE, mode = "list")
        names(model.tte) <- endpoint.TTE

        txt.modelTTE <- paste0("prodlim::Hist(",endpoint.TTE,",",censoring,") ~ ",treatment)
        if(!is.null(strata)){
            txt.modelTTE <- paste0(txt.modelTTE," + .allStrata")
        }
        
        for(iEndpoint.TTE in 1:D.TTE){ ## iEndpoint.TTE <- 1
            model.tte[[iEndpoint.TTE]] <- prodlim::prodlim(as.formula(txt.modelTTE[iEndpoint.TTE]),
                                                           data = data)
        }
    }else{ ## convert treatment to numeric
        for(iEndpoint.TTE in 1:D.TTE){ ## iEndpoint.TTE <- 1
            model.tte[[iEndpoint.TTE]]$X[[treatment]] <- as.numeric(factor(model.tte[[iEndpoint.TTE]]$X[[treatment]], levels = level.treatment))-1
        }
    }

    ## ** predict individual survival
    ## *** dataset
    colnames.obs <- c("time",
                      paste("SurvivalC",c("-threshold","_0","+threshold"),sep=""),
                      paste("SurvivalT",c("-threshold","_0","+threshold"),sep="")) 

    colnames.jump <- c("time","survival","dSurvival") 

    list.survTimeC <- vector(mode = "list", length = D.TTE) # list of matrix of control in the treatment arm at each observation time
    list.survTimeT <- vector(mode = "list", length = D.TTE) # list of matrix of survival in the treatment arm at each observation time
    list.survJumpC <- vector(mode = "list", length = D.TTE)  # list of matrix of survival in the control arm at each jump time
    list.survJumpT <- vector(mode = "list", length = D.TTE)  # list of matrix of survival in the treatment arm at each jump time
    list.lastSurv <- vector(mode = "list", length = D.TTE)

    ## *** fill
    for(iEndpoint.TTE in 1:D.TTE){ ## iEndpoint.TTE <- 1

        iModelStrata <- data.table(model.tte[[iEndpoint.TTE]]$X,
                                   start = model.tte[[iEndpoint.TTE]]$first.strata,
                                   stop = model.tte[[iEndpoint.TTE]]$first.strata + model.tte[[iEndpoint.TTE]]$size.strata - 1)

        if(is.null(strata)){
            iModelStrata[,".allStrata" := 1]
        }
        setkeyv(iModelStrata, c(treatment,".allStrata"))
        n.strata <- length(unique(iModelStrata$.allStrata))
        
        ## initialization
        list.survTimeC[[iEndpoint.TTE]] <- vector(mode = "list", length = n.strata)
        list.survTimeT[[iEndpoint.TTE]] <- vector(mode = "list", length = n.strata)
        list.survJumpC[[iEndpoint.TTE]] <- vector(mode = "list", length = n.strata)
        list.survJumpT[[iEndpoint.TTE]] <- vector(mode = "list", length = n.strata)

        list.lastSurv[[iEndpoint.TTE]] <- matrix(NA, ncol = 2, nrow = n.strata,
                                                 dimnames = list(NULL,c("Control","Treatment")))
            
        for(iStrata in 1:n.strata){ ## iStrata <- 1

            ## **** identify jump times and model
            for(iGroup in 0:1){ ## iGroup <- 0

                iStart <- iModelStrata[list(iGroup,iStrata),.SD$start]
                iStop <- iModelStrata[list(iGroup,iStrata),.SD$stop]
                iAllTime <- model.tte[[iEndpoint.TTE]]$time[iStart:iStop]
                iAllSurv <- model.tte[[iEndpoint.TTE]]$surv[iStart:iStop]
                i0LastSurv <- (utils::tail(iAllSurv,1)==0)
                iOrder <- order(iAllTime)

                iIndex.jump <- intersect(iStart:iStop,which(model.tte[[iEndpoint.TTE]]$hazard>0))
                if(length(iIndex.jump)==0){ ## i.e. no event
                    iIndex.jump <- iStart
                }
                iJump <- model.tte[[iEndpoint.TTE]]$time[iIndex.jump]
                
                if(iGroup == 0){
                    ## 0 at infinity if last event is a death
                    predSurvC <- stats::approxfun(x = iAllTime[iOrder],
                                                  y = iAllSurv[iOrder],
                                                  yleft=1, yright=switch(as.character(i0LastSurv),
                                                                         "TRUE" = 0,
                                                                         "FALSE" = NA),
                                                  f=0,
                                                  method = "constant")
                    jumpC <- iJump
                    list.lastSurv[[iEndpoint.TTE]][iStrata,"Control"] <- tail(iAllSurv,1)
                }else if(iGroup == 1){
                    ## 0 at infinity if last event is a death
                    predSurvT <- stats::approxfun(x = iAllTime[iOrder],
                                                  y = iAllSurv[iOrder],
                                                  yleft=1, yright=switch(as.character(i0LastSurv),
                                                                         "TRUE" = 0,
                                                                         "FALSE" = NA),
                                                  f=0,
                                                  method = "constant")

                    jumpT <- iJump
                    list.lastSurv[[iEndpoint.TTE]][iStrata,"Treatment"] <- tail(iAllSurv,1)
                }
            }

            ## **** survival at jump times
            list.survJumpC[[iEndpoint.TTE]][[iStrata]] <- cbind(time = jumpC,
                                                                survival = predSurvT(jumpC + threshold.TTE[iEndpoint.TTE]),
                                                                dSurvival = predSurvC(jumpC) - predSurvC(jumpC - 1e-10))
            
            list.survJumpT[[iEndpoint.TTE]][[iStrata]] <- cbind(time = jumpT,
                                                                survival = predSurvC(jumpT + threshold.TTE[iEndpoint.TTE]),
                                                                dSurvival = predSurvT(jumpT) - predSurvT(jumpT-1e-10))
            
            ## **** survival at observation time (+/- threshold)
            list.survTimeC[[iEndpoint.TTE]][[iStrata]] <- matrix(NA,
                                                                 nrow = length(index.strataC[[iStrata]]), ncol = length(colnames.obs),
                                                                 dimnames = list(NULL, colnames.obs))
            list.survTimeT[[iEndpoint.TTE]][[iStrata]] <- matrix(NA,
                                                                 nrow = length(index.strataT[[iStrata]]), ncol = length(colnames.obs),
                                                                 dimnames = list(NULL, colnames.obs))

            iControl.time <- dataC[index.strataC[[iStrata]],.SD[[endpoint.TTE[iEndpoint.TTE]]]]
            iTreatment.time <- dataT[index.strataT[[iStrata]],.SD[[endpoint.TTE[iEndpoint.TTE]]]]

            for(iGroup in 0:1){ ## iGroup <- 0

                if(iGroup==0){
                    iPredSurv <- predSurvC
                    iCol <- c("SurvivalC-threshold","SurvivalC_0","SurvivalC+threshold")
                }else if(iGroup==1){
                    iPredSurv <- predSurvT
                    iCol <- c("SurvivalT-threshold","SurvivalT_0","SurvivalT+threshold")
                }
                list.survTimeC[[iEndpoint.TTE]][[iStrata]][,"time"] <- iControl.time
                list.survTimeC[[iEndpoint.TTE]][[iStrata]][,iCol[1]] <- iPredSurv(iControl.time - threshold.TTE[iEndpoint.TTE]) # survival at t - tau
                list.survTimeC[[iEndpoint.TTE]][[iStrata]][,iCol[2]] <- iPredSurv(iControl.time) # survival at t
                list.survTimeC[[iEndpoint.TTE]][[iStrata]][,iCol[3]] <- iPredSurv(iControl.time + threshold.TTE[iEndpoint.TTE]) # survival at t + tau

                list.survTimeT[[iEndpoint.TTE]][[iStrata]][,"time"] <- iTreatment.time
                list.survTimeT[[iEndpoint.TTE]][[iStrata]][,iCol[1]] <- iPredSurv(iTreatment.time - threshold.TTE[iEndpoint.TTE]) # survival at t - tau
                list.survTimeT[[iEndpoint.TTE]][[iStrata]][,iCol[2]] <- iPredSurv(iTreatment.time) # survival at t
                list.survTimeT[[iEndpoint.TTE]][[iStrata]][,iCol[3]] <- iPredSurv(iTreatment.time + threshold.TTE[iEndpoint.TTE]) # survival at t + tau

            }

        }
    }

    return(list(survTimeC = list.survTimeC,
                survTimeT = list.survTimeT,
                survJumpC = list.survJumpC,
                survJumpT = list.survJumpT,
                lastSurv = list.lastSurv))
    
}






