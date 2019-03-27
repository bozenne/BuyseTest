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
                           censoring,
                           correction.uninf = NULL,
                           cpus = NULL,
                           data,
                           endpoint,
                           formula,
                           hierarchical = NULL,
                           keep.pairScore = NULL,
                           method.inference = NULL,
                           method.tte = NULL,
                           model.tte,
                           n.resampling = NULL,
                           name.call,
                           neutral.as.uninf = NULL,
                           operator,
                           option,
                           seed = NULL,
                           strata,
                           threshold,
                           trace = NULL,
                           treatment,
                           type,
                           weight){

    ## ** apply default options
    if(is.null(cpus)){ cpus <- option$cpus }
    if(is.null(keep.pairScore)){ keep.pairScore <- option$keep.pairScore }
    if(is.null(method.tte)){ method.tte <- option$method.tte }
    if(is.null(hierarchical)){ hierarchical <- option$hierarchical }
    if(is.null(correction.uninf)){ correction.uninf <- option$correction.uninf }
    if(is.null(method.inference)){ method.inference <- option$method.inference }
    if(is.null(n.resampling)){ n.resampling <- option$n.resampling }
    if(is.null(neutral.as.uninf)){ neutral.as.uninf <- option$neutral.as.uninf }
    if(is.null(trace)){ trace <- option$trace }
    if(is.null(alternative)){ alternative <- option$alternative }

    ## ** convert formula into separate arguments
    if(!missing(formula)){
        ## the missing is for BuysePower where the arguments are not necessarily specified
        test.null <- c(censoring = !missing(censoring) && !is.null(censoring),
                       endpoint = !missing(endpoint) && !is.null(endpoint),
                       operator = !missing(operator) && !is.null(operator),
                       strata = !missing(strata) && !is.null(strata),
                       threshold = !missing(threshold) && !is.null(threshold),
                       treatment = !missing(treatment) && !is.null(treatment),
                       type = !missing(type) && !is.null(type),
                       weight = !missing(weight) && !is.null(weight)
                       )
        if(any(test.null)){
            txt <- names(test.null)[test.null]
            warning("Argument",if(sum(test.null)>1){"s"}," \'",paste(txt, collpase="\' \'"),if(sum(test.null)>1){" are "}else{" is "}," ignored when argument \'formula\' has been specified\n")
        }
        
        resFormula <- initializeFormula(formula)

        treatment <- resFormula$treatment
        type <- resFormula$type
        endpoint <- resFormula$endpoint
        threshold <- resFormula$threshold
        censoring <- resFormula$censoring
        weight <- resFormula$weight
        operator <- resFormula$operator
        strata <- resFormula$strata
    }else{
        if(is.null(operator)){
            operator <- rep(">0",length(endpoint))
        }
        formula <- NULL
        if(is.null(weight)){
            weight <- rep(1,length(endpoint))
        }
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
    attr(method.inference,"permutation") <- grepl("permutation",method.inference)
    attr(method.inference,"bootstrap") <- grepl("bootstrap",method.inference)
    attr(method.inference,"studentized") <- grepl("studentized",method.inference)
    attr(method.inference,"stratified") <- grepl("stratified",method.inference)
    attr(method.inference,"ustatistic") <- grepl("asymptotic",method.inference)

    iid <- any(c(attr(method.inference,"studentized"), method.inference == "asymptotic"))
    if(is.null(strata) && length(grep("stratified ",method.inference))>0){ ## remove stratified if no strata variable
        method.inference <- gsub("stratified ","",method.inference)
        attr(method.inference,"stratified") <- FALSE
    }

    ## ** censoring
    if(D.TTE==0){
        censoring <- rep("..NA..", D)
    }else if(length(censoring) == D.TTE){
        censoring.save <- censoring
        censoring <- rep("..NA..", D)
        censoring[type==3] <- censoring.save 
    }
    ## from now, censoring contains for each endpoint the name of variable indicating censoring (0) or event (1) or NA
    
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
        iid = iid,
        keep.pairScore = keep.pairScore,
        keep.survival = option$keep.survival,
        method.tte = method.tte,
        model.tte = model.tte,
        method.inference = method.inference,
        n.resampling = n.resampling,
        hierarchical = hierarchical,
        neutral.as.uninf = neutral.as.uninf,
        operator = operator,
        seed = seed,
        strata = strata,
        threshold = threshold,
        trace = trace,
        treatment = treatment,
        type = type,
        weight = weight
    ))
}

## * initializeData
#' @rdname internal-initialization
initializeData <- function(data, type, endpoint, method.tte, censoring, operator, strata, treatment, copy){

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }else if(copy){
        data <- data.table::copy(data)
    }

    ## ** convert character/factor to numeric for binary endpoints
    name.bin <- endpoint[which(type %in% 1)]
    if(length(name.bin)>0){
        data.class <- sapply(data,class)
        test.num <- (data.class %in% c("numeric","integer"))
        if(any(test.num==FALSE)){
            endpoint.char <- names(data.class)[test.num==FALSE]
            for(iE in endpoint.char){
                data[, c(iE) := as.double(as.factor(.SD[[1]]))-1.0, .SDcols = iE]
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

    ## ** n.obs
    n.obs <- data[,.N]

    ## ** strata
    if(!is.null(strata)){  
    
        data[ , c("..strata..") := interaction(.SD, drop = TRUE, lex.order = FALSE, sep = "."), .SDcols = strata]
        level.strata <- levels(data[["..strata.."]])        
        data[ , c("..strata..") := as.numeric(.SD[["..strata.."]])] # convert to numeric
        data.table::setkeyv(data, cols = c("..strata..", treatment)) ## important to first sort by strata when doing the resampling
        
        n.obsStrata <- data[,.N, by = "..strata.."][,setNames(.SD[[1]],.SD[[2]]),.SD = c("N","..strata..")]
    }else{
        data.table::setkeyv(data, cols = treatment)
        
        data[ , c("..strata..") := 1]
        n.obsStrata <- n.obs
        level.strata <- 1
    }

    n.strata <- length(level.strata)

    ## ** convert treatment to binary indicator
    level.treatment <- levels(as.factor(data[[treatment]]))
    trt2bin <- setNames(0:1,level.treatment)
    data[ , c(treatment) := trt2bin[as.character(.SD[[1]])], .SDcols = treatment]

    ## ** rowIndex
    data[,c("..rowIndex..") := 1:.N]

    ## ** unique censoring
    Ucensoring <- unique(censoring)
    if(any(censoring == "..NA..")){
        data[,c("..NA..") := as.numeric(NA)]
    }

    ## ** unique endpoint
    Uendpoint <- unique(endpoint)

    ## ** scoring method for each endpoint
    method.score <- 1 + (type==3) ## 1 binary/continuous and 2 Gehan
    if(method.tte > 0){ ## if Peron
        test.CR <- sapply(Ucensoring, function(iC){max(data[[iC]])>1})[censoring]
        method.score[type == 3] <- method.score[type == 3] + 1 + test.CR[type == 3]
        ## 3 Peron survival
        ## 4 Peron CR
    }
    
    ## ** export
    return(list(data = data[,.SD,.SDcols =  c(treatment,"..strata..")],
                M.endpoint = as.matrix(data[, .SD, .SDcols = Uendpoint]),
                M.censoring = as.matrix(data[, .SD, .SDcols = Ucensoring]),
                index.C = which(data[[treatment]] == 0),
                index.T = which(data[[treatment]] == 1),
                index.strata = tapply(data[["..rowIndex.."]], data[["..strata.."]], list),
                index.endpoint = match(endpoint, Uendpoint) - 1,
                index.censoring = match(censoring, Ucensoring) - 1,
                level.treatment = level.treatment,
                level.strata = level.strata,
                method.score = method.score,
                n.strata = n.strata,
                n.obs = n.obs,
                n.obsStrata = n.obsStrata,
                cumn.obsStrata = cumsum(c(0,n.obsStrata))[1:n.strata]
                ))
}

## * buildWscheme
buildWscheme <- function(method.tte, endpoint, D.TTE, D, n.strata,
                         type, threshold){

    Wscheme <- matrix(0,nrow=D,ncol=D) # design matrix indicating to combine the weights obtained at differents endpoints
    rownames(Wscheme) <- paste("weigth of ",endpoint,"(",threshold,")",sep="")
    colnames(Wscheme) <- paste("for ",endpoint,"(",threshold,")",sep="")
    Wscheme[upper.tri(Wscheme)] <- 1 ## only previous endpoint can contribute to the current weights
    Wscheme[lower.tri(Wscheme)] <- NA ## do not look at future endpoint 
        
    ## take care of repeated survival endpoints
    if(D.TTE>1 && method.tte > 0){

        index.TTE <- which(type == 3)
        indexDuplicated.TTE <- which(duplicated(endpoint[index.TTE]))

        for(iEndpoint in indexDuplicated.TTE){    ## iEndpoint <- indexDuplicated.endpoint.TTE[1]
            iEndpoint2 <- index.TTE[iEndpoint] ## position of the current endpoint relative to all endpoint
            iEndpoint2_M1 <- which(endpoint[1:(iEndpoint2-1)] == endpoint[iEndpoint2])  ## position of the previous endpoint(s) relative to all endpoints
            Wscheme[iEndpoint2_M1,iEndpoint2] <- 0
        }
    }

    ## skeleton
    skeleton <- lapply(1:D, function(iE){
        lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})
    })

    ## unique tte endpoint
    endpoint.UTTE <- unique(endpoint[type==3])
    
    ## export
    return(list(Wscheme = Wscheme,
                endpoint.UTTE = endpoint.UTTE,
                index.UTTE = match(endpoint, endpoint.UTTE, nomatch = 0) - 1,
                D.UTTE = length(endpoint.UTTE),
                reanalyzed = rev(duplicated(rev(endpoint))),
                outSurv = list(survTimeC = skeleton,
                               survTimeT = skeleton,
                               survJumpC = skeleton,
                               survJumpT = skeleton,
                               lastSurv = lapply(1:D, function(iS){matrix(nrow = n.strata, ncol = 4)}) ## 4 for competing risk setting, 2 is enough for survival
                               ))
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
    if(n.endpoint==0){
        stop("initFormula: x must contain endpoints \n",
             "nothing of the form type(endpoint,threshold,censoring) found in the formula \n")
    }

    ## ** extract endpoints and additional arguments 
    threshold <- rep(NA, n.endpoint)
    censoring <- rep("..NA..", n.endpoint)
    endpoint <- rep(NA, n.endpoint)
    operator <- rep(">0", n.endpoint)
    weight <- rep(1, n.endpoint)
    validArgs <- c("endpoint","threshold","censoring","operator","weight")

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
            thresholdTempo <- eval(expr = parse(text = iArg[iName=="threshold"]))

            if(inherits(thresholdTempo, "function")){
                packageTempo <- environmentName(environment(thresholdTempo))
                if(nchar(packageTempo)>0){
                    txt <- paste0("(package ",packageTempo,")")
                }else{
                    txt <- ""
                }
                stop(iArg[iName=="threshold"]," is already defined as a function ",txt,"\n",
                     "cannot be used to specify the threshold \n")
            }
            
            threshold[iE] <- as.numeric(thresholdTempo)
        }
        if("censoring" %in% iName){
            censoring[iE] <- gsub("\"","",iArg[iName=="censoring"])
        }
        if("operator" %in% iName){
            operator[iE] <- gsub("\"","",iArg[iName=="operator"])
        }
        if("weight" %in% iName){
            weight[iE] <- as.numeric(eval(expr = parse(text = iArg[iName=="weight"])))
        }
    }

    ## ** export
    return(list(treatment = treatment,
                type = type,
                endpoint = endpoint,
                threshold = threshold,
                censoring = censoring,
                operator = operator,
                weight = weight,
                strata = strata))
}

## * initializePeron
#' @rdname internal-initialization
initializePeron <- function(data,
                            model.tte,
                            method.score,
                            treatment,
                            level.treatment,
                            endpoint,
                            endpoint.UTTE,
                            censoring,
                            D.TTE,
                            D.UTTE,
                            type,
                            strata,
                            threshold,
                            n.strata,
                            out){

    . <- NULL ## for CRAN check
        
    ## ** prepare
    if(n.strata == 1){
        ls.indexC <- list(which(data[[treatment]]==0))
        ls.indexT <- list(which(data[[treatment]]==1))
    }else{
        indexC <- which(data[[treatment]]==0)
        indexT <- which(data[[treatment]]==1)
        ls.indexC <- vector(mode = "list", length = n.strata)
        ls.indexT <- vector(mode = "list", length = n.strata)
        for(iStrata in 1:n.strata){
            iIndex.strata <- which(data[["..strata.."]]==iStrata)
            ls.indexC[[iStrata]] <- intersect(indexC,iIndex.strata)
            ls.indexT[[iStrata]] <- intersect(indexT,iIndex.strata)
        }
    }            
    
    zeroPlus <- 1e-12

    ## ** unique tte endpoints
    endpoint.TTE <- endpoint[type==3]
    censoring.TTE <- censoring[type==3]
    
    index.endpoint.UTTE <- which(!duplicated(endpoint.TTE))
    censoring.UTTE <- censoring.TTE[index.endpoint.UTTE]

    ## ** estimate cumulative incidence function (survival case or competing risk case)
    if(is.null(model.tte)){
        model.tte <- vector(length = D.UTTE, mode = "list")
        names(model.tte) <- endpoint.UTTE

        txt.modelUTTE <- paste0("prodlim::Hist(",endpoint.UTTE,",",censoring.UTTE,") ~ ",treatment," + ..strata..")

        for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.UTTE <- 1
            model.tte[[iEndpoint.UTTE]] <- prodlim::prodlim(as.formula(txt.modelUTTE[iEndpoint.UTTE]),
                                                            data = data)
        }
        
    }else{ 
        for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.TTE <- 1
            ## convert treatment to numeric
            model.tte[[iEndpoint.UTTE]]$X[[treatment]] <- as.numeric(factor(model.tte[[iEndpoint.UTTE]]$X[[treatment]], levels = level.treatment))-1
            p <- NCOL(model.tte[[iEndpoint.UTTE]]$X)

            ## create ..strata..
            if(p==1){
                model.tte[[iEndpoint.UTTE]]$X <- cbind(model.tte[[iEndpoint.UTTE]]$X,
                                                       "..strata.." = 1)
            }else{
                col.strata <- setdiff(1:p,which(colnames(model.tte[[iEndpoint.UTTE]]$X)==treatment))
                value.strata <- apply(model.tte[[iEndpoint.UTTE]]$X[,col.strata],1,paste0,collapse="")
                model.tte[[iEndpoint.UTTE]]$X <- cbind(model.tte[[iEndpoint.UTTE]]$X,
                                                       "..strata.." = as.numeric(as.factor(value.strata)))

            }
        }
    }

    ## ** predict individual survival
    ## *** fill
    for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.TTE <- 1
        iEndpoint.UTTE.name <- endpoint.UTTE[iEndpoint.UTTE]
        iIndex.associatedEndpoint <- which(endpoint == iEndpoint.UTTE.name)
        iTest.CR <- method.score[iIndex.associatedEndpoint[1]]==4

        if(iTest.CR){
            index.jump1 <- which(model.tte[[iEndpoint.UTTE]]$cause.hazard[[1]]>0)
            index.jump2 <- which(model.tte[[iEndpoint.UTTE]]$cause.hazard[[2]]>0)
        }else{
            index.jump <- which(model.tte[[iEndpoint.UTTE]]$hazard>0)
        }
        
        indexX.C <- which(model.tte[[iEndpoint.UTTE]]$X[[treatment]]==0)
        indexX.T <- which(model.tte[[iEndpoint.UTTE]]$X[[treatment]]==1)

        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            iNcontrol <- length(ls.indexC[[iStrata]])
            iNtreatment <- length(ls.indexT[[iStrata]])

            if("..strata.." %in% colnames(model.tte[[iEndpoint.UTTE]]$X)){
                indexX.strata <- which(model.tte[[iEndpoint.UTTE]]$X[["..strata.."]]==iStrata)
                indexX.strataC <- intersect(indexX.C,indexX.strata)
                indexX.strataT <- intersect(indexX.T,indexX.strata)
            }

            iIndex.startC <- model.tte[[iEndpoint.UTTE]]$first.strata[indexX.strataC]
            iIndex.startT <- model.tte[[iEndpoint.UTTE]]$first.strata[indexX.strataT]

            iIndex.stopC <- iIndex.startC + model.tte[[iEndpoint.UTTE]]$size.strata[indexX.strataC] - 1
            iIndex.stopT <- iIndex.startT + model.tte[[iEndpoint.UTTE]]$size.strata[indexX.strataT] - 1

            iTimeC <- data[ls.indexC[[iStrata]],.SD[[iEndpoint.UTTE.name]]]
            iTimeT <- data[ls.indexT[[iStrata]],.SD[[iEndpoint.UTTE.name]]]

            if(iTest.CR){
                iJump1C <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump1,iIndex.startC:iIndex.stopC)]
                iJump1T <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump1,iIndex.startT:iIndex.stopT)]
                iJump2C <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump2,iIndex.startC:iIndex.stopC)]
                iJump2T <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump2,iIndex.startT:iIndex.stopT)]

                iLast.cif1C <- model.tte[[iEndpoint.UTTE]]$cuminc[[1]][iIndex.stopC]
                iLast.cif1T <- model.tte[[iEndpoint.UTTE]]$cuminc[[1]][iIndex.stopT]
                iLast.cif2C <- model.tte[[iEndpoint.UTTE]]$cuminc[[2]][iIndex.stopC]
                iLast.cif2T <- model.tte[[iEndpoint.UTTE]]$cuminc[[2]][iIndex.stopT]

                sumCifC = iLast.cif1C + iLast.cif2C
                sumCifT = iLast.cif1T + iLast.cif2T
            
                iPredCif1C <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startC:iIndex.stopC], 
                                               y = model.tte[[iEndpoint.UTTE]]$cuminc[[1]][iIndex.startC:iIndex.stopC],
                                               yleft = 0, yright = switch(as.character(sumCifC == 1),
                                                                          "TRUE" = iLast.cif1C,
                                                                          "FALSE" = NA), f = 0, method = "constant")
            
                iPredCif1T <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startT:iIndex.stopT], 
                                               y = model.tte[[iEndpoint.UTTE]]$cuminc[[1]][iIndex.startT:iIndex.stopT],
                                               yleft = 0, yright = switch(as.character(sumCifT == 1),
                                                                          "TRUE" = iLast.cif1T,
                                                                          "FALSE" = NA), f = 0, method = "constant")
            
                iPredCif2C <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startC:iIndex.stopC], 
                                               y = model.tte[[iEndpoint.UTTE]]$cuminc[[2]][iIndex.startC:iIndex.stopC],
                                               yleft = 0, yright = switch(as.character(sumCifC == 1),
                                                                          "TRUE" = iLast.cif2C,
                                                                          "FALSE" = NA), f = 0, method = "constant")
            
                iPredCif2T <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startT:iIndex.stopT], 
                                               y = model.tte[[iEndpoint.UTTE]]$cuminc[[2]][iIndex.startT:iIndex.stopT],
                                               yleft = 0, yright = switch(as.character(sumCifT == 1),
                                                                          "TRUE" = iLast.cif2T,
                                                                          "FALSE" = NA), f = 0, method = "constant")

                ## independent of the threshold i.e. of the priority
                ## avoid repeated calculation when the same endpoint is used several times with different thresholds
                iDCif1C.jumpC <- iPredCif1C(iJump1C) - iPredCif1C(iJump1C - zeroPlus)
                iCif1C.timeC <- iPredCif1C(iTimeC)
                iCif1C.timeT <- iPredCif1C(iTimeT)
                iCif1T.timeC <- iPredCif1T(iTimeC)
                iCif1T.timeT <- iPredCif1T(iTimeT)
                iCif2C.timeC <- iPredCif2C(iTimeC)
                                        #iCif2C.timeT <- iPredCif2C(iTimeT)
                                        #iCif2T.timeC <- iPredCif2T(iTimeC)
                iCif2T.timeT <- iPredCif2T(iTimeT)
            
            for(iEndpoint in iIndex.associatedEndpoint){ ## iEndpoint <- 1
              iThreshold <- threshold[iEndpoint] ## iThreshold = 1
              
              ## **** last survival
              out$lastSurv[[iEndpoint]][iStrata,] <- cbind(iLast.cif1C, iLast.cif1T, iLast.cif2C, iLast.cif2T)
              
              ## **** survival at jump times
              out$survJumpC[[iEndpoint]][[iStrata]] <- cbind(time = iJump1C,
                                                            "CIF1T-threshold" = iPredCif1T(iJump1C - iThreshold),
                                                            "CIF1T+threshold" = iPredCif1T(iJump1C + iThreshold),
                                                            dCIF = iDCif1C.jumpC)
              
              ## **** survival at observation time (+/- threshold)
              out$survTimeC[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeC,
                                                            "CIF1C-threshold" = iPredCif1C(iTimeC - iThreshold),
                                                            "CIF1C_0" = iCif1C.timeC,
                                                            "CIF1C+threshold" = iPredCif1C(iTimeC + iThreshold),
                                                            "CIF1T-threshold" = iPredCif1T(iTimeC - iThreshold),
                                                            "CIF1T_0" = iCif1T.timeC,
                                                            "CIF1T+threshold" = iPredCif1T(iTimeC + iThreshold),
                                                            "CIF2C_0" = iCif2C.timeC)#,
                                                            #"CIF2T_0" = iCif2T.timeC)
              
              out$survTimeT[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeT,
                                                             "CIF1C-threshold" = iPredCif1C(iTimeT - iThreshold),
                                                             "CIF1C_0" = iCif1C.timeT,
                                                             "CIF1C+threshold" = iPredCif1C(iTimeT + iThreshold),
                                                             "CIF1T-threshold" = iPredCif1T(iTimeT - iThreshold),
                                                             "CIF1T_0" = iCif1T.timeT,
                                                             "CIF1T+threshold" = iPredCif1T(iTimeT + iThreshold),
                                        #"CIF2C_0" = iCif2C.timeT,
                                                             "CIF2T_0" = iCif2T.timeT)
            }

            }else{
                iJumpC <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump,iIndex.startC:iIndex.stopC)]
                iJumpT <- model.tte[[iEndpoint.UTTE]]$time[intersect(index.jump,iIndex.startT:iIndex.stopT)]

                iLast.survC <- model.tte[[iEndpoint.UTTE]]$surv[iIndex.stopC]
                iLast.survT <- model.tte[[iEndpoint.UTTE]]$surv[iIndex.stopT]

                iPredSurvC <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startC:iIndex.stopC],
                                               y = model.tte[[iEndpoint.UTTE]]$surv[iIndex.startC:iIndex.stopC],
                                               yleft=1, yright=switch(as.character(iLast.survC<zeroPlus),
                                                                      "TRUE" = 0,
                                                                      "FALSE" = NA),
                                               f=0,
                                               method = "constant")
            
                iPredSurvT <- stats::approxfun(x = model.tte[[iEndpoint.UTTE]]$time[iIndex.startT:iIndex.stopT],
                                               y = model.tte[[iEndpoint.UTTE]]$surv[iIndex.startT:iIndex.stopT],
                                               yleft=1, yright=switch(as.character(iLast.survT<zeroPlus),
                                                                      "TRUE" = 0,
                                                                      "FALSE" = NA),
                                               f=0,
                                               method = "constant")


                ## independent of the threshold i.e. of the priority
                ## avoid repeated calculation when the same endpoint is used several times with different thresholds
                iDSurvivalC.jumpC <- iPredSurvC(iJumpC) - iPredSurvC(iJumpC - zeroPlus)
                iDSurvivalT.jumpT <- iPredSurvT(iJumpT) - iPredSurvT(iJumpT - zeroPlus)
                iSurvivalC.timeC <- iPredSurvC(iTimeC)
                iSurvivalC.timeT <- iPredSurvC(iTimeT)
                iSurvivalT.timeC <- iPredSurvT(iTimeC)
                iSurvivalT.timeT <- iPredSurvT(iTimeT)
            
            for(iEndpoint in iIndex.associatedEndpoint){ ## iEndpoint <- 1
                iThreshold <- threshold[iEndpoint]

                ## **** last survival
                out$lastSurv[[iEndpoint]][iStrata,1:2] <- c(iLast.survC, iLast.survT)

                ## **** survival at jump times
                out$survJumpC[[iEndpoint]][[iStrata]] <- cbind(time = iJumpC,
                                                               survival = iPredSurvT(iJumpC + iThreshold),
                                                               dSurvival = iDSurvivalC.jumpC)
            
                out$survJumpT[[iEndpoint]][[iStrata]] <- cbind(time = iJumpT,
                                                               survival = iPredSurvC(iJumpT + iThreshold),
                                                               dSurvival = iDSurvivalT.jumpT)

                ## **** survival at observation time (+/- threshold)
                out$survTimeC[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeC,
                                                               "SurvivalC-threshold" = iPredSurvC(iTimeC - iThreshold),
                                                               "SurvivalC_0" = iSurvivalC.timeC,
                                                               "SurvivalC+threshold" = iPredSurvC(iTimeC + iThreshold),
                                                               "SurvivalT-threshold" = iPredSurvT(iTimeC - iThreshold),
                                                               "SurvivalT_0" = iSurvivalT.timeC,
                                                               "SurvivalT+threshold" = iPredSurvT(iTimeC + iThreshold)
                                                               )

                out$survTimeT[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeT,
                                                               "SurvivalC-threshold" = iPredSurvC(iTimeT - iThreshold),
                                                               "SurvivalC_0" = iSurvivalC.timeT,
                                                               "SurvivalC+threshold" = iPredSurvC(iTimeT + iThreshold),
                                                               "SurvivalT-threshold" = iPredSurvT(iTimeT - iThreshold),
                                                               "SurvivalT_0" = iSurvivalT.timeT,
                                                               "SurvivalT+threshold" = iPredSurvT(iTimeT + iThreshold)
                                                               )
            }

            }
            
            


        }
    }
    ## export
    return(out)
    
}



