## * Documentation initialization functions called by BuyseTest

#' @title internal functions for BuyseTest - initialization
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' @name internal-initialization
#' 
#' @details
#' 
#' \code{initializeArgs}: Normalize the argument 
#' \itemize{
#' \item scoring.rule, neutral.as.uninf, keep.pairScore, n.resampling, seed, cpus, trace: set to default value when not specified.
#' \item formula: call \code{initializeFormula} to extract arguments.
#' \item type: convert to numeric.
#' \item censoring: only keep censoring relative to TTE endpoint. Set to \code{NULL} if no TTE endpoint.
#' \item threshold: set default threshold to 1e-12 expect for binary variable where it is set to 1/2.
#' the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12.
#' When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' \item data: convert to data.table object.
#' \item scoring.rule: convert to numeric.
#' }
#'
#' \code{initializeFormula}:  extract \code{treatment}, \code{type}, \code{endpoint}, \code{threshold}, \code{censoring}, \code{operator}, and \code{strata}
#' from the formula. \cr \cr
#'
#' \code{initializeData}: Divide the dataset into two, one relative to the treatment group and the other relative to the control group.
#' Merge the strata into one with the interaction variable.
#' Extract for each strata the index of the observations within each group.
#'
#' \code{initializePeron}: Compute the survival via KM.
#' 
#' @keywords function internal BuyseTest

## * initializeArgs
#' @rdname internal-initialization
initializeArgs <- function(censoring,
                           correction.uninf = NULL,
                           cpus = NULL,
                           data,
                           endpoint,
                           formula,
                           hierarchical = NULL,
                           keep.pairScore = NULL,
                           method.inference = NULL,
                           scoring.rule = NULL,
                           model.tte,
                           n.resampling = NULL,
                           strata.resampling = NULL,
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
    if(is.null(scoring.rule)){ scoring.rule <- option$scoring.rule }
    if(is.null(hierarchical)){ hierarchical <- option$hierarchical }
    if(is.null(correction.uninf)){ correction.uninf <- option$correction.uninf }
    if(is.null(method.inference)){ method.inference <- option$method.inference }
    if(is.null(n.resampling)){ n.resampling <- option$n.resampling }
    if(is.null(strata.resampling)){ strata.resampling <- option$strata.resampling }
    if(is.null(neutral.as.uninf)){ neutral.as.uninf <- option$neutral.as.uninf }
    if(is.null(trace)){ trace <- option$trace }
    alternative <- option$alternative
    
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
    
    ## ** type
    if(!is.numeric(type)){
        validType1 <- c("b","bin","binary")
        validType2 <- c("c","cont","continuous")
        validType3 <- c("t","tte","time","timetoevent") ## [if modified, remember to change the corresponding vector in initFormula]
        type <- tolower(type)

        type[grep(paste(validType1,collapse="|"), type)] <- "1" 
        type[grep(paste(validType2,collapse="|"), type)] <- "2" 
        type[grep(paste(validType3,collapse="|"), type)] <- "3"
        type <- sapply(type, function(iType){
            switch(iType,
                   "1" = 1, ## binary endpoint
                   "2" = 2, ## continuous endpoint
                   "3" = 3, ## time to event endpoint
                   NA)})
    }
 
    ## ** endpoint
    index.type3 <- which(type==3)
    
    D <- length(endpoint)
    Uendpoint <- unique(endpoint)
    
    ## time to event endpoints 
    endpoint.TTE <- endpoint[index.type3]
    threshold.TTE <- threshold[index.type3]
    D.TTE <- length(endpoint.TTE) 

    ## ** censoring
    if(D.TTE==0){
        censoring <- rep("..NA..", D)
    }else if(length(censoring) == D.TTE){
        censoring.save <- censoring
        censoring <- rep("..NA..", D)
        censoring[index.type3] <- censoring.save 
    }
    Ucensoring <- unique(censoring)
    censoring.TTE <- censoring[index.type3]
    ## from now, censoring contains for each endpoint the name of variable indicating censoring (0) or event (1) or NA
    
    
    ## ** scoring.rule
    ## WARNING: choices must be lower cases
    ##          remember to update check scoring.rule (in BuyseTest-check.R)
    scoring.rule <- tolower(scoring.rule)
    scoring.rule <- switch(scoring.rule,
                           "gehan" = 0,
                           "peron" = 1,
                           NA
                           )

    if (D.TTE == 0) {
        scoring.rule <- 0
        if ("scoring.rule" %in% name.call && trace > 0) {
            message("NOTE : there is no survival endpoint, \'scoring.rule\' argument is ignored \n")
        }
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

    ## ** method.inference
    method.inference <- tolower(method.inference)    
    attr(method.inference,"permutation") <- grepl("permutation",method.inference)
    attr(method.inference,"bootstrap") <- grepl("bootstrap",method.inference)
    attr(method.inference,"studentized") <- grepl("studentized",method.inference)
    attr(method.inference,"ustatistic") <- grepl("u-statistic",method.inference)
    if(identical(strata.resampling,"treatment")){
        attr(method.inference,"resampling-strata:treatment") <- TRUE
        attr(method.inference,"resampling-strata:strata") <- FALSE
    }else if(identical(strata.resampling,"strata")){
        attr(method.inference,"resampling-strata:treatment") <- FALSE
        attr(method.inference,"resampling-strata:strata") <- TRUE
    }else if(is.na(strata.resampling) || length(strata.resampling)== 0){
        attr(method.inference,"resampling-strata:treatment") <- FALSE
        attr(method.inference,"resampling-strata:strata") <- FALSE
    }else{
        attr(method.inference,"resampling-strata:treatment") <- NA
        attr(method.inference,"resampling-strata:strata") <- NA
    }
    
    ## ** correction.uninf
    correction.uninf <- as.numeric(correction.uninf)

    ## ** model.tte
    if(identical(scoring.rule,1)){
        if((!is.null(model.tte)) && (length(unique(endpoint.TTE)) == 1) && inherits(model.tte, "prodlim")){
            attr.save <- attr(model.tte,"iidNuisance")
            
            model.tte <- list(model.tte)
            names(model.tte) <- unique(endpoint.TTE)
            attr(model.tte,"iidNuisance") <- attr.save
        }
    }else{
        model.tte <- NULL
    }
    
    ## ** iid
    iid <- attr(method.inference,"studentized") || (method.inference == "u-statistic")
    if(iid){
        attr(method.inference,"hprojection") <- option$order.Hprojection
        if(attr(method.inference,"hprojection")==2 & identical(scoring.rule,1)){
            keep.pairScore <- TRUE ## need the detail of the score to perform the 2nd order projection
        }
    }else{
        attr(method.inference,"hprojection") <- NA
    }
    iidNuisance <- iid && identical(scoring.rule,1) && (is.null(model.tte) || identical(attr(model.tte,"iidNuisance"),TRUE))

    ## ** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }

    ## ** export
    return(list(
        name.call = name.call,
        censoring = censoring,
        censoring.TTE = censoring.TTE,
        correction.uninf = correction.uninf,
        cpus = cpus,
        D = D,
        D.TTE = D.TTE,
        data = data,
        endpoint = endpoint,
        endpoint.TTE = endpoint.TTE,
        formula = formula,
        iid = iid,
        iidNuisance = iidNuisance,
        index.endpoint = match(endpoint, Uendpoint) - 1,
        index.censoring = match(censoring, Ucensoring) - 1,
        keep.pairScore = keep.pairScore,
        keep.survival = option$keep.survival,
        scoring.rule = scoring.rule,
        model.tte = model.tte,
        method.inference = method.inference,
        n.resampling = n.resampling,
        hierarchical = hierarchical,
        neutral.as.uninf = neutral.as.uninf,
        operator = operator,
        order.Hprojection = option$order.Hprojection,
        seed = seed,
        strata = strata,
        threshold = threshold,
        trace = trace,
        treatment = treatment,
        type = type,
        Uendpoint = Uendpoint,
        Ucensoring = Ucensoring,
        weight = weight,
        debug = option$debug
    ))
}

## * initializeData
#' @rdname internal-initialization
initializeData <- function(data, type, endpoint, Uendpoint, D, scoring.rule, censoring, Ucensoring, method.inference, operator, strata, treatment, hierarchical, copy,
                           endpoint.TTE, censoring.TTE, iidNuisance){

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
        
        n.obsStrata <- data[,.N, by = "..strata.."][,setNames(.SD[[1]],.SD[[2]]),.SD = c("N","..strata..")]
    }else{
        
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
    if(any(censoring == "..NA..")){
        data[,c("..NA..") := as.numeric(NA)]
    }

    ## ** TTE with censoring
    if(scoring.rule>0){
        test.censoring <- sapply(censoring.TTE, function(iC){any(data[[iC]]==0)})
        if(all(test.censoring==FALSE)){
            scoring.rule <- 0
            iidNuisance <- FALSE
        }        

        ## distinct time to event endpoints
        endpoint.UTTE <- unique(endpoint.TTE[test.censoring])
        censoring.UTTE <- unique(censoring.TTE[test.censoring])
        D.UTTE <- length(endpoint.UTTE)

        ## correspondance endpoint, TTE endpoint (non TTEe endpoint are set to -100)
        index.UTTE = match(endpoint, endpoint.UTTE, nomatch = -99) - 1
    }else{
        endpoint.UTTE <- numeric(0)
        censoring.UTTE <- numeric(0)
        D.UTTE <- 0
        index.UTTE <- rep(-100, D)
    }
    
    ## ** scoring method for each endpoint
    ## check if censoring
    method.score <- 1 + (type==3) ## 1 binary/continuous and 2 Gehan
    if(scoring.rule > 0){ ## if Peron
        test.CR <- sapply(Ucensoring, function(iC){max(data[[iC]])>1})[censoring]
        method.score[type == 3] <- method.score[type == 3] + test.censoring * (1 + test.CR[type == 3])
        ## 3 Peron survival
        ## 4 Peron CR
    }


    ## ** previously analyzed distinct TTE endpoints
    if((scoring.rule==1) && hierarchical){ ## only relevant when using Peron scoring rule with hierarchical GPC
        ## number of distinct, previously analyzed, TTE endpoints
        nUTTE.analyzedPeron_M1 <- sapply(1:D, function(iE){
            if(iE>1){
                sum(endpoint.UTTE %in% endpoint[1:(iE-1)])
            }else{
                return(0)
            }
        })
    }else{
        nUTTE.analyzedPeron_M1 <- rep(0,D)
    }

    ## ** number of observations per strata used when resampling
    index.C <- which(data[[treatment]] == 0)
    index.T <- which(data[[treatment]] == 1)
    
    if(attr(method.inference,"resampling-strata:treatment")){
        n.obsStrataResampling <- c(length(index.C), length(index.T))
    }else if(attr(method.inference,"resampling-strata:strata")){
        n.obsStrataResampling <- n.obsStrata
    }else{
        n.obsStrataResampling <- n.obs
    }
    
    ## ** skeleton for survival proba (only relevant for Peron scoring rule)
    skeletonPeron <- list(survTimeC = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survTimeT = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survJumpC = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survJumpT = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          lastSurv = lapply(1:D, function(iS){matrix(nrow = n.strata, ncol = 4)}), ## 4 for competing risk setting, 2 is enough for survival
                          p.C = matrix(NA, nrow = n.strata, ncol = D),
                          p.T = matrix(NA, nrow = n.strata, ncol = D)
                          )


    ## ** export
    return(list(data = data[,.SD,.SDcols =  c(treatment,"..strata..")],
                M.endpoint = as.matrix(data[, .SD, .SDcols = Uendpoint]),
                M.censoring = as.matrix(data[, .SD, .SDcols = Ucensoring]),
                index.C = index.C,
                index.T = index.T,
                index.strata = tapply(data[["..rowIndex.."]], data[["..strata.."]], list),
                level.treatment = level.treatment,
                level.strata = level.strata,
                method.score = method.score,
                n.strata = n.strata,
                n.obs = n.obs,
                n.obsStrata = n.obsStrata,
                n.obsStrataResampling = n.obsStrataResampling,
                skeletonPeron = skeletonPeron,
                scoring.rule = scoring.rule,
                iidNuisance = iidNuisance,
                nUTTE.analyzedPeron_M1 = nUTTE.analyzedPeron_M1,
                endpoint.UTTE = endpoint.UTTE,
                censoring.UTTE = censoring.UTTE,
                D.UTTE = D.UTTE,
                index.UTTE = index.UTTE
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
  
    ## remove all blanks
    x.rhs <- gsub("[[:blank:]]", "", x.rhs)

    ## find endpoints
    ## https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
    ## (*SKIP)(*FAIL): ignore
    ## \\( \\): inside brackets
    ## [^()]*: anything but ()
    magic.formula <- "\\([^()]*\\)(*SKIP)(*FAIL)|\\h*\\+\\h*"
    vec.x.rhs <- unlist(strsplit(x.rhs, split = magic.formula, perl = TRUE))
    ## find all element in the vector corresponding to endpoints (i.e. ...(...) )
    ## \\w* any letter/number
    ## [[:print:]]* any letter/number/punctuation/space
    index.endpoint <- grep("\\w*\\([[:print:]]*\\)$", vec.x.rhs)
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
    validArgs <- c("endpoint","censoring","threshold","operator","weight")

    ## split around parentheses
    ls.x.endpoint <- strsplit(vec.x.endpoint, split = "(", fixed = TRUE)

    type <- character(length = n.endpoint)
    for(iE in 1:n.endpoint){
        ## extract type
        type[iE] <- tolower(ls.x.endpoint[[iE]][1])
        if(type[iE] %in% c("b","bin","binary")){
            iValidArgs <- setdiff(validArgs,c("censoring","threshold"))
        }else if(type[iE] %in% c("c","cont","continuous")){
            iValidArgs <- setdiff(validArgs,"censoring")
        }else{ ## if(type[iE] %in% c("t","tte","time","timetoevent"))
            iValidArgs <- validArgs
        }
        
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
            if(any(iiName %in% iValidArgs == FALSE)){
                stop("initFormula: invalid formula \n",
                     vec.x.rhs[iE]," contains arguments that are not \"",paste0(iValidArgs,sep = "\" \""),"\" \n")
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
            iName[setdiff(1:n.args,iIndex.name)] <- setdiff(iValidArgs,iiName)[1:n.missingNames]
        }

        ## extract arguments
        endpoint[iE] <- gsub("\"","",iArg[iName=="endpoint"])
        if("threshold" %in% iName){
            thresholdTempo <- try(eval(expr = parse(text = iArg[iName=="threshold"])), silent = TRUE)
            if(inherits(thresholdTempo,"try-error")){
                stop(iArg[iName=="threshold"]," does not refer to a valid threshold \n",
                     "Should be numeric or the name of a variable in the global workspace \n")
            }
                
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
                            endpoint.TTE,
                            endpoint.UTTE,
                            censoring,
                            censoring.TTE,
                            censoring.UTTE,
                            D.TTE,
                            D.UTTE,
                            type,
                            strata,
                            threshold,
                            n.strata,
                            iidNuisance,
                            out){

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

    ## ** estimate cumulative incidence function (survival case or competing risk case)
    if(is.null(model.tte)){        
        model.tte <- vector(length = D.UTTE, mode = "list")
        names(model.tte) <- endpoint.UTTE

        txt.modelUTTE <- paste0("prodlim::Hist(",endpoint.UTTE,",",censoring.UTTE,") ~ ",treatment," + ..strata..")

        for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.UTTE <- 1
            model.tte[[iEndpoint.UTTE]] <- prodlim::prodlim(as.formula(txt.modelUTTE[iEndpoint.UTTE]),
                                                            data = data)
            model.tte[[iEndpoint.UTTE]]$XX <- model.tte[[iEndpoint.UTTE]]$X
        }
        
    }else{
        for(iEndpoint.UTTE in 1:D.UTTE){ ## iEndpoint.TTE <- 1
            ## convert treatment to numeric
            model.tte[[iEndpoint.UTTE]]$XX <- model.tte[[iEndpoint.UTTE]]$X
            model.tte[[iEndpoint.UTTE]]$XX[[treatment]] <- as.numeric(factor(model.tte[[iEndpoint.UTTE]]$XX[[treatment]], levels = level.treatment))-1
            p <- NCOL(model.tte[[iEndpoint.UTTE]]$XX)

            ## create ..strata..
            if(p==1){
                model.tte[[iEndpoint.UTTE]]$XX <- cbind(model.tte[[iEndpoint.UTTE]]$XX,
                                                        "..strata.." = 1)
            }else{
                col.strata <- setdiff(1:p,which(colnames(model.tte[[iEndpoint.UTTE]]$XX)==treatment))
                value.strata <- apply(model.tte[[iEndpoint.UTTE]]$XX[,col.strata],1,paste0,collapse="")
                model.tte[[iEndpoint.UTTE]]$XX <- cbind(model.tte[[iEndpoint.UTTE]]$XX,
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
        
        indexX.C <- which(model.tte[[iEndpoint.UTTE]]$XX[[treatment]]==0)
        indexX.T <- which(model.tte[[iEndpoint.UTTE]]$XX[[treatment]]==1)

        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            iNcontrol <- length(ls.indexC[[iStrata]])
            iNtreatment <- length(ls.indexT[[iStrata]])

            if("..strata.." %in% colnames(model.tte[[iEndpoint.UTTE]]$XX)){
                indexX.strata <- which(model.tte[[iEndpoint.UTTE]]$XX[["..strata.."]]==iStrata)
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
                ## jump times
                iIndexJumpC <- intersect(index.jump,iIndex.startC:iIndex.stopC)
                iIndexJumpT <- intersect(index.jump,iIndex.startT:iIndex.stopT)
                    
                iJumpC <- model.tte[[iEndpoint.UTTE]]$time[iIndexJumpC]
                iJumpT <- model.tte[[iEndpoint.UTTE]]$time[iIndexJumpT]

                ## last survival times
                iLast.survC <- model.tte[[iEndpoint.UTTE]]$surv[iIndex.stopC]
                iLast.survT <- model.tte[[iEndpoint.UTTE]]$surv[iIndex.stopT]

                ## survival at each jump
                iSurvTimeC <- c(-1e12,iJumpC)
                iSurvC <- c(1,model.tte[[iEndpoint.UTTE]]$surv[iIndexJumpC])
                if(iLast.survC!=0){ ## just after last event is unknown when the survival curve does not ends at 0
                    iSurvTimeC <- c(iSurvTimeC, model.tte[[iEndpoint.UTTE]]$time[iIndex.stopC] + 1e-12)
                    iSurvC <- c(iSurvC,NA)
                }

                iSurvTimeT <- c(-1e12,iJumpT)
                iSurvT <- c(1,model.tte[[iEndpoint.UTTE]]$surv[iIndexJumpT])
                if(iLast.survT!=0){ ## just after last event is unknown when the survival curve does not ends at 0
                    iSurvTimeT <- c(iSurvTimeT, model.tte[[iEndpoint.UTTE]]$time[iIndex.stopT] + 1e-12)
                    iSurvT <- c(iSurvT, NA)
                }

                ## dSurvival at each jump
                if(length(iJumpC)>0){
                    iIndexSurvivalC.JumpCm <- prodlim::sindex(iSurvTimeC, iJumpC - 1e-12)
                    iIndexSurvivalC.JumpCp <- prodlim::sindex(iSurvTimeC, iJumpC + 1e-12)
                    iDSurvC <- iSurvC[iIndexSurvivalC.JumpCp] - iSurvC[iIndexSurvivalC.JumpCm]
                }
                
                if(length(iJumpT)>0){
                    iIndexSurvivalT.JumpTm <- prodlim::sindex(iSurvTimeT, iJumpT - 1e-12)
                    iIndexSurvivalT.JumpTp <- prodlim::sindex(iSurvTimeT, iJumpT + 1e-12)
                    iDSurvT <- iSurvT[iIndexSurvivalT.JumpTp] - iSurvT[iIndexSurvivalT.JumpTm]
                }
                
                ## independent of the threshold i.e. of the priority
                ## avoid repeated calculation when the same endpoint is used several times with different thresholds
                iIndexSurvivalC.timeC <- prodlim::sindex(iSurvTimeC, iTimeC)
                iIndexSurvivalC.timeT <- prodlim::sindex(iSurvTimeC, iTimeT)
                iIndexSurvivalT.timeC <- prodlim::sindex(iSurvTimeT, iTimeC)
                iIndexSurvivalT.timeT <- prodlim::sindex(iSurvTimeT, iTimeT)
                
                iSurvivalC.timeC <- iSurvC[iIndexSurvivalC.timeC]
                iSurvivalC.timeT <- iSurvC[iIndexSurvivalC.timeT]
                iSurvivalT.timeC <- iSurvT[iIndexSurvivalT.timeC]
                iSurvivalT.timeT <- iSurvT[iIndexSurvivalT.timeT]
                
            for(iEndpoint in iIndex.associatedEndpoint){ ## iEndpoint <- 1
                iThreshold <- threshold[iEndpoint]

                ## **** last survival
                out$lastSurv[[iEndpoint]][iStrata,1:2] <- c(iLast.survC, iLast.survT)

                ## **** survival at jump times
                if(length(iJumpC)>0){                    
                    iIndexSurvivalT.JumpCpTau <- prodlim::sindex(iSurvTimeT, iJumpC + iThreshold)
                    out$survJumpC[[iEndpoint]][[iStrata]] <- cbind(time = iJumpC,
                                                                   survival = iSurvT[iIndexSurvivalT.JumpCpTau],
                                                                   dSurvival = iDSurvC)

                    if(iidNuisance){
                        out$survJumpC[[iEndpoint]][[iStrata]] <- cbind(out$survJumpC[[iEndpoint]][[iStrata]],
                                                                       index.survival = iIndexSurvivalT.JumpCpTau - 1,
                                                                       index.dSurvival1 = iIndexSurvivalC.JumpCm - 1,
                                                                       index.dSurvival2 = iIndexSurvivalC.JumpCp - 1)
                        ## iSurvT[iIndexSurvivalT.JumpCpTau]
                        ## iSurvC[iIndexSurvivalC.JumpCm]
                    }
                }else{
                    out$survJumpC[[iEndpoint]][[iStrata]] <- matrix(nrow = 0, ncol = 3,
                                                                    dimnames = list(NULL, c("time","surival","dSurvival")))
                }
                
                if(length(iJumpT)>0){                    
                    iIndexSurvivalC.JumpTpTau <- prodlim::sindex(iSurvTimeC, iJumpT + iThreshold)                
                    out$survJumpT[[iEndpoint]][[iStrata]] <- cbind(time = iJumpT,
                                                                   survival = iSurvC[iIndexSurvivalC.JumpTpTau],
                                                                   dSurvival = iDSurvT)
                    if(iidNuisance){
                        out$survJumpT[[iEndpoint]][[iStrata]] <- cbind(out$survJumpT[[iEndpoint]][[iStrata]],
                                                                       index.survival = iIndexSurvivalC.JumpTpTau - 1,
                                                                       index.dSurvival1 = iIndexSurvivalT.JumpTm - 1,
                                                                       index.dSurvival2 = iIndexSurvivalT.JumpTp - 1)
                    }
                }else{
                    out$survJumpT[[iEndpoint]][[iStrata]] <- matrix(nrow = 0, ncol = 3,
                                                                    dimnames = list(NULL, c("time","surival","dSurvival")))
                }

                ## **** survival at observation time (+/- threshold)
                iIndexSurvivalC.timeCmTau <- prodlim::sindex(iSurvTimeC, iTimeC - iThreshold)
                iIndexSurvivalC.timeCpTau <- prodlim::sindex(iSurvTimeC, iTimeC + iThreshold)
                iIndexSurvivalT.timeCmTau <- prodlim::sindex(iSurvTimeT, iTimeC - iThreshold)
                iIndexSurvivalT.timeCpTau <- prodlim::sindex(iSurvTimeT, iTimeC + iThreshold)

                out$survTimeC[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeC,
                                                               "SurvivalC-threshold" = iSurvC[iIndexSurvivalC.timeCmTau],
                                                               "SurvivalC_0" = iSurvivalC.timeC,
                                                               "SurvivalC+threshold" = iSurvC[iIndexSurvivalC.timeCpTau],
                                                               "SurvivalT-threshold" = iSurvT[iIndexSurvivalT.timeCmTau],
                                                               "SurvivalT_0" = iSurvivalT.timeC,
                                                               "SurvivalT+threshold" = iSurvT[iIndexSurvivalT.timeCpTau])
                if(iidNuisance){                    
                    out$survTimeC[[iEndpoint]][[iStrata]] <- cbind(out$survTimeC[[iEndpoint]][[iStrata]],
                                                                   "index.SurvivalC-threshold" = iIndexSurvivalC.timeCmTau - 1,
                                                                   "index.SurvivalC_0" = iIndexSurvivalC.timeC - 1,
                                                                   "index.SurvivalC+threshold" = iIndexSurvivalC.timeCpTau - 1,
                                                                   "index.SurvivalT-threshold" = iIndexSurvivalT.timeCmTau - 1,
                                                                   "index.SurvivalT_0" = iIndexSurvivalT.timeC - 1,
                                                                   "index.SurvivalT+threshold" = iIndexSurvivalT.timeCpTau - 1
                                                                   )
                }

                iIndexSurvivalC.timeTmTau <- prodlim::sindex(iSurvTimeC, iTimeT - iThreshold)
                iIndexSurvivalC.timeTpTau <- prodlim::sindex(iSurvTimeC, iTimeT + iThreshold)
                iIndexSurvivalT.timeTmTau <- prodlim::sindex(iSurvTimeT, iTimeT - iThreshold)
                iIndexSurvivalT.timeTpTau <- prodlim::sindex(iSurvTimeT, iTimeT + iThreshold)

                out$survTimeT[[iEndpoint]][[iStrata]] <- cbind("time" = iTimeT,
                                                               "SurvivalC-threshold" = iSurvC[iIndexSurvivalC.timeTmTau],
                                                               "SurvivalC_0" = iSurvivalC.timeT,
                                                               "SurvivalC+threshold" = iSurvC[iIndexSurvivalC.timeTpTau],
                                                               "SurvivalT-threshold" = iSurvT[iIndexSurvivalT.timeTmTau],
                                                               "SurvivalT_0" = iSurvivalT.timeT,
                                                               "SurvivalT+threshold" = iSurvT[iIndexSurvivalT.timeTpTau]
                                                               )
                if(iidNuisance){
                    out$survTimeT[[iEndpoint]][[iStrata]] <- cbind(out$survTimeT[[iEndpoint]][[iStrata]],
                                                                   "index.SurvivalC-threshold" = iIndexSurvivalC.timeTmTau - 1,
                                                                   "index.SurvivalC_0" = iIndexSurvivalC.timeT - 1,
                                                                   "index.SurvivalC+threshold" = iIndexSurvivalC.timeTpTau - 1,
                                                                   "index.SurvivalT-threshold" = iIndexSurvivalT.timeTmTau - 1,
                                                                   "index.SurvivalT_0" = iIndexSurvivalT.timeT - 1,
                                                                   "index.SurvivalT+threshold" = iIndexSurvivalT.timeTpTau - 1
                                                                   )
                }
            }

            }
            
            


        }
    }

    ## ** prepare influence function
    out$iid <- vector(mode = "list", length = 4)
    template <- lapply(1:D.UTTE, function(IE){
        lapply(1:n.strata, matrix, nrow = 0, ncol = 0)
        })
    out$iid <- lapply(out$iid, function(x){template})
    names(out$iid) <- c("survJumpC","dSurvJumpC","survJumpT","dSurvJumpT")

    if(iidNuisance){
        iid.model.tte <- lapply(model.tte, function(iModel){ ## iModel <- model.tte[[1]]
            iOut <- iidProdlim(iModel, add0 = TRUE)
            iOut$IFsurvival.control <- iOut$IFsurvival[which(iOut$X[,treatment]==0)]
            iOut$IFsurvival.treatment <- iOut$IFsurvival[which(iOut$X[,treatment]==1)]
            return(iOut)
        })
        for(iEndpoint in 1:D.UTTE){  ## iEndpoint <- 1
            iEndpoint.UTTE.name <- endpoint.UTTE[iEndpoint.UTTE]
            iIndex.associatedEndpoint <- which(endpoint == iEndpoint.UTTE.name)

            for(iStrata in 1:n.strata){  ## iStrata <- 1
                iIID.control <- iid.model.tte[[iEndpoint]]$IFsurvival.control[[iStrata]]
                iIID.treatment <- iid.model.tte[[iEndpoint]]$IFsurvival.treatment[[iStrata]]

                ## iid.model.tte[[iEndpoint]]$time
                out$iid$survJumpC[[iEndpoint]][[iStrata]] <- iIID.control
                if(NCOL(iIID.control)>1){
                    out$iid$dSurvJumpC[[iEndpoint]][[iStrata]] <- iIID.control - cbind(0,iIID.control[,1:(NCOL(iIID.control)-1),drop=FALSE])
                }else{
                    out$iid$dSurvJumpC[[iEndpoint]][[iStrata]] <- iIID.control
                }
                
                out$iid$survJumpT[[iEndpoint]][[iStrata]] <- iIID.treatment
                if(NCOL(iIID.treatment)>1){
                    out$iid$dSurvJumpT[[iEndpoint]][[iStrata]] <-  iIID.treatment - cbind(0,iIID.treatment[,1:(NCOL(iIID.treatment)-1),drop=FALSE])
                }else{
                    out$iid$dSurvJumpT[[iEndpoint]][[iStrata]] <-  iIID.treatment
                }
                out$p.C[iStrata, iIndex.associatedEndpoint] <- NCOL(iIID.control)
                out$p.T[iStrata, iIndex.associatedEndpoint] <- NCOL(iIID.treatment)
            }
        }
        
    }
    
    ## ** export
    return(out)
    
}



