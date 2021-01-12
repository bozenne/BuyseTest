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
#' \item status: only keep status relative to TTE endpoint. Set to \code{NULL} if no TTE endpoint.
#' \item threshold: set default threshold to 1e-12 expect for binary variable where it is set to 1/2.
#' the rational being we consider a pair favorable if X>Y ie X>=Y+1e-12.
#' When using a threshold e.g. 5 we want X>=Y+5 and not X>Y+5, especially when the measurement is discrete. \cr
#' \item data: convert to data.table object.
#' \item scoring.rule: convert to numeric.
#' }
#'
#' \code{initializeFormula}:  extract \code{treatment}, \code{type}, \code{endpoint}, \code{threshold}, \code{status}, \code{operator}, and \code{strata}
#' from the formula. \cr \cr
#'
#' \code{initializeData}: Divide the dataset into two, one relative to the treatment group and the other relative to the control group.
#' Merge the strata into one with the interaction variable.
#' Extract for each strata the index of the observations within each group.
#'
#' @keywords function internal BuyseTest
#' @author Brice Ozenne

## * initializeArgs
#' @rdname internal-initialization
initializeArgs <- function(status,
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
                           operator = NULL,
                           censoring,
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
    fitter.model.tte <- option$fitter.model.tte
    engine <- option$engine
    alternative <- option$alternative
    precompute <- option$precompute
    
    ## ** convert formula into separate arguments
    if(!missing(formula)){
        ## the missing is for BuysePower where the arguments are not necessarily specified
        test.null <- c(status = !missing(status) && !is.null(status),
                       endpoint = !missing(endpoint) && !is.null(endpoint),
                       operator = !missing(operator) && !is.null(operator),
                       censoring = !missing(censoring) && !is.null(censoring),
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
        status <- resFormula$status
        weight <- resFormula$weight
        operator <- resFormula$operator
        censoring <- resFormula$censoring
        strata <- resFormula$strata
    }else{
        formula <- NULL
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
    endpoint.TTE <- endpoint[index.type3]
    threshold.TTE <- threshold[index.type3]

    D <- length(endpoint)
    D.TTE <- length(endpoint.TTE)
    
    Uendpoint <- unique(endpoint) 
    Uendpoint.TTE <- unique(endpoint.TTE) 
    
    ## ** default values 
    if(is.null(formula)){

        if(is.null(operator)){
            operator <- rep(">0",D)
        }
        if(is.null(weight)){
            weight <- rep(1,D)
        }
        if(is.null(status)){
            status <- rep("..NA..",D)
        }else if(length(status) != D && length(status) == D.TTE){
            status.save <- status
            status <- rep("..NA..", D)
            status[index.type3] <- status.save             
        }
        
        if(is.null(censoring)){
            censoring <- rep("right",D)
        }else if(length(status) != D && length(status) == D.TTE){
            censoring.save <- status
            censoring <- rep("right", D)
            censoring[index.type3] <- status.save             
        }
    }

    ## ** status
    Ustatus <- unique(status)
    status.TTE <- status[index.type3]
    ## from now, status contains for each endpoint the name of variable indicating status (0) or event (1) or NA

    ## ** censoring
    ## if(any(type %in% 1:2)){
    ##     censoring[type %in% 1:2] <- as.character(NA)
    ## }
    censoring.save <- censoring
    censoring <- sapply(unname(censoring),function(iC){
        if(identical(iC,"NA")){
            return(0)
        }else if(identical(iC,"right")){
            return(1)
        }else if(identical(iC,"left")){
            return(2)
        }else{
            return(NA)
        }
    })
    attr(censoring,"original") <- censoring.save
    
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
    if(is.na(strata.resampling) || length(strata.resampling)== 0){
        attr(method.inference,"resampling-strata") <- as.character(NA)
    }else{
        attr(method.inference,"resampling-strata") <- strata.resampling
    }

    ## ** neutral.as.uninf
    if(length(neutral.as.uninf)==1 && D>1){
        neutral.as.uninf <- rep(neutral.as.uninf,D)
    }
    
    ## ** correction.uninf
    correction.uninf <- as.numeric(correction.uninf)
    if(correction.uninf>0){
        engine <- "GPC_cpp"
    }
    
    ## ** model.tte
    if(identical(scoring.rule,1)){
        if((!is.null(model.tte)) && (length(unique(endpoint.TTE)) == 1) && !inherits(model.tte, "list")){
            attr.save <- attr(model.tte,"iidNuisance")
            
            model.tte <- list(model.tte)
            names(model.tte) <- unique(endpoint.TTE)
            attr(model.tte,"iidNuisance") <- attr.save
        }
    }else{
        model.tte <- NULL
    }
    if(!is.null(model.tte)){
        fitter.model.tte <- unlist(lapply(model.tte, class))
    }else{
        fitter.model.tte <- setNames(rep(fitter.model.tte, length(Uendpoint.TTE)), Uendpoint.TTE)
    }

    ## ** iid
    iid <- attr(method.inference,"studentized") || (method.inference == "u-statistic")
    if(iid){
        attr(method.inference,"hprojection") <- option$order.Hprojection
    }else{
        attr(method.inference,"hprojection") <- NA
    }
    iidNuisance <- iid && identical(scoring.rule,1) && (is.null(model.tte) || identical(attr(model.tte,"iidNuisance"),TRUE))

    ## ** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }

    ## ** trace
    if(is.logical(trace)){
        trace <- as.numeric(trace)
    }

    ## ** operator
    operator <- sapply(operator, switch, ">0"=1, "<0"=-1, NA)

    ## ** export
    return(list(
        name.call = name.call,
        status = status,
        status.TTE = status.TTE,
        correction.uninf = correction.uninf,
        cpus = cpus,
        D = D,
        D.TTE = D.TTE,
        data = data,
        endpoint = endpoint,
        endpoint.TTE = endpoint.TTE,
        engine = engine,
        fitter.model.tte = fitter.model.tte,
        formula = formula,
        iid = iid,
        iidNuisance = iidNuisance,
        index.endpoint = match(endpoint, Uendpoint) - 1,
        index.status = match(status, Ustatus) - 1,
        keep.pairScore = keep.pairScore,
        keep.survival = option$keep.survival,
        scoring.rule = scoring.rule,
        model.tte = model.tte,
        method.inference = method.inference,
        n.resampling = n.resampling,
        hierarchical = hierarchical,
        neutral.as.uninf = neutral.as.uninf,
        operator = operator,
        censoring = censoring,
        order.Hprojection = option$order.Hprojection,
        precompute = precompute,
        seed = seed,
        strata = strata,
        threshold = threshold,
        trace = trace,
        treatment = treatment,
        type = type,
        Uendpoint = Uendpoint,
        Ustatus = Ustatus,
        weight = weight,
        debug = option$debug
    ))
}

## * initializeData
#' @rdname internal-initialization
initializeData <- function(data, type, endpoint, Uendpoint, D, scoring.rule, status, Ustatus, method.inference, censoring, strata, treatment, hierarchical, copy,
                           keep.pairScore, endpoint.TTE, status.TTE, iidNuisance){

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
            endpoint.char <- setdiff(names(data.class)[test.num==FALSE],treatment)
            for(iE in endpoint.char){
                data[, c(iE) := as.double(as.factor(.SD[[1]]))-1.0, .SDcols = iE]
            }
        }
    }


    ## ** n.obs
    n.obs <- data[,.N]

    ## ** strata
    if(!is.null(strata)){  
    
        data[ , c("..strata..") := interaction(.SD, drop = TRUE, lex.order = FALSE, sep = "."), .SDcols = strata]
        level.strata <- levels(data[["..strata.."]])        
        data[ , c("..strata..") := as.numeric(.SD[["..strata.."]])] # convert to numeric
        
        n.obsStrata <- data[,.N, by = "..strata.."][,stats::setNames(.SD[[1]],.SD[[2]]),.SD = c("N","..strata..")]
    }else{
        
        data[ , c("..strata..") := 1]
        n.obsStrata <- n.obs
        level.strata <- 1
    }

    n.strata <- length(level.strata)

    ## ** convert treatment to binary indicator
    level.treatment <- levels(as.factor(data[[treatment]]))
    trt2bin <- stats::setNames(0:1,level.treatment)
    data[ , c(treatment) := trt2bin[as.character(.SD[[1]])], .SDcols = treatment]

    ## ** rowIndex
    data[,c("..rowIndex..") := 1:.N]

    ## ** unique status
    if(any(status == "..NA..")){
        data[,c("..NA..") := -100]
    }

    ## ** TTE with status
    if(scoring.rule>0){
        test.status <- sapply(status.TTE, function(iC){any(data[[iC]]==0)})
        if(all(test.status==FALSE)){
            scoring.rule <- 0
            iidNuisance <- FALSE            
        }else if(identical(attr(method.inference,"hprojection"),2)){
            keep.pairScore <- TRUE ## need the detail of the score to perform the 2nd order projection
        }
        
        ## distinct time to event endpoints
        endpoint.UTTE <- unique(endpoint.TTE[test.status])
        status.UTTE <- unique(status.TTE[test.status])
        D.UTTE <- length(endpoint.UTTE)

        ## correspondance endpoint, TTE endpoint (non TTEe endpoint are set to -100)
        index.UTTE <- match(endpoint, endpoint.UTTE, nomatch = -99) - 1
    }else{
        endpoint.UTTE <- numeric(0)
        status.UTTE <- numeric(0)
        D.UTTE <- 0
        index.UTTE <- rep(-100, D)
    }
    
    ## ** scoring method for each endpoint
    ## check if status
    n.CR <- sapply(status, function(iC){max(data[[iC]])})
    test.CR <- n.CR[status]>1
    test.censoring <- sapply(Ustatus, function(iC){any(data[[iC]]==0)})[status]

    method.score <- sapply(1:D, function(iE){ ## iE <- 1
        if(type[iE] %in% 1:2 || (test.censoring[iE]==FALSE && test.CR[iE]==FALSE)){
            return(1) ## 1 binary/continuous
        }else if(scoring.rule == 0){ ## 2/3 Gehan (right/left censoring)
            return(1 + censoring[iE])
        }else{
            return(4 + test.CR[iE]) ## 4/5 Peron (survival/competing risks)
        }
    })
    attr(method.score,"test.censoring") <- test.censoring
    attr(method.score,"test.CR") <- test.CR
    
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
    if(!is.na(attr(method.inference,"resampling-strata"))){
        n.obsStrataResampling <- table(data[,interaction(.SD), .SDcols = attr(method.inference,"resampling-strata")])
    }else{
        n.obsStrataResampling <- n.obs
    }
    
    ## ** skeleton for survival proba (only relevant for Peron scoring rule)
    skeletonPeron <- list(survTimeC = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survTimeT = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survJumpC = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          survJumpT = lapply(1:D, function(iE){lapply(1:n.strata, function(iS){matrix(nrow=0,ncol=0)})}),
                          lastSurv = lapply(1:D, function(iS){matrix(nrow = n.strata, ncol = 2*max(1,n.CR[iS]))}), ## 4 for competing risk setting, 2 is enough for survival
                          p.C = matrix(-100, nrow = n.strata, ncol = D),
                          p.T = matrix(-100, nrow = n.strata, ncol = D),
                          iid = list(survJumpC = lapply(1:D.UTTE, function(IE){lapply(1:n.strata, matrix, nrow = 0, ncol = 0)}),
                                     survJumpT = lapply(1:D.UTTE, function(IE){lapply(1:n.strata, matrix, nrow = 0, ncol = 0)})
                                     )
                          )

    ## ** export
    keep.cols <- union(c(treatment, "..strata.."),
                       na.omit(attr(method.inference,"resampling-strata")))

    return(list(data = data[,.SD,.SDcols = keep.cols],
                M.endpoint = as.matrix(data[, .SD, .SDcols = Uendpoint]),
                M.status = as.matrix(data[, .SD, .SDcols = Ustatus]),
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
                cumn.obsStrataResampling = c(0,cumsum(n.obsStrataResampling)),
                skeletonPeron = skeletonPeron,
                scoring.rule = scoring.rule,
                iidNuisance = iidNuisance,
                nUTTE.analyzedPeron_M1 = nUTTE.analyzedPeron_M1,
                endpoint.UTTE = endpoint.UTTE,
                status.UTTE = status.UTTE,
                D.UTTE = D.UTTE,
                index.UTTE = index.UTTE,
                keep.pairScore = keep.pairScore
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
             "nothing of the form type(endpoint,threshold,status) found in the formula \n")
    }

    ## ** extract endpoints and additional arguments 
    threshold <- rep(NA, n.endpoint)
    status <- rep("..NA..", n.endpoint)
    endpoint <- rep(NA, n.endpoint)
    operator <- rep(">0", n.endpoint)
    censoring <- rep("right", n.endpoint)
    weight <- rep(1, n.endpoint)
    validArgs <- c("endpoint","status","threshold","operator","weight","censoring")

    ## split around parentheses
    ls.x.endpoint <- strsplit(vec.x.endpoint, split = "(", fixed = TRUE)

    type <- character(length = n.endpoint)
    for(iE in 1:n.endpoint){
        ## extract type
        type[iE] <- tolower(ls.x.endpoint[[iE]][1])
        if(type[iE] %in% c("b","bin","binary")){
            iValidArgs <- setdiff(validArgs,c("status","threshold"))
        }else if(type[iE] %in% c("c","cont","continuous")){
            iValidArgs <- setdiff(validArgs,"status")
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
                 x[iE]," has too many arguments (maximum 4: endpoint, threshold, status variable, operator) \n")
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
        if("status" %in% iName){
            status[iE] <- gsub("\"","",iArg[iName=="status"])
        }
        if("operator" %in% iName){
            operator[iE] <- gsub("\"","",iArg[iName=="operator"])
        }
        if("weight" %in% iName){
            weight[iE] <- as.numeric(eval(expr = parse(text = iArg[iName=="weight"])))
        }
        if("censoring" %in% iName){
            censoring[iE] <- gsub("\"","",iArg[iName=="censoring"])
        }
    }

    ## ** export
    return(list(treatment = treatment,
                type = type,
                endpoint = endpoint,
                threshold = threshold,
                status = status,
                operator = operator,
                weight = weight,
                censoring = censoring,
                strata = strata))
}




