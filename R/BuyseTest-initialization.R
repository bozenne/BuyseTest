## * Documentation initialization functions called by BuyseTest

#' @title internal functions for BuyseTest - initialization
#' @name BuyseTest-initialization
#' @description Functions called by \code{\link{BuyseTest}} to initialize the arguments.
#' @noRd
#' 
#' @details
#' 
#' \code{initializeArgs}: Normalize the argument 
#' \itemize{
#' \item scoring.rule, pool.strata, neutral.as.uninf, add.halfNeutral, keep.pairScore, n.resampling, seed, cpus, trace: set to default value when not specified.
#' \item formula: call \code{initializeFormula} to extract arguments.
#' \item type: convert to numeric.
#' \item status: only keep status relative to TTE endpoint. Set to \code{NULL} if no TTE endpoint.
#' \item threshold: set default threshold to 1e-12.
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
#' @author Brice Ozenne

## * initializeArgs
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
                           pool.strata = NULL,
                           model.tte,
                           n.resampling = NULL,
                           strata.resampling = NULL,
                           call,
                           neutral.as.uninf = NULL,
                           add.halfNeutral = NULL,
                           operator = NULL,
                           censoring,
                           restriction,
                           option,
                           seed = NULL,
                           strata,
                           threshold,
                           trace = NULL,
                           treatment,
                           type,
                           weightEndpoint = NULL,
                           weightObs = NULL,
                           envir){

    name.call <- names(call)

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
    if(is.null(add.halfNeutral)){ add.halfNeutral <- option$add.halfNeutral }
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
                       restriction = !missing(restriction) && !is.null(restriction),
                       strata = !missing(strata) && !is.null(strata),
                       threshold = !missing(threshold) && !is.null(threshold),
                       treatment = !missing(treatment) && !is.null(treatment),
                       type = !missing(type) && !is.null(type),
                       weightEndpoint = !missing(weightEndpoint) && !is.null(weightEndpoint)
                       )
        if(any(test.null)){
            txt <- names(test.null)[test.null]
            warning("Argument",if(sum(test.null)>1){"s"}," \'",paste(txt, collpase="\' \'"),if(sum(test.null)>1){" are "}else{" is "}," ignored when argument \'formula\' has been specified\n")
        }
        
        resFormula <- initializeFormula(formula, hierarchical = hierarchical, envir = envir)

        treatment <- resFormula$treatment
        type <- resFormula$type
        endpoint <- resFormula$endpoint
        threshold <- resFormula$threshold
        status <- resFormula$status
        weightEndpoint <- resFormula$weightEndpoint
        operator <- resFormula$operator
        censoring <- resFormula$censoring
        restriction <- resFormula$restriction
        strata <- resFormula$strata
        if(!is.null(strata)){
            attr(strata,"match") <- resFormula$match
        }
    }else{
        formula <- NULL
    }

    ## ** type
    validType1 <- paste0("^",c("b","bin","binary"),"$")
    validType2 <- paste0("^",c("c","cont","continuous"),"$")
    validType3 <- paste0("^",c("t","tte","time","timetoevent"),"$") ## [if modified, remember to change the corresponding vector in initFormula]
    validType4 <- paste0("^",c("g","gaus","gaussian"),"$") ## [if modified, remember to change the corresponding vector in initFormula]
    type <- tolower(type)

    type[grep(paste(validType1,collapse="|"), type)] <- "bin" 
    type[grep(paste(validType2,collapse="|"), type)] <- "cont" 
    type[grep(paste(validType3,collapse="|"), type)] <- "tte"
    type[grep(paste(validType4,collapse="|"), type)] <- "gaus"

    ## ** endpoint
    index.typeTTE <- which(type=="tte")    
    endpoint.TTE <- endpoint[index.typeTTE]
    threshold.TTE <- threshold[index.typeTTE]

    D <- length(endpoint)
    D.TTE <- length(endpoint.TTE)
    
    Uendpoint <- unique(endpoint) 
    Uendpoint.TTE <- unique(endpoint.TTE) 
    
    ## ** default values
    if(is.null(formula)){
        if(is.null(threshold)){
            threshold <- rep(10^{-12},D)  # if no treshold is proposed all threshold are by default set to 10^{-12}
            attr(threshold, "original") <- ifelse(type=="bin",as.numeric(NA),0)
        }else{
            attr(threshold, "original") <- threshold
        }

        if(is.null(restriction)){
            restriction <- rep(as.numeric(NA),D)  
        }

        if(is.null(operator)){
            operator <- rep(">0",D)
        }
        if(is.null(weightEndpoint)){
            if(hierarchical){
                weightEndpoint <- rep(1,D)
            }else{
                weightEndpoint <- rep(1/D,D)
            }
        }
        if(is.null(status)){
            status <- rep("..NA..",D)
        }else if(length(status) != D && length(status) == D.TTE){
            status.save <- status
            status <- rep("..NA..", D)
            status[index.typeTTE] <- status.save             
        }
        
        if(is.null(censoring)){
            censoring <- rep("right",D)
        }else if(length(status) != D && length(status) == D.TTE){
            censoring.save <- status
            censoring <- rep("right", D)
            censoring[index.typeTTE] <- status.save             
        }
    }else{
        attr(threshold, "original") <- threshold
    }

    ## ** status
    Ustatus <- unique(status)
    status.TTE <- status[index.typeTTE]
    ## from now, status contains for each endpoint the name of variable indicating status (0) or event (1) or NA

    ## ** censoring
    ## ## if(any(type %in% 1:2)){
    ## ##     censoring[type %in% 1:2] <- as.character(NA)
    ## ## }
    ## if(!is.numeric(censoring)){
    ##     censoring.save <- censoring
    ##     censoring <- sapply(unname(censoring),function(iC){
    ##         if(identical(iC,"NA")){
    ##             return(0)
    ##         }else if(identical(iC,"right")){
    ##             return(1)
    ##         }else if(identical(iC,"left")){
    ##             return(2)
    ##         }else{
    ##             return(NA)
    ##         }
    ##     })
    ##     attr(censoring,"original") <- censoring.save
    ## }

    ## ** scoring.rule
    ## WARNING: choices must be lower cases
    ##          remember to update check scoring.rule (in BuyseTest-check.R)
    if(is.character(scoring.rule)){
        scoring.rule <- switch(tolower(scoring.rule),
                               "gehan" = 0,
                               "peron" = 1,
                               "efron" = 2,
                               NA
                               )
    }

    if (D.TTE == 0) {
        scoring.rule <- 0
        if ("scoring.rule" %in% name.call && trace > 0) {
            message("NOTE : there is no survival endpoint, \'scoring.rule\' argument is ignored \n")
        }
    }

    ## ** strata
    if(!is.null(strata) && is.null(attr(strata,"match"))){
        attr(strata,"match") <- FALSE
    }
    
    ## ** pool.strata
    if(is.null(strata)){
        pool.strata <- 0
        attr(pool.strata,"type") <- "none"
        attr(pool.strata,"original") <- NA
    }else if(is.null(pool.strata)){
        pool.strata <- switch(tolower(option$pool.strata),
                              "buyse" = 0,
                              "cmh" = 1,
                              "equal" = 2,
                              "standardisation" = 3,
                              "standardization" = 3,
                              "var-favorable" = 4.1,
                              "var-unfavorable" = 4.2,
                              "var-netbenefit" = 4.3,
                              "var-winratio" = 4.4,
                              NA
                              )
        if(tolower(option$pool.strata)=="standardisation"){
            attr(pool.strata,"type") <- "standardization"
        }else{
            attr(pool.strata,"type") <- tolower(option$pool.strata)
        }
        attr(pool.strata,"original") <- NA   
        
    }else if(is.character(pool.strata)){
        pool.strata_save <- tolower(pool.strata)
        pool.strata <- switch(pool.strata_save,
                              "buyse" = 0,
                              "cmh" = 1,
                              "equal" = 2,
                              "standardisation" = 3,
                              "standardization" = 3,
                              "var-favorable" = 4.1,
                              "var-unfavorable" = 4.2,
                              "var-netbenefit" = 4.3,
                              "var-winratio" = 4.4,
                              NA
                              )
        if(!is.na(pool.strata_save) && pool.strata_save=="standardisation"){
            attr(pool.strata,"type") <- "standardization"
        }else{
            attr(pool.strata,"type") <- pool.strata_save
        }
        attr(pool.strata,"original") <- pool.strata_save
    }else if(is.numeric(pool.strata)){
        pool.strata_save <- switch(as.character(pool.strata),
                                   "0" = "buyse",
                                   "1" = "cmh",
                                   "2" = "equal",
                                   "3" = "standardization",
                                   "4.1" = "var-favorable",
                                   "4.2" = "var-unfavorable",
                                   "4.3" = "var-netbenefit",
                                   "4.4" = "var-winratio",
                                   NA
                                   )
        attr(pool.strata,"type") <- pool.strata_save
        attr(pool.strata,"original") <- pool.strata_save
    }else{
        pool.strata <- NA
    }
    
    ## ** threshold
    if(any(is.na(threshold))){
        threshold[which(is.na(threshold))] <- 10^{-12}
    }
    if(any(abs(threshold)<10^{-12})){
        threshold[which(abs(threshold)<10^{-12})] <- 10^{-12}
    }

    ## ** method.inference
    method.inference <- gsub("-"," ",tolower(method.inference),fixed = TRUE)
    attr(method.inference,"permutation") <- grepl("permutation",method.inference)
    attr(method.inference,"bootstrap") <- grepl("bootstrap",method.inference)
    attr(method.inference,"studentized") <- grepl("studentized",method.inference)
    attr(method.inference,"ustatistic") <- grepl("u statistic",method.inference)
    if(all(is.na(strata.resampling)) || length(strata.resampling)== 0){
        attr(method.inference,"resampling-strata") <- as.character(NA)
    }else{
        attr(method.inference,"resampling-strata") <- strata.resampling
    }
    if(method.inference == "varexact permutation"){
        n.resampling <- Inf
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
    if(scoring.rule>0){
        if((!is.null(model.tte))){
            if((length(unique(endpoint.TTE)) == 1) && !inherits(model.tte, "list")){
                attr.save <- attr(model.tte,"iidNuisance")
            
                model.tte <- list(model.tte)
                names(model.tte) <- unique(endpoint.TTE)
                attr(model.tte,"iidNuisance") <- attr.save
            }
            attr(data,"model.tte_regressor") <- unique(unlist(lapply(model.tte, function(iM){
                all.vars(stats::delete.response(terms(formula(iM))))
            })))
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
    iid <- attr(method.inference,"studentized") || (method.inference == "u statistic")
    if(iid){
        attr(method.inference,"hprojection") <- option$order.Hprojection
    }else{
        attr(method.inference,"hprojection") <- NA
    }
    if(iid && scoring.rule>0){ ## Peron/Efron scoring rule
        if(is.null(model.tte)){
            iidNuisance <- TRUE
        }else if(!is.null(attr(model.tte,"iidNuisance"))){
            iidNuisance <- attr(model.tte,"iidNuisance")
        }else if(all(deparse(expr = call$data) == sapply(model.tte, function(iModel){deparse(expr=iModel$call$data)}))){
            iidNuisance <- TRUE
        }else{
            iidNuisance <- FALSE
            message("Uncertainty related to the estimation of the survival probabilities is ignored. \n",
                    "Consider adding an attribute \"iidNuisance\" to the argument \'model.tte\' taking value TRUE to change this default behavior. \n")
        }
    }else{
        iidNuisance <- FALSE
    }
    

    ## ** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }

    ## ** trace
    if(is.logical(trace)){
        trace <- as.numeric(trace)
    }

    ## ** seed
    if(!is.null(seed)){
        attr(seed,"max") <- 10^(floor(log10(.Machine$integer.max))-1)        
    }

    ## ** operator
    if(!is.numeric(operator)){
        operator <- sapply(operator, switch, ">0"=1, "<0"=-1, NA)
    }

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
        pool.strata = pool.strata,
        model.tte = model.tte,
        method.inference = method.inference,
        n.resampling = n.resampling,
        hierarchical = hierarchical,
        neutral.as.uninf = neutral.as.uninf,
        add.halfNeutral = add.halfNeutral,
        operator = operator,
        censoring = censoring,
        restriction = restriction,
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
        weightEndpoint = weightEndpoint,
        weightObs = weightObs,
        debug = option$debug
    ))
}

## * initializeData
initializeData <- function(data, type, endpoint, Uendpoint, D, scoring.rule, status, Ustatus, method.inference, censoring, strata, pool.strata, treatment, hierarchical, copy,
                           keep.pairScore, endpoint.TTE, status.TTE, iidNuisance, weightEndpoint, weightObs){

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }else if(copy){
        data <- data.table::copy(data)
    }

    ## ** convert character/factor to numeric for binary endpoints
    name.bin <- endpoint[which(type == "bin")]
    if(length(name.bin)>0){
        data.class <- sapply(data,class)
        test.num <- (data.class %in% c("numeric","integer"))
        if(any(test.num==FALSE)){
            endpoint.char <- setdiff(names(data.class)[test.num==FALSE],c(treatment,strata))
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

    nlevel.strata <- length(level.strata)
    if(pool.strata==3){
        grid.strata <- as.matrix(expand.grid(0:(nlevel.strata-1), 0:(nlevel.strata-1)))
        rownames(grid.strata) <- ifelse(level.strata[grid.strata[,1]+1]==level.strata[grid.strata[,2]+1],
                                        level.strata[grid.strata[,1]+1],
                                        paste(level.strata[grid.strata[,1]+1],level.strata[grid.strata[,2]+1],sep="."))
    }else{
        grid.strata <- cbind(0:(nlevel.strata-1), 0:(nlevel.strata-1))
        rownames(grid.strata) <- level.strata
    }
    n.strata <- NROW(grid.strata)
    
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
        if(type[iE] %in% c("bin","cont")){
            return("continuous") 
        }else if(type[iE] == "gaus"){
            return("gaussian")
        }else if(type[iE] == "tte"){
            if(test.censoring[iE]==FALSE && test.CR[iE]==FALSE){
                return("continuous")
            }else if(scoring.rule == 0){ ## 3/4 Gehan (right/left censoring)
                return(switch(censoring[iE],
                              "left" = "TTEgehan2",
                              "right" = "TTEgehan"))
            }else if(scoring.rule>0){
                return(switch(as.character(test.CR[iE]),
                              "FALSE" = "SurvPeron",
                              "TRUE" = "CRPeron"))
            }
        }
    })
    attr(method.score,"test.censoring") <- test.censoring
    attr(method.score,"test.CR") <- test.CR
    
    ## ** previously analyzed distinct TTE endpoints
    if(scoring.rule>0 && hierarchical){ ## only relevant when using Peron scoring rule with hierarchical GPC
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
    if(any(!is.na(attr(method.inference,"resampling-strata")))){
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

    ## ** iid for gaussian endpoints
    n.endpoint <- length(endpoint)
    index.gaussiid <- which(type == "gaus")
    if(length(index.gaussiid)>0 && any(!is.na(censoring[index.gaussiid]))){
        index.gaussiid2 <- intersect(which(!is.na(censoring)),index.gaussiid)
        for(iE in index.gaussiid2){ ## iE <- 1
            skeletonPeron$survTimeC[[iE]] <- data[index.C,list(list(do.call(cbind,.SD[[1]]))), .SDcols = censoring[iE], by = "..strata.."][[2]]
            skeletonPeron$survTimeT[[iE]] <- data[index.T,list(list(do.call(cbind,.SD[[1]]))), .SDcols = censoring[iE], by = "..strata.."][[2]]
        }
    }

    ## ** weightEndpoint
    if(missing(weightObs) || is.null(weightObs)){
        data$..weight.. <- 1
        weightObs <- data$..weight..
    }else if( (length(weightObs)==1) && (weightObs %in% names(data)) ){
        names(data)[names(data)==weightObs] <- "..weight.."
        weightObs <- data$..weight..
    }
    
    ## ** pool.strata
    ## set default pool.strata to Buyse for match data
    ## otherwise pooling will do something strange
    if(!is.null(strata) && attr(strata,"match") && pool.strata !=0){
        if(is.na(attr(pool.strata,"original"))){
            pool.strata[] <- 0
        }else{
            warning("Weights from the \"buyse\" pooling scheme (argument \'pool.strata\') are recommended for matched data. \n")
        }
    }

    ## ** keep.pairScore
    if(identical(attr(method.inference,"hprojection"),2) && scoring.rule>0){
        ## need the detail of the score to perform the 2nd order projection
        keep.pairScore <- TRUE 
    }else if(identical(attr(method.inference,"hprojection"),2) && pool.strata == 3){
        table2x2.nobs <- table(strata = data[["..strata.."]],data[[treatment]])
        dfStrata2 <- expand.grid(strata0 = 1:length(level.strata),
                                 strata1 = 1:length(level.strata))
        dfStrata2$N <- NROW(data)
        dfStrata2$m <- length(index.C)
        dfStrata2$n <- length(index.T)
        dfStrata2$m.X0 <- table2x2.nobs[dfStrata2$strata0,1]
        dfStrata2$n.X1 <- table2x2.nobs[dfStrata2$strata1,2]
        dfStrata2$X0 <- table(data[["..strata.."]])[dfStrata2$strata0]
        dfStrata2$X1 <- table(data[["..strata.."]])[dfStrata2$strata1]
        dfStrata2$weight <- (dfStrata2$X0*dfStrata2$X1/dfStrata2$N^2)/(dfStrata2$m.X0*dfStrata2$n.X1/(dfStrata2$m*dfStrata2$n))
        if(any(dfStrata2$weight %% 1 != 0)){
            ## if decimal weigths need the detail of the score to perform the 2nd order projection
            keep.pairScore <- TRUE
        }
    }
    
    ## ** export
    keep.cols <- union(union(c(treatment, strata, "..strata.."),
                             na.omit(attr(method.inference,"resampling-strata"))),
                       attr(data,"model.tte_regressor")) ## add regressor from survival models in case they do not match GPC strata variable (user-specific survival)

    return(list(data = data[,.SD,.SDcols = keep.cols],
                M.endpoint = as.matrix(data[, .SD, .SDcols = Uendpoint]),
                M.status = as.matrix(data[, .SD, .SDcols = Ustatus]),
                index.C = index.C,
                index.T = index.T,
                weightObs = weightObs,
                index.strata = tapply(data[["..rowIndex.."]], data[["..strata.."]], list),
                level.treatment = level.treatment,
                level.strata = level.strata, pool.strata = pool.strata, ## distinct strata levels (e.g. M, F)
                method.score = method.score, 
                grid.strata = grid.strata, ## strata (e.g. M, F, M.F, F.M) - different from level.strata when using standardisation
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
initializeFormula <- function(x, hierarchical, envir){

    ## ** check format
    validClass(x, valid.class = "formula")
    if(length(as.character(x))!=3){
        stop("Argument \'formula\' has unexpected length, as.character(x) should have length 3\n",
             "  length founded: ",length(as.character(x)),"\n")
    }

    ## ** categorize elements in the formula
    x.var <- all.vars(x)
    x.rhs <- stats::delete.response(terms(x))
    x.lhs <- terms(stats::as.formula(paste("~", deparse(x[[2]]))))

    config.rhs <- data.frame(label = attr(x.rhs,"term.labels"), var = as.character(NA), operator = as.character(NA), arguments = as.character(NA))
    config.rhs[config.rhs$label %in% x.var,"var"] <- config.rhs$label[config.rhs$label %in% x.var]
    split.rhs <- strsplit(config.rhs[config.rhs$label %in% x.var == FALSE,"label"], split = "(", fixed = TRUE)
    config.rhs[config.rhs$label %in% x.var == FALSE,"operator"] <- tolower(sapply(split.rhs,"[",1))
    config.rhs[config.rhs$label %in% x.var == FALSE,"arguments"] <- gsub(pattern = "\\)$", replacement = "", sapply(split.rhs,"[",2))
    
    config.lhs <- data.frame(label = attr(x.lhs,"term.labels"), var = as.character(NA), operator = as.character(NA), arguments = as.character(NA))
    config.lhs[config.lhs$label %in% x.var,"var"] <- config.lhs$label[config.lhs$label %in% x.var]
    split.lhs <- strsplit(config.lhs[config.lhs$label %in% x.var == FALSE,"label"], split = "(", fixed = TRUE)
    config.lhs[config.lhs$label %in% x.var == FALSE,"operator"] <- tolower(sapply(split.lhs,"[",1))
    config.lhs[config.lhs$label %in% x.var == FALSE,"arguments"] <- gsub(pattern = "\\)$", replacement = "", sapply(split.lhs,"[",2))

    ## normalize operator
    config.rhs$operator[config.rhs$operator %in% c("b","bin","binary")] <- "binary"
    config.rhs$operator[config.rhs$operator %in% c("c","cont","continuous")] <- "continuous"
    config.rhs$operator[config.rhs$operator %in% c("t","tte","time","timetoevent")] <- "timetoevent"
    config.rhs$operator[config.rhs$operator %in% c("g","gaus","gaussian")] <- "gaussian"
    config.rhs$operator[config.rhs$operator %in% c("s","strat","strata")] <- "strata"
    
    config.lhs$operator[config.lhs$operator %in% c("b","bin","binary")] <- "binary"
    config.lhs$operator[config.lhs$operator %in% c("c","cont","continuous")] <- "continuous"
    config.lhs$operator[config.lhs$operator %in% c("t","tte","time","timetoevent")] <- "timetoevent"
    config.lhs$operator[config.lhs$operator %in% c("g","gaus","gaussian")] <- "gaussian"
    config.lhs$operator[config.lhs$operator %in% c("s","strat","strata")] <- "strata"

    valid.operator <- c("binary","continuous","timetoevent","gaussian","strata")
    if(any(stats::na.omit(config.lhs$operator) %in% valid.operator == FALSE)){
        pb.lhs <- config.lhs[!is.na(config.lhs$operator) & config.lhs$operator %in% valid.operator == FALSE,]
        stop("Unknown operator \"",paste(pb.lhs$operator, collapse="\", \""),"\" on the right hand side of the \'formula\' argument.\n",
             "  It is part of the element(s) \"",paste(pb.lhs$label,collapse="\", \""),"\" \n", sep = "")
    }else if(any(stats::na.omit(config.rhs$operator) %in% valid.operator == FALSE)){
        pb.rhs <- config.rhs[!is.na(config.rhs$operator) & config.rhs$operator %in% valid.operator == FALSE,]
        stop("Unknown operator \"",paste(pb.rhs$operator, collapse="\", \""),"\" on the right hand side of the \'formula\' argument.\n",
             "  It is part of the element(s) \"",paste(pb.rhs$label,collapse="\", \""),"\" \n", sep = "")
    }

    ## ** extract number of endpoint variables
    endpoint.operator <- valid.operator[1:4]
    if(any(config.lhs$operator %in% endpoint.operator) && any(config.rhs$operator %in% endpoint.operator)){
        stop("Endpoint variables should not appear on both sides of argument \'formula\': \n",
             "  Detected endpoints (left): \"",paste(config.lhs[config.lhs$operator %in% endpoint.operator,"label"], collapse = "\", \""),"\". \n",
             "  Detected endpoints (right): \"",paste(config.rhs[config.rhs$operator %in% endpoint.operator,"label"], collapse = "\", \""),"\". \n"
             )
    }else if(all(config.lhs$operator %in% endpoint.operator == FALSE) && all(config.rhs$operator %in% endpoint.operator == FALSE)){
        stop("No endpoint detected in argument \'formula\'. \n",
             "  The formula should contain terms of the form \'type(endpoint,threshold)\' where type is one of \"",paste(endpoint.operator, collapse ="\", \""),"\". \n")
    }else if(config.lhs$operator %in% endpoint.operator){
        config.endpoint <- config.lhs ## may contain strata
        config.treatment <- config.rhs ## may contain strata

        label.endpoint <- config.lhs[config.lhs$operator %in% endpoint.operator,"label"]
        operator.endpoint <- config.lhs[config.lhs$operator %in% endpoint.operator,"operator"]
        arg.endpoint <- config.lhs[config.lhs$operator %in% endpoint.operator,"arguments"]
    }else{
        config.endpoint <- config.rhs ## may contain strata
        config.treatment <- config.lhs ## may contain strata

        label.endpoint <- config.rhs[config.rhs$operator %in% endpoint.operator,"label"]
        operator.endpoint <- config.rhs[config.rhs$operator %in% endpoint.operator,"operator"]
        arg.endpoint <- config.rhs[config.rhs$operator %in% endpoint.operator,"arguments"]
    }

    ## ** extract treatment variable
    if(all(is.na(config.treatment$var))){
        if(config.lhs$operator %in% endpoint.operator){
            stop("No treatment variable detected on the right hande side of argument \'formula\'. \n",
                 "  Should be something like: ",paste(deparse(x[[2]]),"~ treatment")," \n")
        }else{
            stop("No treatment variable detected on the left hande side of argument \'formula\'. \n",
                 "  Should be something like: ",paste("treatment ~", deparse(x[[3]]))," \n")
        }
    }

    treatment <- config.treatment[is.na(config.treatment$operator),"var"]
    if(length(treatment)>1){
        if(config.lhs$operator %in% endpoint.operator){
            stop("There should be exactly one treatment variable on right hand side of argument \'formula\'. \n",
                 "  Detected treatment variables: \"",paste(treatment, collapse = "\", \""),"\". \n")
        }else{
            stop("There should be exactly one treatment variable on left hand side of argument \'formula\'. \n",
                 "  Detected treatment variables: \"",paste(treatment, collapse = "\", \""),"\". \n")
        }
    }
  
    ## ** extract strata terms
    if(any(config.treatment$operator %in% "strata") && (any(config.endpoint$operator %in% "strata") || any(is.na(config.endpoint$operator))) ){
        if(config.lhs$operator %in% endpoint.operator){
            stop("Strata variables should not appear on both sides of argument \'formula\': \n",
                 "  Detected strata (left): \"",paste(config.endpoint[config.endpoint$operator %in% "strata" | is.na(config.endpoint$operator),"label"], collapse = "\", \""),"\". \n",
                 "  Detected strata (right): \"",paste(config.treatment[config.treatment$operator %in% "strata","label"], collapse = "\", \""),"\". \n"
                 )
        }else{
            stop("Strata variables should not appear on both sides of argument \'formula\': \n",
                 "  Detected strata (left): \"",paste(config.treatment[config.treatment$operator %in% "strata","label"], collapse = "\", \""),"\". \n",
                 "  Detected strata (right): \"",paste(config.endpoint[config.endpoint$operator %in% "strata" | is.na(config.endpoint$operator),"label"], collapse = "\", \""),"\". \n"
                 )
        }
    }else  if(any(is.na(config.endpoint$operator))){

        if(any(config.endpoint$operator %in% "strata")){ ## two ways to encode the strata variable
            config.strata <- config.endpoint[(is.na(config.endpoint$operator) | config.endpoint$operator %in% "strata"),,drop=FALSE]
            config.strata$label.fix <- config.strata$label
            config.strata[is.na(config.strata$operator),"label.fix"] <- paste0("strata(",config.strata[is.na(config.strata$operator),"label.fix"],")")
            stop("Using ",paste(config.strata$label,collapse ="+")," in argument \'formula\' is confusing \n",
                 "                           as it does not make consistent use of the strata operator. \n",
                 "  Consider using ", paste(config.strata$label.fix,collapse =" + "), " instead. \n",sep="")
        }

        strata <- config.endpoint[is.na(config.endpoint$operator),"var"]
        match <- FALSE
    }else if(any(config.treatment$operator %in% "strata") || any(config.endpoint$operator %in% "strata")){

        if(any(config.treatment$operator %in% "strata")){
            config.strata <- config.treatment[config.treatment$operator %in% "strata",,drop=FALSE]
        }else if(any(config.endpoint$operator %in% "strata")){
            config.strata <- config.endpoint[config.endpoint$operator %in% "strata",,drop=FALSE]
        }

        ## get each argument
        lsStrata.args <- strsplit(config.strata$arguments, split = ",", fixed = TRUE)
        
        ## identify each argument
        dfStrata.args <- do.call(rbind,lapply(1:NROW(config.strata), function(iS){
            catchArgument(lsStrata.args[[iS]], valid.arguments = c("variable","match"),
                          label = config.strata$label, name.operator = "strata", keep.all = TRUE)
        }))
        strata <- dfStrata.args$variable
        if(any(is.na(strata))){
            stop("Missing argument \'variable\' for strata operator in argument \'formula\'. \n",
                 "Detected in: ",paste(config.strata[is.na(strata),"label"], collapse = ", "),"\n")
        }
        if(all(is.na(dfStrata.args$match))){
            match <- FALSE
        }else{
            test <- try(match <- eval(expr = parse(text = unique(dfStrata.args$match))), silent = TRUE)
            if(inherits(test,"try-error")){
                stop(test,
                     "(occurred when evaluating argument match of ",paste(config.strata$label, collapse = "+")," in argument \'formula\') \n")
            }
            if(length(match)>1){
                stop("Argument \'match\' for strata operator should take the same value for all strata in argument \'formula\'. \n",
                     "Detected values: ",paste(match, collapse = ", "),"\n")
            }
        }
    }else{
        strata <- NULL
        match <- FALSE
    }


    ## ** extract endpoints and additional arguments
    n.endpoint <- length(label.endpoint)
    validArgs <- list(binary = c("endpoint","operator","weight"),
                      continuous = c("endpoint","threshold","restriction","operator","weight"),
                      timetoevent = c("endpoint","status","threshold","restriction","censoring","operator","weight"),
                      gaussian = c("mean","std","iid","threshold","restriction","operator","weight"))
    default.censoring <- c(binary = NA, continuous = NA, timetoevent = "right", gaussian = NA)
    df.args <- data.frame(type = operator.endpoint,
                          endpoint = rep(as.numeric(NA), n.endpoint),
                          status = rep("..NA..", n.endpoint),
                          threshold = rep(as.numeric(NA), n.endpoint),
                          censoring = default.censoring[operator.endpoint],
                          restriction = rep(as.numeric(NA), n.endpoint),
                          operator = rep(">0", n.endpoint),
                          weight = rep(ifelse(hierarchical,1,NA), n.endpoint))
    rownames(df.args) <- label.endpoint

    split.endpoint2args <- strsplit(arg.endpoint, split = ",", fixed = TRUE)
    
    for(iE in 1:n.endpoint){ ## iE <- 1

        ## identify each argument
        iArgs <- catchArgument(split.endpoint2args[[iE]], valid.arguments = validArgs[[operator.endpoint[iE]]],
                               label = label.endpoint[iE], name.operator = operator.endpoint[iE], keep.all = FALSE)
    
        ## trick to normalize argument names of Gaussian operator
        if("mean" %in% names(iArgs)){
            names(iArgs)[names(iArgs)=="mean"] <- "endpoint"
        }
        if("std" %in% names(iArgs)){
            names(iArgs)[names(iArgs)=="std"] <- "status"
        }
        if("iid" %in% names(iArgs)){
            names(iArgs)[names(iArgs)=="iid"] <- "censoring"
        }
        iArgs_save <- iArgs ## useful when errors are triggered

        ## process each argument
        if("endpoint" %in% names(iArgs)){
            iArgs[["endpoint"]] <- gsub("\"","",iArgs[["endpoint"]])
        }
        if("status" %in% names(iArgs)){
            iArgs[["status"]] <- gsub("\"","",iArgs[["status"]])
        }
        for(iName in intersect(c("threshold","restriction","weight"), names(iArgs))){
            ## handle the case where the numeric value is not defined directly but via a variable

            iTest <- try(iArgs[[iName]] <- eval(expr = parse(text = iArgs[[iName]])), silent = TRUE)

            if(inherits(iTest,"try-error") || (!is.numeric(iTest) && !is.logical(iTest) && !is.integer(iTest))){
                iTest <- try(iArgs[[iName]] <- eval(expr = parse(text = iArgs_save[[iName]]), envir = envir), silent = TRUE)
                if(inherits(iTest,"try-error")){
                    stop(iTest,
                         "(occurred when evaluating argument \'",iName,"\' of ",label.endpoint[iE]," in argument \'formula\') \n")
                }
            }

            if(inherits(iArgs[[iName]], "function")){
                packageTempo <- environmentName(environment(iArgs[[iName]]))
                txt <- ifelse(nchar(packageTempo)>0,paste0(" (package ",packageTempo,")"),"")
                stop("Argument \'",iName,"\' of ",label.endpoint[iE]," in argument \'formula\' is a function. \n",
                     "  Function ",iArgs_save[[iName]],txt," \n",
                     "  It should refer to a numeric value or a variable defining a numeric value instead. \n")
            }
            
        }

        if("censoring" %in% names(iArgs)){
            iArgs["censoring"] <- gsub("\"","",iArgs["censoring"])
        }
        if("operator" %in% names(iArgs)){
            iArgs["operator"] <- gsub("\"","",iArgs["operator"])
        }

        ## store
        df.args[iE,names(iArgs)] <- iArgs
        
    }

    ## ** export
    if(all(is.na(df.args$weight))){
        df.args$weight <- rep(1/n.endpoint, n.endpoint)
    }else if(sum(df.args$weight, na.rm = TRUE)<1){
        df.args$weight[!is.na(df.args$weight)] <- (1-sum(df.args$weight, na.rm = TRUE))/sum(is.na(df.args$weight))
    }

    out <- list(treatment = treatment,
                type = df.args$type,
                endpoint = df.args$endpoint,
                threshold = df.args$threshold,
                status = df.args$status,
                operator = df.args$operator,
                weightEndpoint = df.args$weight,
                censoring = df.args$censoring,
                restriction = df.args$restriction,
                strata = strata,
                match = match)
    return(out)
}

## * catchArgument
##' @title Associate input to argument
##' @noRd
##' @examples
##' ## check error
##' ## catchArgument(c("time","event","restriction=5=5"), valid.argument = c("endpoint","status","threshold","restriction"), name.operator = "tte")
##' ## catchArgument(c("time","event","restriction=5","restriction=5"), valid.argument = c("endpoint","status","threshold","restriction"), name.operator = "tte")
##' ## catchArgument(c("time","event","x","y","restriction=5"), valid.argument = c("endpoint","status","threshold","restriction"), name.operator = "tte")
##' ## catchArgument("", valid.argument = c("endpoint","status","threshold","restriction"), name.operator = "tte")
##'
##' ## normal
##' catchArgument(c("time","event","restriction=5"), valid.argument = c("endpoint","status","threshold","restriction"), name.operator = "tte")
catchArgument <- function(object, valid.arguments, label = NULL, name.operator = NULL, keep.all = FALSE){

    nvalid.arguments <- length(valid.arguments)

    ## prepare output
    out <- stats::setNames(rep(as.character(NA), nvalid.arguments), valid.arguments)
    keep.out <- stats::setNames(rep(FALSE, nvalid.arguments), valid.arguments)

    ## identify equal signs
    if(all(sapply(object,nchar)==0)){
        stop(name.operator," operator in argument \'formula\' must contain a name of variable between the parentheses \n",
             "This is not the case in: ",label,"\"\n")
    }
    object.split <- lapply(strsplit(trimws(object, which = "both"), split = "=", fixed = TRUE), trimws, which = "both")
    if(length(object.split)>nvalid.arguments){
        stop("Too many arguments for ",name.operator," operator in argument \'formula\'. \n",
             "  It can handle at most ",nvalid.arguments," arguments (\"",paste(valid.arguments, collapse="\", \""),"\") \n",
             "  ",length(object.split)," arguments detected in ",label,". \n")
    }
    object.length <- lengths(object.split)
    if(any(object.length>2)){
        stop(name.operator," operator in argument \'formula\' should not contain elements with multiple equal signs. \n",
             "  Multiple equal signs detected (\"",object[object.length>2][1],"\") in ",label,". \n",sep="")
    }

    ## match arguments
    if(any(object.length==2)){ ## named by the users
        object2.arg <- sapply(object.split[object.length==2],"[",1)
        object2.value <- sapply(object.split[object.length==2],"[",2)
        if(any(duplicated(object2.arg))){
            stop("Each argument can only appear once for each ",name.operator," operator in argument \'formula\'. \n",
                 "  Argument(s) \"",paste(unique(object2.arg[duplicated(object2.arg)]), collapse = "\", \""),"\" are duplicated in ",label,". \n")
        }
        if(any(object2.arg %in% valid.arguments == FALSE)){
            stop(name.operator," operator in argument \'formula\' only recognize arguments \"",paste(valid.arguments, collapse = "\", \""),"\". \n",
                 "Invalid arguments: \"",paste(setdiff(object2.arg, valid.arguments), collapse = "\", \""),"\" in ",label," \n",
                 sep="")
        }
        keep.out[object2.arg] <- TRUE
        out[object2.arg] <- object2.value
    }

    if(any(object.length==1)){ ## unname
        remain.arguments <- valid.arguments[is.na(out)]
        keep.out[remain.arguments[1:sum(object.length==1)]] <- TRUE
        out[remain.arguments[1:sum(object.length==1)]] <- sapply(object.split[object.length==1],"[",1)
    }    
    ## export
    if(keep.all){
        return(as.data.frame(as.list(out)))
    }else{
        return(as.data.frame(as.list(out[keep.out])))
    }
    
}



