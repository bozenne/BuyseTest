## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Generalized Pairwise Comparisons (GPC)
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event outcomes.
#' @param formula [formula] a symbolic description of the model to be fitted. The response variable should be a binary variable defining the treatment arms. 
#' The rest of the formula should indicate the strata variables (if any) and the endpoints by order of priority. 
#' @param data [data.frame] dataset.
#' @param treatment [character] the name of the treatment variable identifying the control and the experimental group.
#' Disregarded if the argument \code{formula} is defined.
#' @param endpoint [character vector] the name of the endpoint variable(s).
#' Disregarded if the argument \code{formula} is defined.
#' @param operator [character vector] the sign defining a favorable endpoint:
#' ">0" indicates that higher values are favorable while "<0" indicates the opposite.
#' Disregarded if the argument \code{formula} is defined.
#' @param threshold [numeric vector] critical values used to compare the pairs.
#' There must be one threshold for each endpoint variable.
#' Disregarded if the argument \code{formula} is defined.
#' @param strata [numeric vector] if not \code{NULL}, the GPC will be applied within each group of patient defined by the strata variable(s).
#' Disregarded if the argument \code{formula} is defined.
#' @param censoring [character vector] the name of the binary variable(s) indicating whether the endpoint was observed or censored.
#' There must be one threshold for each endpoint variable.
#' Must value \code{NA} when the endpoint is not a time to event.
#' Disregarded if the argument \code{formula} is defined.
#' @param type [character vector] the type of each endpoint: \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.
#' @param method [character] defines the method used to handle pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations (see details).
#' Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}.
#' Only relevant when there is one or more time-to-event endpoints.
#' Default value read from \code{BuyseTest.options()}.
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analysed using endpoints of lower priority.
#' Default value read from \code{BuyseTest.options()}.
#' @param n.permutation [integer] the number of permutations used for computing the confidence interval and the p.values. See details.
#' Default value read from \code{BuyseTest.options()}.
#' @param prob.alloc [0<double<1] the resampling probability for assignement to the experimental group in the permutation test.
#' Can also be \code{NULL} to use the proportion of patients in the experimental group.
#' @param stratified [logical] should the assignement in the permutation test be performed within strata.
#' This means that the \code{prob.alloc} will be satisfyied within strata not only globally.
#' @param keep.comparison [logical] should the result of each pairwise comparison be kept?
#' @param alternative [character] the alternative hypothesis.
#' Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. 
#' @param seed [integer, >0] the seed to consider for the permutation test.
#' Default value read from \code{BuyseTest.options()}.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' Default value read from \code{BuyseTest.options()}.
#' @param trace [integer] should the execution of the function be traced ? See details.
#' Default value read from \code{BuyseTest.options()}.
#' 
#' @details 
#' \bold{treatment:} The variable corresponding to \code{treatment} in data must have only two levels (e.g. \code{0} and \code{1}). \cr
#' \bold{endpoint, threshold, censoring, operator, and type:}  they must have the same length. \cr
#' \code{threshold} must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. \cr
#' \code{censoring} must be \code{NA} for binary or continuous endpoints and indicate a variable in data for time to event endpoints. 
#' Short forms for endpoint \code{type} are \code{"bin"} (binary endpoint), \code{"cont"} (continuous endpoint), \
#' code{"TTE"} (time-to-event endpoint). 
#' \bold{operator:} when the operator is set to \code{"<0"} the corresponding column in the dataset is multiplied by \code{-1}.
#' 
#' \bold{n.permutation:} The number of permutation replications must be specified to enable the computation of the confidence intervals and the p.value. 
#' A large number of permutations (e.g. \code{n.permutation=10000}) are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 
#' 
#' \bold{trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the permutation test, and \code{0} remains silent.
#' 
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus. The parallelization relies on the \emph{snowfall} package (function \emph{sfClusterApplyLB}). The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevend from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.
#' 
#' \bold{method:} Pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}, Gehan 1965). 
#' Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions. 
#' \code{method="Peto"} estimate these probabilities using the common Kaplan-Meier estimator of the survival function for treated and control patients (Peto et al. 1972). 
#' \code{method="Efron"}, and \code{method="Peron"} estimate these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. 
#' When the largest observation is censored, it is not possible to estimate the survival probability by the Kaplan-Meier estimator beyond this time point.  
#' \code{method="Efron"} treats the largest observations in each patient group as if it were uncensored (Efron 1967). 
#' \code{method="Peron"} treats the probability of survival beyond the last observation as NA, resulting in a non null probability that the pair is uninformative    
#' See Peron et al. (2016) for more details.
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references 
#' Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clincal trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2016). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. Statistical Methods in Medical Research. \cr
#' Efron B (1967). \bold{The two sample problem with censored data}. \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} 4:831-583 \cr
#' Peto R, Peto J (1972). \bold{Asymptotically efficient rank invariant test procedures}. \emph{Journal of the Royal Statistical Society - series A} 135(2):185-198 \cr
#' Gehan EA (1965). \bold{A generalized two-sample Wilcoxon test for doubly censored data}. \emph{Biometrika}  52(3):650-653 \cr
#'
#' @seealso 
#' \code{\link{BuyseRes-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' 
#' @example
#' R/examples/EX_BuyseTest.R
#'     
#' @keywords function BuyseTest
#'

## * Function - BuyseTest
#' @rdname BuyseTest
#' @export
BuyseTest <- function(formula,
                      data, 
                      treatment = NULL,
                      endpoint = NULL,
                      type = NULL,
                      threshold = NULL,
                      censoring = NULL,
                      operator = NULL,
                      strata = NULL, 
                      method = NULL,
                      neutral.as.uninf = NULL,
                      n.permutation = NULL,
                      prob.alloc = NULL,
                      stratified = FALSE,
                      keep.comparison = NULL,
                      alternative = "two.sided", 
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL){
    
### ** normalize arguments
    BuyseCall <- match.call()
    option <- BuyseTest.options()
    if(is.null(method)){ method <- option$method }
    if(is.null(neutral.as.uninf)){ neutral.as.uninf <- option$neutral.as.uninf }
    if(is.null(keep.comparison)){ keep.comparison <- option$keep.comparison }
    if(is.null(n.permutation)){ n.permutation <- option$n.permutation }
    if(is.null(seed)){ seed <- option$seed }
    if(is.null(cpus)){ cpus <- option$cpus }
    if(is.null(trace)){ trace <- option$trace }

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    
    ## *** formula interface
    if(!missing(formula)){
        argnames <- c("treatment", "endpoint", "type", "threshold", "censoring", "strata")
        if(any(names(BuyseCall) %in% argnames)){
            txt <- paste(names(BuyseCall)[names(BuyseCall) %in% argnames], collapse = "\' \'")
            warning("BuyseTest : argument",if(length(txt)>1){"s"}," \'",txt,"\' ha",if(length(txt)>1){"ve"}else{"s"}," been ignored \n",
                    "when specified, only argument \'formula\' is used \n")
        }
        
        resFormula <- initFormula(formula)
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
    }
    
    ## *** Treatment: extract the 2 levels
    validCharacter(treatment,
                   valid.length = 1,
                   method = "BuyseTest")
    validNames(data,
               required.values = treatment,
               valid.length = NULL,
               method = "BuyseTest")
    levels.treatment <- levels(as.factor(data[[treatment]])) # extraction of the levels of the treatment variable
    if (length(levels.treatment) != 2) {
        stop("BuyseTest : wrong specification of \'treatment\' \n",
             "the corresponding column in \'data\' must have exactly 2 levels \n",
             "proposed levels : ",paste(levels.treatment,collapse = " "),"\n")
    }
    
    ## *** endpoint
    validNames(data,
               required.values = endpoint,
               valid.length = NULL,
               method = "BuyseTest")
    D <- length(endpoint) # number of endpoints
    
    ## *** type: convert type to numeric and count the number of endpoints
    validType1 <- c("b","bin","binary")
    validType2 <- c("c","cont","continuous")
    validType3 <- c("t","tte","time","timetoevent")
    type <- tolower(type)

    validCharacter(type,
                   valid.values = c(validType1,validType2,validType3),
                   valid.length = D,
                   method = "BuyseTest")

    type[grep(paste(validType1,collapse="|"), type)] <- "1" 
    type[grep(paste(validType2,collapse="|"), type)] <- "2" 
    type[grep(paste(validType3,collapse="|"), type)] <- "3" 
    type <- as.numeric(type) # type is an integer equal to 1 (binary endpoint), 2 (continuous endpoint) or 3 (time to event endpoint)
    
    D.TTE <- sum(type == 3) # number of time to event endpoints

    n.typePerEndpoint <- tapply(type,endpoint, function(x){length(unique(x))})
    if(any(n.typePerEndpoint>1)){
        message <- paste0("several types have been specified for endpoint(s) ",
                          paste0(unique(endpoint)[n.typePerEndpoint>1],collapse = ""),
                          "\n")        
        stop("BuyseTest: wrong specification of \'endpoint\' or \'type\' \n",message)
    }

    ##  convert character/factor to numeric for binary endpoints
    name.bin <- endpoint[which(type %in% 1)]
    if(length(name.bin)>0){
        data.class <- sapply(data,class)
        
        for(iBin in name.bin){
            if(data.class[iBin] %in% c("numeric","integer") == FALSE){
                data[[iBin]] <- as.numeric(as.factor(data[[iBin]])) - 1
            }
        }
    }

    ## *** censoring
    censoring <- initCensoring(censoring = censoring,
                               endpoint = endpoint,
                               type = type,
                               D = D,
                               D.TTE = D.TTE,
                               treatment = treatment,
                               strata = strata)
    
    validNames(data,
               name1 = "data",
               required.values = censoring,
               valid.length = NULL,
               refuse.NULL = FALSE,
               method = "BuyseTest")

    ## *** operator
    data <- applyOperator(data, operator = operator,
                          type = type, endpoint = endpoint, D = D)
    
    ## *** data: split the data according to the two levels
    if (!is.null(strata)) {
        validNames(data,
                   name1 = "strata",
                   required.values = strata,
                   valid.length = NULL,
                   method = "BuyseTest")
    }
    
    indexT <- which(data[[treatment]] == levels.treatment[2])
    indexC <- which(data[[treatment]] == levels.treatment[1])
    dataT <- data[indexT, c(endpoint, strata, censoring), with = FALSE]
    dataC <- data[indexC, c(endpoint, strata, censoring), with = FALSE]

    n.Treatment <- NROW(dataT) # number of patient in the treatment arm
    n.Control <- NROW(dataC) # number of patient in the control arm
    
    ## *** threshold
    threshold <- initThreshold(threshold = threshold, type = type, D = D,
                               endpoint = endpoint)
    
    ## *** strata
    res <- initStrata(strata = strata,
                      dataT = dataT, dataC = dataC,
                      n.Treatment = n.Treatment, n.Control = n.Control,
                      endpoint = endpoint, censoring = censoring)
    
    index.strataT <- res$index.strataT 
    index.strataC <- res$index.strataC 
    n.strata <- res$n.strata 
    levels.strata <- res$levels.strata 
    
    ## *** method
    validCharacter(method,
                   valid.length = 1,
                   valid.values = c("Gehan","Peto","Efron","Peron"),
                   method = "BuyseTest")
    
    method <- switch(method,
                     "Gehan" = 0,
                     "Peto" = 1,
                     "Efron" = 2,
                     "Peron" = 3)
    
    if (D.TTE == 0) {
        method <- 0
        if ("method" %in% names(BuyseCall) && trace > 0) {
            message("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
        }
    }
    
    ## *** alternative
    validCharacter(alternative,
                   valid.length = 1,
                   valid.values = c("two.sided", "less", "greater"),
                   method = "BuyseTest")
    
    ## *** n.permutation
    validInteger(n.permutation,
                 valid.length = 1,
                 min = 0,
                 method = "BuyseTest")
    
    ## *** proba
    if (is.null(prob.alloc)) { # if prob.alloc is not set by the user
        prob.alloc <- n.Treatment/(n.Treatment + n.Control)  # it is set to the proportion of patients in the treatment arm.
    }else{
        validNumeric(prob.alloc,
                     valid.length = 1,
                     min = 0,
                     max = 1,
                     method = "BuyseTest")
    }
    
                                        # stratified
    validLogical(stratified,
                 valid.length = 1,
                 method = "BuyseTest")
    
    ## *** seed
    validInteger(seed,
                 valid.length = 1,
                 refuse.NULL = FALSE,
                 min = 1,
                 method = "BuyseTest")
    
    ## *** cpu
    if (cpus == "all") { 
        cpus <- parallel::detectCores() # this function detect the number of CPU cores 
    }else{# if several cpus are intended to be used, check this correspond to a valid number of CPU cores
        validInteger(cpus,
                     valid.length = 1,
                     valid.values = 1:parallel::detectCores(),
                     method = "BuyseTest")
    }
    
    ## *** data
    res <- initData(dataT = dataT,
                    dataC = dataC,
                    type = type,
                    endpoint = endpoint,
                    D = D,
                    operator = operator,
                    censoring = censoring,
                    index.strataT = index.strataT,
                    index.strataC = index.strataC,
                    n.strata = n.strata,                  
                    method = method,
                    D.TTE = D.TTE,
                    threshold = threshold,
                    Wscheme = NULL,
                    test = TRUE,
                    trace = trace)
    
    M.Treatment <- res$M.Treatment
    M.Control <- res$M.Control 
    M.delta_Treatment <- res$M.delta_Treatment
    M.delta_Control <- res$M.delta_Control
    Wscheme <- res$Wscheme
    threshold_TTEM1 <- res$threshold_TTEM1
    index_survivalM1 <- res$index_survivalM1
    list_survivalT <- res$list_survivalT
    list_survivalC <- res$list_survivalC
    
### ** 2- General display
    if (trace > 1) {
        printGeneral(levels.treatment = levels.treatment,
                     levels.strata = levels.strata,
                     n.strata = n.strata,
                     endpoint = endpoint,
                     threshold = threshold,
                     censoring = censoring,
                     operator = operator,
                     type = type,
                     D = D,
                     D.TTE = D.TTE,
                     method = method,
                     neutral.as.uninf = neutral.as.uninf,
                     Wscheme = if (D>1 && method %in% 1:3) {Wscheme} else {NULL}, 
                     threshold_TTEM1 = if (D>1 && method %in% 1:3) {threshold_TTEM1} else {NULL})
    }
    
### ** 3- Punctual estimation
    if (trace > 1) {cat("Punctual estimation \n")}

    time <- system.time({
        resPonctual <- GPC_cpp(Treatment = M.Treatment,
                               Control = M.Control,
                               threshold = threshold,
                               survEndpoint = (type == 3),
                               delta_Treatment = M.delta_Treatment,
                               delta_Control = M.delta_Control,
                               D = D,
                               returnIndex = TRUE,
                               strataT = index.strataT,
                               strataC = index.strataC,
                               n_strata = n.strata,
                               n_TTE = D.TTE,
                               Wscheme = Wscheme,
                               index_survivalM1 = index_survivalM1,
                               threshold_TTEM1 = threshold_TTEM1, 
                               list_survivalT = list_survivalT,
                               list_survivalC = list_survivalC,
                               methodTTE = method,
                               neutralAsUninf = neutral.as.uninf,
                               keepComparison = keep.comparison
                               )
    })
    if (trace > 1) {cat("   > done \n")}
    
### ** 4- Transfomration into BuyseRes object
    if (trace > 1) {cat("Conversion to BuyseRes object \n")}
    BuyseRes.object <- BuyseRes(
        count_favorable = resPonctual$count_favorable,      
        count_unfavorable = resPonctual$count_unfavorable,
        count_neutral = resPonctual$count_neutral,    
        count_uninf = resPonctual$count_uninf,
        delta = list(netChance = resPonctual$delta_netChance,
                     winRatio = resPonctual$delta_winRatio),
        Delta = list(netChance = resPonctual$Delta_netChance,
                     winRatio = resPonctual$Delta_winRatio),
        endpoint = endpoint,
        index_neutralT = resPonctual$index_neutralT,
        index_neutralC = resPonctual$index_neutralC,
        index_uninfT = resPonctual$index_uninfT,
        index_uninfC = resPonctual$index_uninfC,
        levels.treatment = levels.treatment,
        n_pairs = resPonctual$n_pairs,
        strata = levels.strata,
        threshold = threshold,
        conf.level = as.numeric(NA),
        tableComparison = resPonctual$tableComparison,
        args = list(indexT = indexT, indexC = indexC)
    )
    if (trace > 1) {cat("   > done \n")}
    
### ** 5- Permutation test 
    if (n.permutation > 0) {
        if (trace > 1) {
            cat("\n")
            printPermutation(prob.alloc, n.permutation, stratified, cpus, time, seed)
        }
        
        n.eachStrataT <- unlist(lapply(index.strataT, length)) # those variables are used during the permutation test and passed to the wrapper through the environment
        n.eachStrataC <- unlist(lapply(index.strataC, length))
        nCumSum.strataControl <- cumsum(c(1,n.eachStrataC))
        nCumSum.strataTreatment <- cumsum(c(1,n.eachStrataT))
        
        delta_permutation <- list(netChance = array(NA, dim = c(n.strata, D, n.permutation)),
                                 winRatio = array(NA, dim = c(n.strata, D, n.permutation))
                                 )
        Delta_permutation <- list(netChance = matrix(NA, nrow = D, ncol = n.permutation),
                                 winRatio = matrix(NA, nrow = D, ncol = n.permutation)
                                 )
        
        ## *** computation
        calcPermutation(environment())
        
        ## *** post treatment
        if (trace > 1) {cat("Post-Treatment and update of the BuyseRes object \n")}
        resCI <- calcCI(Delta = BuyseRes.object@Delta,
                        Delta_permutation = Delta_permutation, 
                        endpoint = endpoint,
                        D = D,
                        alternative = alternative,
                        alpha = 1 - option$conf.level,
                        n.permutation = n.permutation,
                        cpus = cpus,
                        trace = trace)
        
        ## *** update BuyseRes object 
        if(option$keep.permutation){
            BuyseRes.object@delta_permutation <- delta_permutation
        }
        
        BuyseRes.object@Delta_quantile <- resCI$Delta_quantile
        BuyseRes.object@n_permutation <- resCI$n_permutation
        BuyseRes.object@p.value <- resCI$p.value
        BuyseRes.object@conf.level <- option$conf.level
        validObject(BuyseRes.object)
        if (trace > 1) {cat("   > done \n")}
    }
    
### ** 6- export
    return(BuyseRes.object)
}
