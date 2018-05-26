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
#' @param method.tte [character] defines the method used to handle pairs
#' which can not be decidedly classified as favorable, unfavorable, or neutral because of censored observations (see details).
#' Can be \code{"Gehan"}, \code{"Gehan corrected"}, \code{"Peron"}, or \code{"Peron corrected"}.
#' Only relevant when there is one or more time-to-event endpoints.
#' Default value read from \code{BuyseTest.options()}.
#' @param method.inference [character] should a permutation test (\code{"permutation"} or \code{"stratified permutation"}),
#' or bootstrap resampling (\code{"bootstrap"} or \code{"stratified boostrap"})
#' be used to compute p-values and confidence intervals.
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analysed using endpoints of lower priority.
#' Default value read from \code{BuyseTest.options()}.
#' @param n.resampling [integer] the number of simulations used for computing the confidence interval and the p.values. See details.
#' Default value read from \code{BuyseTest.options()}.
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
#' @param method Obsolete. Alias for 'method.tte'.
#' @param n.bootstrap Obsolete. Alias for 'n.resampling'.
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
#' \bold{n.resampling:} The number of permutation replications must be specified to enable the computation of the confidence intervals and the p.value. 
#' A large number of permutations (e.g. \code{n.resampling=10000}) are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 
#' 
#' \bold{trace:} \code{2} reports all messages, \code{1} reports only the percentage of advancement of the permutation test, and \code{0} remains silent.
#' 
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus.
#' The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevent from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.
#' 
#' \bold{method:} Pairs which can not be decidedly classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}, Gehan 1965). 
#' Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions.
#' \code{method="Peron"} estimates these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. 
#' Probabilities of survival beyond the last observation are set NA, resulting in a non null probability that the pair is informative.
#' See Peron et al. (2016) for more details. \cr
#' Due to the presence of uninformative pairs, the proportion of favorable, unfavorable, or neutral pairs is underestimated. 
#' \code{method="Gehan corrected"} and \code{method="Peron corrected"} aim at correcting this bias
#' by multiplying the contribution of each pair by the inverse of the total number of pairs minus the number of uninformative pairs
#' and setting the number of uninformative pairs to 0.
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references 
#' Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
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

## * BuyseTest
##' @rdname BuyseTest
##' @export
BuyseTest <- function(formula,
                      data, 
                      treatment = NULL,
                      endpoint = NULL,
                      type = NULL,
                      threshold = NULL,
                      censoring = NULL,
                      operator = NULL,
                      strata = NULL, 
                      method.tte = NULL,
                      neutral.as.uninf = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      keep.comparison = NULL,
                      alternative = "two.sided", 
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL,
                      n.bootstrap,
                      method){

    name.call <- names(match.call())
    option <- BuyseTest.options()
    
    ## ** compatibility with previous version
    if(!missing(n.bootstrap)){
        stop("Argument \'n.bootstrap\' is obsolete. \n",
             "It has been replaced by the argument \'n.resampling\' \n")
    }
    if(!missing(method)){
        stop("Argument \'method\' is obsolete. \n",
             "It has been replaced by the argument \'method.tte\' \n")
    }
    
    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(alternative = alternative,
                              name.call = name.call,
                              censoring = censoring,
                              cpus = cpus,
                              data = data,
                              endpoint = endpoint,
                              formula = formula,
                              keep.comparison = keep.comparison,
                              method.tte = method.tte,
                              n.resampling = n.resampling,
                              neutral.as.uninf = neutral.as.uninf,
                              operator = operator,
                              option = option,
                              seed = seed,
                              strata = strata,
                              threshold = threshold,
                              trace = trace,
                              treatment = treatment,
                              type = type,
                              method.inference = method.inference)

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = outArgs)
    }

    ## ** initialization data
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    outArgs[c("data","level.treatment","level.strata","n.strata","allstrata")] <- initializeData(data = outArgs$data,
                                                                                                 type = outArgs$type,
                                                                                                 endpoint = outArgs$endpoint,
                                                                                                 operator = outArgs$operator,
                                                                                                 strata = outArgs$strata,
                                                                                                 treatment = outArgs$treatment)

    ## ** create weights matrix for survival endpoints
    if(outArgs$D.TTE>0 && outArgs$D>1){
        ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
        outArgs[c("Wscheme","index.survivalM1","threshold.TTEM1")] <- buildWscheme(endpoint = outArgs$endpoint,
                                                                                   D = outArgs$D,
                                                                                   D.TTE = outArgs$D.TTE,
                                                                                   type = outArgs$type,
                                                                                   threshold = outArgs$threshold)
    }else{ #  factice arguments. Will be sent to the C++ arguments to fill the argument but not used by the function.
        outArgs$Wscheme <- matrix(nrow=0,ncol=0)
        outArgs$index.survivalM1 <- numeric(0)
        outArgs$threshold.TTEM1 <- numeric(0)
    }
    
    ## ** Displayd
    if (outArgs$trace > 1) {
        outPrint <- do.call(printGeneral, args = outArgs)
    }
    
    ## ** Punctual estimation
    if (outArgs$trace > 1) {cat("Punctual estimation ")}

    time <- system.time({
        outPunctual <- .BuyseTest(data = outArgs$data,
                                  censoring = outArgs$censoring,
                                  correction.tte = outArgs$correction.tte,
                                  D = outArgs$D,
                                  D.TTE = outArgs$D.TTE,
                                  endpoint = outArgs$endpoint,
                                  index.survivalM1 = outArgs$index.survivalM1,                       
                                  keep.comparison = outArgs$keep.comparison,
                                  level.treatment = outArgs$level.treatment,
                                  method.inference = "none",
                                  method.tte = outArgs$method.tte,
                                  neutral.as.uninf = outArgs$neutral.as.uninf,
                                  n.strata = outArgs$n.strata,
                                  returnIndex = option$returnIndex,
                                  strata = outArgs$allstrata,
                                  threshold = outArgs$threshold,
                                  treatment = outArgs$treatment,
                                  threshold.TTEM1 = outArgs$threshold.TTEM1,
                                  type = outArgs$type,
                                  Wscheme = outArgs$Wscheme)
    })
    if (outArgs$trace > 1) {cat("(done) \n")}
    
    ## ** Permutation test
    if (outArgs$method.inference %in% c("permutation","bootstrap","stratified permutation", "stratified bootstrap")) {

        ## ** define environment
        envirBT <- environment()
        envirBT$.BuyseTest <- .BuyseTest
        envirBT$initializeData <- initializeData
        envirBT$initializeSurvival_Peto <- initializeSurvival_Peto
        envirBT$initializeSurvival_Peron <- initializeSurvival_Peron
       
        ## ** run
        outResampling <- inferenceResampling(envirBT)
        
    }else if(outArgs$method.inference %in% c("asymptotic")){
        stop("Not implemented yet \n")
        outResampling <- list()
    }else{
        outResampling <- list()
    }
    
    ## ** Gather results into a BuyseRes object
    method.tte <- c("Gehan","Peto","Efron","Peron")[outArgs$method.tte+1]
    type <- c("Binary","Continuous","TimeToEvent")[outArgs$type]

    BuyseRes.object <- BuyseRes(
        count.favorable = outPunctual$count_favorable,      
        count.unfavorable = outPunctual$count_unfavorable,
        count.neutral = outPunctual$count_neutral,    
        count.uninf = outPunctual$count_uninf,
        n.pairs = outPunctual$n_pairs,
        delta.netChance = outPunctual$delta_netChance,
        delta.winRatio = outPunctual$delta_winRatio,
        Delta.netChance = outPunctual$Delta_netChance,
        Delta.winRatio = outPunctual$Delta_winRatio,
        index.neutralT = outPunctual$index_neutralT,
        index.neutralC = outPunctual$index_neutralC,
        index.uninfT = outPunctual$index_uninfT,
        index.uninfC = outPunctual$index_uninfC,
        type = type,
        endpoint = outArgs$endpoint,
        level.treatment = outArgs$level.treatment,
        method.tte = data.frame(method = method.tte, correction = outArgs$correction.tte, stringsAsFactors = FALSE),
        method.inference = outArgs$method.inference,
        strata = outArgs$strata,
        level.strata = outArgs$level.strata,
        threshold = outArgs$threshold,
        n.resampling = outArgs$n.resampling,
        deltaResampling.netChance = outResampling$deltaResampling_netChance,
        deltaResampling.winRatio = outResampling$deltaResampling_winRatio,
        DeltaResampling.netChance = outResampling$DeltaResampling_netChance,
        DeltaResampling.winRatio = outResampling$DeltaResampling_winRatio,
        tableComparison = outPunctual$tableComparison,
        args = list(indexC = which(outArgs$data[[outArgs$treatment]]==outArgs$level.treatment[1]),
                    indexT = which(outArgs$data[[outArgs$treatment]]==outArgs$level.treatment[2])
                    )
    )

    ## ** export
    return(BuyseRes.object)
}

## * .BuyseTest
.BuyseTest <- function(data,
                       censoring,
                       correction.tte,
                       D,
                       D.TTE,
                       endpoint,
                       index.survivalM1,                       
                       keep.comparison,
                       level.treatment,
                       method.inference,
                       method.tte,
                       neutral.as.uninf,
                       returnIndex,
                       n.strata,
                       strata,
                       threshold,
                       treatment,
                       threshold.TTEM1,
                       type,
                       Wscheme){

    ## ** Resampling
    data <- data.table::copy(data)
    if(method.inference == "none"){
        ## do nothing
    }else if(method.inference == "permutation"){
        data[, c(treatment) := .SD[[1]][sample.int(.N, size = .N, replace = FALSE)], .SDcols = treatment]
    }else if(method.inference == "bootstrap"){
        data <- data[sample.int(.N, size = .N, replace = TRUE)]    
    }else if(method.inference == "stratified permutation"){
        data[, c(treatment) := .SD[[1]][sample.int(.N, size = .N, replace = FALSE)], .SDcols = treatment, by = strata]
    }else if(method.inference == "stratified bootstrap"){
        data <- data[,.SD[sample.int(.N, size = .N, replace = TRUE)], by = strata]
    }

    ## ** Check valid resampling
    if(is.null(strata)){
        n.groups <- length(unique(data[[treatment]]))
    }else{
        n.groups <- data[,length(unique(.SD[[1]])), by = strata, .SDcols = treatment][[2]]
    }
    if (any(n.groups!=2) || length(n.groups) != n.strata) { ## failure of the resampling
        ##        return(matrix(NA, nrow = n.strata + 1, ncol = 2*D))
        return(NULL)
    }
    
    ## ** Initialize data    

    ## *** data: split the data according to the two levels
    indexT <- which(data[[treatment]] == level.treatment[2])
    indexC <- which(data[[treatment]] == level.treatment[1])

    ## *** data: extract endpoint 
    M.Treatment <- as.matrix(data[indexT,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
    M.Control <- as.matrix(data[indexC,endpoint,with=FALSE]) # matrix of endpoints for the control arm

    ## *** strata
    if(!is.null(strata)){
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
         index.strataT <- lapply(1:n.strata,function(iS){
            which(data[indexT,strata[1],with=FALSE] == iS) - 1
        })
        index.strataC <- lapply(1:n.strata,function(iS){
            which(data[indexC,strata[1],with=FALSE] == iS) - 1
        })
                                        
    }else{ # if there is no strata variable the same strata is used for all patient
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
        index.strataT <- list(0:(length(indexT)-1))
        index.strataC <- list(0:(length(indexC)-1))
    }

    ## *** data: censoring
    if(!is.null(censoring)){
        M.delta.Treatment <- as.matrix(data[indexT,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
        M.delta.Control <- as.matrix(data[indexC,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    }else{ # if the is no time to event variables
        M.delta.Treatment <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
        M.delta.Control <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }

    ## *** efron correction
    if(method.tte==2){ # "Efron"
        for(iStrata in 1:n.strata){ ## iStrata <- 1

                index.typeTTE <- which(type==3)
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

    ## *** Update survival
    if(method.tte == 0){ ## Gehan
        outSurv <- list(list.survivalT = lapply(1:D.TTE, matrix),
                        list.survivalC = lapply(1:D.TTE, matrix))        
    }else if(method.tte == 1){ ## Peto
        outSurv <- initializeSurvival_Peto(M.Treatment = M.Treatment,
                                           M.Control = M.Control,
                                           M.delta.Treatment = M.delta.Treatment,
                                           M.delta.Control = M.delta.Control,
                                           endpoint = endpoint,
                                           D.TTE = D.TTE,
                                           type = type,
                                           threshold = threshold,
                                           index.strataT = index.strataT,
                                           index.strataC = index.strataC,
                                           n.strata = n.strata)
    }else if(method.tte %in% 2:3){
        outSurv <- initializeSurvival_Peron(M.Treatment = M.Treatment,
                                            M.Control = M.Control,
                                            M.delta.Treatment = M.delta.Treatment,
                                            M.delta.Control = M.delta.Control,
                                            endpoint = endpoint,
                                            D.TTE = D.TTE,
                                            type = type,
                                            threshold = threshold,
                                            index.strataT = index.strataT,
                                            index.strataC = index.strataC,
                                            n.strata = n.strata)
    }
    
    ## ** Computation
    resBT <-   GPC_cpp(Treatment = M.Treatment,
                       Control = M.Control,
                       threshold = threshold,
                       survEndpoint = (type == 3),
                       delta_Treatment = M.delta.Treatment,
                       delta_Control = M.delta.Control,
                       D = D,
                       returnIndex = returnIndex,
                       strataT = index.strataT,
                       strataC = index.strataC,
                       n_strata = n.strata,
                       n_TTE = D.TTE,
                       Wscheme = Wscheme,
                       index_survivalM1 = index.survivalM1,
                       threshold_TTEM1 = threshold.TTEM1,
                       list_survivalT = outSurv$list.survivalT,
                       list_survivalC = outSurv$list.survivalC,
                       methodTTE = method.tte,
                       correctionTTE = correction.tte,
                       neutralAsUninf = neutral.as.uninf,
                       keepComparison = keep.comparison
                       )
    
    ## ** export
    if(method.inference == "none"){
        return(resBT)
    }else{
        Mout <- cbind(rbind(resBT$delta_netChance, resBT$Delta_netChance),
                      rbind(resBT$delta_winRatio, resBT$Delta_winRatio))
        dimnames(Mout) <- list(c(paste0("delta.",1:n.strata),"Delta"),
                               c(paste0("netChance.",1:D),paste0("winRatio.",1:D))
                               )
        return(Mout)
    }
}







