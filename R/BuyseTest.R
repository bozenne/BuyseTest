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
#' @param method [character] defines the method used to handle pairs which can not be decidedly classified as favorable, unfavorable, or neutral because of censored observations (see details).
#' Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}.
#' Only relevant when there is one or more time-to-event endpoints.
#' Default value read from \code{BuyseTest.options()}.
#' @param correctionTTE [logical] bias correction for non-informative pairs. Experimental.
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
#' \bold{trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the permutation test, and \code{0} remains silent.
#' 
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus.
#' The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevent from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.
#' 
#' \bold{method:} Pairs which can not be decidedly classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}, Gehan 1965). 
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
                      method = NULL,
                      correctionTTE = FALSE,
                      neutral.as.uninf = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      keep.comparison = NULL,
                      alternative = "two.sided", 
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL){

    name.call <- names(match.call())
    option <- BuyseTest.options()

    ## ** initialize arguments
    outArgs <- initializeArgs(alternative = alternative,
                              name.call = name.call,
                              censoring = censoring,
                              correctionTTE = correctionTTE,
                              cpus = cpus,
                              data = data,
                              endpoint = endpoint,
                              formula = formula,
                              keep.comparison = keep.comparison,
                              method = method,
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
    trace <- outArgs$trace

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = outArgs)
    }
    
    ## ** Display
    if (trace > 1) {
            outPrint <- do.call(printGeneral, args = outArgs)
    }
    
    ## ** Computations
    outArgs$name.call <- NULL
    outArgs$formula <- NULL
    level.strata <- outArgs$level.strata
    outArgs$level.strata <- NULL
    level.treatment <- outArgs$level.treatment
    outArgs$level.treatment <- NULL
    outArgs$formula <- NULL

    outBT <- do.call(.BuyseTest, args = outArgs)
    ## ** Gather results into a BuyseRes object
    method <- switch(as.character(outArgs$method),
                     "0" = "Gehan",
                     "1" = "Peto" ,
                     "2" = "Efron",
                     "3" = "Peron"
                     )

    BuyseRes.object <- BuyseRes(
        count.favorable = outBT$count_favorable,      
        count.unfavorable = outBT$count_unfavorable,
        count.neutral = outBT$count_neutral,    
        count.uninf = outBT$count_uninf,
        n.pairs = outBT$n_pairs,
        delta.netChance = outBT$delta_netChance,
        delta.winRatio = outBT$delta_winRatio,
        Delta.netChance = outBT$Delta_netChance,
        Delta.winRatio = outBT$Delta_winRatio,
        index.neutralT = outBT$index_neutralT,
        index.neutralC = outBT$index_neutralC,
        index.uninfT = outBT$index_uninfT,
        index.uninfC = outBT$index_uninfC,
        endpoint = outArgs$endpoint,
        level.treatment = level.treatment,
        method = method,
        method.inference = outArgs$method.inference,
        strata = outArgs$strata,
        level.strata = level.strata,
        threshold = outArgs$threshold,
        n.resampling = outArgs$n.resampling,
        deltaResampling.netChance = outBT$deltaResampling_netChance,
        deltaResampling.winRatio = outBT$deltaResampling_winRatio,
        DeltaResampling.netChance = outBT$DeltaResampling_netChance,
        DeltaResampling.winRatio = outBT$DeltaResampling_winRatio,
        tableComparison = outBT$tableComparison,
        args = list(indexC = which(outArgs$data[[outArgs$treatment]]==level.treatment[1]),
                    indexT = which(outArgs$data[[outArgs$treatment]]==level.treatment[2])
                    )
    )

    ## ** export
    return(BuyseRes.object)
}

## * .BuyseTest
.BuyseTest <- function(alternative, 
                       censoring,
                       correctionTTE,
                       cpus,
                       data, 
                       endpoint,
                       index.survivalM1,
                       keep.comparison,                       
                       method,
                       method.inference,
                       n.resampling,
                       neutral.as.uninf,
                       operator,
                       seed,
                       strata, 
                       threshold,
                       threshold.TTEM1,
                       trace,
                       treatment,
                       type,
                       Wscheme){
    
    ## ** Initialize dataset
    outData <- initializeData(data = data,
                              treatment = treatment,
                              endpoint = endpoint,
                              censoring = censoring,
                              type = type,
                              operator = operator,
                              method = method,
                              strata = strata)
    
    index.strataT <- outData$index.strataT
    index.strataC <- outData$index.strataC
    n.strata <- outData$n.strata
    level.strata <- outData$level.strata
    M.Treatment <- outData$M.Treatment
    M.Control <- outData$M.Control
    M.delta.Treatment <- outData$M.delta.Treatment
    M.delta.Control <- outData$M.delta.Control
    D <- length(endpoint)
    D.TTE <- sum(type==3)
        
    ## ** Initialize survival
    if(method==0){
        ## factice lists. Will be sent to the C++ arguments to fill the argument but not used by the function.
        list.survivalT <- list()
        list.survivalC <- list()
    }else if(method==1){
        outSurv <- initializeSurvival_Peto(M.Treatment=M.Treatment,
                                          M.Control=M.Control,
                                          M.delta.Treatment=M.delta.Treatment,
                                          M.delta.Control=M.delta.Control,
                                          endpoint=endpoint,
                                          D.TTE=D.TTE,
                                          type=type,
                                          threshold=threshold,
                                          index.strataT=index.strataT,
                                          index.strataC=index.strataC,
                                          n.strata=n.strata)
        list.survivalT <- outSurv$list.survivalT
        list.survivalC <- outSurv$list.survivalC
    }else if(method %in% 2:3){
        outSurv <- initializeSurvival_Peron(M.Treatment=M.Treatment,
                                           M.Control=M.Control,
                                           M.delta.Treatment=M.delta.Treatment,
                                           M.delta.Control=M.delta.Control,
                                           endpoint=endpoint,
                                           D.TTE=D.TTE,
                                           type=type,
                                           threshold=threshold,
                                           index.strataT=index.strataT,
                                           index.strataC=index.strataC,
                                           n.strata=n.strata)
        list.survivalT <- outSurv$list.survivalT
        list.survivalC <- outSurv$list.survivalC
    }

    ## ** Punctual estimation
    if (trace > 1) {cat("Punctual estimation ")}

    time <- system.time({
        outPunctual <- GPC_cpp(Treatment = M.Treatment,
                               Control = M.Control,
                               threshold = threshold,
                               survEndpoint = (type == 3),
                               delta_Treatment = M.delta.Treatment,
                               delta_Control = M.delta.Control,
                               D = D,
                               returnIndex = TRUE,
                               strataT = index.strataT,
                               strataC = index.strataC,
                               n_strata = n.strata,
                               n_TTE = D.TTE,
                               Wscheme = Wscheme,
                               index_survivalM1 = index.survivalM1,
                               threshold_TTEM1 = threshold.TTEM1, 
                               list_survivalT = list.survivalT,
                               list_survivalC = list.survivalC,
                               methodTTE = method,
                               correctionTTE = correctionTTE,
                               neutralAsUninf = neutral.as.uninf,
                               keepComparison = keep.comparison
                               )
    })
    if (trace > 1) {cat("(done) \n")}

    ## ** Permutation test
    
    if (method.inference %in% c("permutation","bootstrap","stratified permutation", "stratified bootstrap")) {
    
        ## *** display
        if (trace > 1) {
            if (time[3] == 0) {
                time.punctual <- "<0.001 s"
                time.permutation <- paste0("<",signif(0.001*n.resampling/cpus,4)," s")
            }else{
                time.punctual <- paste(time[3],"s")
                time.permutation <- paste(signif(time[3]*n.resampling/cpus,4),"s")
            }
            txt.type <- switch(method.inference,
                               "bootstrap" = "bootstrap resampling",
                               "stratified bootstrap" = "stratified bootstrap resampling",
                               "permutation" = "permutation test",
                               "stratified permutation" = "stratified permutation test")
            
            cat("Settings (",txt.type,"): \n",
                "   > requested time for one sample: ", time.punctual, "\n",
                "   > estimated time for ", n.resampling, " samples with ", cpus, " core", if (cpus > 1) {"s"}, ": ", time.permutation, "\n", sep = "")
            if (!is.null(seed)) {
                cat("   > seed", if (cpus > 1) {"s"}, ": ",paste(seq(seed,seed + cpus - 1), collapse = " "), sep = "")       
            }
            cat("\n")         
        }

        ## *** computations
        n.eachStrataT <- unlist(lapply(index.strataT, length)) 
        n.eachStrataC <- unlist(lapply(index.strataC, length))
        n.T <- sum(n.eachStrataT)
        n.C <- sum(n.eachStrataC)
        nCumSum.strataControl <- cumsum(c(1,n.eachStrataC))
        nCumSum.strataTreatment <- cumsum(c(1,n.eachStrataT))
        envirBT <- environment()
        envirBT$.BuyseTest <- .BuyseTest
        envirBT$warperResampling <- warperResampling
        envirBT$initializeSurvival_Peto <- initializeSurvival_Peto
        envirBT$initializeSurvival_Peron <- initializeSurvival_Peron

        outResampling <- inferenceResampling(envirBT)
        outPunctual <- c(outPunctual,
                         outResampling)
    }else if(method.inference %in% c("asymptotic")){
        stop("Not implemented yet \n")
    }
    ## ** Export
    return(outPunctual)
}


