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
#' Can be \code{"Gehan"}, \code{"Gehan corrected"}, \code{"Gehan IPCW"}, \code{"Peron"}, \code{"Peron corrected"}, or \code{"Peron IPCW"}.
#' Only relevant when there is one or more time-to-event endpoints.
#' Default value read from \code{BuyseTest.options()}.
#' @param model.tte [list] optionnal survival models relative to each time to each time to event outcome.
#' Models must \code{prodlim} objects and stratified on the treatment and strata variable.
#' @param method.inference [character] should a permutation test (\code{"permutation"} or \code{"stratified permutation"}),
#' or bootstrap resampling (\code{"bootstrap"} or \code{"stratified boostrap"})
#' be used to compute p-values and confidence intervals.
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analysed using endpoints of lower priority.
#' Default value read from \code{BuyseTest.options()}.
#' @param n.resampling [integer] the number of simulations used for computing the confidence interval and the p.values. See details.
#' Default value read from \code{BuyseTest.options()}.
#' @param keep.pairScore [logical] should the result of each pairwise comparison be kept?
#' @param alternative [character] the alternative hypothesis.
#' Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. 
#' @param seed [integer, >0] the seed to consider for the permutation test.
#' Default value read from \code{BuyseTest.options()}.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' Default value read from \code{BuyseTest.options()}.
#' @param trace [integer] should the execution of the function be traced ? See details.
#' Default value read from \code{BuyseTest.options()}.
#' @param keep.comparison Obsolete. Alias for 'keep.pairScore'.
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
#' \bold{method.tte:} Pairs which can not be decidedly classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method.tte="Gehan"}, Gehan 1965). 
#' Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions.
#' \code{method.tte="Peron"} estimates these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. 
#' Probabilities of survival beyond the last observation are set NA, resulting in a non null probability that the pair is informative.
#' See Peron et al. (2016) for more details. \cr
#' Due to the presence of uninformative pairs, the proportion of favorable, unfavorable, or neutral pairs is underestimated. 
#' \code{method.tte="Gehan corrected"} and \code{method.tte="Peron corrected"} aim at correcting this bias
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
#' Gehan EA (1965). \bold{A generalized two-sample Wilcoxon test for doubly censored data}. \emph{Biometrika}  52(3):650-653 \cr
#'
#' @seealso 
#' \code{\link{BuyseRes-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' @keywords function BuyseTest

## * BuyseTest (example)
#' @rdname BuyseTest
#' @examples
#' # reset the default value of the number of permuation sample
#' BuyseTest.options(method.inference = "none") # no permutation test
#'
#' #### simulate some data ####
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#'                                        # display 
#' if(require(survival)){
#'    resKM_tempo <- survfit(Surv(eventtime,status)~Treatment, data = df.data)
#'    plot(resKM_tempo)
#' }
#'
#' #### one time to event endpoint ####
#' BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data)
#' 
#' summary(BT) # net chance in favor of treatment
#' summary(BT, percentage = FALSE)  
#' summary(BT, statistic = "winRatio") # win Ratio
#' 
#' ## bootstrap to compute the CI
#' \dontrun{
#'     BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3)
#' }
#' \dontshow{
#'     BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e1, trace = 0)
#' }
#' summary(BT, statistic = "netChance") ## default
#' summary(BT, statistic = "winRatio") 
#' 
#' ## parallel boostrap
#' \dontrun{
#'     BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3, cpus = 2)
#'     summary(BT)
#' }
#' 
#' ## method Gehan is much faster but does not optimally handle censored observations
#' BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                 method.tte = "Gehan", trace = 0)
#' summary(BT)
#' 
#' #### one time to event endpoint: only differences in survival over 1 unit ####
#' BT <- BuyseTest(Treatment ~ TTE(eventtime, threshold = 1, censoring = status), data=df.data)
#' summary(BT)
#' 
#' #### one time to event endpoint with a strata variable
#' BT <- BuyseTest(Treatment ~ strata + TTE(eventtime, censoring = status), data=df.data)
#' summary(BT)
#' 
#' #### several endpoints with a strata variable
#' f <- Treatment ~ strata + T(eventtime, 1, status) + B(toxicity) 
#' f <- update(f, 
#'             ~. + T(eventtime, 0.5, status) + C(score, 1) + T(eventtime, 0.25, status))
#' 
#' BT <- BuyseTest(f, data=df.data)
#' summary(BT)
#' 
#' #### real example : Veteran dataset of the survival package ####
#' #### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#' #### method.tte = "Gehan"
#' 
#' if(require(survival)){
#' \dontrun{
#'   data(veteran,package="survival")
#'  
#'   ## method.tte = "Gehan"
#'   BT_Gehan <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
#'                         data=veteran, method.tte="Gehan",
#'                         method.inference = "permutation", n.resampling = 1e3)
#'   
#'   summary_Gehan <- summary(BT_Gehan)
#'   summary_Gehan <- summary(BT_Gehan, statistic = "winRatio")
#'   
#'   ## method.tte = "Peron"
#'   BT_Peron <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
#'                         data=veteran, method.tte="Peron",
#'                         method.inference = "permutation", n.resampling = 1e3)
#' 
#'   class(BT_Peron)
#'   summary(BT_Peron)
#' }
#' }

## * BuyseTest (code)
##' @rdname BuyseTest
##' @export
BuyseTest <- function(formula,
                      data,
                      method.tte = NULL,
                      model.tte = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      neutral.as.uninf = NULL,
                      keep.pairScore = NULL,
                      treatment = NULL,
                      endpoint = NULL,
                      type = NULL,
                      threshold = NULL,
                      censoring = NULL,
                      operator = NULL,
                      strata = NULL, 
                      alternative = "two.sided", 
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL,
                      keep.comparison){

    name.call <- names(match.call())
    option <- BuyseTest.options()
    
    ## ** compatibility with previous version
    if(!missing(keep.comparison)){
        stop("Argument \'keep.comparison\' is obsolete. \n",
             "It has been replaced by the argument \'keep.pairScore\' \n")
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
                              keep.pairScore = keep.pairScore,
                              method.tte = method.tte,
                              model.tte = model.tte,
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

    if(outArgs$keep.pairScore || outArgs$method.inference == "asymptotic"){
        keep.pairScore <- TRUE
    }else{
        keep.pairScore <- FALSE
    }

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = outArgs)
    }

    ## ** initialization data
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    outArgs[c("data","level.treatment","level.strata","n.strata","n.obs","n.obsStrata","allstrata")] <- initializeData(data = outArgs$data,
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

    ## ** Display
    if (outArgs$trace > 1) {
        outPrint <- do.call(printGeneral, args = outArgs)
    }

    ## ** define environment
    envirBT <- environment()
    envirBT$.BuyseTest <- .BuyseTest
    envirBT$initializeData <- initializeData
    envirBT$initializeSurvival_Peron <- initializeSurvival_Peron
    
    ## ** Point estimation
    if (outArgs$trace > 1) {cat("Point estimation ")}
    time <- system.time({
        outPoint <- .BuyseTest(envir = envirBT,
                               keep.pairScore = keep.pairScore,
                               method.inference = "none")
    })
    ## convert from a list of vector (output of C++) to a list of data.table
    if(keep.pairScore){
        ## needed for inference
        envirBT$indexT <- which(outArgs$data[[outArgs$treatment]]==1)
        envirBT$indexC <- which(outArgs$data[[outArgs$treatment]]==0)

        outPoint$tablePairScore <- pairScore2dt(outPoint$tableScore,
                                                            correction.tte = outArgs$correction.tte,
                                                            level.treatment = outArgs$level.treatment,
                                                            level.strata = outArgs$level.strata,
                                                            n.strata = outArgs$n.strata,
                                                            endpoint = outArgs$endpoint,
                                                            threshold = outArgs$threshold,
                                                            indexT = envirBT$indexT,
                                                            indexC = envirBT$indexC)
    }
    
    if (outArgs$trace > 1) {cat("(done) \n")}
    
    ## ** Permutation test
    if (outArgs$method.inference %in% c("permutation","bootstrap","stratified permutation", "stratified bootstrap")) {
          
        ## ** run
        outResampling <- inferenceResampling(envirBT)
        outCovariance <- matrix(nrow = 0, ncol = 0)
        
    }else if(outArgs$method.inference %in% c("asymptotic")){

        outCovariance <- inferenceUstatistic(envirBT)

        outResampling <- list(deltaResampling.netChance = array(dim=c(0,0,0)),
                              deltaResampling.winRatio = array(dim=c(0,0,0)),
                              DeltaResampling.netChance = matrix(NA, nrow = 0, ncol = 0),
                              DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
                              n.resampling = as.double(NA))
    }else{
        outResampling <- list(deltaResampling.netChance = array(dim=c(0,0,0)),
                              deltaResampling.winRatio = array(dim=c(0,0,0)),
                              DeltaResampling.netChance = matrix(NA, nrow = 0, ncol = 0),
                              DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
                              n.resampling = as.double(NA))
        outCovariance <- matrix(nrow = 0, ncol = 0)
    }
    
    ## ** Gather results into a BuyseRes object
    method.tte <- c("Gehan","Peron")[outArgs$method.tte+1]
    type <- c("Binary","Continuous","TimeToEvent")[outArgs$type]

    BuyseRes.object <- BuyseRes(
        count.favorable = outPoint$count_favorable,      
        count.unfavorable = outPoint$count_unfavorable,
        count.neutral = outPoint$count_neutral,    
        count.uninf = outPoint$count_uninf,
        n.pairs = outPoint$n_pairs,
        delta.netChance = outPoint$delta_netChance,
        delta.winRatio = outPoint$delta_winRatio,
        Delta.netChance = outPoint$Delta_netChance,
        Delta.winRatio = outPoint$Delta_winRatio,
        type = type,
        endpoint = outArgs$endpoint,
        level.treatment = outArgs$level.treatment,
        method.tte = data.frame(method = method.tte, correction = outArgs$correction.tte, stringsAsFactors = FALSE),
        method.inference = outArgs$method.inference,
        strata = outArgs$strata,
        level.strata = outArgs$level.strata,
        threshold = outArgs$threshold,
        n.resampling = outArgs$n.resampling,
        deltaResampling.netChance = outResampling$deltaResampling.netChance,
        deltaResampling.winRatio = outResampling$deltaResampling.winRatio,
        DeltaResampling.netChance = outResampling$DeltaResampling.netChance,
        DeltaResampling.winRatio = outResampling$DeltaResampling.winRatio,
        covariance = outCovariance,
        tablePairScore = if(outArgs$keep.pairScore){outPoint$tablePairScore}else{list()},
        tableSurvival = if(outArgs$keep.survival){outPoint$tableSurvival}else{list()}
    )

    ## ** export
    return(BuyseRes.object)
}

## * .BuyseTest (code)
.BuyseTest <- function(envir,
                       keep.pairScore,
                       method.inference){
   
    strata <- envir$outArgs$strata ## for by in data.table otherwise it cant find what it is
    n.strata <- envir$outArgs$n.strata ## to simplify code
    treatment <- envir$outArgs$treatment ## to simplify code
    censoring <- envir$outArgs$censoring ## to simplify code
    endpoint <- envir$outArgs$endpoint ## to simplify code
    method.tte <- envir$outArgs$method.tte ## to simplify code
    type <- envir$outArgs$type ## to simplify code
    D.TTE <- envir$outArgs$D.TTE ## to simplify code
    D <- envir$outArgs$D ## to simplify code
    
    ## ** Resampling
    if(method.inference == "none"){
        data <- envir$outArgs$data
    }else {
        if(method.inference == "permutation"){
            ## permute the treatment variable over all strata
            data <- data.table::copy(envir$outArgs$data)
            data[[treatment]] <- data[[treatment]][sample.int(envir$outArgs$n.obs, replace = FALSE)]
        }else if(method.inference == "bootstrap"){
            ## randomly pick observations over all strata
            data <- envir$outArgs$data[sample.int(envir$outArgs$n.obs, replace = TRUE)]    
        }else if(method.inference == "stratified permutation"){
            ## permute the treatment variable within each strata
            data <- data.table::copy(envir$outArgs$data)
            data[, c(treatment) := .SD[[1]][sample.int(.N, size = .N, replace = FALSE)], .SDcols = treatment, by = strata]
        }else if(method.inference == "stratified bootstrap"){
            ## randomly pick observations within each
            data <- envir$outArgs$data[,.SD[sample.int(.N, size = .N, replace = TRUE)], by = strata]
        }
    }

    ## ** Check valid resampling
    if(is.null(strata)){
        n.groups <- length(unique(data[[treatment]]))
    }else{
        n.groups <- sum(data[,.N, by = c(".allStrata",treatment)][["N"]]>0)
    }
    if (n.groups != 2*n.strata) { ## failure of the resampling
        ##        return(matrix(NA, nrow = n.strata + 1, ncol = 2*D))
        return(NULL)
    }
    
    ## ** Initialize data    

    ## *** data: split the data according to the two levels
    indexT <- which(data[[treatment]] == 1)
    dataT <- data[indexT]
    dataC <- data[setdiff(1:.N,indexT)]
    
    ## *** data: extract endpoint 
    M.Treatment <- as.matrix(dataT[,endpoint,with=FALSE]) # matrix of endpoints for the treatment arm 
    M.Control <- as.matrix(dataC[,endpoint,with=FALSE]) # matrix of endpoints for the control arm

    ## *** strata
    if(!is.null(strata)){
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
        index.strataT <- lapply(1:n.strata,function(iS){
            which(dataT[[".allStrata"]] == iS) - 1
        })
        index.strataC <- lapply(1:n.strata,function(iS){
            which(dataC[[".allStrata"]] == iS) - 1
        })
        
    }else{ # if there is no strata variable the same strata is used for all patient
        ## For each strata, the index of the patients belonging to each strata, by treatment arm
        ## Index begins at 0. This is compulsory for C++.
        index.strataT <- list(0:(NROW(dataT)-1))
        index.strataC <- list(0:(NROW(dataC)-1))
    }

    ## *** data: censoring
    if(!is.null(censoring)){
        M.delta.Treatment <- as.matrix(dataT[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
        M.delta.Control <- as.matrix(dataC[,censoring,with=FALSE]) # matrix of censoring variables for the treatment arm : censored (0) event time (1)
    }else{ # if the is no time to event variables
        M.delta.Treatment <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
        M.delta.Control <- matrix(nrow=0,ncol=0) # factice censoring matrix. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }

    ## *** Update survival
    if(method.tte == 0){ ## Gehan
        outSurv <- list(survTimeC = list(list(matrix())),
                        survTimeT = list(list(matrix())),
                        survJumpC = list(list(matrix())),
                        survJumpT = list(list(matrix()))
                        )        
    }else{ ## Peron
        outSurv <- initializeSurvival_Peron(data =  data, dataT = dataT, dataC = dataC,
                                            model.tte = envir$outArgs$model.tte,
                                            n.T = NROW(M.Treatment), n.C = NROW(M.Control),
                                            treatment = treatment,
                                            level.treatment = envir$outArgs$level.treatment,
                                            endpoint = endpoint,
                                            censoring = censoring,
                                            D.TTE = D.TTE,
                                            type = type,
                                            strata = strata,
                                            threshold = envir$outArgs$threshold,
                                            index.strataT = index.strataT,
                                            index.strataC = index.strataC,
                                            n.strata = n.strata)
    }

    ## ** Computation
    resBT <- GPC_cpp(Control = M.Control,
                     Treatment = M.Treatment,
                     threshold = envir$outArgs$threshold,
                     survEndpoint = (type == 3),
                     delta_Control = M.delta.Control,
                     delta_Treatment = M.delta.Treatment,
                     D = D,
                     strataC = index.strataC,
                     strataT = index.strataT,
                     n_strata = n.strata,
                     n_TTE = D.TTE,
                     Wscheme = envir$outArgs$Wscheme,
                     index_survivalM1 = envir$outArgs$index.survivalM1,
                     threshold_TTEM1 = envir$outArgs$threshold.TTEM1,
                     list_survTimeC = outSurv$survTimeC,
                     list_survTimeT = outSurv$survTimeT,
                     list_survJumpC = outSurv$survJumpC,
                     list_survJumpT = outSurv$survJumpT,
                     list_lastSurv = outSurv$lastSurv,
                     methodTTE = method.tte,
                     correctionTTE = envir$outArgs$correction.tte,
                     neutralAsUninf = envir$outArgs$neutral.as.uninf,
                     keepScore = keep.pairScore
                     )
    
    ## ** export
    if(method.inference == "none"){
        if(envir$outArgs$keep.survival){ ## useful to test initSurvival 
            resBT$tableSurvival <- outSurv
        }
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







