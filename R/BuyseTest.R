## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Generalized Pairwise Comparisons (GPC)
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event endpoints.
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
#' @param threshold [numeric vector] critical values used to compare the pairs (threshold of minimal important difference).
#' There must be one threshold for each endpoint variable.
#' Disregarded if the argument \code{formula} is defined.
#' @param strata [numeric vector] if not \code{NULL}, the GPC will be applied within each group of patient defined by the strata variable(s).
#' Disregarded if the argument \code{formula} is defined.
#' @param censoring [character vector] the name of the binary variable(s) indicating whether the endpoint was observed or censored.
#' There must be one threshold for each endpoint variable.
#' Must value \code{NA} when the endpoint is not a time to event.
#' Disregarded if the argument \code{formula} is defined.
#' @param type [character vector] the type of each endpoint: \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.
#' @param method.tte [character] defines the method used to compare the observations of a pair in presence of right censoring (i.e. \code{"timeToEvent"} endpoints).
#' Can be \code{"Gehan"} or \code{"Peron"}. \code{"Gehan"} only scores pairs that can be decidedly classified as favorable, unfavorable, or neutral.
#' \code{"Peron"} uses the empirical survival curves of each group to also score the pairs that cannot be decidedly classified (see Peron et al. for more details).
#' Default value read from \code{BuyseTest.options()}.
#' @param correction.uninf [integer] should a correction be applied to remove the bias due to the presence of uninformative pairs?
#' 0 indicates no correction, 1 impute the average score of the informative pair, and 2 performs inverse probability of censoring weights.
#' Default value read from \code{BuyseTest.options()}.
#' @param model.tte [list] optional survival models relative to each time to each time to event endpoint.
#' Models must \code{prodlim} objects and stratified on the treatment and strata variable.
#' @param method.inference [character] should the asymptotic theory (\code{"asymptotic"}),
#' or a permutation test (\code{"permutation"} or \code{"stratified permutation"}),
#' or bootstrap resampling (\code{"bootstrap"} or \code{"stratified boostrap"})
#' be used to compute p-values and confidence intervals.
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analyzed using endpoints of lower priority.
#' Default value read from \code{BuyseTest.options()}.
#' @param n.resampling [integer] the number of simulations used for computing the confidence interval and the p.values. See details.
#' Default value read from \code{BuyseTest.options()}.
#' @param keep.pairScore [logical] should the result of each pairwise comparison be kept?
#' @param alternative [character] the alternative hypothesis.
#' Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
#' Default value read from \code{BuyseTest.options()}.
#' @param seed [integer, >0] the seed to consider for the permutation test.
#' If \code{NULL} no seed is set.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' Default value read from \code{BuyseTest.options()}.
#' @param trace [integer] should the execution of the function be traced ? \code{0} remains silent
#' and \code{1}-\code{3} correspond to a more and more verbose output in the console.
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
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus.
#' The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal endpoint.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevent from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analyzed on the endpoint with priority \code{l-1}.
#' 
#' \bold{method.tte:} the \code{method.tte="Peron"} is recommended in presence of right censored observations since it gives a more efficient estimator than \code{method.tte="Gehan"}.
#' 
#' \bold{method.inference:} the \code{method.inference="asymptotic"} estimate the distribution of the net benefit or win ratio statistics
#' based on the asymptotic theory of the U-statistics (see formula 2.2 in Bebu et al. 2016).
#' Its current implementation is not valid in small sample or when using \code{method.tte="Peron"}, or  \code{correction.uninf=1}, or \code{correction.uninf=2}.
#'
#' \bold{correction.uninf:} in presence of uninformative pairs, the proportion of favorable, unfavorable, or neutral pairs is underestimated.
#' Inverse probability of censoring weights (\code{correction.uninf=2}) is only recommanded when the analysis is stopped after the first endpoint with uninformative pairs.
#' Imputing the average score of the informative pairs (\code{correction.uninf=1}) gives equivalent results at the first endpoint but better behaves at latter endpoints.
#' Note that both corrections will convert the whole proportion of uninformative pairs of a given endpoint into favorable, unfavorable, or neutral pairs.
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references
#' 
#' J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2018). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. \emph{Statistical Methods in Medical Research} 27: 1230-1239  \cr 
#' D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' I. Bebu, J. M. Lachin. \bold{Large sample inference for a win ratio analysis of a composite outcome based on prioritized components}. \emph{Biostatistics} 17(1):178-187 \cr
#' Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
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
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#'                                        # display 
#' if(require(prodlim)){
#'    resKM_tempo <- prodlim(Hist(eventtime,status)~Treatment, data = df.data)
#'    plot(resKM_tempo)
#' }
#'
#' #### one time to event endpoint ####
#' BT <- BuyseTest(Treatment ~ TTE(eventtime, censoring = status), data= df.data)
#'
#' summary(BT) # net benefit
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
#' summary(BT, statistic = "netBenefit") ## default
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
                      correction.uninf = NULL,
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
                      alternative = NULL, 
                      seed = 10,
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
                              correction.uninf = correction.uninf,
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
    outArgs[c("data","M.endpoint","M.censoring","level.treatment","level.strata","n.strata","n.obs","n.obsStrata")] <- initializeData(data = outArgs$data,
                                                                                                                                      type = outArgs$type,
                                                                                                                                      method.tte = outArgs$method.tte,
                                                                                                                                      endpoint = outArgs$endpoint,
                                                                                                                                      censoring = outArgs$censoring,
                                                                                                                                      operator = outArgs$operator,
                                                                                                                                      strata = outArgs$strata,
                                                                                                                                      treatment = outArgs$treatment)
    
    ## ** create weights matrix for survival endpoints
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    outArgs[c("Wscheme","index.survival_M1","threshold.TTE_M1","outSurv")] <- buildWscheme(method.tte = outArgs$method.tte,
                                                                                           endpoint = outArgs$endpoint,
                                                                                           D = outArgs$D,
                                                                                           D.TTE = outArgs$D.TTE,
                                                                                           n.strata = outArgs$n.strata,
                                                                                           type = outArgs$type,
                                                                                           threshold = outArgs$threshold)
    
    ## ** Display
    if (outArgs$trace > 1) {
        cat("\n         Generalized Pairwise Comparisons\n\n")
        do.call(printGeneral, args = outArgs)
        cat("\n")
    }

    ## ** define environment
    envirBT <- environment()
    envirBT$.BuyseTest <- .BuyseTest
    envirBT$initializeData <- initializeData
    envirBT$initializeSurvival_Peron <- initializeSurvival_Peron
    
    ## ** Point estimation
    if (outArgs$trace > 1) {
        cat("Point estimation")
    }
    outPoint <- .BuyseTest(envir = envirBT,
                           keep.pairScore = keep.pairScore,
                           method.inference = "none")
    
    if (outArgs$trace > 1) {
        cat("\n\n")
    }
    
    ## convert from a list of vector (output of C++) to a list of data.table
    if(keep.pairScore){
        ## needed for inference
        envirBT$indexT <- which(outArgs$data[[outArgs$treatment]]==1)
        envirBT$indexC <- which(outArgs$data[[outArgs$treatment]]==0)

        outPoint$tablePairScore <- pairScore2dt(outPoint$tableScore,
                                                level.treatment = outArgs$level.treatment,
                                                level.strata = outArgs$level.strata,
                                                n.strata = outArgs$n.strata,
                                                endpoint = outArgs$endpoint,
                                                threshold = outArgs$threshold,
                                                indexT = envirBT$indexT,
                                                indexC = envirBT$indexC)
    }
    
    ## ** Inference
    if (outArgs$method.inference!="none") {
        if(outArgs$trace > 1){
            do.call(printInference, args = outArgs)
        }

        if(outArgs$method.inference %in% c("asymptotic")){

            if(outArgs$method.tte > 0){
                warning("The current implementation of the asymptotic distribution is not valid for method.tte=\"Peron\" \n",
                        "Standard errors / confidence intervals / p-values should not be trusted \n")
            }
            if(outArgs$correction.uninf > 0){
                warning("The current implementation of the asymptotic distribution is not valid when a correction is used \n",
                        "Standard errors / confidence intervals / p-values should not be trusted \n")
            }
            
            outCovariance <- inferenceUstatistic(tablePairScore = outPoint$tablePairScore,
                                                 count.favorable = outPoint$count_favorable, count.unfavorable = outPoint$count_unfavorable,
                                                 n.pairs = outPoint$n_pairs, n.C = length(envirBT$indexC), n.T = length(envirBT$indexC),                                
                                                 n.strata = outArgs$n.strata, endpoint = outArgs$endpoint)$Sigma

            outResampling <- list(deltaResampling.netBenefit = array(dim=c(0,0,0)),
                                  deltaResampling.winRatio = array(dim=c(0,0,0)),
                                  DeltaResampling.netBenefit = matrix(NA, nrow = 0, ncol = 0),
                                  DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
                                  n.resampling = as.double(NA))
        }else{
            outResampling <- inferenceResampling(envirBT)
            outCovariance <- matrix(nrow = 0, ncol = 0)
        }
        if(outArgs$trace > 1){
            cat("\n")
        }
        
    }else{
        outResampling <- list(deltaResampling.netBenefit = array(dim=c(0,0,0)),
                              deltaResampling.winRatio = array(dim=c(0,0,0)),
                              DeltaResampling.netBenefit = matrix(NA, nrow = 0, ncol = 0),
                              DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
                              n.resampling = as.double(NA))
        outCovariance <- matrix(nrow = 0, ncol = 0)
    }
    
    ## ** Gather results into a BuyseRes object
    if(outArgs$trace > 1){
        cat("Gather the results in a BuyseRes object \n")
    }
    method.tte <- c("Gehan","Peron")[outArgs$method.tte+1]
    type <- c("Binary","Continuous","TimeToEvent")[outArgs$type]

    BuyseRes.object <- BuyseRes(
        count.favorable = outPoint$count_favorable,      
        count.unfavorable = outPoint$count_unfavorable,
        count.neutral = outPoint$count_neutral,    
        count.uninf = outPoint$count_uninf,
        n.pairs = outPoint$n_pairs,
        delta.netBenefit = outPoint$delta_netBenefit,
        delta.winRatio = outPoint$delta_winRatio,
        Delta.netBenefit = outPoint$Delta_netBenefit,
        Delta.winRatio = outPoint$Delta_winRatio,
        type = type,
        endpoint = outArgs$endpoint,
        level.treatment = outArgs$level.treatment,
        method.tte = method.tte,
        correction.uninf = outArgs$correction.uninf,
        method.inference = outArgs$method.inference,
        strata = outArgs$strata,
        level.strata = outArgs$level.strata,
        threshold = outArgs$threshold,
        n.resampling = outArgs$n.resampling,
        deltaResampling.netBenefit = outResampling$deltaResampling.netBenefit,
        deltaResampling.winRatio = outResampling$deltaResampling.winRatio,
        DeltaResampling.netBenefit = outResampling$DeltaResampling.netBenefit,
        DeltaResampling.winRatio = outResampling$DeltaResampling.winRatio,
        covariance = outCovariance,
        tablePairScore = if(outArgs$keep.pairScore){outPoint$tablePairScore}else{list()},
        tableSurvival = if(outArgs$keep.survival){outPoint$tableSurvival}else{list()}
    )
    if(outArgs$trace > 1){
        cat("\n")
    }

    ## ** export
    return(BuyseRes.object)
}

## * .BuyseTest (code)
.BuyseTest <- function(envir,
                       keep.pairScore,
                       method.inference){
   
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
            data <- data.table::copy(envir$outArgs$data)
            data[ ,c(treatment) := .SD[[1]][sample.int(.N, replace = FALSE)], .SDcols = treatment]            
        }else if(method.inference == "bootstrap"){
            ## randomly pick observations over all strata
            data <- envir$outArgs$data[sample.int(envir$outArgs$n.obs, replace = TRUE)]            
        }else if(method.inference == "stratified permutation"){
            ## permute the treatment variable within each strata
            data <- data.table::copy(envir$outArgs$data)
            data[, c(treatment) := .SD[[1]][sample.int(.N, replace = FALSE)], by = "..strata..", .SDcols = treatment]
        }else if(method.inference == "stratified bootstrap"){
            ## randomly pick observations within each strata
            data <- envir$outArgs$data[,.SD[sample.int(.N, replace = TRUE)], by = "..strata.."]
        }        
    }
    
    ## ** split data
    if(n.strata==1){        
        ls.indexC <- list(which(data[[treatment]]==0)-1)
        ls.indexT <- list(which(data[[treatment]]==1)-1)
    }else{        
        indexC <- which(data[[treatment]]==0)-1
        indexT <- which(data[[treatment]]==1)-1
        indexStrata <- lapply(1:n.strata, function(iS){
            which(data[["..strata.."]]==iS)-1
        })
        
        ls.indexC <- vector(mode = "list", length = n.strata)
        ls.indexT <- vector(mode = "list", length = n.strata)
        for(iStrata in 1:n.strata){ ## iStrata <- 1  
            ls.indexC[[iStrata]] <- intersect(indexC, indexStrata[[iStrata]])
            ls.indexT[[iStrata]] <- intersect(indexT, indexStrata[[iStrata]])
        }
    }
    ## table(duplicated(c(unlist(ls.indexC),unlist(ls.indexT))))
    ## sort(c(unlist(ls.indexC),unlist(ls.indexT)))
    ## data

    ## Check valid resampling
    if (any(c(sapply(ls.indexC,length),sapply(ls.indexT,length))==0)) {
        return(NULL)
    }

    ## *** Update survival
    if(method.tte == 0){ ## Gehan
        outSurv <- envir$outArgs$outSurv

    }else{ ## Peron
        outSurv <- initializeSurvival_Peron(data = data, ls.indexC = ls.indexC, ls.indexT = ls.indexT,
                                            model.tte = envir$outArgs$model.tte,
                                            index.survival_M1 = envir$outArgs$index.survival_M1,
                                            treatment = treatment,
                                            level.treatment = envir$outArgs$level.treatment,
                                            endpoint = endpoint,
                                            censoring = censoring,
                                            D.TTE = D.TTE,
                                            type = type,
                                            threshold = envir$outArgs$threshold,
                                            n.strata = n.strata,
                                            strata = envir$outArgs$strata,
                                            out = envir$outArgs$outSurv)
    }

    ## ** Computation
    resBT <- GPC_cpp(endpoint = envir$outArgs$M.endpoint[data[["..rowIndex.."]],,drop=FALSE],
                     censoring = envir$outArgs$M.censoring[data[["..rowIndex.."]],,drop=FALSE],
                     indexC = ls.indexC,
                     indexT = ls.indexT,                     
                     threshold = envir$outArgs$threshold,
                     method = envir$outArgs$method.score,
                     D = D,
                     n_strata = n.strata,
                     n_TTE = D.TTE,
                     Wscheme = envir$outArgs$Wscheme,
                     index_survival_M1 = envir$outArgs$index.survival_M1,
                     threshold_M1 = envir$outArgs$threshold.TTE_M1,
                     list_survTimeC = outSurv$survTimeC,
                     list_survTimeT = outSurv$survTimeT,
                     list_survJumpC = outSurv$survJumpC,
                     list_survJumpT = outSurv$survJumpT,
                     list_lastSurv = outSurv$lastSurv,
                     correctionUninf = envir$outArgs$correction.uninf,
                     neutralAsUninf = envir$outArgs$neutral.as.uninf,
                     keepScore = keep.pairScore,
                     reserve = TRUE
                     )

    ## ** export
    if(method.inference == "none"){
        if(envir$outArgs$keep.survival){ ## useful to test initSurvival 
            resBT$tableSurvival <- outSurv
        }
        return(resBT)
    }else{
        return(list(delta_netBenefit = resBT$delta_netBenefit,
                    Delta_netBenefit = resBT$Delta_netBenefit,
                    delta_winRatio = resBT$delta_winRatio,
                    Delta_winRatio = resBT$Delta_winRatio))
}
}







