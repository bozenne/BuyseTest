## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Generalized Pairwise Comparisons (GPC)
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event endpoints.
#' @param formula [formula] a symbolic description of the GPC model,
#' typically \code{treatment ~ type1(endpoint1) + type2(endpoint2, threshold2) + strata}.
#' See Details, section "Specification of the GPC model".
#' @param treatment,endpoint,type,threshold,censoring,operator,strata Alternative to \code{formula} for describing the GPC model.
#' See Details, section "Specification of the GPC model".
#' @param data [data.frame] dataset.
#' @param scoring.rule [character] method used to compare the observations of a pair in presence of right censoring (i.e. \code{"timeToEvent"} endpoints).
#' Can be \code{"Gehan"} or \code{"Peron"}.
#' See Details, section "Handling missing values".
#' @param correction.uninf [integer] should a correction be applied to remove the bias due to the presence of uninformative pairs?
#' 0 indicates no correction, 1 impute the average score of the informative pairs, and 2 performs IPCW.
#' See Details, section "Handling missing values".
#' @param model.tte [list] optional survival models relative to each time to each time to event endpoint.
#' Models must \code{prodlim} objects and stratified on the treatment and strata variable. When used, the uncertainty from the estimates of these survival models is ignored.
#' @param method.inference [character] method used to compute confidence intervals and p-values.
#' Can be \code{"none"}, \code{"u-statistic"}, \code{"permutation"}, \code{"stratified permutation"}, \code{"bootstrap"}, \code{"stratified bootstrap"}, \code{"studentized bootstrap"}, or \code{"studentized stratified bootstrap"}.
#' See Details, section "Statistical inference".
#' @param n.resampling [integer] the number of permutations/samples used for computing the confidence intervals and the p.values. 
#' See Details, section "Statistical inference".
#' @param hierarchical [logical] should only the uninformative pairs be analyzed at the lower priority endpoints (hierarchical GPC)?
#' Otherwise all pairs will be compaired for all endpoint (full GPC).
#' @param weight [numeric vector] weights used to cumulating the pairwise scores over the endpoints.
#' Only used when \code{hierarchical=FALSE}. Disregarded if the argument \code{formula} is defined.
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analyzed using endpoints of lower priority (as it is done for uninformative pairs).
#' See Details, section "Handling missing values".
#' @param keep.pairScore [logical] should the result of each pairwise comparison be kept?
#' @param seed [integer, >0] the seed to consider for the permutation test.
#' If \code{NULL} no seed is set.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' See Details, section "Statistical inference".
#' @param trace [integer] should the execution of the function be traced ? \code{0} remains silent
#' and \code{1}-\code{3} correspond to a more and more verbose output in the console.
#' @param keep.comparison Obsolete. Alias for 'keep.pairScore'.
#' @param method.tte Obsolete. Alias for 'scoring.rule'.
#' 
#' @details
#'
#' \bold{Specification of the GPC model}: \cr
#' There are two way to specify the GPC model in \code{BuyseTest}.
#' A \emph{Formula interface} via the argument \code{formula} where the response variable should be a binary variable defining the treatment arms. 
#' The rest of the formula should indicate the endpoints by order of priority and the strata variables (if any).
#' A \emph{Vector interface} using  the following arguments \itemize{
#'   \item \code{treatment}: [character] name of the treatment variable identifying the control and the experimental group.
#' Must have only two levels (e.g. \code{0} and \code{1}).
#'   \item \code{endpoint}: [character vector] the name of the endpoint variable(s).
#'   \item \code{threshold}: [numeric vector] critical values used to compare the pairs (threshold of minimal important difference).
#' There must be one threshold for each endpoint variable; it must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. 
#'   \item \code{censoring}: [character vector] the name of the binary variable(s) indicating whether the endpoint was observed or censored.
#' Must value \code{NA} when the endpoint is not a time to event.
#'   \item \code{operator}: [character vector] the sign defining a favorable endpoint.
#' \code{">0"} indicates that higher values are favorable while "<0" indicates the opposite.
#' When the operator is set to \code{"<0"} the corresponding column in the dataset is internally multiplied by \code{-1}.#' 
#'   \item \code{type}: [character vector] indicates whether it is
#' a binary outcome  (\code{"b"}, \code{"bin"}, or \code{"binary"}),
#' a continuous outcome  (\code{"c"}, \code{"cont"}, or \code{"continuous"}),
#' or a time to event outcome  (\code{"t"}, \code{"tte"}, \code{"time"}, or \code{"timetoevent"})
#'   \item \code{strata}: [character vector] if not \code{NULL}, the GPC will be applied within each group of patient defined by the strata variable(s).
#' }
#' The formula interface can be more concise, especially when considering few outcomes, but may be more difficult to apprehend for new users.
#' Note that arguments \code{endpoint}, \code{threshold}, \code{censoring}, \code{operator}, and \code{type} must have the same length. \cr \cr \cr
#'
#' 
#' \bold{GPC procedure} \cr
#' The GPC procedure form all pairs of observations, one belonging to the experimental group and the other to the control group, and class them in 4 categories: \itemize{
#'  \item \emph{Favorable pair}: the endpoint is better for the observation in the experimental group.
#'  \item \emph{Unfavorable pair}: the endpoint is better for the observation in the control group.
#'  \item \emph{Neutral pair}: the difference between the endpoints of the two observations is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal endpoint. Lower-priority outcomes (if any) are then used to classified the pair into favorable/unfavorable.
#'  \item \emph{Uninformative pair}: censoring/missingness prevents from classifying into favorable, unfavorable or neutral.
#' }
#' With complete data, pairs can be decidely classified as favorable/unfavorable/neutral.
#' In presence of missing values, the GPC procedure uses the scoring rule (argument \code{scoring.rule}) and the correction for uninformative pairs (argument \code{correction.uninf}) to classify the pairs.
#' The classification may not be 0,1, e.g. the probability that the pair is favorable/unfavorable/neutral with the Peron's scoring rule.
#' To export the classification of each pair set the argument code{keep.pairScore} to \code{TRUE} and call the function \code{getPairScore} on the result of the \code{BuyseTest} function. \cr \cr \cr
#' 
#' 
#' \bold{Handling missing values}
#' \itemize{
#'   \item \code{scoring.rule}: indicates how to handle right-censoring in time to event endpoints.
#' The Gehan's scoring rule (argument \code{scoring.rule="Gehan"}) only scores pairs that can be decidedly classified as favorable, unfavorable, or neutral
#' while the "Peron"'s scoring rule (argument \code{scoring.rule="Peron"}) uses the empirical survival curves of each group to also score the pairs that cannot be decidedly classified.
#' The Peron's scoring rule is the recommanded scoring rule.
#'   \item \code{correction.uninf}: indicates how to handle missing values that could not be classified by the scoring rule. \code{0} treat them as uninformative:
#' if \code{neutral.as.uninf=FALSE}  - this is an equivalent to complete case analysis -
#' while for \code{neutral.as.uninf=TRUE} uninformative pairs are treated as neutral, i.e., analyzed at the following endpoint (if any).
#' However both will (in general) lead to biased estimates for the proportion of favorable, unfavorable, or neutral pairs.
#' Inverse probability of censoring weights (IPCW, \code{correction.uninf=2}) is only recommanded when the analysis is stopped after the first endpoint with uninformative pairs.
#' Imputing the average score of the informative pairs (\code{correction.uninf=1}) is the recommanded approach.
#' Note that both corrections will convert the whole proportion of uninformative pairs of a given endpoint into favorable, unfavorable, or neutral pairs. \cr \cr
#' }
#' 
#'
#' \bold{Statistical inference} \cr
#' The argument \code{method.inference} defines how to approximate the distribution of the GPC estimators and so how standard errors, confidence intervals, and p-values are computed.
#' Available methods are:
#' \itemize{
#'   \item argument \code{method.inference="none"}: only the point estimate is computed which makes the execution of the \code{BuyseTest} faster than with the other methods.
#'   \item argument \code{method.inference="u-statistic"}: uses a Gaussian approximation to obtain the distribution of the GPC estimators.
#' The U-statistic theory indicates that this approximation is asymptotically exact.
#' The variance is computed using a H-projection of order 1 (default option), which is a consistent but downward biased estimator.
#' An unbiased estimator can be obtained using a H-projection of order 2 (only available for the uncorrected Gehan's scoring rule, see \code{BuyseTest.options}).
#' \bold{WARNING}: the current implementation of the H-projection is not valid when using corrections for uninformative pairs (\code{correction.uninf=1}, or \code{correction.uninf=2}) and the Peron scoring rule (\code{scoring.rule="Peron"}).
#'   \item argument \code{method.inference="permutation"} or \code{method.inference="stratified permutation"}:
#'   \item argument \code{method.inference="bootstrap"} or \code{method.inference="stratified bootstrap"}:
#' }
#' Additional arguments for permutation and bootstrap resampling:
#' \itemize{
#'    \item when stratified, permutation and bootstrap are performed separately in each strata level (and not each treatment group).
#' This is therefore only relevant for stratified analyses.
#'    \item \code{n.resampling} set the number of permutations/samples used.
#' A large number of permutations (e.g. \code{n.resampling=10000}) are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details.
#'    \item \code{cpus} indicates whether the resampling procedure can be splitted on several cpus to save time. Can be set to \code{"all"} to use all available cpus.
#' The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package. \cr \cr
#' }
#'
#'
#' \bold{Default values} \cr
#' The default of the arguments
#' \code{scoring.rule}, \code{correction.uninf}, \code{method.inference}, \code{n.resampling},
#' \code{hierarchical}, \code{neutral.as.uninf}, \code{keep.pairScore}, \code{n.resampling},
#' \code{cpus}, \code{trace} is read from \code{BuyseTest.options()}. \cr
#' Additional (hidden) arguments are \itemize{
#'  \item \code{alternative} [character] the alternative hypothesis. Must be one of "two.sided", "greater" or "less" (used by \code{confint}).
#'  \item \code{conf.level} [numeric] level for the confidence intervals (used by \code{confint}).
#'  \item \code{keep.survival} [logical] export the survival values used by the Peron's scoring rule.
#'  \item \code{order.Hprojection} [1 or 2] the order of the H-projection used to compute the variance when \code{method.inference="u-statistic"}. 
#' }
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Peron's scoring rule: J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2018). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. \emph{Statistical Methods in Medical Research} 27: 1230-1239  \cr 
#' On the Gehan's scoring rule: Gehan EA (1965). \bold{A generalized two-sample Wilcoxon test for doubly censored data}. \emph{Biometrika}  52(3):650-653 \cr
#' On inference in GPC using the U-statistic theory: I. Bebu, J. M. Lachin (2015). \bold{Large sample inference for a win ratio analysis of a composite outcome based on prioritized components}. \emph{Biostatistics} 17(1):178-187 \cr
#'
#' @seealso 
#' \code{\link{BuyseRes-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' @keywords function BuyseTest

## * BuyseTest (example)
#' @rdname BuyseTest
#' @examples
#' library(data.table)
#' 
#' # reset the default value of the number of permuation sample
#' BuyseTest.options(method.inference = "none") # no permutation test
#'
#' #### simulate some data ####
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#'                                        # display 
#' if(require(prodlim)){
#'    resKM_tempo <- prodlim(Hist(eventtime,status)~treatment, data = df.data)
#'    plot(resKM_tempo)
#' }
#'
#' #### one time to event endpoint ####
#' BT <- BuyseTest(treatment ~ TTE(eventtime, censoring = status), data= df.data)
#'
#' summary(BT) # net benefit
#' summary(BT, percentage = FALSE)  
#' summary(BT, statistic = "winRatio") # win Ratio
#' 
#' ## bootstrap to compute the CI
#' \dontrun{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3)
#' }
#' \dontshow{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e1, trace = 0)
#' }
#' summary(BT, statistic = "netBenefit") ## default
#' summary(BT, statistic = "winRatio") 
#' 
#' ## parallel bootstrap
#' \dontrun{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3, cpus = 2)
#'     summary(BT)
#' }
#' 
#' ## method Gehan is much faster but does not optimally handle censored observations
#' BT <- BuyseTest(treatment ~ TTE(eventtime, censoring = status), data=df.data,
#'                 scoring.rule = "Gehan", trace = 0)
#' summary(BT)
#' 
#' #### one time to event endpoint: only differences in survival over 1 unit ####
#' BT <- BuyseTest(treatment ~ TTE(eventtime, threshold = 1, censoring = status), data=df.data)
#' summary(BT)
#' 
#' #### one time to event endpoint with a strata variable
#' BT <- BuyseTest(treatment ~ strata + TTE(eventtime, censoring = status), data=df.data)
#' summary(BT)
#' 
#' #### several endpoints with a strata variable
#' f <- treatment ~ strata + T(eventtime, status, 1) + B(toxicity) 
#' f <- update(f, 
#'             ~. + T(eventtime, status, 0.5) + C(score, 1) + T(eventtime, status, 0.25))
#' 
#' BT <- BuyseTest(f, data=df.data)
#' summary(BT)
#' 
#' #### real example : Veteran dataset of the survival package ####
#' #### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#' #### scoring.rule = "Gehan"
#' 
#' if(require(survival)){
#' \dontrun{
#'   data(veteran,package="survival")
#'  
#'   ## scoring.rule = "Gehan"
#'   BT_Gehan <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
#'                         data=veteran, scoring.rule="Gehan",
#'                         method.inference = "permutation", n.resampling = 1e3)
#'   
#'   summary_Gehan <- summary(BT_Gehan)
#'   summary_Gehan <- summary(BT_Gehan, statistic = "winRatio")
#'   
#'   ## scoring.rule = "Peron"
#'   BT_Peron <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,censoring=status), 
#'                         data=veteran, scoring.rule="Peron",
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
                      scoring.rule = NULL,
                      correction.uninf = NULL,
                      model.tte = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      hierarchical = NULL,
                      weight = NULL,
                      neutral.as.uninf = NULL,
                      keep.pairScore = NULL,
                      seed = 10,
                      cpus = NULL,
                      trace = NULL,
                      treatment = NULL,
                      endpoint = NULL,
                      type = NULL,
                      threshold = NULL,                      
                      censoring = NULL,
                      operator = NULL,
                      strata = NULL, 
                      keep.comparison,
                      method.tte){

    name.call <- names(match.call())
    option <- BuyseTest.options()

    ## ** compatibility with previous version
    if(!missing(keep.comparison)){
        stop("Argument \'keep.comparison\' is obsolete. \n",
             "It has been replaced by the argument \'keep.pairScore\' \n")
    }
    if(!missing(method.tte)){
        stop("Argument \'method.tte\' is obsolete. \n",
             "It has been replaced by the argument \'scoring.rule\' \n")
    }
    if(!is.null(method.inference) && (method.inference=="asymptotic")){
        stop("Value \"asymptotic\" for argument \'method.inference\' is obsolete. \n",
             "Use \"u-statistic\" instead \n")
    }

    ## ** initialize arguments (all expect data that is just converted to data.table)
    ## initialized arguments are stored in outArgs
    outArgs <- initializeArgs(censoring = censoring,
                              correction.uninf = correction.uninf,
                              cpus = cpus,
                              data = data,
                              endpoint = endpoint,
                              formula = formula,
                              hierarchical = hierarchical,
                              keep.pairScore = keep.pairScore,
                              method.inference = method.inference,
                              scoring.rule = scoring.rule,
                              model.tte = model.tte,
                              n.resampling = n.resampling,
                              name.call = name.call,
                              neutral.as.uninf = neutral.as.uninf,
                              operator = operator,
                              option = option,
                              seed = seed,
                              strata = strata,
                              threshold = threshold,
                              trace = trace,
                              treatment = treatment,
                              type = type,
                              weight = weight)

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = outArgs)
        if((attr(outArgs$method.inference,"ustatistic") || attr(outArgs$method.inference,"studentized")) & outArgs$scoring.rule==1 & length(outArgs$endpoint)>1){
            warning("Computation of the standard error with the Peron scoring rule for several endpoints is still under development\n",
                    "Will probably fail or give incorrect results\n",
                    "Wait future version of the package (aim: december 2020)\n.")
        }
    }

    ## ** initialization data
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    out.name <- c("data","M.endpoint","M.censoring",
                  "index.C","index.T","index.strata",
                  "level.treatment","level.strata", "method.score",
                  "n.strata","n.obs","n.obsStrata","cumn.obsStrata","skeletonPeron")
    outArgs[out.name] <- initializeData(data = outArgs$data,
                                        type = outArgs$type,
                                        endpoint = outArgs$endpoint,
                                        Uendpoint = outArgs$Uendpoint,
                                        D = outArgs$D,
                                        scoring.rule = outArgs$scoring.rule,
                                        censoring = outArgs$censoring,
                                        Ucensoring = outArgs$Ucensoring,
                                        operator = outArgs$operator,
                                        strata = outArgs$strata,
                                        treatment = outArgs$treatment,
                                        copy = TRUE)
    
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
    envirBT$initializePeron <- initializePeron
    
    ## ** Point estimation
    if (outArgs$trace > 1) {
        cat("Point estimation")
    }
    outPoint <- .BuyseTest(envir = envirBT,
                           iid = outArgs$iid,
                           method.inference = "none",
                           pointEstimation = TRUE)

    if (outArgs$trace > 1) {
        cat("\n\n")
    }

    ## check number of pairs
    if(option$check){
        vec.nPair <- (outPoint$count_favorable + outPoint$count_unfavorable + outPoint$count_neutral + outPoint$count_uninf )[,1]
        if(any(abs(outPoint$n_pairs - vec.nPair) > 0.01)){
            warning("Incorrect estimation of the number of pairs \n",
                    "Something probably went wrong - contact the package maintainer\n")
        }
    }

    
    ## convert from a list of vector (output of C++) to a list of data.table
    if(outArgs$keep.pairScore){
        ## needed for inference with bebu
        outPoint$tablePairScore <- pairScore2dt(outPoint$tableScore,
                                                level.treatment = outArgs$level.treatment,
                                                level.strata = outArgs$level.strata,
                                                n.strata = outArgs$n.strata,
                                                endpoint = outArgs$endpoint,
                                                threshold = outArgs$threshold)
    }
    
    ## ** Inference
    if((outArgs$method.inference != "none") && (outArgs$trace > 1)){
        do.call(printInference, args = outArgs)
    }

    outResampling <- list(deltaResampling.netBenefit = array(dim=c(0,0,0)),
                          deltaResampling.winRatio = array(dim=c(0,0,0)),
                          DeltaResampling.netBenefit = matrix(NA, nrow = 0, ncol = 0),
                          DeltaResampling.winRatio = matrix(NA, nrow = 0, ncol = 0),
                          covariance = array(NA, dim = c(0,0,0)),
                          n.resampling = as.double(NA))

    if(outArgs$method.inference == "none"){
        outPoint$Mvar <- matrix(nrow = 0, ncol = 0)
        outPoint$iid_favorable <- NULL
        outPoint$iid_unfavorable <- NULL
    }else if(outArgs$method.inference == "u-statistic"){
        ## done in the C++ code
        ## outCovariance <- inferenceUstatistic(tablePairScore = outPoint$tablePairScore, order = option$order.Hprojection,
        ##                                      count.favorable = colSums(outPoint$count_favorable), count.unfavorable = colSums(outPoint$count_unfavorable),
        ##                                      n.pairs = sum(outPoint$n_pairs), n.C = length(envirBT$outArgs$index.C), n.T = length(envirBT$outArgs$index.T),
        ##                                      level.strata = outArgs$level.strata, n.strata = outArgs$n.strata, endpoint = outArgs$endpoint)
    }else if(outArgs$method.inference == "u-statistic-bebu"){
        if(outArgs$keep.pairScore == FALSE){
            stop("Argument \'keep.pairScore\' needs to be TRUE when argument \'method.inference\' is \"u-statistic-bebu\" \n")
        }

        ## direct computation of the variance
        outCovariance <- inferenceUstatisticBebu(tablePairScore = outPoint$tablePairScore, order = option$order.Hprojection,
                                                 weight = outArgs$weight,
                                                 count.favorable = colSums(outPoint$count_favorable), count.unfavorable = colSums(outPoint$count_unfavorable),
                                                 n.pairs = outPoint$n_pairs, n.C = length(envirBT$outArgs$index.C), n.T = length(envirBT$outArgs$index.T),                                                                                   level.strata = outArgs$level.strata, n.strata = outArgs$n.strata, endpoint = outArgs$endpoint)

        outPoint$Mvar <- outCovariance$Sigma
        outPoint$iid_favorable <- NULL
        outPoint$iid_unfavorable <- NULL
        attr(outArgs$method.inference,"Hprojection") <- option$order.Hprojection
    }else if(grepl("bootstrap|permutation",outArgs$method.inference)){
        outResampling <- inferenceResampling(envirBT)
        if(outArgs$iid==FALSE){
            outPoint$Mvar <- matrix(nrow = 0, ncol = 0)
            outPoint$iid_favorable <- NULL
            outPoint$iid_unfavorable <- NULL
        }
    }
    if((outArgs$method.inference != "none") && (outArgs$trace > 1)){
        cat("\n")
    }

    ## ** Gather results into a BuyseRes object
    if(outArgs$trace > 1){
        cat("Gather the results in a BuyseRes object \n")
    }
    scoring.rule <- c("Gehan","Peron")[outArgs$scoring.rule+1]
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
        scoring.rule = scoring.rule,
        hierarchical = outArgs$hierarchical,
        neutral.as.uninf = outArgs$neutral.as.uninf,
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
        covarianceResampling = outResampling$covariance,
        covariance = outPoint$Mvar,
        weight = outArgs$weight,
        iid_favorable = outPoint$iid_favorable,
        iid_unfavorable = outPoint$iid_unfavorable,
        iidNuisance_favorable = outPoint$iidNuisance_favorable,
        iidNuisance_unfavorable = outPoint$iidNuisance_unfavorable,
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
                       iid,
                       method.inference,
                       pointEstimation){
   
    n.strata <- envir$outArgs$n.strata ## to simplify code
    treatment <- envir$outArgs$treatment ## to simplify code
    censoring <- envir$outArgs$censoring ## to simplify code
    endpoint <- envir$outArgs$endpoint ## to simplify code
    scoring.rule <- envir$outArgs$scoring.rule ## to simplify code
    type <- envir$outArgs$type ## to simplify code
    D.TTE <- envir$outArgs$D.TTE ## to simplify code
    D <- envir$outArgs$D ## to simplify code


    ## ** Resampling
    ls.indexC <- vector(mode = "list", length = n.strata)
    ls.indexT <- vector(mode = "list", length = n.strata)
   
    if(method.inference == "none"){

        ## find groups
        if(n.strata==1){        
            ls.indexC[[1]] <- envir$outArgs$index.C - 1
            ls.indexT[[1]] <- envir$outArgs$index.T - 1
        }else{        
            for(iStrata in 1:n.strata){ ## iStrata <- 1  
                ls.indexC[[iStrata]] <- intersect(envir$outArgs$index.C, envir$outArgs$index.strata[[iStrata]]) - 1
                ls.indexT[[iStrata]] <- intersect(envir$outArgs$index.T, envir$outArgs$index.strata[[iStrata]]) - 1
            }
        }

        ## rebuild dataset
        if(scoring.rule>0){
            data <- data.table::data.table(envir$outArgs$data,envir$outArgs$M.endpoint,envir$outArgs$M.censoring)
        }
        
    }else if(attr(method.inference, "permutation")){

        ## permute
        if(attr(method.inference, "stratified")){
            index.resampling <- NULL
            for(iStrata in 1:n.strata){ ## iStrata <- 1  
                index.resampling <- c(index.resampling,envir$outArgs$cumn.obsStrata[iStrata] + sample.int(envir$outArgs$n.obsStrata[iStrata], replace = FALSE))
            }
        }else{
            index.resampling <- sample.int(envir$outArgs$n.obs, replace = FALSE)

        }
        ## find groups
        if(n.strata==1){        
            ls.indexC[[1]] <- which(index.resampling %in% envir$outArgs$index.C) - 1
            ls.indexT[[1]] <- which(index.resampling %in% envir$outArgs$index.T) - 1            
        }else{
            index.C <- which(index.resampling %in% envir$outArgs$index.C)
            index.T <- which(index.resampling %in% envir$outArgs$index.T)
            for(iStrata in 1:n.strata){ ## iStrata <- 1  
                ls.indexC[[iStrata]] <- intersect(index.C, envir$outArgs$index.strata[[iStrata]]) - 1
                ls.indexT[[iStrata]] <- intersect(index.T, envir$outArgs$index.strata[[iStrata]]) - 1
            }
            ## ls.indexC[[1]]
        }

        ## rebuild dataset
        if(scoring.rule>0){
            data <- data.table::data.table(envir$outArgs$data[[treatment]][index.resampling],
                                           "..strata.." = envir$outArgs$data[["..strata.."]],
                                           envir$outArgs$M.endpoint,envir$outArgs$M.censoring)
            data.table::setnames(data, old = names(data)[1], new = treatment)
        }
        
    }else if(attr(method.inference, "bootstrap")){ 

        ## bootstrap
        if(attr(method.inference, "stratified")){
            index.resampling <- NULL
            for(iStrata in 1:n.strata){ ## iStrata <- 1  
                index.resampling <- c(index.resampling,envir$outArgs$cumn.obsStrata[iStrata] + sample.int(envir$outArgs$n.obsStrata[iStrata], replace = TRUE))
            }            
        }else{
            index.resampling <- sample.int(envir$outArgs$n.obs, replace = TRUE)
        }
            
        ## find groups
        if(n.strata==1){
            ls.indexC[[1]] <- index.resampling[index.resampling %in% envir$outArgs$index.C] - 1
            ls.indexT[[1]] <- index.resampling[index.resampling %in% envir$outArgs$index.T] - 1
        }else{
            index.C <- index.resampling[index.resampling %in% envir$outArgs$index.C]
            index.T <- index.resampling[index.resampling %in% envir$outArgs$index.T]
            ## update strata!!!!
            for(iStrata in 1:n.strata){ ## iStrata <- 1
                ## do not use intersect since it removes duplicated elements
                iIndex.strata <- index.resampling[index.resampling %in% envir$outArgs$index.strata[[iStrata]]]                
                ls.indexC[[iStrata]] <- index.C[index.C %in% iIndex.strata] - 1
                ls.indexT[[iStrata]] <- index.T[index.T %in% iIndex.strata] - 1
            }
        }

        ## lapply(ls.indexC,length)
        ## lapply(ls.indexT,length)
        ## table(data[["..strata.."]],data[[treatment]])
        
        ## rebuild dataset
        if(scoring.rule>0){
            data <- data.table::data.table(envir$outArgs$data,
                                           envir$outArgs$M.endpoint,envir$outArgs$M.censoring)[index.resampling]
        }
    }

    ## new ordering of the observations in the dataset
    if(method.inference == "none"){
        ls.posC <- ls.indexC
        ls.posT <- ls.indexT
    }else if(attr(method.inference,"bootstrap") || attr(method.inference,"permutation")){
        ## new position in the dataset
        new.cumn <- c(0,cumsum(sapply(ls.indexC,length))+cumsum(sapply(ls.indexT,length)))
        ls.posC <- vector(mode = "list", length = n.strata)
        ls.posT <- vector(mode = "list", length = n.strata)
        index.tempo <- 0
        for(iStrata in 1:n.strata){
            ls.posC[[iStrata]] <- new.cumn[iStrata] + 0:(length(ls.indexC[[iStrata]])-1)
            ls.posT[[iStrata]] <- new.cumn[iStrata] + length(ls.indexC[[iStrata]]) + 0:(length(ls.indexT[[iStrata]])-1)
        }
        ## Check valid resampling
        if (any(c(sapply(ls.indexC,length),sapply(ls.indexT,length))==0)) {
            return(NULL)
        }
    }

    ## *** Update survival
    if(scoring.rule == 0){ ## Gehan
        outSurv <- envir$outArgs$skeletonPeron
    }else{ ## Peron 
        outSurv <- initializePeron(data = data,
                                   model.tte = envir$outArgs$model.tte,
                                   method.score = envir$outArgs$method.score,
                                   treatment = treatment,
                                   level.treatment = envir$outArgs$level.treatment,
                                   endpoint = endpoint,
                                   endpoint.TTE = envir$outArgs$endpoint.TTE,
                                   endpoint.UTTE = envir$outArgs$endpoint.UTTE,
                                   censoring = censoring,
                                   D.TTE = D.TTE,
                                   D.UTTE = envir$outArgs$D.UTTE,
                                   type = type,
                                   threshold = envir$outArgs$threshold,
                                   n.strata = n.strata,
                                   strata = envir$outArgs$strata,
                                   iid = iid,
                                   out = envir$outArgs$outSurv)
    }
    
    ## ** Computation
    resBT <- GPC_cpp(endpoint = envir$outArgs$M.endpoint,
                     censoring = envir$outArgs$M.censoring,
                     indexC = ls.indexC,
                     posC = ls.posC,
                     indexT = ls.indexT,                     
                     posT = ls.posT,                     
                     threshold = envir$outArgs$threshold,
                     weight = envir$outArgs$weight,
                     method = envir$outArgs$method.score,
                     D = D,
                     D_UTTE = envir$outArgs$D.UTTE,
                     n_strata = n.strata,
                     nUTTE_analyzed_M1 = envir$outArgs$nUTTE.analyzedPeron_M1,
                     index_endpoint = envir$outArgs$index.endpoint,
                     index_censoring = envir$outArgs$index.censoring,
                     index_UTTE = envir$outArgs$index.UTTE,
                     list_survTimeC = outSurv$survTimeC,
                     list_survTimeT = outSurv$survTimeT,
                     list_survJumpC = outSurv$survJumpC,
                     list_survJumpT = outSurv$survJumpT,
                     list_lastSurv = outSurv$lastSurv,
                     p_C = outSurv$p.C,
                     p_T = outSurv$p.T,
                     iid_survJumpC = outSurv$iid$survJumpC,
                     iid_survJumpT = outSurv$iid$survJumpT,
                     zeroPlus = 1e-8,
                     correctionUninf = envir$outArgs$correction.uninf,
                     hierarchical = envir$outArgs$hierarchical,
                     hprojection = envir$outArgs$order.Hprojection,
                     neutralAsUninf = envir$outArgs$neutral.as.uninf,
                     keepScore = (pointEstimation && envir$outArgs$keep.pairScore),
                     reserve = 1e4,
                     returnIID = iid + iid*iidNuisance
                     )

    ## ** export
    if(pointEstimation){
        if(envir$outArgs$keep.survival){ ## useful to test initSurvival 
            resBT$tableSurvival <- outSurv
        }
        return(resBT)
    }else{
        return(list(delta_netBenefit = resBT$delta_netBenefit,
                    Delta_netBenefit = resBT$Delta_netBenefit,
                    delta_winRatio = resBT$delta_winRatio,
                    Delta_winRatio = resBT$Delta_winRatio,
                    Mvar = resBT$Mvar))
    }
}







