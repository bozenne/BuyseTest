## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Generalized Pairwise Comparisons (GPC)
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event endpoints.
#' 
#' @param formula [formula] a symbolic description of the GPC model,
#' typically \code{treatment ~ type1(endpoint1) + type2(endpoint2, threshold2) + strata}.
#' See Details, section "Specification of the GPC model".
#' @param treatment,endpoint,type,threshold,status,operator,censoring,strata Alternative to \code{formula} for describing the GPC model.
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
#' Can be \code{"none"}, \code{"u-statistic"}, \code{"permutation"}, \code{"studentized permutation"}, \code{"bootstrap"}, \code{"studentized bootstrap"}.
#' See Details, section "Statistical inference".
#' @param n.resampling [integer] the number of permutations/samples used for computing the confidence intervals and the p.values. 
#' See Details, section "Statistical inference".
#' @param strata.resampling [character] the variable on which the permutation/sampling should be stratified. 
#' See Details, section "Statistical inference".
#' @param hierarchical [logical] should only the uninformative pairs be analyzed at the lower priority endpoints (hierarchical GPC)?
#' Otherwise all pairs will be compaired for all endpoint (full GPC).
#' @param weight [numeric vector] weights used to cumulating the pairwise scores over the endpoints.
#' Only used when \code{hierarchical=FALSE}. Disregarded if the argument \code{formula} is defined.
#' @param neutral.as.uninf [logical vector] should paired classified as neutral be re-analyzed using endpoints of lower priority (as it is done for uninformative pairs).
#' See Details, section "Handling missing values".
#' @param keep.pairScore [logical] should the result of each pairwise comparison be kept?
#' @param seed [integer, >0] the seed to consider when performing resampling.
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
#' A pair will be classified as neutral if the difference in endpoint is strictly below this threshold.
#' There must be one threshold for each endpoint variable; it must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. 
#'   \item \code{status}: [character vector] the name of the binary variable(s) indicating whether the endpoint was observed or censored.
#' Must value \code{NA} when the endpoint is not a time to event.
#'   \item \code{operator}: [character vector] the sign defining a favorable endpoint.
#' \code{">0"} indicates that higher values are favorable while "<0" indicates the opposite.
#'   \item \code{type}: [character vector] indicates whether it is
#' a binary outcome  (\code{"b"}, \code{"bin"}, or \code{"binary"}),
#' a continuous outcome  (\code{"c"}, \code{"cont"}, or \code{"continuous"}),
#' or a time to event outcome  (\code{"t"}, \code{"tte"}, \code{"time"}, or \code{"timetoevent"})
#'   \item \code{censoring}: [character vector] is the endpoint subject to right or left censoring (\code{"left"} or \code{"right"}). The default is right-censoring.
#'   \item \code{strata}: [character vector] if not \code{NULL}, the GPC will be applied within each group of patient defined by the strata variable(s).
#' }
#' The formula interface can be more concise, especially when considering few outcomes, but may be more difficult to apprehend for new users.
#' Note that arguments \code{endpoint}, \code{threshold}, \code{status}, \code{operator},  \code{type}, and \code{censoring} must have the same length. \cr \cr \cr
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
#' The Peron's scoring rule is the recommanded scoring rule but only handles right-censoring.
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
#' \bold{WARNING}: the current implementation of the H-projection is not valid when using corrections for uninformative pairs (\code{correction.uninf=1}, or \code{correction.uninf=2}).
#'   \item argument \code{method.inference="permutation"}: perform a permutation test, estimating in each sample the summary statistics (net benefit, win ratio).
#'   \item argument \code{method.inference="studentized permutation"}: perform a permutation test, estimating in each sample the summary statistics (net benefit, win ratio) and the variance-covariance matrix of the estimate.
#'   \item argument \code{method.inference="bootstrap"}: perform a non-parametric boostrap, estimating in each sample the summary statistics (net benefit, win ratio).
#'   \item argument \code{method.inference=" studentized bootstrap"}: perform a non-parametric boostrap, estimating in each sample the summary statistics (net benefit, win ratio) and the variance-covariance matrix of the estimator.
#' }
#' Additional arguments for permutation and bootstrap resampling:
#' \itemize{
#'    \item \code{strata.resampling} If \code{NA} or of length 0, the permutation/non-parametric boostrap will be performed by resampling in the whole sample.
#' Otherwise, the permutation/non-parametric boostrap will be performed separately for each level that the variable defined in \code{strata.resampling} take.
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
#' \code{hierarchical}, \code{neutral.as.uninf}, \code{keep.pairScore}, \code{strata.resampling},
#' \code{cpus}, \code{trace} is read from \code{BuyseTest.options()}. \cr
#' Additional (hidden) arguments are \itemize{
#'  \item \code{alternative} [character] the alternative hypothesis. Must be one of "two.sided", "greater" or "less" (used by \code{confint}).
#'  \item \code{conf.level} [numeric] level for the confidence intervals (used by \code{confint}).
#'  \item \code{keep.survival} [logical] export the survival values used by the Peron's scoring rule.
#'  \item \code{order.Hprojection} [1 or 2] the order of the H-projection used to compute the variance when \code{method.inference="u-statistic"}. 
#' }
#' 
#' @return An \R object of class \code{\linkS4class{S4BuyseTest}}.
#' 
#' @references 
#' On the GPC procedure: Marc Buyse (2010). \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257 \cr
#' On the win ratio: D. Wang, S. Pocock (2016). \bold{A win ratio approach to comparing continuous non-normal outcomes in clinical trials}. \emph{Pharmaceutical Statistics} 15:238-245 \cr
#' On the Peron's scoring rule: J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2018). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. \emph{Statistical Methods in Medical Research} 27: 1230-1239  \cr 
#' On the Gehan's scoring rule: Gehan EA (1965). \bold{A generalized two-sample Wilcoxon test for doubly censored data}. \emph{Biometrika}  52(3):650-653 \cr
#' On inference in GPC using the U-statistic theory: I. Bebu, J. M. Lachin (2015). \bold{Large sample inference for a win ratio analysis of a composite outcome based on prioritized components}. \emph{Biostatistics} 17(1):178-187 \cr
#'
#' @seealso 
#' \code{\link{S4BuyseTest-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-class}} for a presentation of the \code{S4BuyseTest} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' @keywords function BuyseTest
#' @author Brice Ozenne

## * BuyseTest (example)
#' @rdname BuyseTest
#' @examples
#' library(data.table)
#' 
#' #### simulate some data ####
#' set.seed(10)
#' df.data <- simBuyseTest(1e2, n.strata = 2)
#'
#' ## display 
#' if(require(prodlim)){
#'    resKM_tempo <- prodlim(Hist(eventtime,status)~treatment, data = df.data)
#'    plot(resKM_tempo)
#' }
#'
#' #### one time to event endpoint ####
#' BT <- BuyseTest(treatment ~ TTE(eventtime, status = status), data= df.data)
#'
#' summary(BT) # net benefit
#' summary(BT, percentage = FALSE)  
#' summary(BT, statistic = "winRatio") # win Ratio
#' 
#' ## permutation instead of asymptotics to compute the p-value
#' \dontrun{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3)
#' }
#' \dontshow{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e1, trace = 0)
#' }
#' summary(BT, statistic = "netBenefit") ## default
#' summary(BT, statistic = "winRatio") 
#' 
#' ## parallel permutation
#' \dontrun{
#'     BT <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3, cpus = 2)
#'     summary(BT)
#' }
#' 
#' ## method Gehan is much faster but does not optimally handle censored observations
#' BT <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                 scoring.rule = "Gehan", trace = 0)
#' summary(BT)
#' 
#' #### one time to event endpoint: only differences in survival over 1 unit ####
#' BT <- BuyseTest(treatment ~ TTE(eventtime, threshold = 1, status = status), data=df.data)
#' summary(BT)
#' 
#' #### one time to event endpoint with a strata variable
#' BT <- BuyseTest(treatment ~ strata + TTE(eventtime, status = status), data=df.data)
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
#' #### real example : veteran dataset of the survival package ####
#' ## Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#' ## scoring.rule = "Gehan"
#' 
#' if(require(survival)){
#' \dontrun{
#'   library(survival) ## import veteran
#'  
#'   ## scoring.rule = "Gehan"
#'   BT_Gehan <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,status=status), 
#'                         data=veteran, scoring.rule="Gehan")
#'   
#'   summary_Gehan <- summary(BT_Gehan)
#'   summary_Gehan <- summary(BT_Gehan, statistic = "winRatio")
#'   
#'   ## scoring.rule = "Peron"
#'   BT_Peron <- BuyseTest(trt ~ celltype + TTE(time,threshold=0,status=status), 
#'                         data=veteran, scoring.rule="Peron")
#' 
#'   class(BT_Peron)
#'   summary(BT_Peron)
#' }
#' }

## * BuyseTest (code)
##' @export
BuyseTest <- function(formula,
                      data,
                      scoring.rule = NULL,
                      correction.uninf = NULL,
                      model.tte = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      strata.resampling = NULL,
                      hierarchical = NULL,
                      weight = NULL,
                      neutral.as.uninf = NULL,
                      keep.pairScore = NULL,
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL,
                      treatment = NULL,
                      endpoint = NULL,
                      type = NULL,
                      threshold = NULL,                      
                      status = NULL,
                      operator = NULL,
                      censoring = NULL,
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
    outArgs <- initializeArgs(status = status,
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
                              strata.resampling = strata.resampling,
                              name.call = name.call,
                              neutral.as.uninf = neutral.as.uninf,
                              operator = operator,
                              censoring = censoring,
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
    }

    ## ** initialization data
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    out.name <- c("data","M.endpoint","M.status",
                  "index.C","index.T","index.strata",
                  "level.treatment","level.strata", "method.score",
                  "n.strata","n.obs","n.obsStrata","n.obsStrataResampling","cumn.obsStrataResampling","skeletonPeron",
                  "scoring.rule", "iidNuisance", "nUTTE.analyzedPeron_M1", "endpoint.UTTE", "status.UTTE", "D.UTTE","index.UTTE","keep.pairScore")
    outArgs[out.name] <- initializeData(data = outArgs$data,
                                        type = outArgs$type,
                                        endpoint = outArgs$endpoint,
                                        Uendpoint = outArgs$Uendpoint,
                                        D = outArgs$D,
                                        scoring.rule = outArgs$scoring.rule,
                                        status = outArgs$status,
                                        Ustatus = outArgs$Ustatus,
                                        method.inference = outArgs$method.inference,
                                        censoring = outArgs$censoring,
                                        strata = outArgs$strata,
                                        treatment = outArgs$treatment,
                                        hierarchical = outArgs$hierarchical,
                                        copy = TRUE,
                                        keep.pairScore = outArgs$keep.pairScore,
                                        endpoint.TTE = outArgs$endpoint.TTE,
                                        status.TTE = outArgs$status.TTE,
                                        iidNuisance = outArgs$iidNuisance)
    
    if(option$check){
        if(outArgs$iidNuisance && any(outArgs$method.score == 5)){
            warning("Inference via the asymptotic theory  for competing risks when using the Peron's scoring rule has not been validating \n",
                    "Consider setting \'method.inference\' to \"none\", \"bootstrap\", or \"permutation\" \n")
        }
        ## if(outArgs$precompute && any(outArgs$method.score == 5)){
        ##     stop("Option \'precompute\' is not available for the Peron scoring rule in the competing risk case \n")
        ## }
        if(outArgs$precompute && any(outArgs$method.score == 5)){
            outArgs$precompute <- FALSE
        }
    }

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
    envirBT$calcPeron <- calcPeron
    envirBT$outArgs$args.model.tte <- option$args.model.tte

    ## ** Point estimation
    if (outArgs$trace > 1) {
        if(outArgs$iid){
            cat("Point estimation and calculation of the iid decomposition")
        }else{
            cat("Point estimation")
        }
        
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
        outPoint$tableScore <- pairScore2dt(outPoint$tableScore,
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

    outResampling <- NULL
    if(outArgs$method.inference == "u-statistic"){
        ## done in the C++ code
        ## outCovariance <- inferenceUstatistic(tablePairScore = outPoint$tableScore, order = option$order.Hprojection,
        ##                                      count.favorable = colSums(outPoint$count_favorable), count.unfavorable = colSums(outPoint$count_unfavorable),
        ##                                      n.pairs = sum(outPoint$n_pairs), n.C = length(envirBT$outArgs$index.C), n.T = length(envirBT$outArgs$index.T),
        ##                                      level.strata = outArgs$level.strata, n.strata = outArgs$n.strata, endpoint = outArgs$endpoint)
    }else if(outArgs$method.inference == "u-statistic-bebu"){
        if(outArgs$keep.pairScore == FALSE){
            stop("Argument \'keep.pairScore\' needs to be TRUE when argument \'method.inference\' is \"u-statistic-bebu\" \n")
        }

        ## direct computation of the variance
        outCovariance <- inferenceUstatisticBebu(tablePairScore = outPoint$tableScore,
                                                 order = option$order.Hprojection,
                                                 weight = outArgs$weight,
                                                 count.favorable = colSums(outPoint$count_favorable),
                                                 count.unfavorable = colSums(outPoint$count_unfavorable),
                                                 n.pairs = outPoint$n_pairs,
                                                 n.C = length(envirBT$outArgs$index.C),
                                                 n.T = length(envirBT$outArgs$index.T),
                                                 level.strata = outArgs$level.strata,
                                                 n.strata = outArgs$n.strata,
                                                 endpoint = outArgs$endpoint)

        outPoint$covariance <- outCovariance$Sigma
        attr(outArgs$method.inference,"Hprojection") <- option$order.Hprojection
    }else if(grepl("bootstrap|permutation",outArgs$method.inference)){
        outResampling <- inferenceResampling(envirBT)
    }
    if((outArgs$method.inference != "none") && (outArgs$trace > 1)){
        cat("\n")
    }

    ## ** Gather results into a S4BuyseTest object
    if(outArgs$trace > 1){
        cat("Gather the results in a S4BuyseTest object \n")
    }
    keep.args <- c("index.T", "index.C", "type","endpoint","level.strata","level.treatment","scoring.rule","hierarchical","neutral.as.uninf",
                   "correction.uninf","method.inference","method.score","strata","threshold","weight","n.resampling")
    BuyseTest.object <- do.call("S4BuyseTest", args = c(outPoint, outArgs[keep.args], outResampling))
    
    if(outArgs$trace > 1){
        cat("\n")
    }
    
    ## ** export
    return(BuyseTest.object)
}

## * .BuyseTest (code)
.BuyseTest <- function(envir,
                       iid,
                       method.inference,
                       pointEstimation){

    ## ** Resampling
    outSample <- calcSample(envir = envir, method.inference = method.inference)
    if(is.null(outSample)){return(NULL)}

    ## ** Estimate survival curves with its iid
    if(envir$outArgs$scoring.rule == 0){ ## Gehan
        outSurv <- envir$outArgs$skeletonPeron
    }else{ ## Peron
        outSurv <- calcPeron(data = outSample$data,
                             model.tte = envir$outArgs$model.tte,                             
                             method.score = envir$outArgs$method.score,
                             treatment = envir$outArgs$treatment,
                             level.treatment = envir$outArgs$level.treatment,
                             endpoint = envir$outArgs$endpoint,
                             endpoint.TTE = envir$outArgs$endpoint.TTE,
                             endpoint.UTTE = envir$outArgs$endpoint.UTTE,
                             status = envir$outArgs$status,
                             status.TTE = envir$outArgs$status.TTE,
                             status.UTTE = envir$outArgs$status.UTTE,
                             D.TTE = envir$outArgs$D.TTE,
                             D.UTTE = envir$outArgs$D.UTTE,
                             type = envir$outArgs$type,
                             threshold = envir$outArgs$threshold,
                             level.strata = envir$outArgs$level.strata,
                             n.strata = envir$outArgs$n.strata,
                             strata = envir$outArgs$strata,
                             precompute = envir$outArgs$precompute,
                             iidNuisance = envir$outArgs$iidNuisance * iid,
                             out = envir$outArgs$skeletonPeron,
                             fitter = envir$outArgs$fitter.model.tte,
                             args = envir$outArgs$args.model.tte)
    }

    ## ** Perform GPC
    resBT <- do.call(envir$outArgs$engine,
                     args = list(endpoint = envir$outArgs$M.endpoint,
                                 status = envir$outArgs$M.status,
                                 indexC = outSample$ls.indexC,
                                 posC = outSample$ls.posC,
                                 indexT = outSample$ls.indexT,                     
                                 posT = outSample$ls.posT,                     
                                 threshold = envir$outArgs$threshold,
                                 weight = envir$outArgs$weight,
                                 method = envir$outArgs$method.score,
                                 op = envir$outArgs$operator,
                                 D = envir$outArgs$D,
                                 D_UTTE = envir$outArgs$D.UTTE,
                                 n_strata = envir$outArgs$n.strata,
                                 nUTTE_analyzedPeron_M1 = envir$outArgs$nUTTE.analyzedPeron_M1,
                                 index_endpoint = envir$outArgs$index.endpoint,
                                 index_status = envir$outArgs$index.status,
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
                                 precompute = envir$outArgs$precompute,
                                 returnIID = iid + iid*envir$outArgs$iidNuisance,
                                 debug = envir$outArgs$debug
                                 ))
    
    ## print(range(resBT$iidNuisance_favorable))
    ## print(range(resBT$iidNuisance_unfavorable))
    
    ## ** export
    if(pointEstimation){
        if(envir$outArgs$keep.survival){ ## useful to test initSurvival 
            resBT$tableSurvival <- outSurv
        }
        return(resBT)
    }else{
        return(list(delta = resBT$delta,
                    Delta = resBT$Delta,
                    covariance = resBT$covariance))
    }
}

## * calcSample
#' @rdname internal-initialization
calcSample <- function(envir, method.inference){

    ## ** initialization
    out <- list(## rows in M.endpoint/M.status corresponding to observations from the control/treatment group (not unique when boostraping)
        ls.indexC = vector(mode = "list", length = envir$outArgs$n.strata), 
        ls.indexT = vector(mode = "list", length = envir$outArgs$n.strata),
        ## identifier for each observation from the control/treatment group (unique even when boostrap)
        ls.posC = vector(mode = "list", length = envir$outArgs$n.strata),
        ls.posT = vector(mode = "list", length = envir$outArgs$n.strata),
        data = data.table::data.table()
    )

    if(method.inference == "none"){

        ## ** no resampling
        if(envir$outArgs$n.strata==1){        
            out$ls.indexC[[1]] <- envir$outArgs$index.C - 1
            out$ls.indexT[[1]] <- envir$outArgs$index.T - 1
        }else{        
            for(iStrata in 1:envir$outArgs$n.strata){ ## iStrata <- 1  
                out$ls.indexC[[iStrata]] <- intersect(envir$outArgs$index.C, envir$outArgs$index.strata[[iStrata]]) - 1
                out$ls.indexT[[iStrata]] <- intersect(envir$outArgs$index.T, envir$outArgs$index.strata[[iStrata]]) - 1
            }
        }
        out$ls.posC <- out$ls.indexC
        out$ls.posT <- out$ls.indexT

        if(envir$outArgs$scoring.rule>0){
            out$data <- data.table::data.table(envir$outArgs$data,envir$outArgs$M.endpoint,envir$outArgs$M.status)
        }
    }else{

        ## ** stratified resampling
        n.strataResampling <- length(envir$outArgs$n.obsStrataResampling)
        index.resampling <- NULL
        for (iSR in 1:n.strataResampling) {
            index.resampling <- c(index.resampling,
                                  envir$outArgs$cumn.obsStrataResampling[iSR] + sample.int(envir$outArgs$n.obsStrataResampling[iSR], replace = grepl("bootstrap",method.inference)))
        }

        ## ** reconstruct groups
        ## index: index of the new observations in the old dataset by treatment group
        ## pos: unique identifier for each observation
        if(envir$outArgs$n.strata==1){ ## no strata
            
            if(grepl("permutation",method.inference)){
                out$ls.indexC[[1]] <- which(index.resampling %in% envir$outArgs$index.C) - 1
                out$ls.indexT[[1]] <- which(index.resampling %in% envir$outArgs$index.T) - 1
                out$ls.posC[[1]] <- out$ls.indexC[[1]]
                out$ls.posT[[1]] <- out$ls.indexT[[1]]
            }else if(grepl("bootstrap",method.inference)){
                out$ls.posC[[1]] <- which(index.resampling %in% envir$outArgs$index.C) - 1
                out$ls.posT[[1]] <- which(index.resampling %in% envir$outArgs$index.T) - 1
                out$ls.indexC[[1]] <- index.resampling[out$ls.posC[[1]] + 1] - 1
                out$ls.indexT[[1]] <- index.resampling[out$ls.posT[[1]] + 1] - 1
            }
            ## check that each group has at least one observation
            if(length(out$ls.indexC[[1]])==0 || length(out$ls.indexT[[1]])==0){return(NULL)}
            ## out$data[treatment == 0,eventtime1] - envir$outArgs$M.endpoint[out$ls.indexC[[1]]+1,1]
            ## out$data[treatment == 1,eventtime1] - envir$outArgs$M.endpoint[out$ls.indexT[[1]]+1,1]
            
        }else{ ## strata

            if (grepl("permutation",method.inference)) {
                index.C <- which(index.resampling %in% envir$outArgs$index.C)
                index.T <- which(index.resampling %in% envir$outArgs$index.T)
            }
            
            for(iStrata in 1:envir$outArgs$n.strata){ ## iStrata <- 1  
                ## index of the new observation in the old dataset by treatment group
                if(grepl("permutation",method.inference)){
                    out$ls.indexC[[iStrata]] <- intersect(index.C, envir$outArgs$index.strata[[iStrata]]) - 1
                    out$ls.indexT[[iStrata]] <- intersect(index.T, envir$outArgs$index.strata[[iStrata]]) - 1
                    out$ls.posC[[iStrata]] <- out$ls.indexC[[iStrata]]
                    out$ls.posT[[iStrata]] <- out$ls.indexT[[iStrata]]
                }else if(grepl("bootstrap",method.inference)){
                    out$ls.posC[[iStrata]] <- which(index.resampling %in% intersect(envir$outArgs$index.C, envir$outArgs$index.strata[[iStrata]])) - 1
                    out$ls.posT[[iStrata]] <- which(index.resampling %in% intersect(envir$outArgs$index.T, envir$outArgs$index.strata[[iStrata]])) - 1
                    out$ls.indexC[[iStrata]] <- index.resampling[out$ls.posC[[iStrata]] + 1] - 1
                    out$ls.indexT[[iStrata]] <- index.resampling[out$ls.posT[[iStrata]] + 1] - 1
                }
                ## check that each group has at least one observation
                if(length(out$ls.indexC[[iStrata]])==0 || length(out$ls.indexT[[iStrata]])==0){return(NULL)} 
            }
            
        }

        ## ** rebuild dataset
        if(envir$outArgs$scoring.rule>0){
            if(grepl("permutation",method.inference)){
                out$data <- data.table::data.table(envir$outArgs$data[[envir$outArgs$treatment]][index.resampling],
                                                   "..strata.." = envir$outArgs$data[["..strata.."]],
                                                   envir$outArgs$M.endpoint,envir$outArgs$M.status)
                data.table::setnames(out$data, old = names(out$data)[1], new = envir$outArgs$treatment)
            }else{
                out$data <- data.table::data.table(envir$outArgs$data[,.SD,.SDcols = c(envir$outArgs$treatment,"..strata..")],
                                                   envir$outArgs$M.endpoint,
                                                   envir$outArgs$M.status)[index.resampling]
            }
        }

    }
    return(out)
}









