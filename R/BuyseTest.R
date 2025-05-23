## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Two-group GPC
#' 
#' @description Performs Generalized Pairwise Comparisons (GPC) between two groups.
#' Can handle one or several binary, continuous and time-to-event endpoints.
#' 
#' @param formula [formula] a symbolic description of the GPC model,
#' typically \code{type1(endpoint1) + type2(endpoint2, threshold2) ~ treatment + strata(gender)}.
#' See Details, section "Specification of the GPC model".
#' @param treatment,endpoint,type,threshold,status,operator,censoring,restriction,strata Alternative to \code{formula} for describing the GPC model.
#' See Details, section "Specification of the GPC model".
#' @param data [data.frame] dataset.
#' @param scoring.rule [character] method used to compare the observations of a pair in presence of right censoring (i.e. \code{"timeToEvent"} endpoints).
#' Can be \code{"Gehan"}, \code{"Peron"}, or \code{"Efron"}.
#' See Details, section "Handling missing values".
#' @param pool.strata [character] weights used to combine estimates across strata. Can be
#' \code{"Buyse"} to weight proportionally to the number of pairs in the strata,
#' \code{"CMH"} to weight proportionally to the ratio between the number of pairs in the strata and the number of observations in the strata.
#' \code{"equal"} to weight equally each strata,
#' \code{"standardisation"} to recover a marginal comparison,
#' or \code{"var-netBenefit"} to weight each strata proportionally to the precision of its estimated net benefit (similar syntax for the win ratio: \code{"var-winRatio"})
#' @param correction.uninf [integer] should a correction be applied to remove the bias due to the presence of uninformative pairs?
#' 0 indicates no correction, 1 impute the average score of the informative pairs, and 2 performs IPCW.
#' See Details, section "Handling missing values".
#' @param model.tte [list] optional survival models relative to each time to each time to event endpoint.
#' Models must \code{prodlim} objects and stratified on the treatment and strata variable. When used, the uncertainty from the estimates of these survival models is ignored.
#' @param method.inference [character] method used to compute confidence intervals and p-values.
#' Can be \code{"none"}, \code{"u-statistic"}, \code{"permutation"}, \code{"studentized permutation"}, \code{"bootstrap"}, \code{"studentized bootstrap"}, \code{"varExact permutation"}.
#' See Details, section "Statistical inference".
#' @param n.resampling [integer] the number of permutations/samples used for computing the confidence intervals and the p.values. 
#' See Details, section "Statistical inference".
#' @param strata.resampling [character] the variable on which the permutation/sampling should be stratified. 
#' See Details, section "Statistical inference".
#' @param hierarchical [logical] should only the uninformative pairs be analyzed at the lower priority endpoints (hierarchical GPC)?
#' Otherwise all pairs will be compaired for all endpoint (full GPC).
#' @param weightEndpoint [numeric vector] weights used to cumulating the pairwise scores over the endpoints.
#' Only used when \code{hierarchical=FALSE}. Disregarded if the argument \code{formula} is defined.
#' @param weightObs [character or numeric vector] weights or variable in the dataset containing the weight associated to each observation.
#' These weights are only considered when performing GPC (but not when fitting surival models).
#' @param neutral.as.uninf [logical vector] should pairs classified as neutral be re-analyzed using endpoints of lower priority (as it is done for uninformative pairs).
#' See Details, section "Handling missing values".
#' @param add.halfNeutral [logical] should half of the neutral score be added to the favorable and unfavorable scores?
#' @param keep.pairScore [logical] should the result of each pairwise comparison be kept?
#' @param seed [integer, >0] Random number generator (RNG) state used when starting resampling.
#' If \code{NULL} no state is set.
#' @param cpus [integer, >0] the number of CPU to use.
#' Only the permutation test can use parallel computation.
#' See Details, section "Statistical inference".
#' @param trace [integer] should the execution of the function be traced ? \code{0} remains silent
#' and \code{1}-\code{3} correspond to a more and more verbose output in the console.
#' 
#' @details
#'
#' \bold{Specification of the GPC model} \cr
#' There are two way to specify the GPC model in \code{BuyseTest}. 
#' A \emph{Formula interface} via the argument \code{formula} with:
#' \itemize{
#' \item on one side of the \code{~} symbol the endpoints by order of priority, surrounded by parentheses and a character string indicating the type of variable (binary,continuous,timetoevent). Additional argument may be included in the parenthesis: threshold of clinical relevance (\code{threshold}), restriction time for time to event endpoints (\code{restriction}), ...
#' \item in the other side a binary variable defining the two treatment arms.
#' \item the right-hand side of the formula may also contain a strata variable(s) surrounded by parenthese and the character strata. The argument \code{match} should be set to \code{TRUE} when strata are independent but not observations (e.g. cross-over trial with strata as patients).
#' }
#'
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
#'   \item \code{censoring}: [character vector] is the endpoint subject to right or left censoring (\code{"left"} or \code{"right"}). The default is right-censoring and left-censoring is only implemented with the Gehan's scoring rule.
#'   \item \code{restriction}: [numeric vector] value above which any difference is classified as neutral.
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
#' The classification may not be 0,1, e.g. the probability that the pair is favorable/unfavorable/neutral with the Peron's/Efron's scoring rule.
#' To export the classification of each pair set the argument \code{keep.pairScore} to \code{TRUE} and call the function \code{getPairScore} on the result of the \code{BuyseTest} function. \cr \cr \cr
#' 
#' 
#' \bold{Handling missing values}: the recommended default approach is to use the Peron's scoring rule with a restriction time if a non-neglectable part of the survival is unknown and otherwise analyse uniformative pairs using the following endpoint(s) if any. 
#' \itemize{
#'   \item \code{scoring.rule}: indicates how to handle right-censoring in time to event endpoints using information from the survival curves. \describe{
#'    \item{\code{scoring.rule="Gehan"}}{When no observations is censored or only the pair with the largest timepoint is censored, the pair is decidedly classified as favorable, unfavorable, or neutral. Otherwise pairs are classified as uninformative.}
#'    \item{\code{scoring.rule="Peron"}}{Score pairs involving censored observations using the (group-specific) survival curves. It may still lead to uninformative pairs when the survival curve is only partially known.}
#'    \item{\code{scoring.rule="Efron"}}{Same as the Peron's scoring rule except that the survival curve is extrapolated to 0 when its tail is unknown. Only relevant when using a (stratified) Kaplan-Meier estimator and no competing risks.}
#' }
#' 
#'   \item \code{correction.uninf}: indicates how to handle pairs that were classified as uninformative by the scoring rule.  \describe{
#'    \item{\code{correction.uninf=0}}{ treat them as uninformative: this is an equivalent to complete case analysis when \code{neutral.as.uninf=FALSE}, while when \code{neutral.as.uninf=TRUE}, uninformative pairs are treated as neutral, i.e., analyzed at the following endpoint (if any). This approach will (generally) lead to biased estimates for the proportion of favorable, unfavorable, or neutral pairs.}
#'    \item{\code{correction.uninf=1}}{ imputes to the uninformative pairs the average score of the informative pairs, i.e. assumes that uninformative pairs would on average behave like informative pairs. This is therefore the recommanded approach when this assumption is resonnable, typically when the the tail of the survival function estimated by the Kaplan–Meier method is close to 0.}
#'    \item{\code{correction.uninf=2}}{ up-weight informative pairs to represent uninformative pairs. It also assumes that uninformative pairs would on average behave like informative pairs and is only recommanded when the analysis is stopped after the first endpoint with uninformative pairs.}
#' }
#' }
#' 
#'
#' \bold{Statistical inference} \cr
#' The argument \code{method.inference} defines how to approximate the distribution of the GPC estimators and so how standard errors, confidence intervals, and p-values are computed.
#' Available methods are:
#' \itemize{
#'   \item argument \code{method.inference="none"}: only the point estimate is computed which makes the execution of the \code{BuyseTest} faster than with the other methods.
#'   \item argument \code{method.inference="u-statistic"}: compute the variance of the estimate using a H-projection of order 1 (default option) or 2 (see \code{BuyseTest.options}). The first order is downward biased but consistent. When considering the Gehan scoring rule, no transformation nor correction, the second order is unbiased and equivalent to the variance of the (stratified) bootstrap distribution. P-values and confidence intervals are then evaluated assuming that the estimates follow a Gaussian distribution. In presence of strata, the weights used to combine estimates across strata are assumed known. This is the case the type of weights used is pre-defined and the covariate distribution is fixed by design as in a blocked randomized trial.
#' \bold{WARNING}: the current implementation of the H-projection is not valid when using corrections for uninformative pairs (\code{correction.uninf=1}, or \code{correction.uninf=2}).
#'   \item argument \code{method.inference="permutation"}: perform a permutation test, estimating in each sample the summary statistics (net benefit, win ratio).
#'   \item argument \code{method.inference="studentized permutation"}: perform a permutation test, estimating in each sample the summary statistics (net benefit, win ratio) and the variance-covariance matrix of the estimate.
#'   \item argument \code{method.inference="varExact permutation"}: compute the variance of the permutation distribution using a closed-form formula (Anderson and Verbeeck 2023). P-values and confidence intervals are then evaluated assuming that the estimates follow a Gaussian distribution.
#' \bold{WARNING}: the current implementation of the variance estimator for the permutation distribution is not valid when using the Peron scoring rule or corrections for uninformative pairs.
#'   \item argument \code{method.inference="bootstrap"}: perform a non-parametric boostrap, estimating in each sample the summary statistics (net benefit, win ratio).
#'   \item argument \code{method.inference="studentized bootstrap"}: perform a non-parametric boostrap, estimating in each sample the summary statistics (net benefit, win ratio) and the variance-covariance matrix of the estimator.
#' }
#' Additional arguments for permutation and bootstrap resampling:
#' \itemize{
#'    \item \code{strata.resampling} If \code{NA} or of length 0, the permutation/non-parametric boostrap will be performed by resampling in the whole sample.
#' Otherwise, the permutation/non-parametric boostrap will be performed separately for each level that the variable defined in \code{strata.resampling} take.
#'    \item \code{n.resampling} set the number of permutations/samples used.
#' A large number of permutations (e.g. \code{n.resampling=10000}) are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details.
#'    \item \code{seed}: the seed is used to generate one seed per sample. These seeds are the same whether one or several CPUs are used.
#'    \item \code{cpus} indicates whether the resampling procedure can be splitted on several cpus to save time. Can be set to \code{"all"} to use all available cpus.
#' The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package. \cr \cr
#' }
#'
#' \bold{Pooling results across strata} \cr Consider \eqn{K} strata and denote by \eqn{m_k} and \eqn{n_k} the sample size in the control and active arm (respectively) for strata \eqn{k}. Let \eqn{\sigma_k} be the standard error of the strata-specific summary statistic (e.g. net benefit). The strata specific weights, \eqn{w_k}, are given by:
#' \itemize{
#' \item \code{"CMH"}: \eqn{w_k=\frac{\frac{m_k \times n_k}{m_k + n_k}}{\sum_{l=1}^K \frac{m_l \times n_l}{m_l + n_l}}}. Optimal if the if the odds ratios are constant across strata.
#' \item \code{"equal"}:  \eqn{w_k=\frac{1}{K}}
#' \item \code{"Buyse"}:  \eqn{w_k=\frac{m_k \times n_k}{\sum_{l=1}^K m_l \times n_l}}. Optimal if the risk difference is constant across strata
#' \item \code{"var-*"} (e.g. \code{"var-netBenefit"}): . \eqn{w_k=\frac{1/\sigma^2_k}{\sum_{l=1}^K 1/\sigma^2_k}}
#' }
#' Only when using \code{"var-winRatio"}, the pooled Win Ratio is computed by pooling the strata-specific win-ratios. Otherwise the pooled Win Ratio is obtained by dividing the pooled number of favorable pairs divided by the pooled number of unfavorable pairs, possibly adding half the pooled neutral pairs, according to formula (1) in Dong et al. (2018). \cr \cr
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
#' On the stratified win ratio: G. Dong et al. (2018). \bold{The stratified win ratio}. \emph{Journal of biopharmaceutical statistics}. 28(4):778-796 \cr
#' On the Peron's scoring rule: J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2018). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. \emph{Statistical Methods in Medical Research} 27: 1230-1239. \cr
#' On the Gehan's scoring rule: Gehan EA (1965). \bold{A generalized two-sample Wilcoxon test for doubly censored data}. \emph{Biometrika}  52(3):650-653 \cr
#' On inference in GPC using the U-statistic theory: Ozenne B, Budtz-Jorgensen E, Peron J (2021). \bold{The asymptotic distribution of the Net Benefit estimator in presence of right-censoring}. \emph{Statistical Methods in Medical Research} 2021 doi:10.1177/09622802211037067 \cr
#' On how using a restriction time: Piffoux M, Ozenne B, De Backer M, Buyse M, Chiem JC, Péron J (2024). \bold{Restricted Net Treatment Benefit in oncology}. \emph{Journal of Clinical Epidemiology}. Jun;170:111340. \cr
#' On the argument \code{correction.uninf}: J. Peron, M. Idlhaj, D. Maucort-Boulch, et al. (2021) \bold{Correcting the bias of the net benefit estimator due to right-censored observations}. \emph{Biometrical Journal} 63: 893–906. \cr
#' On closed-form formula for permutation variance:  W.N. Anderson and J. Verbeeck (2023). \bold{Exact Permutation and Bootstrap Distribution of Generalized Pairwise Comparisons Statistics}. \emph{Mathematics} , 11, 1502. doi:10.3390/math11061502.
#'
#' @seealso 
#' \code{\link{S4BuyseTest-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{S4BuyseTest-confint}} for exporting estimates with confidence intervals and p-values. \cr
#' \code{\link{S4BuyseTest-model.tables}} for exporting the number or percentage of favorable/unfavorable/neutral/uninformative pairs. \cr
#' \code{\link{sensitivity}} for performing a sensitivity analysis on the choice of the threshold(s). \cr
#' \code{\link{S4BuyseTest-plot}} for graphical display of the pairs across endpoints. \cr
#' \code{\link{getIid}} for exporting the first order H-decomposition. \cr
#' \code{\link{getPairScore}} for exporting the scoring of each pair. 
#' @keywords models
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
#' summary(BT) ## net benefit
#' model.tables(BT) ## export the table at the end of summary
#' summary(BT, percentage = FALSE)  
#' summary(BT, statistic = "winRatio") ## win Ratio
#' 
#' ## permutation instead of asymptotics to compute the p-value
#' \dontrun{
#'     BTperm <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3)
#' }
#' \dontshow{
#'     BTperm <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e1, trace = 0)
#' }
#' summary(BTperm)
#' summary(BTperm, statistic = "winRatio") 
#' 
#' ## same with parallel calculations
#' \dontrun{
#'     BTperm <- BuyseTest(treatment ~ TTE(eventtime, status = status), data=df.data,
#'                     method.inference = "permutation", n.resampling = 1e3, cpus = 8)
#'     summary(BTperm)
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
#' BTS <- BuyseTest(treatment ~ strata + TTE(eventtime, status = status), data=df.data)
#' summary(BTS)
#' 
#' #### several endpoints with a strata variable
#' ff <- treatment ~ strata + T(eventtime, status, 1) + B(toxicity) 
#' ff <- update(ff, 
#'             ~. + T(eventtime, status, 0.5) + C(score, 1) + T(eventtime, status, 0.25))
#' 
#' BTM <- BuyseTest(ff, data=df.data)
#' summary(BTM)
#' plot(BTM)
#' 
#' #### real example : veteran dataset of the survival package ####
#' ## Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#' ## scoring.rule = "Gehan"
#' 
#' if(require(survival)){
#' \dontrun{
#'   data(cancer, package = "survival") ## import veteran
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
#'   summary(BT_Peron)
#' }
#' }

## * BuyseTest (code)
##' @export
BuyseTest <- function(formula,
                      data,
                      scoring.rule = NULL,
                      pool.strata = NULL,
                      correction.uninf = NULL,
                      model.tte = NULL,
                      method.inference = NULL,
                      n.resampling = NULL,
                      strata.resampling = NULL,
                      hierarchical = NULL,
                      weightEndpoint = NULL,
                      weightObs = NULL,
                      neutral.as.uninf = NULL,
                      add.halfNeutral = NULL,
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
                      restriction = NULL,
                      strata = NULL){

    mycall <- match.call()
    option <- BuyseTest.options()

    ## ** compatibility with previous version
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
                              pool.strata = pool.strata,
                              model.tte = model.tte,
                              n.resampling = n.resampling,
                              strata.resampling = strata.resampling,
                              call = mycall,
                              neutral.as.uninf = neutral.as.uninf,
                              add.halfNeutral = add.halfNeutral,
                              operator = operator,
                              censoring = censoring,
                              option = option,
                              seed = seed,
                              strata = strata,
                              threshold = threshold,
                              restriction = restriction,
                              trace = trace,
                              treatment = treatment,
                              type = type,
                              weightEndpoint = weightEndpoint,
                              weightObs = weightObs,
                              envir = parent.frame())

    ## ** test arguments
    if(option$check){
        outTest <- do.call(testArgs, args = outArgs)        
    }
    ## ** initialization data
    ## WARNING when updating code: names in the c() must precisely match output of initializeData, in the same order
    out.name <- c("data","M.endpoint","M.status",
                  "index.C","index.T","weightObs","index.strata",
                  "level.treatment","level.strata", "pool.strata", "method.score", 
                  "grid.strata","n.obs","n.obsStrata","n.obsStrataResampling","cumn.obsStrataResampling","skeletonPeron",
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
                                        pool.strata = outArgs$pool.strata,
                                        treatment = outArgs$treatment,
                                        hierarchical = outArgs$hierarchical,
                                        copy = TRUE,
                                        keep.pairScore = outArgs$keep.pairScore,
                                        endpoint.TTE = outArgs$endpoint.TTE,
                                        status.TTE = outArgs$status.TTE,
                                        iidNuisance = outArgs$iidNuisance,
                                        weightEndpoint = outArgs$weightEndpoint,
                                        weightObs = outArgs$weightObs)

    if(option$check){
        if(outArgs$iidNuisance && any(outArgs$method.score == "CRPeron")){
            warning("Inference via the asymptotic theory  for competing risks when using the Peron's scoring rule has not been validating \n",
                    "Consider setting \'method.inference\' to \"none\", \"bootstrap\", or \"permutation\" \n")
        }
        ## if(outArgs$precompute && any(outArgs$method.score == "CRPeron")){
        ##     stop("Option \'precompute\' is not available for the Peron scoring rule in the competing risk case \n")
        ## }
        if(outArgs$precompute && any(outArgs$method.score == "CRPeron")){
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
                           method.inference = "none", ## do not use outArgs$method.inference as when it is equal to "bootstrap" or "permutation" we need the point estimate first.
                           pointEstimation = TRUE)

    if (outArgs$trace > 1) {
        cat("\n\n")
    }

    ## check number of pairs
    if(option$check){
        vec.nPair <- (outPoint$count_favorable + outPoint$count_unfavorable + outPoint$count_neutral + outPoint$count_uninf)[,1]
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
                                            level.strata = rownames(outArgs$grid.strata),
                                            n.strata = NROW(outArgs$grid.strata),
                                            endpoint = outArgs$endpoint,
                                            threshold = outArgs$threshold,
                                            restriction = outArgs$restriction)
    }
    
    ## ** Inference
    if((outArgs$method.inference != "none") && (outArgs$trace > 1)){
        do.call(printInference, args = outArgs)
    }

    outResampling <- NULL
    if(outArgs$method.inference == "u statistic"){
        ## done in the C++ code        
    }else if(outArgs$method.inference == "u statistic bebu"){
        if(outArgs$keep.pairScore == FALSE){
            stop("Argument \'keep.pairScore\' needs to be TRUE when argument \'method.inference\' is \"u-statistic-bebu\" \n")
        }

        ## direct computation of the variance
        outCovariance <- inferenceUstatisticBebu(tablePairScore = outPoint$tableScore,
                                                 order = option$order.Hprojection,
                                                 weightEndpoint = outArgs$weightEndpoint,
                                                 n.pairs = outPoint$n_pairs,
                                                 n.C = length(envirBT$outArgs$index.C),
                                                 n.T = length(envirBT$outArgs$index.T),
                                                 level.strata = outArgs$level.strata,
                                                 n.strata = NROW(outArgs$grid.strata),
                                                 endpoint = outArgs$endpoint)

        outPoint$covariance <- outCovariance$Sigma
        attr(outArgs$method.inference,"Hprojection") <- option$order.Hprojection
    }else if(outArgs$method.inference == "varexact permutation"){

        if(!is.null(outArgs$weightObs) && any(abs(outArgs$weightObs-1)>1e-10)){
            warning("Argument \'weightObs\' is being ignored when computing the exact permutation variance. \n")
        }

        outVariance <- inferenceVarPermutation(data = as.data.frame(data),
                                               treatment = outArgs$treatment,
                                               level.treatment = outArgs$level.treatment,
                                               weightStrata = outPoint$weightStrata,
                                               scoring.rule = outArgs$scoring.rule,
                                               pool.strata = outArgs$pool.strata,
                                               correction.uninf = outArgs$correction.uninf,
                                               model.tte = outArgs$model.tte,
                                               hierarchical = outArgs$hierarchical,
                                               weightEndpoint = outArgs$weightEndpoint,
                                               neutral.as.uninf = outArgs$neutral.as.uninf,
                                               add.halfNeutral = outArgs$add.halfNeutral,
                                               endpoint = outArgs$endpoint,
                                               type = outArgs$type,
                                               threshold = outArgs$threshold,                      
                                               status = outArgs$status,
                                               operator = outArgs$operator,
                                               censoring = outArgs$censoring,
                                               restriction = outArgs$restriction,
                                               strata = outArgs$strata)

        outPoint$covariance <- outVariance
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
    keep.args <- c("index.T", "index.C", "index.strata", "type","endpoint","level.strata","level.treatment","scoring.rule","hierarchical","neutral.as.uninf","add.halfNeutral",
                   "correction.uninf","method.inference","method.score","strata","threshold","restriction","weightObs","weightEndpoint","pool.strata","grid.strata","n.resampling")
    mycall2 <- setNames(as.list(mycall),names(mycall))
    if(!missing(formula)){
        mycall2$formula <- formula ## change name of the variable into actual value
    }
    mycall2$data <- data ## change name of the variable into actual value
    BuyseTest.object <- do.call("S4BuyseTest", args = c(list(call = mycall2),
                                                        outPoint, outArgs[keep.args], outResampling))
    
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
                             threshold = envir$outArgs$threshold,
                             restriction = envir$outArgs$restriction,
                             level.strata = envir$outArgs$level.strata,
                             grid.strata = envir$outArgs$grid.strata,
                             strata = envir$outArgs$strata,
                             precompute = envir$outArgs$precompute,
                             iidNuisance = envir$outArgs$iidNuisance * iid,
                             out = envir$outArgs$skeletonPeron,
                             fitter = envir$outArgs$fitter.model.tte,
                             efron = envir$outArgs$scoring.rule==2,
                             args = envir$outArgs$args.model.tte)

        index.test <- which(envir$outArgs$method.score == "SurvPeron")
        if(!grepl("permutation|bootstrap",method.inference) && envir$outArgs$correction.uninf>0 && length(index.test)>0 && all(is.na(envir$outArgs$restriction))){
            maxLastSurv <- setNames(sapply(outSurv$lastSurv[index.test],max),envir$outArgs$endpoint[index.test])[!duplicated(envir$outArgs$endpoint[index.test])]
            Wtau <- BuyseTest.options("warning.correction")
            if(any(maxLastSurv>Wtau)){
                warning("Some of the survival curves for endpoint(s) \"",paste(names(which(maxLastSurv>Wtau)),collapse = "\", \""),"\" are unknown beyond a survival of ",Wtau,".\n",
                        "The correction of uninformative pairs assume that uninformative pairs would on average behave like informative pairs. \n",
                        "This can be a strong assumption and have substantial impact when the tail of the survival curve is unknown. \n")
            }
        }
    }

    ## ** Restriction
    if(any(!is.na(envir$outArgs$restriction))){
        ## index restricted endpoint
        index.rendpoint <- setdiff(which(!is.na(envir$outArgs$restriction)), ## non-NA value
                                   which(duplicated(envir$outArgs$index.endpoint))) ## not already visitied
        for(iE in index.rendpoint){ ## iE <- 1
            iRestriction <- envir$outArgs$restriction[iE]
            iStatus <- envir$outArgs$index.status[iE]+1
            if(envir$outArgs$operator[iE]==1){ ## ">0"
                if(envir$outArgs$method.score[iE] %in% c("TTEgehan","SurvPeron","CRPeron")){ ## right censoring
                    envir$outArgs$M.status[envir$outArgs$M.endpoint[,iE]>iRestriction,iStatus] <- 1/2
                    envir$outArgs$M.status[envir$outArgs$M.endpoint[,iE]==iRestriction & envir$outArgs$M.status[,iStatus]==0,iStatus] <- 1/2 ## rm censoring when restriction at the censoring time
                }
                envir$outArgs$M.endpoint[envir$outArgs$M.endpoint[,iE]>iRestriction,iE] <- iRestriction
            }else if(envir$outArgs$operator[iE]==-1){ ## "<0"
                if(envir$outArgs$method.score[iE] %in% c("TTEgehan2")){ ## left censoring
                    envir$outArgs$M.status[envir$outArgs$M.endpoint[,iE]<iRestriction,iStatus] <- 1/2
                    envir$outArgs$M.status[envir$outArgs$M.endpoint[,iE]==iRestriction & envir$outArgs$M.status[,iStatus]==0,iStatus] <- 1/2 ## rm censoring when restriction at the censoring time
                }
                envir$outArgs$M.endpoint[envir$outArgs$M.endpoint[,iE]<iRestriction,iE] <- iRestriction
            }
        }
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
                                 threshold0 = attr(envir$outArgs$threshold,"original")==0,
                                 restriction = envir$outArgs$restriction,
                                 weightEndpoint = envir$outArgs$weightEndpoint,
                                 weightObs = envir$outArgs$weightObs,
                                 method = sapply(envir$outArgs$method.score, switch, "continuous" = 1, "gaussian" = 2, "TTEgehan" = 3, "TTEgehan2" = 4, "SurvPeron" = 5, "CRPeron" = 6),
                                 pool = envir$outArgs$pool.strata,
                                 op = envir$outArgs$operator,
                                 D = envir$outArgs$D,
                                 D_UTTE = envir$outArgs$D.UTTE,
                                 grid_strata = envir$outArgs$grid.strata,
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
                                 addHalfNeutral = envir$outArgs$add.halfNeutral,
                                 keepScore = (pointEstimation && envir$outArgs$keep.pairScore),
                                 precompute = envir$outArgs$precompute,
                                 match = !is.null(envir$outArgs$strata) && attr(envir$outArgs$strata,"match"),
                                 returnIID = c(iid,envir$outArgs$iidNuisance),
                                 debug = envir$outArgs$debug
                                 ))

    
    ## ** export
    if(pointEstimation){
        if(envir$outArgs$keep.survival){ ## useful to test initSurvival 
            resBT$tableSurvival <- outSurv
        }
        return(resBT)
    }else{
        ## index <- 5
        ## resBT$Delta[,index]
        ## sum(resBT$delta[,,index][,1] * resBT$weightStrata)
        return(list(n = rbind(T = lengths(outSample$ls.indexT), C = lengths(outSample$ls.indexC)),
                    delta = resBT$delta,
                    Delta = resBT$Delta,
                    weightStrata = resBT$weightStrata,
                    covariance = resBT$covariance))
    }
}

## * calcSample
calcSample <- function(envir, method.inference){

    ## ** initialization
    nlevel.strata <- length(envir$outArgs$level.strata)
    out <- list(## rows in M.endpoint/M.status corresponding to observations from the control/treatment group (not unique when boostraping)
        ls.indexC = vector(mode = "list", length = nlevel.strata), 
        ls.indexT = vector(mode = "list", length = nlevel.strata),
        ## identifier for each observation from the control/treatment group (unique even when boostrap)
        ls.posC = vector(mode = "list", length = nlevel.strata),
        ls.posT = vector(mode = "list", length = nlevel.strata),
        ## dataset
        data = data.table::data.table()
    )

    if(method.inference %in% c("none","u statistic")){

        ## ** no resampling
        if(nlevel.strata==1){        
            out$ls.indexC[[1]] <- envir$outArgs$index.C - 1
            out$ls.indexT[[1]] <- envir$outArgs$index.T - 1
        }else{        
            for(iStrata in 1:nlevel.strata){ ## iStrata <- 1  
                out$ls.indexC[[iStrata]] <- intersect(envir$outArgs$index.C, envir$outArgs$index.strata[[iStrata]]) - 1
                out$ls.indexT[[iStrata]] <- intersect(envir$outArgs$index.T, envir$outArgs$index.strata[[iStrata]]) - 1
            }
        }
        out$ls.posC <- out$ls.indexC
        out$ls.posT <- out$ls.indexT

        if(envir$outArgs$scoring.rule>0){
            out$data <- data.table::data.table(envir$outArgs$data,
                                               envir$outArgs$M.endpoint,
                                               envir$outArgs$M.status)
        }
    }else{

        ## ** stratified resampling
        n.strataResampling <- length(envir$outArgs$n.obsStrataResampling)
        index.resampling <- NULL
        for (iSR in 1:n.strataResampling) { ## iSR <- 1
            index.resampling <- c(index.resampling,
                                  envir$outArgs$cumn.obsStrataResampling[iSR] + sample.int(envir$outArgs$n.obsStrataResampling[iSR], replace = attr(method.inference, "bootstrap")))
        }

        ## ** reconstruct groups
        ## index: index of the new observations in the old dataset by treatment group
        ## pos: unique identifier for each observation
        if(nlevel.strata==1){ ## no strata
            
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

            for(iStrata in 1:nlevel.strata){ ## iStrata <- 1  
                ## index of the new observations in the old dataset by treatment group
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
            ## re-create original strata variables                
            if(nlevel.strata>1 && length(envir$outArgs$strata)==1 && envir$outArgs$strata %in% names(out$data) == FALSE){ 
                out$data[[envir$outArgs$strata]] <- envir$outArgs$level.strata[out$data[["..strata.."]]]
                if(is.factor(envir$outArgs$data[[envir$outArgs$strata]])){
                    out$data[[envir$outArgs$strata]] <- factor(out$data[[envir$outArgs$strata]], levels = envir$outArgs$level.strata)
                }
            }else if(nlevel.strata>1 && length(envir$outArgs$strata)>1 && all(envir$outArgs$strata %in% names(out$data) == FALSE)){
                grid.strata <- unique(envir$outArgs$data[,.SD,.SDcols = c("..strata..",envir$outArgs$strata)])

                out$data <- cbind(out$data,grid.strata[match(out$data[["..strata.."]], grid.strata[["..strata.."]]),.SD,.SDcols = envir$outArgs$strata] )
            }
        }

    }
    return(out)
}









