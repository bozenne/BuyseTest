## * Documentation - BuyseTest
#' @name BuyseTest
#' @title Generalized Pairwise Comparisons
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event outcomes.
#' @param formula a symbolic description of the model to be fitted. The response variable should be a binary variable defining the treatment arms. 
#' The rest of the formula should indicate the strata variables (if any) and the endpoints by order of priority. \emph{formula}.
#' @param data A \code{data.frame} containing the variables.
#' @param treatment the name of the treatment variable identifying the control and the experimental group. \emph{character}.
#' @param endpoint the name of the endpoint variable(s). \emph{character vector}.
#' @param threshold the thresholds, one for each endpoint variable. \emph{numeric vector}.
#' @param strata the name of the strata variable(s). \emph{numeric vector}. 
#' @param censoring the name of the censoring variable(s), one for each endpoint. \emph{character vector}.
#' @param type the type of each endpoint. \emph{character vector}. Can be \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.
#' @param method Is defined when at least one time-to-event outcome is analyzed. Defines the method used to handle pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations.  Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}. See details. 
#' @param neutralAsUninf Should paired classified as neutral be re-analysed using endpoints of lower priority. \emph{logical}.
#' @param n.bootstrap the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}. 
#' @param prob.alloc the resampling probability for assignement to the experimental group in the bootstrap samples. \emph{double}. Can also be \code{NULL} to use the proportion of patients in the experimental group.
#' @param stratified Should the boostrap be a stratified bootstrap? \emph{logical}. 
#' @param alternative a \emph{character} specifying the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. 
#' @param seed the seed to consider for the bootstrap. \emph{integer}. 
#' @param cpus the number of CPU to use. \emph{integer}. 
#' @param trace Should the execution of the function be traced ? \emph{integer}. 
#' 
#' @details 
#' \bold{treatment:} The variable corresponding to \code{treatment} in data must have only two levels (e.g. \code{0} and \code{1}). \cr
#' \bold{endpoint, threshold, censoring, and type:} Arguments \code{endpoint}, \code{threshold}, \code{censoring}  and \code{type} must have the same length. \cr
#' \code{threshold} must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. \cr
#' \code{censoring} must be \code{NA} for binary or continuous endpoints and indicate a variable in data for time to event endpoints. 
#' Short forms for endpoint \code{type} are \code{"bin"} (binary endpoint), \code{"cont"} (continuous endpoint), \code{"TTE"} (time-to-event endpoint). 
#' 
#' \bold{Bootstrap:} The number of bootstrap replications (argument \code{n.bootstrap}) must be specified to enable the computation of the confidence intervals and the p.value. 
#' A large number of bootstrap samples (e.g. \code{n.bootstrap=10000})  are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 
#' 
#' \bold{trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the bootstrap,  and \code{0} remains silent.
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
BuyseTest <- function(formula, data, 
                      treatment = NULL, endpoint = NULL, type = NULL, threshold = NULL, censoring = NULL, strata = NULL, 
                      method = BuyseTest.options()$method, neutralAsUninf = BuyseTest.options()$neutralAsUninf,
                      n.bootstrap = BuyseTest.options()$n.bootstrap, prob.alloc = NULL, stratified = FALSE, alternative = "two.sided", 
                      seed = BuyseTest.options()$seed, cpus = BuyseTest.options()$cpus, trace = BuyseTest.options()$trace){
  
  
  ## ** 1- data management + tests
  Buysecall <- match.call()
  
  ## 
  if(!missing(formula)){
    argnames <- c("treatment", "endpoint", "type", "threshold", "censoring", "strata")
    if(any(names(Buysecall) %in% argnames)){
      warning("BuyseTest : arguments \'",paste(names(Buysecall)[names(Buysecall) %in% argnames], collapse = " ")," have been ignored \n",
              "when specified, only argument \'formula\' is used \n")
    }
    
    resFormula <- initFormula(formula)
    treatment <- resFormula$treatment
    type <- resFormula$type
    endpoint <- resFormula$endpoint
    
    if(any(unlist(lapply(resFormula$threshold, length))>1)){
      indexF <- which(unlist(lapply(resFormula$threshold, length))>1)
      stop("BuyseTest: several thresholds found for endpoint",if(length(indexF)>1){"s"}," number ",paste(indexF, collapse = " "),"\n",
           "only one threshold can be used per priority, consider adding an additional endpoint using + \n")
    }
    
    threshold <- unlist(resFormula$threshold)
    censoring <- resFormula$censoring
    strata <- resFormula$strata
  }
  
  ## *** Treatment: extract the 2 levels
  validCharacter(treatment, validLength = 1, method = "BuyseTest")
  validNames(data, requiredValues = treatment, validLength = NULL, method = "BuyseTest")
  levels.treatment <- levels(as.factor(data[[treatment]])) # extraction of the levels of the treatment variable
  if (length(levels.treatment) != 2) {
    stop("BuyseTest : wrong specification of \'treatment\' \n",
         "the corresponding column in \'data\' must have exactly 2 levels \n",
         "proposed levels : ",paste(levels.treatment,collapse = " "),"\n")
  }
  
  ## *** endpoint
  validNames(data, requiredValues = endpoint, validLength = NULL, method = "BuyseTest")
  D <- length(endpoint) # number of endpoints
  
  ## *** type: convert type to numeric and count the number of endpoints
  validCharacter(type, validValues = c(1:3,"bin","binary","cont","continuous","TTE","timeToEvent"), validLength = D, method = "BuyseTest")
  type[type %in% c("binary","bin")] <- "1" 
  type[type %in% c("continuous","cont")] <- "2"
  type[type %in% c("timeToEvent","TTE")] <- "3"
  type <- as.numeric(type) # type is an integer equal to 1 (binary endpoint), 2 (continuous endpoint) or 3 (time to event endpoint)
  
  D.TTE <- sum(type == 3) # number of time to event endpoints
  
  ## *** censoring: (see .Rd internal-intilisation, section details for details)
  censoring <- initCensoring(censoring = censoring, endpoint = endpoint, type = type, D = D, D.TTE = D.TTE,
                             treatment = treatment, strata = strata)
  
  validNames(data, name1 = "data", requiredValues = censoring, validLength = NULL, refuse.NULL = FALSE, method = "BuyseTest")
  
  ## *** data: split the data according to the two levels
  if (!is.null(strata)) {
    validNames(data, name1 = "strata", requiredValues = strata, validLength = NULL, method = "BuyseTest")
  }
  
  if (data.table::is.data.table(data)) {
    dataT <- data[data[[treatment]] == levels.treatment[2], c(endpoint, strata, censoring), with = FALSE]
    dataC <- data[data[[treatment]] == levels.treatment[1], c(endpoint, strata, censoring), with = FALSE]
  }else if (is.data.frame(data)) {
    dataT <- data.table::as.data.table(data[data[[treatment]] == levels.treatment[2], c(endpoint, strata, censoring), drop = FALSE])
    dataC <- data.table::as.data.table(data[data[[treatment]] == levels.treatment[1], c(endpoint, strata, censoring), drop = FALSE])
  }else{
    stop("BuyseTest : data must be a data.frame or a data.table \n")
  }
  n.Treatment <- NROW(dataT) # number of patient in the treatment arm
  n.Control <- NROW(dataC) # number of patient in the control arm
  
  ## *** threshold: (see .Rd internal-intilisation, section details for details)
  threshold <- initThreshold(threshold = threshold, type = type, D = D,
                             method = method, endpoint = endpoint)
  
  ## *** strata: (see .Rd internal-intilisation, section details for details)
  res <- initStrata(strata = strata,
                    dataT = dataT, dataC = dataC, n.Treatment = n.Treatment, n.Control = n.Control,
                    endpoint = endpoint, censoring = censoring)
  
  index.strataT <- res$index.strataT 
  index.strataC <- res$index.strataC 
  n.strata <- res$n.strata 
  levels.strata <- res$levels.strata 
  
  ## *** method
  validCharacter(method, validLength = 1, validValues = c("Gehan","Peto","Efron","Peron"), method = "BuyseTest")
  
  method <- switch(method,
                   "Gehan" = 0,
                   "Peto" = 1,
                   "Efron" = 2,
                   "Peron" = 3)
  
  if (D.TTE == 0) {
    method <- 0
    if ("method" %in% names(Buysecall) && trace > 0) {
      message("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
    }
  }
  
  ## *** alternative
  validCharacter(alternative, validLength = 1, validValues = c("two.sided", "less", "greater"), method = "BuyseTest")
  
  ## *** n.bootstrap
  validInteger(n.bootstrap, validLength = 1, min = 0, method = "BuyseTest")
  
  ## *** proba
  if (is.null(prob.alloc)) { # if prob.alloc is not set by the user
    prob.alloc <- n.Treatment/(n.Treatment + n.Control)  # it is set to the proportion of patients in the treatment arm.
  }else{
    validNumeric(prob.alloc, validLength = 1, min = 0, max = 1, method = "BuyseTest")
  }
  
  # stratified
  validLogical(stratified, validLength = 1, method = "BuyseTest")
  
  ## *** seed
  validInteger(seed, validLength = 1, refuse.NULL = FALSE, min = 1, method = "BuyseTest")
  
  ## *** cpu
  if (cpus == "all") { 
    cpus <- parallel::detectCores() # this function detect the number of CPU cores 
  }else{# if several cpus are intended to be used, check this correspond to a valid number of CPU cores
    validInteger(cpus, validLength = 1, validValues = 1:parallel::detectCores(), method = "BuyseTest")
  }
  
  ## *** data (again)
  res <- initData(dataT = dataT, dataC = dataC,type = type,endpoint = endpoint, D = D, censoring = censoring,
                  index.strataT = index.strataT, index.strataC = index.strataC, n.strata = n.strata,                  
                  method = method,D.TTE = D.TTE, threshold = threshold, Wscheme = NULL,
                  test = TRUE, trace = trace)
  
  M.Treatment <- res$M.Treatment
  M.Control <- res$M.Control 
  M.delta_Treatment <- res$M.delta_Treatment
  M.delta_Control <- res$M.delta_Control
  Wscheme <- res$Wscheme
  threshold_TTEM1 <- res$threshold_TTEM1
  index_survivalM1 <- res$index_survivalM1
  list_survivalT <- res$list_survivalT
  list_survivalC <- res$list_survivalC
  
  ## ** 2- General display
  if (trace > 1) {
    printGeneral(levels.treatment = levels.treatment,
                 levels.strata = levels.strata, n.strata = n.strata,
                 endpoint = endpoint, threshold = threshold, censoring = censoring, type = type, D = D, D.TTE = D.TTE,
                 method = method, neutralAsUninf = neutralAsUninf,
                 Wscheme = if (D>1 && method %in% 1:3) {Wscheme} else {NULL}, 
                 threshold_TTEM1 = if (D>1 && method %in% 1:3) {threshold_TTEM1} else {NULL})
  }
  
  ## ** 3- Punctual estimation
  if (trace > 1) {cat("Punctual estimation \n")}
  
  time <- system.time({
    resPonctual <- GPC_cpp(Treatment = M.Treatment, Control = M.Control, threshold = threshold, survEndpoint = (type == 3),
                           delta_Treatment = M.delta_Treatment, delta_Control = M.delta_Control,
                           D = D, returnIndex = TRUE,
                           strataT = index.strataT, strataC = index.strataC, n_strata = n.strata, n_TTE = D.TTE,
                           Wscheme = Wscheme,index_survivalM1 = index_survivalM1, threshold_TTEM1 = threshold_TTEM1, 
                           list_survivalT = list_survivalT, list_survivalC = list_survivalC,
                           methodTTE = method, neutralAsUninf = neutralAsUninf
    )
  })
  if (trace > 1) {cat("   # done \n")}
  
  ## ** 4- Transfomration into BuyseRes object
  BuyseRes.object <- BuyseRes(
    count_favorable = resPonctual$count_favorable,      
    count_unfavorable = resPonctual$count_unfavorable,
    count_neutral = resPonctual$count_neutral,    
    count_uninf = resPonctual$count_uninf,
    delta = list(netChance = resPonctual$delta_netChance, winRatio = resPonctual$delta_winRatio),
    Delta = list(netChance = resPonctual$Delta_netChance, winRatio = resPonctual$Delta_winRatio),
    endpoint = endpoint,
    index_neutralT = resPonctual$index_neutralT,
    index_neutralC = resPonctual$index_neutralC,
    index_uninfT = resPonctual$index_uninfT,
    index_uninfC = resPonctual$index_uninfC,
    levels.treatment = levels.treatment,
    n_pairs = resPonctual$n_pairs,
    strata = levels.strata,
    threshold = threshold
  )
  
  ## ** 5- Bootstrap 
  if (n.bootstrap > 0) {
    if (trace > 1) {
      printBoostrap(prob.alloc, n.bootstrap, stratified, cpus, time, seed)
    }
    
    n.eachStrataT <- unlist(lapply(index.strataT, length)) # those variables are used during the boostrap and passed to the wrapper through the environment
    n.eachStrataC <- unlist(lapply(index.strataC, length))
    nCumSum.strataControl <- cumsum(c(1,n.eachStrataC))
    nCumSum.strataTreatment <- cumsum(c(1,n.eachStrataT))
       
    delta_boot <- list(netChance = array(NA, dim = c(n.strata, D, n.bootstrap)),
                       winRatio = array(NA, dim = c(n.strata, D, n.bootstrap))
    )
    Delta_boot <- list(netChance = matrix(NA, nrow = D, ncol = n.bootstrap),
                       winRatio = matrix(NA, nrow = D, ncol = n.bootstrap)
    )
    
    ## *** computation
    calcBootstrap(environment())
    
    ## *** post treatment
    if (trace > 1) {cat("Post-Treatment \n")}
    resCI <- calcCI(Delta = BuyseRes.object@Delta, Delta_boot = Delta_boot, 
                    endpoint = endpoint, D = D, alternative = alternative, alpha = 1 - BuyseTest.options()$conf.level,
                    n.bootstrap = n.bootstrap, cpus = cpus, trace = trace)
    
    ## *** update BuyseRes object 
    if(BuyseTest.options()$keep.bootstrap){
      BuyseRes.object@delta_boot <- delta_boot
    }
    
    BuyseRes.object@Delta_quantile <- resCI$Delta_quantile
    BuyseRes.object@n_bootstrap <- resCI$n_bootstrap
    BuyseRes.object@p.value <- resCI$p.value
    validObject(BuyseRes.object)
  }
  
  ## ** 6- export
  return(BuyseRes.object)
}
