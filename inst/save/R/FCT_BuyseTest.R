#' @name BuyseTest
#' @title Generalized Pairwise Comparisons
#' @aliases BuyseTest
#' 
#' @description Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event outcomes.
#' @param data A \code{data.frame} containing the variables.
#' @param treatment the name of the treatment variable identifying the control and the experimental group. \emph{character}.
#' @param endpoint the name of the endpoint variable(s). \emph{character vector}.
#' @param threshold the thresholds, one for each endpoint variable. \emph{numeric vector}. Default is \code{NULL} indicating no threshold.
#' @param strata the name of the strata variable(s). \emph{numeric vector}. Default is \code{NULL} indicating only one strata.
#' @param censoring the name of the censoring variable(s), one for each endpoint. \emph{character vector}. Default is \code{NULL}.
#' @param type the type of each endpoint. \emph{character vector}. Can be \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.
#' @param method Is defined when at least one time-to-event outcome is analyzed. Defines the method used to handle pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations.  Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}. See details. 
#' @param n.bootstrap the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}. Default is \code{0} meaning no bootstrap (and thus only ponctual estimation).
#' @param prob.alloc the resampling probability for assignement to the experimental group in the bootstrap samples. \emph{double}. Default is \code{NULL} indicating to use the proportion of patients in the experimental group.
#' @param stratified Should the boostrap be a stratified bootstrap? \emph{logical}. Default is \code{FALSE}.
#' @param alternative a \emph{character} specifying the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. Default is \code{"two.sided"}.
#' @param seed the seed to consider for the bootstrap. \emph{integer}. Default is \code{10}.
#' @param cpus the number of CPU to use. \emph{integer}. Default is \code{1}.
#' @param trace Should the execution of the function be traced ? \emph{integer}. Default is \code{3}.
#' 
#' @details 
#' \bold{Treatment:} The variable corresponding to \code{treatment} in data must have only two levels (e.g. \code{0} and \code{1}). \cr
#' \bold{Endpoint, threshold, censoring, and type:} Arguments \code{endpoint}, \code{threshold}, \code{censoring}  and \code{type} must have the same length. \cr
#' \code{threshold} must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. \cr
#' \code{censoring} must be \code{NA} for binary or continuous endpoints and indicate a variable in data for time to event endpoints. 
#' Short forms for endpoint \code{type} are \code{"bin"} (binary endpoint), \code{"cont"} (continuous endpoint), \code{"TTE"} (time-to-event endpoint). 
#' 
#' \bold{Bootstrap:} The number of bootstrap replications (argument \code{n.bootstrap}) must be specified to enable the computation of the confidence intervals and the p.value. 
#' A large number of bootstrap samples (e.g. \code{n.bootstrap=10000})  are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 
#' 
#' \bold{Trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the bootstrap,  and \code{0} remains silent.
#' 
#' \bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus. The parallelization relies on the \emph{snowfall} package (function \emph{sfClusterApplyLB}). The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .
#' 
#' \bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
#' Uninformative pairs correspond to pairs for which the censoring prevend from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.
#' 
#' \bold{Method:} Pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}). Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions. \code{method="Peto"} estimate these probabilities using the common Kaplan-Meier estimator of the survival function for treated and control patients. \code{method="Efron"}, and \code{method="Peron"} estimate these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. When the largest observation is censored, it is not possible to estimate the survival probability by the Kaplan-Meier estimator beyond this time point.  \code{method="Efron"} treats the largest observations in each patient group as if it were uncensored. \code{method="Peron"} treats the probability of survival beyond the last observation as NA, resulting in a non null probability that the pair is uninformative    
#' 
#' @return An \R object of class \code{\linkS4class{BuyseRes}}.
#' 
#' @references 
#' Marc Buyse (2010) Generalized pairwise comparisons of prioritized endpoints in the two-sample problem. \emph{Statistics in Medicine} \bold{vol. 29} 3245-3257 \cr
#' Efron B (1967) The two sample problem with censored data \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} \bold{vol. 4} 831-583 \cr
#' Peto R, Peto J (1972) Asymptotically efficient rank invariant test procedures \emph{J R Stat Soc A} \bold{vol. 135(2)} 185-198 \cr
#' Gehan EA (1965) A generalized two-sample Wilcoxon test for doubly censored data \emph{Biometrika} \bold{vol. 52(3)} 650-653 \cr
#'
#' @seealso 
#' \code{\link{BuyseRes-summary}} for a summary of the results of generalized pairwise comparison. \cr
#' \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
#' \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
#' 
#' @example
#' examples/EX_BuyseTest.R
#'     
#' @keywords function BuyseTest
#' @export
BuyseTest <- function(data, treatment, endpoint, type, threshold = NULL, strata = NULL, censoring = NULL, method = "Peron",
                      n.bootstrap = 0, prob.alloc = NULL, stratified = FALSE, alternative = "two.sided", seed = 10, cpus = 1, trace = 3){
  
  
  #### 1- data management + tests ####
  Buysecall <- match.call()
  
  ## Treatment: extract the 2 levels
  validCharacter(treatment, validLength = 1, method = "BuyseTest")
  validNames(data, requiredValues = treatment, validLength = NULL, method = "BuyseTest")
  levels.treatment <- levels(as.factor(data[[treatment]])) # extraction of the levels of the treatment variable
  if (length(levels.treatment) != 2) {
    stop("BuyseTest : wrong specification of \'treatment\' \n",
         "the corresponding column in \'data\' must have exactly 2 levels \n",
         "proposed levels : ",paste(levels.treatment,collapse = " "),"\n")
  }
  
  ## endpoint
  validNames(data, requiredValues = endpoint, validLength = NULL, method = "BuyseTest")
  D <- length(endpoint) # number of endpoints
  
  ## type: convert type to numeric and count the number of endpoints
  validCharacter(type, validValues = c("bin","binary","cont","continuous","TTE","timeToEvent"), validLength = D, method = "BuyseTest")
  type[type %in% c("binary","bin")] <- "1" 
  type[type %in% c("continuous","cont")] <- "2"
  type[type %in% c("timeToEvent","TTE")] <- "3"
  type <- as.numeric(type) # type is an integer equal to 1 (binary endpoint), 2 (continuous endpoint) or 3 (time to event endpoint)
  
  D.TTE <- sum(type == 3) # number of time to event endpoints
  
  ## censoring: (see .Rd internal-intilisation, section details for details)
  censoring <- initCensoring(censoring = censoring, endpoint = endpoint, type = type, D = D, D.TTE = D.TTE,
                             treatment = treatment, strata = strata)
  validNames(data, requiredValues = censoring, validLength = NULL, refuse.NULL = FALSE, method = "BuyseTest")
  
  ## data: split the data according to the two levels
  if (!is.null(strata)) {
    validNames(data, requiredValues = strata, validLength = NULL, method = "BuyseTest")
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
  
  ## threshold: (see .Rd internal-intilisation, section details for details)
  threshold <- initThreshold(threshold = threshold, type = type, D = D,
                             method = method, endpoint = endpoint)
  
  ## strata: (see .Rd internal-intilisation, section details for details)
  res <- initStrata(strata = strata,
                    dataT = dataT, dataC = dataC, n.Treatment = n.Treatment, n.Control = n.Control,
                    endpoint = endpoint, censoring = censoring)
  
  index.strataT <- res$index.strataT 
  index.strataC <- res$index.strataC 
  n.strata <- res$n.strata 
  levels.strata <- res$levels.strata 
  
  ## method
  validCharacter(method, validLength = 1, validValues = c("Gehan","Peto","Efron","Peron"), method = "BuyseTest")
  if (D.TTE == 0) {
    method <- "Gehan"
    if ("method" %in% names(Buysecall) && trace > 0) {
      message("NOTE : there is no survival endpoint, \'method\' argument is ignored \n")
    }
  }
  
  ## alternative
  validCharacter(alternative, validLength = 1, validValues = c("two.sided", "less", "greater"), method = "BuyseTest")
  
  ## n.bootstrap
  validInteger(n.bootstrap, validLength = 1, min = 0, method = "BuyseTest")
  
  ## proba
  if (is.null(prob.alloc)) { # if prob.alloc is not set by the user
    prob.alloc <- n.Treatment/(n.Treatment + n.Control)  # it is set to the proportion of patients in the treatment arm.
  }else{
    validNumeric(prob.alloc, validLength = 1, min = 0, max = 1, method = "BuyseTest")
  }
  
  # stratified
  validLogical(stratified, validLength = 1, method = "BuyseTest")
  
  ## seed
  validInteger(seed, validLength = 1, refuse.NULL = FALSE, min = 1, method = "BuyseTest")
  
  ## cpu
  if (cpus == "all") { 
    cpus <- parallel::detectCores() # this function detect the number of CPU cores 
  }else{# if several cpus are intended to be used, check this correspond to a valid number of CPU cores
    validInteger(cpus, validLength = 1, validValues = 1:parallel::detectCores(), method = "BuyseTest")
  }
  
  ## data (again)
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
  
  #### 2- General display #### 
  if (trace > 1) {
    printGeneral(levels.treatment = levels.treatment,
                 levels.strata = levels.strata, n.strata = n.strata,
                 endpoint = endpoint, threshold = threshold, censoring = censoring, type = type, D = D, D.TTE = D.TTE,
                 method = method, 
                 Wscheme = if (method %in% c("Peto","Efron","Peron")) {Wscheme} else {NULL}, 
                 threshold_TTEM1 = if (method %in% c("Peto","Efron","Peron")) {Wscheme} else {NULL})
  }
  
  #### 2- Punctual estimation ####
  if (trace > 1) {cat("Punctual estimation \n")}
  
  if (method %in% c("Peto","Efron","Peron")) {  
    time <- system.time({
      resPonctual <- BuyseTest_PetoEfronPeron_cpp(Treatment = M.Treatment, Control = M.Control, threshold = threshold, type = type,
                                                  delta_Treatment = M.delta_Treatment, delta_Control = M.delta_Control,
                                                  D = D, returnIndex = TRUE,
                                                  strataT = index.strataT, strataC = index.strataC, n_strata = n.strata, n_TTE = D.TTE,
                                                  Wscheme = Wscheme,index_survivalM1 = index_survivalM1, threshold_TTEM1 = threshold_TTEM1, 
                                                  list_survivalT = list_survivalT, list_survivalC = list_survivalC,
                                                  PEP = which(c("Peto", "Efron", "Peron") == method)
      )
    })
  }else if (method == "Gehan") {
    time <- system.time({
      resPonctual <-   BuyseTest_Gehan_cpp(Treatment = M.Treatment, Control = M.Control, threshold = threshold, type = type,
                                           delta_Treatment = M.delta_Treatment, delta_Control = M.delta_Control,
                                           D = D, returnIndex = TRUE,
                                           strataT = index.strataT, strataC = index.strataC, n_strata = n.strata, n_TTE = D.TTE)    
    })
  }
  
  if (trace > 1) {cat("   # done \n")}
 
  #### 3- transfomration into BuyseRes object ####
  BuyseRes.object <- BuyseRes(
    delta = resPonctual$delta, 
    count_favorable = resPonctual$count_favorable,      
    count_unfavorable = resPonctual$count_unfavorable,
    count_neutral = resPonctual$count_neutral,    
    count_uninf = resPonctual$count_uninf,
    index_neutralT = resPonctual$index_neutralT,
    index_neutralC = resPonctual$index_neutralC,
    index_uninfT = resPonctual$index_uninfT,
    index_uninfC = resPonctual$index_uninfC,
    n_pairs = resPonctual$n_pairs,
#     delta_boot = array(NA,dim = c(n.strata,D,1)), 
#     p.value = rep(NA,D),    
#     Delta_quantile = matrix(NA,nrow = 2, ncol = D, dimnames = list(c("2.5%","97.5%"))),
    endpoint = endpoint,
    threshold = threshold,
    strata = levels.strata,
    levels.treatment = levels.treatment
  )
  
  #### 4- Bootstrap ####
  if (n.bootstrap == 0) {
    if (trace > 1) {
      message("*** only ponctual estimation requested ***\n",
              "set \'n.bootstrap\' argument to a strictly postive value for confidence interval and p.value \n"
      )}
    
  }else{
    
    if (trace > 1) {
      printBoostrap(prob.alloc, n.bootstrap, stratified, cpus, time, seed)
    }
    
    n.eachStrataT <- unlist(lapply(index.strataT, length)) # those variables are used during the boostrap and passed to the wrapper through the environment
    n.eachStrataC <- unlist(lapply(index.strataC, length))
    nCumSum.strataControl <- cumsum(c(1,n.eachStrataC))
    nCumSum.strataTreatment <- cumsum(c(1,n.eachStrataT))
       
    delta_boot <- array(NA, dim = c(n.strata, D, n.bootstrap))
    
    #### computation
    calcBootstrap(environment())
    
    #### post treatment
    if (trace > 1) {cat("Post-Treatment \n")}
    res <- calcCI(delta = resPonctual$delta, delta_boot = delta_boot, endpoint = endpoint, D = D, alternative = alternative, alpha = 0.05,
                  n.bootstrap = n.bootstrap, cpus = cpus, trace = TRUE)
    p.value <- res$p.value  
    Delta_quantile <-  res$Delta_quantile
    
    #### update BuyseRes object 
    BuyseRes.object@delta_boot <- delta_boot
    BuyseRes.object@p.value <- p.value
    BuyseRes.object@Delta_quantile <- Delta_quantile
  }
  
  #### export ####
  return(BuyseRes.object)
}
