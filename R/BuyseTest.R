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
#' @param neutral.as.uninf [logical] should paired classified as neutral be re-analysed using endpoints of lower priority.
#' Default value read from \code{BuyseTest.options()}.
#' @param n.permutation [integer] the number of permutations used for computing the confidence interval and the p.values. See details.
#' Default value read from \code{BuyseTest.options()}.
#' @param prob.alloc [0<double<1] the resampling probability for being allocated to the experimental group in the permutation test.
#' Can also be \code{NULL} to use the proportion of patients in the experimental group.
#' @param stratified [logical] should the allocation in the permutation test be performed within strata.
#' This means that the \code{prob.alloc} will be satisfied within strata not only globally.
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
                      n.permutation = NULL,
                      prob.alloc = NULL,
                      stratified = FALSE,
                      keep.comparison = NULL,
                      alternative = "two.sided", 
                      seed = NULL,
                      cpus = NULL,
                      trace = NULL){

    BuyseCall <- match.call()
    option <- BuyseTest.options()

    ## ** initialize arguments
    outArgs <- initializeArgs(alternative = alternative,
                              BuyseCall = BuyseCall,
                              censoring = censoring,
                              cpus = cpus,
                              endpoint = endpoint,
                              formula = formula,
                              keep.comparison = keep.comparison,
                              method = method,
                              n.permutation = n.permutation,
                              neutral.as.uninf = neutral.as.uninf,
                              operator = operator,
                              option = option,
                              seed = seed,
                              threshold = threshold,
                              trace = trace,
                              treatment = treatment,
                              type = type)

    browser()
    alternative  <- outArgs$alternative
    censoring <- outArgs$censoring
    cpus <- outArgs$cpus
    endpoint <- outArgs$endpoint
    keep.comparison <- outArgs$keep.comparison
    method <- outArgs$method
    n.permutation <- outArgs$n.permutation
    neutral.as.uninf <- outArgs$neutral.as.uninf
    operator <- outArgs$operator
    seed <- outArgs$seed
    strata <- outArgs$strata
    threshold <- outArgs$threshold
    trace <- outArgs$trace
    treatment <- outArgs$treatment
    type <- outArgs$type

    ## ** test args
    if(TRUE){
        testArgs(alternative = alternative,
                 BuyseCall = BuyseCall,
                 censoring = censoring,
                 correctionTTE = correctionTTE,
                 cpus = cpus,
                 data = data,
                 endpoint = endpoint,
                 formula = formula,
                 method = method,
                 n.permutation = n.permutation,
                 operator = operator,
                 proba = proba,
                 seed = seed,
                 stratified = stratified,
                 strata = strata,
                 threshold = threshold,
                 treatment = treatment,
                 type = type)

        ## ## ** proba
        ## if(!is.null(prob.alloc)){
        ##     validNumeric(prob.alloc,
        ##                  valid.length = 1,
        ##                  min = 0,
        ##                  max = 1,
        ##                  method = "BuyseTest")
        ## }
    }

    browser()
    ## ** General display
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
    
    ## ** Lauch BuyseTest
    .BuyseTest(data = data, 
               treatment = treatment,
               endpoint = endpoint,
               type = type,
               threshold = threshold,
               censoring = censoring,
               strata = strata, 
               method = method,
               correctionTTE = correctionTTE,
               neutral.as.uninf = neutral.as.uninf,
               n.permutation = n.permutation,
               prob.alloc = prob.alloc,
               stratified = stratified,
               keep.comparison = keep.comparison,
               alternative = alternative, 
               seed = seed,
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
    
    ### ** 4- Transformation into BuyseRes object
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

    
}

## * .BuyseTest
.BuyseTest <- function(data, 
                       treatment = NULL,
                       endpoint = NULL,
                       type = NULL,
                       threshold = NULL,
                       censoring = NULL,
                       strata = NULL, 
                       method = NULL,
                       correctionTTE = FALSE,
                       neutral.as.uninf = NULL,
                       n.permutation = NULL,
                       prob.alloc = NULL,
                       stratified = FALSE,
                       keep.comparison = NULL,
                       alternative = "two.sided", 
                       seed = NULL,
                       cpus = NULL,
                       trace = NULL){

   
    ## ** initialize dataset
    initializeData(BT.envir) ## update of the arguments via environment

    if (is.null(prob.alloc)) { # if prob.alloc is not set by the user
        prob.alloc <- n.Treatment/(n.Treatment + n.Control)  # it is set to the proportion of patients in the treatment arm.
    }
    
    ## ** KM imputation
    if(method %in% 1:3){# c("Peto","Efron","Peron")
    
    endpoint.TTE <- endpoint[type==3] # vector of variable names of the TTE endpoints
    
    ## *** design matrix for the weights
    if(is.null(Wscheme)){
        res_init <- initWscheme(D=D,
                                endpoint=endpoint,
                                endpoint.TTE=endpoint.TTE,
                                D.TTE=D.TTE,
                                threshold=threshold,
                                type=type)
      Wscheme <- res_init$Wscheme  
      index_survivalM1 <- res_init$index_survivalM1
      threshold_TTEM1 <- res_init$threshold_TTEM1       
    }
    ## *** Survival estimate using Kaplan Meier    
    res_init <- initSurvival(M.Treatment=M.Treatment,
                             M.Control=M.Control,
                             M.delta_Treatment=M.delta_Treatment,
                             M.delta_Control=M.delta_Control,
                             endpoint=endpoint,
                             D.TTE=D.TTE,
                             type=type,
                             threshold=threshold,
                             index.strataT=index.strataT,
                             index.strataC=index.strataC,
                             n.strata=n.strata,   
                             method=method)
    
    list_survivalT <- res_init$list_survivalT
    list_survivalC <- res_init$list_survivalC
    
    }else{
        Wscheme <- matrix()  # factice design matrix for the weights. Will be sent to the C++ arguments to fill the argument but not used by the function.
        list_survivalT <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
        list_survivalC <- list() # factice list. Will be sent to the C++ arguments to fill the argument but not used by the function.
        index_survivalM1 <- numeric(0) # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
        threshold_TTEM1 <- numeric(0)  # factice vector. Will be sent to the C++ arguments to fill the argument but not used by the function.
    }


    
    ## ** Punctual estimation
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
                               correctionTTE = correctionTTE,
                               neutralAsUninf = neutral.as.uninf,
                               keepComparison = keep.comparison
                               )
    })
    
    if (trace > 1) {cat("   > done \n")}
    

    ## ** Permutation test 
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
        
       
    }
    
    ## ** Export
    return(BuyseRes.object)
}
