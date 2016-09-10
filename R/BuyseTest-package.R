#' @docType package
#' @title BuyseTest package: Generalized Pairwise Comparisons
#' @name BuyseTest-package
#' 
#' @description Implementation of the Generalized Pairwise Comparisons.
#' \code{\link{BuyseTest}} is the main function of the package. See its documentation for more details or the reference below for a complete description of the method and some examples of application.
#' 
#' @useDynLib BuyseTest
#' @import data.table
#' @importFrom lava categorical coxExponential.lvm distribution eventTime lvm sim
#' @import methods
#' @importFrom parallel detectCores
#' @import Rcpp
#' @import snowfall
#' @importFrom stats as.formula na.omit rbinom setNames
#' @importFrom stats4 summary
#' @importFrom survival survfit Surv
#' @importFrom tcltk tkProgressBar setTkProgressBar
#' @importFrom utils tail
#' @references 
#' Methodological references: \cr
#' Marc Buyse (2010) \bold{Generalized pairwise comparisons of prioritized endpoints in the two-sample problem}. \emph{Statistics in Medicine} 29:3245-3257. \cr
#' J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2016). \bold{An extension of generalized pairwise comparisons for prioritized outcomes in the presence of censoring}. Statistical Methods in Medical Research. \cr
#' 
#' Examples of application: \cr 
#' J. Peron, P Roy, K Ding, W R Parulekar, L Roche, M Buyse (2015). \bold{Assessing the benefit-risk of new treatments using generalised pairwise comparisons: the case of erlotinib in pancreatic cancer}. \emph{British journal of cancer} 112:(6)971-976.  \cr
#' J. Peron, M. Buyse, B. Ozenne, L. Roche and P. Roy (2016). \bold{The Net Chance of a Longer Survival as a Patient-Oriented Measure of Treatment Benefit in Randomized Clinical Trials}. \emph{JAMA Oncology} 2(7):901-5. \cr
NULL


