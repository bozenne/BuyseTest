#' @docType package
#' @title BuyseTest package: Generalized Pairwise Comparisons
#' @name BuyseTest-package
#' 
#' @description Implementation of the Generalized Pairwise Comparisons. 
#' \code{\link{BuyseTest}} is the main function of the package. See the vignette of an overview of the functionalities of the package.
#' Run \code{citation("BuyseTest")} in R for how to cite this package in scientific publications.
#' See the section reference below for examples of application in clinical studies.
#'
#' @details
#' The Generalized Pairwise Comparisons form all possible pairs of observations,
#' one observation being taken from the intervention group and the other is taken from the control group,
#' and compare the difference in endpoints (\eqn{Y-X}) to the threshold of clinical relevance (\eqn{\tau}).
#'
#' For a single endpoint,
#' if the difference is greater or equal than the threshold of clinical relevance (\eqn{Y \ge X + \tau}),
#' the pair is classified as favorable (i.e. win).
#' If the difference is lower or equal than minus the threshold of clinical relevance (\eqn{X \ge Y + \tau}),
#' the pair is classified as unfavorable (i.e. loss).
#' Otherwise the pair is classified as neutral. In presence of censoring, it might not be possible to compare the difference to the threshold. In such cases the pair
#' is classified as uninformative.
#'
#' Simultaneously analysis of several endpoints is performed by prioritizing the endpoints, assigning the highest priority to the endpoint considered the most clinically relevant.
#' The endpoint with highest priority is analyzed first, and neutral and uninformative pair are analyzed regarding endpoint of lower priority.
#'
#' @references Examples of application in clinical studies: \cr 
#' J. Peron, P. Roy, K. Ding, W. R. Parulekar, L. Roche, M. Buyse (2015). \bold{Assessing the benefit-risk of new treatments using generalized pairwise comparisons: the case of erlotinib in pancreatic cancer}. \emph{British journal of cancer} 112:(6)971-976.  \cr
#' J. Peron, P. Roy, T. Conroy, F. Desseigne, M. Ychou, S. Gourgou-Bourgade, T. Stanbury, L. Roche, B. Ozenne, M. Buyse (2016). \bold{An assessment of the benefit-risk balance of FOLFORINOX in metastatic pancreatic adenocarcinoma}. \emph{Oncotarget} 7:82953-60, 2016. \cr
#'
#' Comparison between the net benefit and alternative measures of treatment effect: \cr 
#' J. Peron, P. Roy, B. Ozenne, L. Roche, M. Buyse (2016). \bold{The net chance of a longer survival as a patient-oriented measure of benefit in randomized clinical trials}. \emph{JAMA Oncology} 2:901-5. \cr
#' E. D. Saad , J. R. Zalcberg, J. Peron, E. Coart, T. Burzykowski, M. Buyse (2018). \bold{Understanding and communicating measures of treatment effect on survival: can we do better?}. \emph{J Natl Cancer Inst}.
#' @useDynLib BuyseTest, .registration=TRUE
#' @import data.table
#' @importFrom ggplot2 autoplot
#' @importFrom lava categorical coxExponential.lvm distribution eventTime iid lvm sim vars latent<-
#' @import methods
#' @importFrom parallel detectCores
#' @import Rcpp
#' @importFrom stats as.formula delete.response formula na.omit rbinom predict setNames terms
#' @importFrom stats4 summary
#' @importFrom prodlim prodlim Hist
#' @importFrom utils capture.output tail
NULL


