#' @docType package
#' @name BuyseTest
#' @useDynLib BuyseTest
#' @import data.table
#' @importFrom lava categorical coxExponential.lvm distribution eventTime lvm sim
#' @import methods
#' @importFrom parallel detectCores
#' @import Rcpp
#' @import RcppArmadillo
#' @import snowfall
#' @importFrom stats as.formula na.omit rbinom
#' @importFrom stats4 summary
#' @importFrom survival survfit Surv
#' @importFrom tcltk tkProgressBar setTkProgressBar
#' @importFrom utils tail
NULL


#### #' @importFrom snowfall sfClusterEval sfExport sfInit sfLibrary sfStop
