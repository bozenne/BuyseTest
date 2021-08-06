## * load packages
library(testthat)
library(BuyseTest)

library(lava)
library(data.table)
library(survival)
xx <- capture.output(print(veteran))
library(riskRegression)
library(prodlim)
library(mvtnorm)

## * run tests
## setwd("~/Documents/GitHub/BuyseTest/tests/")
## source("~/Documents/GitHub/BuyseTest/tests/testthat/test-BuyseTest-engine.R")
test_check("BuyseTest")

## * memory check
## /data/gannet/ripley/R/packages/tests-clang-SAN/BuyseTest/src/FCT_buyseTest.cpp:894:20

## ** valgrind
## Information: https://kevinushey.github.io/blog/2015/04/05/debugging-with-valgrind/

## R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes" --vanilla < testthat.R
## R -d "valgrind --dsymutil=yes" -e "Rcpp::sourceCpp('segfault.cpp')"
## R -d valgrind

## ** san
## docker run --rm -ti -v $(pwd):/mydir wch1/r-debug
## wget -O - https://github.com/bozenne/BuyseTest/tarball/master | tar xz
## cd bozenne
## cd test/testthat
## RDcsan
## install.packages(c("Rcpp","RcppArmadillo","devtools","riskRegrssion","cvAUC"))
## devtools::install_github("bozenne/BuyseTest")
## allFiles <- list.files()
## lapply(rev(allFiles), function(x){cat(x,"\n"); source(x)})
## gctorture(TRUE)
##



