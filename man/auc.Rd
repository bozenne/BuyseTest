% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auc.R
\name{auc}
\alias{auc}
\title{Estimation of the Area Under the ROC Curve (EXPERIMENTAL)}
\usage{
auc(
  labels,
  predictions,
  fold = NULL,
  observation = NULL,
  direction = ">",
  add.halfNeutral = TRUE,
  null = 0.5,
  conf.level = 0.95,
  transformation = TRUE,
  order.Hprojection = 2,
  pooling = "mean"
)
}
\arguments{
\item{labels}{[integer/character vector] the disease status (should only take two different values).}

\item{predictions}{[numeric vector] A vector with the same length as \code{labels} containing the biomarker values.}

\item{fold}{[character/integer vector] If using cross validation, the index of the fold. 
Should have the same length as \code{labels}.}

\item{observation}{[integer vector] If using cross validation, the index of the corresponding observation in the original dataset.
Necessary to compute the standard error when using cross validation.}

\item{direction}{[character] \code{">"} lead to estimate P[Y>X],
\code{"<"} to estimate P[Y<X],
and \code{"auto"} to estimate max(P[Y>X],P[Y<X]).}

\item{add.halfNeutral}{[logical] should half of the neutral score be added to the favorable and unfavorable scores?
Useful to match the usual definition of the AUC in presence of ties.}

\item{null}{[numeric, 0-1] the value against which the AUC should be compared when computing the p-value.}

\item{conf.level}{[numeric, 0-1] the confidence level of the confidence intervals.}

\item{transformation}{[logical] should a log-log transformation be used when computing the confidence intervals and the p-value.}

\item{order.Hprojection}{[1,2] the order of the H-projection used to linear the statistic when computing the standard error.
2 is involves more calculations but is more accurate in small samples. Only active when the \code{fold} argument is \code{NULL}.}

\item{pooling}{[character] method used to compute the global AUC from the fold-specific AUC: either an empirical average \code{"mean"}
or a weighted average with weights proportional to the number of pairs of observations in each fold \code{"pairs"}.}
}
\value{
An S3 object of class \code{BuyseTestAUC} that inherits from data.frame.
The last line of the object contains the global AUC value with its standard error.
}
\description{
Estimation of the Area Under the ROC curve, possibly after cross validation,
to assess the discriminant ability of a biomarker regarding a disease status.
}
\details{
The iid decomposition of the AUC is based on a first order decomposition.
So its squared value will not exactly match the square of the standard error estimated with a second order H-projection.
}
\examples{
library(data.table)

n <- 200
set.seed(10)
X <- rnorm(n)
dt <- data.table(Y = as.factor(rbinom(n, size = 1, prob = 1/(1+exp(1/2-X)))),
                 X = X,
                 fold = unlist(lapply(1:10,function(iL){rep(iL,n/10)})))

## compute auc
auc(labels = dt$Y, predictions = dt$X, direction = ">")

## compute auc after 10-fold cross-validation
auc(labels = dt$Y, prediction = dt$X, fold = dt$fold, observation = 1:NROW(dt))

}
\references{
Erin LeDell, Maya Petersen, and Mark van der Laan (2015). \bold{Computationally efficient confidence intervals for cross-validated area under the ROC curve estimates}. \emph{Electron J Stat.} 9(1):1583–1607. \cr
}
\keyword{models}
