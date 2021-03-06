% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discreteRoot.R
\name{boot2pvalue}
\alias{boot2pvalue}
\title{Compute the p.value from the distribution under H1}
\usage{
boot2pvalue(
  x,
  null,
  estimate = NULL,
  alternative = "two.sided",
  FUN.ci = quantileCI,
  checkSign = TRUE,
  tol = .Machine$double.eps^0.5
)
}
\arguments{
\item{x}{[numeric vector] a vector of bootstrap estimates of the statistic.}

\item{null}{[numeric] value of the statistic under the null hypothesis.}

\item{estimate}{[numeric] the estimated statistic.}

\item{alternative}{[character] a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}

\item{FUN.ci}{[function] the function used to compute the confidence interval.
Must take \code{x}, \code{alternative}, \code{conf.level} and \code{sign.estimate} as arguments
and only return the relevant limit (either upper or lower) of the confidence interval.}

\item{checkSign}{[logical] should a warning be output if the sign of the estimate differs
from the sign of the mean bootstrap value?}

\item{tol}{[numeric] the absolute convergence tolerance.}
}
\description{
Compute the p.value associated with the estimated statistic
using a bootstrap sample of its distribution under H1.
}
\details{
For test statistic close to 0, this function returns 1. \cr \cr

For positive test statistic, this function search the quantile alpha such that:
\itemize{
\item \code{quantile(x, probs = alpha)=0} when the argument alternative is set to \code{"greater"}.
\item \code{quantile(x, probs = 0.5*alpha)=0} when the argument alternative is set to \code{"two.sided"}.
}
If the argument alternative is set to \code{"less"}, it returns 1. \cr \cr

For negative test statistic, this function search the quantile alpha such that:
\itemize{
\item \code{quantile(x, probs = 1-alpha=0} when the argument alternative is set to \code{"less"}.
\item \code{quantile(x, probs = 1-0.5*alpha=0} when the argument alternative is set to \code{"two.sided"}.
}
If the argument alternative is set to \code{"greater"}, it returns 1.
}
\examples{
set.seed(10)

#### no effect ####
x <- rnorm(1e3) 
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "two.sided")
## expected value of 1
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "greater")
## expected value of 0.5
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "less")
## expected value of 0.5

#### positive effect ####
x <- rnorm(1e3, mean = 1) 
boot2pvalue(x, null = 0, estimate = 1, alternative = "two.sided")
## expected value of 0.32 = 2*pnorm(q = 0, mean = -1) = 2*mean(x<=0)
boot2pvalue(x, null = 0, estimate = 1, alternative = "greater")  
## expected value of 0.16 = pnorm(q = 0, mean = 1) = mean(x<=0)
boot2pvalue(x, null = 0, estimate = 1, alternative = "less")
## expected value of 0.84 = 1-pnorm(q = 0, mean = 1) = mean(x>=0)

#### negative effect ####
x <- rnorm(1e3, mean = -1) 
boot2pvalue(x, null = 0, estimate = -1, alternative = "two.sided") 
## expected value of 0.32 = 2*(1-pnorm(q = 0, mean = -1)) = 2*mean(x>=0)
boot2pvalue(x, null = 0, estimate = -1, alternative = "greater")
## expected value of 0.84 = pnorm(q = 0, mean = -1) = mean(x<=0)
boot2pvalue(x, null = 0, estimate = -1, alternative = "less") # pnorm(q = 0, mean = -1)
## expected value of 0.16 = 1-pnorm(q = 0, mean = -1) = mean(x>=0)
}
