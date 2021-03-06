% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normexp.R
\name{qnormweibull}
\alias{qnormweibull}
\title{Density of a Gaussian Variable Plus an Weibull Variable}
\usage{
qnormweibull(p, scale, shape, rho)
}
\description{
Density of a Gaussian Variable Plus an Weibull Variable
}
\examples{
\dontrun{
n <- 5e6

c(qnormweibull(0.5, scale = 1/2, shape = 1, rho = 1.5),
 quantile(rnorm(n) + 1.5 * rweibull(n, scale = 1/2, shape = 1), 0.5))
c(qnormweibull(0.95, scale = 10, shape = 1, rho = 1.5),
quantile(rnorm(n) + 1.5 * rweibull(n, scale = 10, shape = 1), 0.95))

c(qnormweibull(0.5, scale = 1/2, shape = 2, rho = 1.5),
quantile(rnorm(n) + 1.5 * rweibull(n, scale = 1/2, shape = 2), 0.5))
c(qnormweibull(0.95, scale = 10, shape = 2, rho = 1.5),
quantile(rnorm(n) + 1.5 * rweibull(n, scale = 10, shape = 2), 0.95))

c(qnormweibull(0.5, scale = 1/2, shape = 1, rho = -1.5),
quantile(rnorm(n) - 1.5 * rweibull(n, scale = 1/2, shape = 1), 0.5))
c(qnormweibull(0.95, scale = 10, shape = 1, rho = -1.5),
quantile(rnorm(n) - 1.5 * rweibull(n, scale = 10, shape = 1), 0.95))

c(qnormweibull(0.5, scale = 1/2, shape = 2, rho = -1.5),
quantile(rnorm(n) - 1.5 * rweibull(n, scale = 1/2, shape = 2), 0.5))
c(qnormweibull(0.95, scale = 10, shape = 2, rho = -1.5),
quantile(rnorm(n) - 1.5 * rweibull(n, scale = 10, shape = 2), 0.95))

}
}
\keyword{internal}
