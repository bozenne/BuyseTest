\name{summarySample-methods}
\title{Summary Method for Class "BuyseSample"}
\alias{summary,BuyseSample-method}
\description{
Summarize the results from the \code{\link{BuysePower}} function.
}
\usage{
\S4method{summary}{BuyseSample}(object,Delta=FALSE,delta=FALSE,censoring.rate=FALSE,digit=3,trace=TRUE)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{BuyseSample}}, i.e.,
  as resulting from \code{\link{BuysePower}}}
  \item{Delta}{should the median and confidence interval  of the cumulative proportion in favor of treatment be printed ? \emph{logical} . Default is \code{FALSE}.}
  \item{delta}{should the median of the proportion in favor of treatment for each outcome be printed ? \emph{logical} . Default is \code{FALSE}.}
  \item{censoring.rate}{should the censoring rate for each time to event outcome be printed ? \emph{logical} . Default is \code{FALSE}.}
  \item{digit}{the number of digit to use for printing the results. \emph{integer}. Default is \code{3}.}
  \item{trace}{should the main results be printed ? \emph{logical}. Default is \code{TRUE}.}
}
\details{
WARNING : the confidence interval is computed using quantiles of the distribution of the cumulative proportion in favor of the treatment under the null hypothesis. It thus may not be valid if this hypothesis is rejected. In particular, if the cumulative proportion in favor of the treatment is close to 1, the upper limit of the confidence interval may exceed 1. 
}
\value{
  A \code{"data.frame"} with in columns the sample size and in line :
  \itemize{
\item\code{"power"} : the power of the test. 
\item\code{"DeltaCI"} : one line for the median cumulative proportion in favor of treatment estimated over simulation, one line for the lower value of the confidence interval and one line for the higher value of the confidence interval.
\item\code{"n.bootstrap"} : the number of bootstrap sample effectively used to compute the confidence interval.
\item\code{"censoring.rate"} : the median percent of censored individuals for each time to event outcome estimated on the samples for each sample size.
}
}
\seealso{
 \code{\link{BuysePower}} for performing a generalized pairwise comparison. \cr
 \code{\link{BuyseSample-class}} for a presentation of the \code{BuyseSample} object.
}
\examples{
hypothesis <- list(c(lambda_t=0.25,lambda_c=0.5,T_inclusion=1,T_followUp=1))

\dontrun{
BuysePower_object <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type="TTE",endpoint="End1",threshold=0)
}
\dontshow{
BuysePower_object <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
                  n.simul=100,n.bootstrap=10,trace=0,
                  type="TTE",endpoint="End1",threshold=0)
}
summary(BuysePower_object,Delta=TRUE,censoring.rate=TRUE)
}
\keyword{methods}
