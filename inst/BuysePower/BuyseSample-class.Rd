\docType{class}
\name{BuyseSample-class}
\title{Class "BuyseSample" (output of BuysePower}
\alias{BuyseSample}
\alias{BuyseSample-class}
\description{
  A \code{\link{BuysePower}} results is represented as a \code{BuysePower} object.
}
\value{
An object with the following attributes :
\itemize{
\item\code{@alpha} :  the risk of first species. \emph{numeric}.
\item\code{@alternative} :  the type of alternative hypothesis. \emph{character}.
\item\code{@KMimputation} :  Whether the Kaplan Meier estimator was used to compare uninformative pairs. \emph{logical}.
\item\code{@n.bootstrap} :  the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}.
\item\code{@type} : the type of each endpoint. \emph{vector}.
\item\code{@threshold} :  the thresholds, one for each endpoint variable. \emph{vector}.
\item\code{@power} : the power of the test for each sample size and each outcome (cumulative) for a risk of first species of \code{@alpha}. \emph{array}.
\item\code{@delta} :  the proportion in favor of the treatment for each simulation, each sample size and each outcome. \emph{array}.
\item\code{@censoring.rate} : the censoring rate  for each simulation, each sample size and each time to event outcome for the treatment and control arm. \emph{array}.
\item\code{@p.value} :  the p.value obtained for each simulation, each sample size and each outcome (cumulative). \emph{array}.
\item\code{@DeltaCI} : the cumulative proportion in favor of treatment and its confidence interval for a risk of first species of \code{@alpha} obtained for each simulation and each sample size (cumulative). \emph{array}.
}
}
\seealso{
  \code{\link{BuysePower}} for performing a generalized pairwise comparison. \cr
  \code{\link{summary,BuyseSample-method}} for a summary of the results of generalized pairwise comparison.
}
\examples{
hypothesis <- list(c(lambda_t=0.25,lambda_c=0.5,T_inclusion=1,T_followUp=1))

\dontrun{
BuysePower_object <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type="TTE",endpoint="End1",threshold=0)
}
\dontshow{
BuysePower_object <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
                  n.simul=10,n.bootstrap=100,trace=0,
                  type="TTE",endpoint="End1",threshold=0)
}
class(BuysePower_object)

}
\keyword{classes}

