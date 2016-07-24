\docType{class}
\name{BuyseRes-class}
\title{Class "BuyseRes" (output of BuyseTest}
\alias{BuyseRes}
\alias{BuyseRes-class}
\description{
  A \code{\link{BuyseTest}} output is reported in a \code{BuyseRes} object.
}
\value{
An object with the following attributes :
\itemize{
\item\code{@delta} :  the chance of a better outcome (net difference between the probability for a random pair to be favorable minus the probability to be unfavorable divided by the total number of pairs) for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
\item\code{@count_favorable} :  the probability for a random pair to be favorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
\item\code{@count_unfavorable.group} : the probability for a random pair to be unfavorable for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
\item\code{@count_neutral} :  the probability for a random pair to be neutral for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
\item\code{@count_uninf} : the probability for a random pair to be uninformative for each strata (in rows) and each endpoint (in columns). \emph{matrix}.
\item\code{@index_neutralT} :  the index in the dataset of the treatment observations from remaining neutral pairs. \emph{integer vector}.
\item\code{@index_neutralC} : the index in the dataset of the control observations from remaining neutral pairs. \emph{integer vector}.
\item\code{@index_uninfT} :  the index in the dataset of the treatment observations from remaining uninformative pairs. \emph{integer vector}.
\item\code{@index_uninfC} : the index in the dataset of the control observations from remaining uninformative pairs. \emph{integer vector}.
\item\code{@n_pairs} :  the total number of pairs. \emph{integer}.
\item\code{@delta_boot} :  the chance of a better outcome within each strata (first dimension), over the endpoint (second dimension) and for each bootstrap dataset (third dimension). \emph{array}.
\item\code{@p.value} : the p.value associated to the chance of a better outcome at each prioritized endpoint. \emph{numeric vector}.
\item\code{@Delta_quantile} :  the randomization test-based 2.5\% and the 97.5\% quantiles of the cumulative chance of a better outcome (in rows) for each endpoint (in columns). \emph{matrix}.
\item\code{@endpoints} : the name of the endpoints. \emph{character vector}.
}
}
\seealso{
  \code{\link{BuyseTest}} for the function computing generalized pairwise comparisons. \cr
  \code{\link{summary,BuyseRes-method}} for the summary of the BuyseTest function results
}
\examples{
n.Treatment_testBin <- 500
n.Control_testBin <- 500
prob.Treatment_testBin <- c(0.5,0.75)
prob.Control_testBin <- c(0.5,0.25)

set.seed(10)
data_testBin <- data.frame(treatment=c(rep(1,n.Treatment_testBin),rep(0,n.Treatment_testBin)))
data_testBin$endpoint1 <- c(rbinom(n.Treatment_testBin,size=1,prob=prob.Treatment_testBin[1]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[1]))
data_testBin$endpoint2 <- c(rbinom(n.Control_testBin,size=1,prob=prob.Treatment_testBin[2]),
                            rbinom(n.Control_testBin,size=1,prob=prob.Control_testBin[2]))
data_testBin$strata <- rbinom(n.Treatment_testBin+n.Control_testBin,size=4,prob=0.5)

#### no strata, n.bootsrap=0
BuyseTest_object <- BuyseTest(data=data_testBin,endpoint=c("endpoint1","endpoint2"),
                              treatment="treatment", type=c("bin","bin"))
								 
class(BuyseTest_object)

}
\keyword{classes}

