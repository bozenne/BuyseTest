\name{BuysePower}
\title{Compute the power associated to a simulation scheme}
\alias{BuysePower}
\description{
Compute the power of the test for a given simulation scheme and for various sample sizes.
}
\usage{
BuysePower(n,hypothesis,type, threshold=NULL,endpoint=NULL,alpha=0.05, method="Gehan",
         n.bootstrap=1000,prob.alloc=0.5,n.simul=1000,
         alternative = "two.sided",seed=10,cpus=1,trace=3)
}
\arguments{
  \item{n}{the sample size to consider. \emph{numeric vector}[3]. See the details section.}
  \item{hypothesis}{a list containing the simulation scheme. \emph{list}. See the details section.}
  \item{type}{the type of each endpoint. \emph{character vector}. Can be \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.}
  \item{threshold}{the thresholds, one for each endpoint variable. \emph{numeric vector}. Default is \code{NULL} indicating no threshold.}
  \item{endpoint}{the name of the endpoint variable(s). \emph{character vector}. Default is \code{NULL} leading to attribute standard endpoint names.}
  \item{alpha}{the risk of first species. \emph{numeric}. Default is \code{0.05}.}
  
  \item{method}{paires with censored data can be either classified as uniformative (\code{"Gehan"}) or compared regarding the predicted survival using a common survival curve for treated and control patients (\code{"Peto"}) or a separate survival curve for treated and control patients (\code{"Efron"}).}  
  \item{n.bootstrap}{the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}. Default is \code{1000}.}
  \item{prob.alloc}{the resampling probability for assignement to the experimental group in the bootstrap samples. \emph{double}. Default is \code{NULL} indicating to use the proportion of patients in the experimental group.}
  \item{n.simul}{the number of simulation used to compute the power. \emph{integer}. Default is \code{1000}.}
 
  \item{alternative}{a \emph{character} specifying the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. Default is \code{"two.sided"}.}
  \item{seed}{the seed to consider when generating the data. \emph{integer}. Default is \code{10}.}
  \item{cpus}{the number of CPU to use. \emph{integer}. Default is \code{1}.}
  \item{trace}{should the execution of the function be traced ? \emph{integer}. Default is \code{3}.}
}
\details{
\code{n} specify the arguments of the \code{\link{seq}} function. The the first element corresponds to the minimum sample size, the second to the maximum sample size and the third to the step between sample sizes (e.g. n=c(100,1000,100)). \cr
Argument can be named (e.g. n = c(from=100,to=1000,by=100)). Equivalent parametrization recognized by \code{seq} can be used (e.g. n = c(from=100,to=1000,length.out=10) or n = c(from=100,ny=100,length.out=10)). \cr

\code{hypothesis} must contains the values of the parameter for each endpoint. It must have the same length as argument \code{type}. \cr
Elements in \code{hypothesis} corresponding to binary endpoints should contains two values corresponding the probability of event in the treatment and control arm (respectively called \code{proba_t} and \code{proba_c}). \cr
Elements in \code{hypothesis} corresponding to continuous endpoints should contains three values corresponding the mean value of the outcome in the treatment and control arm, and its standard deviation (respectively called \code{mu_t}, \code{mu_c} and  \code{sigma}). A normal law is used for the simulation. \cr
Elements in \code{hypothesis} corresponding to time to event endpoints should contains fours values corresponding the hazard in the treatment and control arm and to the inclusion and follow-up time (respectively called \code{lambda_t}, \code{lambda_c}, \code{T_inclusion} and \code{T_followUp}). 
For each patient a inclusion time (here denoted t_i) is drawn from a uniform law between 0 and \code{T_inclusion}. The time to event is simulated using an exponential law and the patient is censored if this time is higher than t_i+\code{T_followUp}.

The expected percent of censoring is computed using the following formula : \cr
Prob[censoring] = exp(-lamda*T_followUp)*[1-exp(-lambda*T_inclusion)]/(lambda*T_inclusion) 

}
\value{
  An \R object of class \code{\linkS4class{BuyseSample}}.
}
\seealso{
   \code{\link{summary,BuyseSample-method}} for a summary of the results of power computation. \cr
   \code{\link{BuyseSample-class}} for a presentation of the \code{BuyseSample} object. \cr
}
\examples{

#### 1- Binary outcome ####
hypothesis <- list(c(proba_t=0.5,proba_c=0.75))

\dontrun{
resBin <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type="bin",endpoint="NoToxicity",threshold=NA)
}
\dontshow{
resBin <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
                     n.simul=10,n.bootstrap=100,
                     type="bin",endpoint="NoToxicity",threshold=NA)
}

summary(resBin)

#### 2- Continuous Outcome ####
hypothesis <- list(c(mu_t=2,mu_c=1,sigma=1))

\dontrun{
resCont <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type="cont",endpoint="Score",threshold=1)
}
\dontshow{
resCont <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
					  n.simul=10,n.bootstrap=50,trace=0,
                      type="cont",endpoint="Score",threshold=1)
}
summary(resCont,Delta=TRUE)

#### 3- Survival Outcome ####
hypothesis <- list(c(lambda_t=0.2,lambda_c=0.5,T_inclusion=1,T_followUp=1))

# no imputation
\dontrun{
resTTE <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type="TTE",endpoint="Survival",threshold=0.5)
}
\dontshow{
resTTE <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
                     n.simul=10,n.bootstrap=50,trace=0,
                     type="TTE",endpoint="Survival",threshold=0.5)
}
summary(resTTE,Delta=TRUE,censoring.rate=TRUE)

# imputation
\dontrun{
resTTEKM <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,method="Peto",
                  type="TTE",endpoint="Survival",threshold=0.5)
}
\dontshow{
resTTEKM <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,method="Peto",
                       n.simul=10,n.bootstrap=50,trace=0,
                       type="TTE",endpoint="Survival",threshold=0.5)
}

summary(resTTEKM,Delta=TRUE,censoring.rate=TRUE)

#### 4- Mixte Outcome ####
hypothesis <- list(c(lambda_t=0.25,lambda_c=0.5,T_inclusion=1,T_followUp=1),
				   c(proba_t=0.5,proba_c=0.75),
				   c(lambda_t=0.25,lambda_c=0.5,T_inclusion=1,T_followUp=1)
				   )

# no imputation
\dontrun{
resMixte <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,
                  type=c("TTE","bin","TTE"),endpoint=c("Survival","NoToxicity","Survival"),
				  threshold=c(0.5,NA,0.25))
}
\dontshow{
resMixte <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,
                  type=c("TTE","bin","TTE"),endpoint=c("Survival","NoToxicity","Survival"),
				  threshold=c(0.5,NA,0.25),
				  n.simul=10,n.bootstrap=50,trace=0)
}
summary(resMixte,Delta=TRUE,censoring.rate=TRUE)

# imputation
\dontrun{
resMixteKM <- BuysePower(n=c(10,110,20),hypothesis=hypothesis,method="Peto",
                  type=c("TTE","bin","TTE"),endpoint=c("Survival","NoToxicity","Survival"),
				  threshold=c(0.5,NA,0.25))
}
\dontshow{
resMixteKM <- BuysePower(n=c(10,20,5),hypothesis=hypothesis,method="Peto",
                  type=c("TTE","bin","TTE"),endpoint=c("Survival","NoToxicity","Survival"),
				  threshold=c(0.5,NA,0.25),
				  n.simul=10,n.bootstrap=50,trace=0)
}
summary(resMixteKM,Delta=TRUE,censoring.rate=TRUE)

}
\keyword{function}
