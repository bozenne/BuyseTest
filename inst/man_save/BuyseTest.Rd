\name{BuyseTest}
\title{Generalized Pairwise Comparisons}
\alias{BuyseTest}
\description{
Performs Generalized Pairwise Comparisons for binary, continuous and time-to-event outcomes.
}
\usage{
BuyseTest(data,treatment,endpoint,threshold=NULL,strata=NULL,censoring=NULL,type, 
         method = "Peron", n.bootstrap=0,prob.alloc=NULL, 
         alternative = "two.sided",seed=10,cpus=1,trace=3)
}
\arguments{
  \item{data}{A \code{data.frame} containing the variables.}
  \item{treatment}{the name of the treatment variable identifying the control and the experimental group. \emph{character}.}
  \item{endpoint}{the name of the endpoint variable(s). \emph{character vector}.}
  \item{threshold}{the thresholds, one for each endpoint variable. \emph{numeric vector}. Default is \code{NULL} indicating no threshold.}
  \item{strata}{the name of the strata variable(s). \emph{numeric vector}. Default is \code{NULL} indicating only one strata.}
  \item{censoring}{the name of the censoring variable(s), one for each endpoint. \emph{character vector}. Default is \code{NULL}.}
  \item{type}{the type of each endpoint. \emph{character vector}. Can be \code{"binary"}, \code{"continuous"} or \code{"timeToEvent"}.}
  \item{method}{Is defined when at least one time-to-event outcome is analyzed. Defines the method used to handle pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations.  Can be \code{"Gehan"}, \code{"Peto"}, \code{"Efron"}, or \code{"Peron"}. See details. }  
  \item{n.bootstrap}{the number of bootstrap samples used for computing the confidence interval and the p.values. \emph{integer}. Default is \code{0} meaning no bootstrap (and thus only ponctual estimation).}
  \item{prob.alloc}{the resampling probability for assignement to the experimental group in the bootstrap samples. \emph{double}. Default is \code{NULL} indicating to use the proportion of patients in the experimental group.}
  \item{alternative}{a \emph{character} specifying the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. Default is \code{"two.sided"}.}
  \item{seed}{the seed to consider for the bootstrap. \emph{integer}. Default is \code{10}.}
  \item{cpus}{the number of CPU to use. \emph{integer}. Default is \code{1}.}
  \item{trace}{Should the execution of the function be traced ? \emph{integer}. Default is \code{3}.}
}
\details{
\bold{Treatment:} The variable corresponding to \code{treatment} in data must have only two levels (e.g. \code{0} and \code{1}).

\bold{Endpoint, threshold, censoring, and type:} Arguments \code{endpoint}, \code{threshold}, \code{censoring}  and \code{type} must have the same length. \cr
\code{threshold} must be \code{NA} for binary endpoints and positive for continuous or time to event endpoints. \cr
\code{censoring} must be \code{NA} for binary or continuous endpoints and indicate a variable in data for time to event endpoints. 
Short forms for endpoint \code{type} are \code{"bin"} (binary endpoint), \code{"cont"} (continuous endpoint), \code{"TTE"} (time-to-event endpoint). 

\bold{Bootstrap:} The number of bootstrap replications (argument \code{n.bootstrap}) must be specified to enable the computation of the confidence intervals and the p.value. 
A large number of bootstrap samples (e.g. \code{n.bootstrap=10000})  are needed to obtain accurate CI and p.value. See (Buyse et al., 2010) for more details. 

\bold{Trace:} \code{3} reports all messages  \code{2} reports all messages except silent parallelization messages, \code{1} reports only the percentage of advancement of the bootstrap,  and \code{0} remains silent.

\bold{cpus parallelization:} Argument \code{cpus} can be set to \code{"all"} to use all available cpus. The parallelization relies on the \emph{snowfall} package (function \emph{sfClusterApplyLB}). The detection of the number of cpus relies on the \code{detectCores} function from the \emph{parallel} package .



\bold{Dealing with neutral or uninformative pairs:} Neutral pairs correspond to pairs for which the difference between the endpoint of the control observation and the endpoint of the treatment observation is (in absolute value) below the threshold. When \code{threshold=0}, neutral pairs correspond to pairs with equal outcome.\cr
Uninformative pairs correspond to pairs for which the censoring prevend from classifying them into favorable, unfavorable or neutral. Neutral or uninformative pairs for an endpoint with priority \code{l} are, when available, analysed on the endpoint with priority \code{l-1}.

\bold{Method:} Pairs which can not be decidely classified as favorable, unfavorable, or neutral because of censored observations can be classified uninformative (\code{method="Gehan"}). Another solution is to estimate the probability for such pair to be classified as favorable, unfavorable, or neutral based on the survival functions. \code{method="Peto"} estimate these probabilities using the common Kaplan-Meier estimator of the survival function for treated and control patients. \code{method="Efron"}, and \code{method="Peron"} estimate these probabilities using separate Kaplan-Meier estimators of the survival functions for the two groups of patients. When the largest observation is censored, it is not possible to estimate the survival probability by the Kaplan-Meier estimator beyond this time point.  \code{method="Efron"} treats the largest observations in each patient group as if it were uncensored. \code{method="Peron"} treats the probability of survival beyond the last observation as NA, resulting in a non null probability that the pair is uninformative}  
\value{
  An \R object of class \code{\linkS4class{BuyseRes}}.
}
\references{
Marc Buyse (2010)
Generalized pairwise comparisons of prioritized endpoints in the two-sample problem
\emph{Statistics in Medicine}
\bold{vol. 29} 3245-3257

Efron B (1967)
The two sample problem with censored data
\emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability}
\bold{vol. 4} 831-583

Peto R, Peto J (1972)
Asymptotically efficient rank invariant test procedures
\emph{J R Stat Soc A}
\bold{vol. 135(2)} 185-198

Gehan EA (1965)
A generalized two-sample Wilcoxon test for doubly censored data
\emph{Biometrika}
\bold{vol. 52(3)} 650-653

}
\seealso{
   \code{\link{summary,BuyseRes-method}} for a summary of the results of generalized pairwise comparison. \cr
   \code{\link{BuyseRes-class}} for a presentation of the \code{BuyseRes} object. \cr
   \code{\link{constStrata}} to create a strata variable from several clinical variables. \cr
}
\examples{

#### real example : Veteran dataset of the survival package ####
#### Only one endpoint. Type = Time-to-event. Thresold = 0. Stratfication by histological subtype
#### method = "Gehan"
\dontrun{
data(veteran,package="survival")
require(BuyseTest)
BuyseTest_veteran_Gehan <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
                     strata="celltype",type="timeToEvent",censoring="status",threshold=0,
					 n.bootstrap=1000,method="Gehan",cpus="all")

summary_veteran_Gehan <- summary(BuyseTest_veteran_Gehan)

#### method = "Peron"

BuyseTest_veteran_Peron <- BuyseTest(data=veteran,endpoint="time",treatment="trt",
                     strata="celltype",type="timeToEvent",censoring="status",threshold=0,
					 n.bootstrap=1000,method="Peron",cpus="all")

summary_veteran_Peron <- summary(BuyseTest_veteran_Peron)



#### Several endpoints :
#######Survival, a time-to-event endpoint
#######Toxicity, a continuous/ordinal endpoint : 6 grades of maximal adverse event 

n.Treatment <- 100
n.Control<- 100
prob.Treatment_TOX <- c(0.5,0.25,0.10,0.075,0.05,0.025)
prob.Control_TOX <- c(0.7,0.15,0.05,0.05,0.025,0.025)

lambda.Treatment_TTE <- 0.6
lambda.Control_TTE <- 1


set.seed(10)
data_test <- data.frame(treatment=c(rep(1,n.Treatment),
                                         rep(0,n.Control)  ))
data_test$toxicity <- c(apply(rmultinom(n.Treatment,size=1,
                        prob=prob.Treatment_TOX)==1,2,which),
                        apply(rmultinom(n.Control,size=1,
						prob=prob.Control_TOX)==1,2,which))

data_test$toxicityInv <-6-data_test$toxicity

data_test$EventTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
                               rexp(n.Control,rate=lambda.Control_TTE))
data_test$CensoringTime <- c(rexp(n.Treatment,rate=lambda.Treatment_TTE),
                                  rexp(n.Control,rate=lambda.Control_TTE))
data_test$CensoringTime[data_test$CensoringTime>4] <- 4

data_test$Survival <- apply(data_test[,c("EventTime","CensoringTime")],1,min)
data_test$event <- as.numeric(apply(data_test[,c("EventTime","CensoringTime")],
                                          1,which.min)==1)

resKM_tempo <- survfit(Surv(data_test[,"Survival"],data_test[,"event"])~data_test$treatment)
plot(resKM_tempo)

#### method = "Gehan". 

BuyseTest_severalendpoint_Gehan <- BuyseTest(data=data_test,method="Gehan",
         endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
         treatment="treatment",
         censoring=c("event",NA,"event1",NA,"event",NA),
         type=c("TTE","cont","TTE","cont","TTE","cont"),
         threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
summary(BuyseTest_severalendpoint_Gehan)

#### method = "Peron". 

BuyseTest_severalendpoint_Peron <- BuyseTest(data=data_test,method="Peron",
         endpoint=c("Survival","toxicityInv","Survival","toxicityInv","Survival","toxicityInv"),
         treatment="treatment",
         censoring=c("event",NA,"event",NA,"event",NA),
         type=c("TTE","cont","TTE","cont","TTE","cont"),
         threshold=c(1.5,3,0.75,2,0.25,1),n.bootstrap=1000,trace=2,cpus="all")
summary(BuyseTest_severalendpoint_Peron)




}

}
\keyword{function}
