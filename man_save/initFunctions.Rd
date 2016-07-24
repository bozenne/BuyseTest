\name{initFunctions}
\title{initialisation functions}
\alias{initCensoring}
\alias{initCPU}
\alias{initData}
\alias{initEndpoint}
\alias{initN}
\alias{initHypothesis}
\alias{initAlloc}
\alias{initSpace}
\alias{initStrata}
\alias{initSurvival}
\alias{initThreshold}
\alias{initTreatment}
\alias{initType}
\alias{initWscheme}
\description{
Functions called by \code{\link{BuyseTest}} %and \code{\link{BuysePower}} 
to initialize the arguments.
}
\usage{
initCensoring(censoring,data,endpoint,type,D,D.TTE,
         treatment,strata,callMethod="initCensoring")
						  
initCPU(cpus,callMethod="initCPU")

initData(data, Ind.Treatment, type, endpoint, D, censoring,
         index.strataT, index.strataC, n.strata, method, D.TTE,
         threshold, Wscheme = NULL, trace, test = TRUE,
         callMethod = "initData")
					 
initEndpoint(endpoint,data,D,
         treatment,censoring,strata,callMethod="initEndpoint")
						 
initN(n,callMethod="initN")
						  
initHypothesis(hypothesis,type,D,callMethod="initParamSimul")

initAlloc(proba,n.Treatment,n.Control,callMethod="initAlloc")

initSpace(nchar)

initStrata(strata, data, Ind.Treatment, n.Treatment, n.Control,
           endpoint, treatment, censoring, callMethod ="initStrata")
	
initSurvival(M.Treatment,M.Control,M.delta_Treatment,M.delta_Control,
             endpoint,D.TTE,type,threshold,
             index.strataT,index.strataC,n.strata,
             method,callMethod="initSurvival")
	
initThreshold(threshold,type,D,method,endpoint,callMethod="initThreshold")

initTreatment(data,treatment,endpoint,censoring,strata,callMethod="initTreatment")

initType(type,callMethod="initType")

initWscheme(endpoint, D, endpoint.TTE,D.TTE,threshold,type)
}
\arguments{
  \item{Ind.Treatment}{The indicator of the treatment observations (1 treatment, 0 control). \code{integer vector}.}
  \item{D}{The number of endpoints. \code{integer}.}
  \item{D.TTE}{The number of event endpoints. \code{integer}.}
  \item{endpoint.TTE}{The time to event endpoints. \code{character vector}.}
  \item{initType}{The name of the function to display in the error message. \code{character}.}
}
\keyword{function}
\keyword{internal}

