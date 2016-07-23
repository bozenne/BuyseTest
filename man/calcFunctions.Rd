\name{calcFunctions}
\title{calculation functions}
\alias{calcBootstrap}
\alias{calcCI}
% \alias{calcSample}
% \alias{calcSimulation}
\description{
Functions called by \code{\link{BuyseTest}} %and \code{\link{BuysePower}} 
to perform calculations.
}
\usage{
calcBootstrap(envir)

calcCI(delta,delta_boot,endpoint,D,alternative,alpha,
         n.bootstrap,cpus,trace)
				   
% calcSample(n.Treatment,n.Control,hypothesis,endpoint,
%          type,D,D.TTE,prob.alloc)
						  
% calcSimulation(envir)
}
\details{
\code{calcBootstrap} performs the bootstrap iterations in sequential or parallel mode. \cr
\code{calcCI} post process the results to obtain the p.value and the confidence interval. \cr
% \code{calcSample} generate a sample according a given parametrization. \cr
% \code{calcSimulation} performs the simulation iterations in sequential or parallel mode. \cr
}
\keyword{function}
\keyword{internal}

