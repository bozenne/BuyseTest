\name{calcOnePairs_cpp}
\alias{calcOnePair_BinaryOutcome_cpp}
\alias{calcOnePair_ContinuousOutcome_cpp}
\alias{calcOnePair_TTEOutcome_Gehan_cpp}
\alias{calcOneProba_TTEOutcome_Peto_cpp}
\alias{calcOneProba_TTEOutcome_Efron_cpp}
\alias{calcIntegralProba_cpp}
\title{
C++ function performing the pairwise comparison over one endpoint. 
}
\description{
\code{calcOnePair_BinaryOutcome_cpp} compairs two binary endpoints.

\code{calcOnePair_ContinuousOutcome_cpp} compairs two continuous endpoints.

\code{calcOnePair_TTEOutcome_Gehan_cpp} compairs two time to event endpoints classifying censored paired as uninformative.

\code{calcOneProba_TTEOutcome_Peto_cpp} computes the probability in favor of treatment, the probability in defavor of treatment for a time to event endpoint using a commun survival estimator between treatment and control patients.

\code{calcOneProba_TTEOutcome_Efron_cpp} computes the probability in favor of treatment, the probability in defavor of treatment for a time to event endpoint using a different survival estimator between treatment and control patients.

\code{calcIntegralProba_cpp} compute an integral over the survival steps.
}
\usage{
calcOnePair_BinaryOutcome_cpp( endpoint_T, endpoint_C, index_T, index_C, 
         Wpair, iter_pair,
         count_favorable, count_unfavorable, count_neutral, count_uninf,
         index_uninfT, index_uninfC, index_neutralT, index_neutralC,
         index_wUninf, index_wNeutral)

calcOnePair_ContinuousOutcome_cpp(endpoint_T, endpoint_C, threshold, index_T, index_C, 
         Wpair, iter_pair, 
         count_favorable, count_unfavorable, count_neutral, count_uninf,
         index_uninfT, index_uninfC, index_neutralT, index_neutralC,
         index_wUninf,  index_wNeutral)

calcOnePair_TTEOutcome_Gehan_cpp(endpoint_T, endpoint_C, delta_T, delta_C, 
         threshold, index_T, index_C, Wpair, iter_pair, 
         count_favorable, count_unfavorable, count_neutral, count_uninf,
         index_uninfT, index_uninfC, index_neutralT, index_neutralC,
         index_wUninf, index_wNeutral)

calcOneProba_TTEOutcome_Peto_cpp(endpoint_T, endpoint_C, delta_T, delta_C, 
         threshold, index_T, index_C,
         survival_T, survival_C)

calcOneProba_TTEOutcome_Efron_cpp(endpoint_T, endpoint_C, delta_T, delta_C,
         threshold, index_T, index_C,
         survival_T, survival_C)

calcIntegralProba_cpp(time_begin, survival, index_col)
}
\keyword{function}
\keyword{internal}

