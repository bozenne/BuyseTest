\name{calcAllPairs_cpp}
\alias{calcAllPairs_BinaryOutcome_cpp}
\alias{calcSubsetPairs_BinaryOutcome_cpp}
\alias{calcSubsetPairs_WeightedBinaryOutcome_cpp}
\alias{calcAllPairs_ContinuousOutcome_cpp}
\alias{calcSubsetPairs_ContinuousOutcome_cpp}
\alias{calcSubsetPairs_WeightedContinuousOutcome_cpp}
\alias{calcAllPairs_TTEOutcome_Gehan_cpp}
\alias{calcSubsetPairs_TTEOutcome_Gehan_cpp}
\alias{calcAllPairs_TTEOutcome_PetoEfron_cpp}
\alias{calcSubsetPairs_TTEOutcome_PetoEfron_cpp}
\title{
C++ function performing the pairwise comparison over one endpoint. 
}
\description{
\code{calcAllPairs_BinaryOutcome_cpp}, \code{calcSubsetPairs_BinaryOutcome_cpp} and \code{calcSubsetPairs_WeightedBinaryOutcome_cpp} functions 
perform pariwise comparison over a binary endpoint. 
\code{calcAllPairs_BinaryOutcome_cpp} compairs all pairs, \code{calcSubsetPairs_BinaryOutcome_cpp} compairs the pred-specified pairs and 
\code{calcSubsetPairs_WeightedBinaryOutcome_cpp} compairs the pred-specified pairs taking into account of the weight associated to each pair.

\code{calcAllPairs_ContinuousOutcome_cpp}, \code{calcSubsetPairs_ContinuousOutcome_cpp} and \code{calcSubsetPairs_WeightedContinuousOutcome_cpp} functions 
perform pariwise comparison over a continuous endpoint.
\code{calcAllPairs_ContinuousOutcome_cpp} compairs all pairs, \code{calcSubsetPairs_ContinuousOutcome_cpp} compairs the pred-specified pairs and 
\code{calcSubsetPairs_WeightedContinuousOutcome_cpp} compairs the pred-specified pairs taking into account of the weight associated to each pair.

\code{calcAllPairs_TTEOutcome_Gehan_cpp}, \code{calcSubsetPairs_TTEOutcome_Gehan_cpp}, \code{calcAllPairs_TTEOutcome_PetoEfron_cpp} and \code{calcSubsetPairs_TTEOutcome_PetoEfron_cpp} functions 
perform pariwise comparison over a time to event endpoint. The first two classify censored paired as uninformative and the last two uses the expected survival to compare the pairs.
\code{calcAllPairs_ContinuousOutcome_cpp} and \code{calcAllPairs_TTEOutcome_PetoEfron_cpp} compairs all pairs, \code{calcSubsetPairs_TTEOutcome_Gehan_cpp} compairs the pred-specified pairs and 
\code{calcSubsetPairs_TTEOutcome_PetoEfron_cpp} compairs the pred-specified pairs taking into account of the weight associated to each pair.
}
\usage{
calcAllPairs_BinaryOutcome_cpp(Treatment, Control)

calcSubsetPairs_BinaryOutcome_cpp(Treatment, Control, 
         index_neutralT, index_neutralC, nNeutral_pairs,
         index_uninfT, index_uninfC, nUninf_pairs)

calcSubsetPairs_WeightedBinaryOutcome_cpp(Treatment, Control, 
         index_neutralT,index_neutralC, nNeutral_pairs,
         index_uninfT, index_uninfC, nUninf_pairs,
         Wpairs) 

calcAllPairs_ContinuousOutcome_cpp(Treatment, Control, threshold);

calcSubsetPairs_ContinuousOutcome_cpp(Treatment, Control, threshold, 
         index_neutralT, index_neutralC, nNeutral_pairs, 
         index_uninfT, index_uninfC, nUninf_pairs);

calcSubsetPairs_WeightedContinuousOutcome_cpp(Treatment, Control, threshold, 
         index_neutralT, index_neutralC, nNeutral_pairs, 
         index_uninfT, index_uninfC, nUninf_pairs,
         Wpairs) ;

calcAllPairs_TTEOutcome_Gehan_cpp(Treatment, Control, threshold, 
         deltaT, deltaC);

calcSubsetPairs_TTEOutcome_Gehan_cpp(Treatment, Control, threshold, deltaT, deltaC, 
         index_neutralT, index_neutralC, nNeutral_pairs,
         index_uninfT, index_uninfC, nUninf_pairs);

calcAllPairs_TTEOutcome_PetoEfron_cpp(Treatment, Control, threshold,
         deltaT, deltaC, matKMT, matKMC,
         Efron)

calcSubsetPairs_TTEOutcome_PetoEfron_cpp(Treatment, Control, threshold, 
         deltaT, deltaC, matKMT, matKMC,
         index_neutralT, index_neutralC, nNeutral_pairs, 
         index_uninfT, index_uninfC, nUninf_pairs,
         Wpairs, threshold_M1, matKMT_M1, matKMC_M1, Efron)

}
\keyword{function}
\keyword{internal}

