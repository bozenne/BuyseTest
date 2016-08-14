\name{BuyseTest_cpp}
\alias{BuyseTest_Gehan_cpp}
\alias{BuyseTest_PetoEfronPeron_cpp}
\title{
C++ function performing the pairwise comparison over several endpoints. 
}
\description{
\code{BuyseTest_Gehan_cpp} and \code{BuyseTest_PetoEfron_cpp} functions calls for each endpoint and each strata the pairwise comparison function suited to the type of endpoint and store the resuts. 
}
\usage{
BuyseTest_Gehan_cpp(Treatment, Control, threshold, type, 
         delta_Treatment, delta_Control,
         D, returnIndex, strataT, strataC, n_strata, n_TTE)
		 
BuyseTest_PetoEfronPeron_cpp(Treatment, Control, threshold, type, 
         delta_Treatment, delta_Control,
         D, returnIndex, strataT, strataC, n_strata, n_TTE, 
         Wscheme, index_survivalM1, threshold_TTEM1, list_survivalT, list_survivalC, PEP)
}
\arguments{
  \item{Treatment}{A matrix containing the values of each endpoint (in columns) for the treatment observations (in rows). \emph{const arma::mat&}. Must have D columns.}
  \item{Control}{A matrix containing the values of each endpoint (in columns) for the control observations (in rows). \emph{const arma::mat&}.}
  \item{threshold}{Store the thresholds associated to each endpoint. \emph{const NumericVector&}. Must have length D. The threshold is ignored for binary endpoints. Must have D columns.}
  \item{type}{The type of each endpoint (1 binary, 2 continuous, 3 TTE). \emph{const IntegerVector&}. Must have length D.}
  \item{delta_Treatment}{A matrix containing in the type of event (0 censoring, 1 event) for each TTE endpoint (in columns) and treatment observations (in rows). \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.}
  \item{delta_Control}{A matrix containing the nature of observations in the control group (in rows) (0 censoring, 1 event) for each TTE endpoint (in columns) . \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.}
  \item{D}{The number of endpoints. Strictly positive \emph{const int}.}
  \item{returnIndex}{Should the indexes of the neutral or uninformative pairs be returned. \emph{const bool}.}
  \item{strataT}{A list containing the indexes of treatment observations belonging for each strata. \emph{List&} of vector containing positive integers. }
  \item{strataC}{A list containing the indexes of control observations belonging for each strata. \emph{List&}  of vector containing positive integers. }
  \item{n_strata}{The number of strata . Strictly positive \emph{const int}.}
  \item{n_TTE}{The number of time-to-event endpoints. Positive \emph{const int}.}

  \item{Wscheme}{The matrix describing the weighting strategy. For each endpoint (except the first) in column, weights of each pair are initialized at 1 and multiplyied by the weight of the endpoints in raws where there is a 1. \emph{const arma::mat&}. Must have n_TTE lines and D-1 columns.}
  \item{index_survivalM1}{The position, among all the survival endpoints, of the last same endpoint (computed with a different threshold). If it is the first time that the TTE endpoint is used it is set to -1. \emph{const IntegerVector}. Must have length n_TTE.}
  \item{threshold_TTEM1}{The previous latest threshold of each TTE endpoint. When it is the first time that the TTE endpoint is used it is set to -1. \emph{const NumericVector}. Must have length n_TTE.}
  \item{list_survivalT}{A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows). \emph{List&}. Must have length n_TTE. Each matrix must have 3 (if method is Peto, only one survival function is computed) or 11 (if method is Efron or Peron, two survival functions are computed) columns. Ignored if method is Gehan.}
  \item{list_survivalC}{A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows). \emph{List&}. Must have length n_TTE. Each matrix must have 3 (if method is Peto) or 11 (if method is Efron or Peron) columns. Ignored if method is Gehan.}
  \item{PEP}{The type of method used to compare censored pairs (1 Peto, 2 Efron, 3 Peron).\emph{const int}.}
}
\keyword{function}
\keyword{internal}
