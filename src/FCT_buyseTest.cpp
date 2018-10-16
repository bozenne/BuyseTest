// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "FCT_calcOnePair.h"
#include "FCT_calcAllPairs.h"
#include "FCT_calcStatistic.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// * Documentation GPC_cpp
//' @title C++ function performing the pairwise comparison over several endpoints. 
//' @description \code{GPC_cpp} call for each endpoint and each strata the pairwise comparison function suited to the type of endpoint and store the results. 
//' @name GPC_cpp
//' 
//' @param endpoint A matrix containing the values of each endpoint (in columns) for each observation (in rows). \emph{const arma::mat&}.
//' @param censoring A matrix containing the values of the censoring variables relative to each endpoint (in columns) for each observation (in rows). \emph{const arma::mat&}.
//' @param indexC A list containing the indexes of control observations belonging for each strata. \emph{List&}  of vector containing positive integers. 
//' @param indexT A list containing the indexes of treatment observations belonging for each strata. \emph{List&} of vector containing positive integers. 
//' @param threshold Store the thresholds associated to each endpoint. \emph{const NumericVector&}. Must have length D. The threshold is ignored for binary endpoints. Must have D columns.
//' @param method The index of the method used to score the pairs. \emph{const IntegerVector&}. Must have length D. 1 for continuous, 2 for Gehan, and 3 for Peron.
//' @param D The number of endpoints. Strictly positive \emph{const int}.
//' @param n_strata The number of strata . Strictly positive \emph{const int}.
//' @param n_TTE The number of time-to-event endpoints. Positive \emph{const int}.
//' @param Wscheme The matrix describing the weighting strategy. For each endpoint (except the first) in column, weights of each pair are initialized at 1 and multiplied by the weight of the endpoints in rows where there is a 1. \emph{const arma::mat&}. Must have D lines and D columns.
//' @param index_survival_M1 The position, among all the survival endpoints, of the last same endpoint (computed with a different threshold). If it is the first time that the TTE endpoint is used it is set to -1. \emph{const IntegerVector}. Must have length n_TTE.
//' @param threshold_M1 The previous latest threshold of each TTE endpoint. When it is the first time that the TTE endpoint is used it is set to -1. \emph{const NumericVector}. Must have length n_TTE.
//' @param list_survTimeC A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows). \emph{List&}.
//' @param list_survTimeT A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows). \emph{List&}. 
//' @param list_survJumpC A list of matrix containing the survival estimates and survival jumps when the survival for the control arm jumps. \emph{List&}.
//' @param list_survJumpT A list of matrix containing the survival estimates and survival jumps when the survival for the treatment arm jumps. \emph{List&}. 
//' @param list_lastSurv A list of matrix containing the last survival estimate in each strata (rows) and treatment group (columns). \emph{List&}. 
//' @param correctionUninf Should the uninformative weight be re-distributed to favorable and unfavorable?
//' @param neutralAsUninf Should paired classified as neutral be re-analyzed using endpoints of lower priority?  \emph{logical}.
//' @param keepScore Should the result of each pairwise comparison be kept? \emph{logical}.
//' @param reserve Should vector storing neutral pairs and uninformative pairs be initialized at their maximum possible length? \emph{logical}.
//' @param returnOnlyDelta Should only the net benefit and win ratio be output? \emph{logical}.
//' @keywords function Cpp BuyseTest

// * Function GPC_cpp
//' @name GPC_cpp
//' @export
// [[Rcpp::export]]
List GPC_cpp(const arma::mat& endpoint,
	     const arma::mat& censoring,
	     const std::vector< arma::uvec >& indexC,
             const std::vector< arma::uvec >& indexT,
	     const NumericVector& threshold,
	     const IntegerVector& method,
	     unsigned int D,
	     int n_strata,
	     int n_TTE, 
             const arma::mat& Wscheme,
	     const IntegerVector& index_survival_M1,
	     const NumericVector& threshold_M1, 
	     const std::vector< std::vector< arma::mat > >& list_survTimeC,
             const std::vector< std::vector< arma::mat > >& list_survTimeT,
             const std::vector< std::vector< arma::mat > >& list_survJumpC,
	     const std::vector< std::vector< arma::mat > >& list_survJumpT,
	     const std::vector< arma::mat >& list_lastSurv,
	     int correctionUninf,
	     bool neutralAsUninf,
	     bool keepScore,
	     bool reserve,
	     bool returnOnlyDelta){

  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.
  
  /// ** initialization
  
  // *** final results
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs by endpoint for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs by endpoint for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs by endpoint for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs by endpoint for each strata
  vector<double> n_pairs(n_strata); // number of pairs sumed over the strats
  
  // *** for a given stata [input]
  // weights of the neutral / uninformative pairs
  arma::mat matWeight;  // for all endpoint up to the current endpoint 
  arma::mat matWeight_M1; // for all endpoint up to the previous endpoint
  
  // *** for a given strata/endpoint
  // index of the pairs
  vector<int> iIndex_control; // in the control arm [current endpoint]
  vector<int> iIndex_treatment;  // in the treatment arm [current endpoint]
  vector<int> iIndex_control_M1; // in the control arm [previous endpoint]
  vector<int> iIndex_treatment_M1;  // in the treatment arm [previous endpoint]

  // weights of the neutral / uninformative pairs
  arma::vec iWeight; // store weights for the current endpoint

  // store indexes of the pairs corresponding to the new weights (iWeight) in the weight vector of the previous endpoint (iCumWeight_M1)
  arma::uvec iIndexWeight_pair; 
  int iSize_weight;
  
  // product of the weights of the neutral / uninformative pairs up to the previous endpoint
  arma::vec iCumWeight_M1;
  // NOTE: it is a special product because weights related to previous survival endpoint are ignored

  // Right-censoring status for TTE endpoint restricted to strata k for a previous endpoint
  arma::mat iSurvTimeC_M1; // at times corresponding to event in the control group
  arma::mat iSurvTimeT_M1; // at times corresponding to event in the treatment group
  arma::mat iSurvJumpC_M1; // at jump times for the survival in the control group
  arma::mat iSurvJumpT_M1; // at jump times for the survival in the treatment group

  // threshold
  double iThreshold_M1;
  
  // *** others
  unsigned int iter_d; // the index of the endpoints
  unsigned int iter_dTTE; // number of time to event endpoints that have been used
  bool iMoreEndpoint;
  arma::uvec iUvec;
  
  // *** keep track of each comparison
  vector<arma::mat> lsScore(D);
  arma::mat iScore;
  arma::mat iMat;
  int iNpairs;

  // ** loop over strata
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
    
    iter_d = 0;
    iter_dTTE = 0;
    
    // *** first endpoint
    // Rcout << endl << "** endpoint 0 **" << endl;
    iMoreEndpoint = (D>1);
    iUvec = arma::uvec({0});
	
    // Rcout << "score" << endl;
    iScore = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec), endpoint.submat(indexT[iter_strata],iUvec), threshold[0],
			  censoring.submat(indexC[iter_strata],iUvec), censoring.submat(indexT[iter_strata],iUvec),
			  list_survTimeC[0][iter_strata], list_survTimeT[0][iter_strata], list_survJumpC[0][iter_strata], list_survJumpT[0][iter_strata],
			  list_lastSurv[0](iter_strata,0), list_lastSurv[0](iter_strata,1), 
			  method[0], correctionUninf,	
			  Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
			  iIndex_control, iIndex_treatment,
			  iWeight,
			  neutralAsUninf, keepScore, iMoreEndpoint, reserve);
    R_CheckUserInterrupt();
	
    if(method[0]>1){ // time to event endpoint
      iter_dTTE++;
    }
	
    // **** add to the total number of pairs the number of pairs founded for this endpoint
    n_pairs[iter_strata] = Mcount_favorable(iter_strata,0) + Mcount_unfavorable(iter_strata,0) + Mcount_neutral(iter_strata,0) + Mcount_uninf(iter_strata,0);
    	
    // **** update weights associated to the remaing pairs
    // Rcout << "> update matWeights" << endl;
    // only relevant if there is one more endpoint
    iSize_weight = iWeight.size();
    if(iMoreEndpoint && (iSize_weight>0) ){
      matWeight.resize(iSize_weight,1);
      matWeight.col(0) = iWeight;
      iWeight.resize(0);
    }

    // Rcout << "> update lsScore" << endl;
    // **** update all Scores
    if(keepScore){
      iNpairs = iScore.n_rows;
      iMat.resize(iNpairs,3);
      iMat.col(0).fill(iter_strata);
      for(int iPair=0; iPair < iNpairs; iPair++){
	iMat(iPair,1) = indexC[iter_strata](iScore(iPair,0));
	iMat(iPair,2) = indexT[iter_strata](iScore(iPair,1));
      }
      // merge with current table and store
      if(iter_strata==0){
	lsScore[0] = arma::join_rows(iMat,iScore);
      }else{
	lsScore[0] = arma::join_cols(lsScore[0],arma::join_rows(iMat,iScore));
      }
    }

    // *** following endpoints
    while((D > iter_d + 1) && (iSize_weight > 0) ){ // loop over the following endpoints

      // while there are remaining endpoints and remaining neutral or uniformative pairs
      iter_d++; // increment the index of the endpoints
      // Rcout << "** endpoint " << iter_d << " (npairs = " << matWeight.n_rows << ") **" << endl;

      // **** copy from previous step
      matWeight_M1 = matWeight; 
      iIndex_control_M1 = iIndex_control;
      iIndex_treatment_M1 = iIndex_treatment;
      iIndex_control.resize(0); iIndex_treatment_M1.resize(0);
	  
      // **** compute the current weights of the pairs
      // Rcout << "> compute cumweight " << endl;
      // matWeight.print("matWeight:");
      // Wscheme.print("Wscheme:");
      // initialize iCumWeight_M1
      iCumWeight_M1.resize(matWeight.n_rows);
      iCumWeight_M1.fill(1.0);
      for(unsigned int iter_endpoint=0 ; iter_endpoint<iter_d ; iter_endpoint++){
	if(Wscheme(iter_endpoint,iter_d)==1){iCumWeight_M1 %= matWeight.col(iter_endpoint);}
      }
      // Rcout << "total weight M1 = " << sum(iCumWeight_M1) << endl;

      // **** computes scores
      iMoreEndpoint = (D>(iter_d+1));
      iUvec = arma::uvec({iter_d});
	
      // Rcout << "> score " << endl;
      // Rcout << method[iter_d] << endl;
      if(method[iter_d]==3 && index_survival_M1[iter_dTTE]>=0){
	// survival endpoint already analyzed with Peron
	iSurvTimeC_M1 = list_survTimeC[index_survival_M1[iter_dTTE]][iter_strata];
	iSurvTimeT_M1 = list_survTimeT[index_survival_M1[iter_dTTE]][iter_strata];
	iSurvJumpC_M1 = list_survJumpC[index_survival_M1[iter_dTTE]][iter_strata];
	iSurvJumpT_M1 = list_survJumpT[index_survival_M1[iter_dTTE]][iter_strata];
	iThreshold_M1 = threshold_M1[iter_dTTE];
      }else{
	iSurvTimeC_M1 = arma::mat(0,0);
	iSurvTimeT_M1 = arma::mat(0,0);
	iSurvJumpC_M1 = arma::mat(0,0);
	iSurvJumpT_M1 = arma::mat(0,0);
	iThreshold_M1 = NA_REAL;
      }
	  
      iScore = calcSubsetPairs(endpoint.submat(indexC[iter_strata],iUvec), endpoint.submat(indexT[iter_strata],iUvec), threshold[iter_d],
			       censoring.submat(indexC[iter_strata],iUvec), censoring.submat(indexT[iter_strata],iUvec),
			       list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata], list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
			       list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1),
			       iIndex_control_M1, iIndex_treatment_M1,
			       iCumWeight_M1, iThreshold_M1,
			       iSurvTimeC_M1, iSurvTimeT_M1, iSurvJumpC_M1, iSurvJumpT_M1,
			       method[iter_d], correctionUninf,	
			       Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
			       iIndex_control, iIndex_treatment,
			       iWeight, iIndexWeight_pair,
			       neutralAsUninf, keepScore, iMoreEndpoint, reserve);
      R_CheckUserInterrupt();
	  
      if(method[iter_d]>1){ // time to event endpoint
	iter_dTTE++;
      }
	  
      // **** update weights/index associated to the remaing pairs
      // Rcout << "> update matWeight " << endl;
      // only relevant if there is one more endpoint and the weights may differ from 1
      iSize_weight = iIndexWeight_pair.size();
      if(iMoreEndpoint && (iSize_weight > 0) ){
	matWeight.set_size(iSize_weight, iter_d+1);
		  
	// store iMweight_M1 in the iMweight restrected to the remaining pairs
	// iUvec.print("indexRow");
	iUvec = arma::regspace<uvec>(0, 1, iter_d - 1);
	matWeight.cols(iUvec) = matWeight_M1.submat(iIndexWeight_pair, iUvec);

	// add the weight relative to the new endpoint
	// matWeight.col(iter_d) = conv_to<colvec>::from(iWeight);
	matWeight.col(iter_d) = iWeight;
		
        iWeight.resize(0); iIndexWeight_pair.resize(0);
      }

      // **** update all Scores
      // Rcout << "> update Scores" << endl;
      if(keepScore){
	iNpairs = iScore.n_rows;
	iMat.resize(iNpairs,3);
	iMat.col(0).fill(iter_strata);
	for(int iPair=0; iPair < iNpairs; iPair++){
	  iMat(iPair,1) = indexC[iter_strata](iScore(iPair,0));
	  iMat(iPair,2) = indexT[iter_strata](iScore(iPair,1));
	}
	// merge with current table and store
	if(iter_strata==0){
	  lsScore[iter_d] = arma::join_rows(iMat,iScore);
	}else{
	  lsScore[iter_d] = arma::join_cols(lsScore[iter_d], arma::join_rows(iMat,iScore));
	}
      }
    } // end endpoint

  }
  
  // ** proportion in favor of treatment 
  arma::mat delta_netBenefit(n_strata,D), delta_winRatio(n_strata,D); // matrix containing for each strata and each endpoint the statistic
  vector<double> Delta_netBenefit(D), Delta_winRatio(D); // vector containing for each endpoint the overall statistic
  
  calcStatistic(delta_netBenefit, delta_winRatio, Delta_netBenefit, Delta_winRatio,
                Mcount_favorable, Mcount_unfavorable, 
                D, n_strata, n_pairs);

  // ** export
  if(returnOnlyDelta){
      return(List::create(
		      Named("delta_netBenefit")  = delta_netBenefit,
		      Named("delta_winRatio")  = delta_winRatio,
		      Named("Delta_netBenefit")  = Delta_netBenefit,
		      Named("Delta_winRatio")  = Delta_winRatio
		      ));
  }else{
    return(List::create(
			Named("count_favorable")  = Mcount_favorable,
			Named("count_unfavorable")  = Mcount_unfavorable,
			Named("count_neutral")  = Mcount_neutral,           
			Named("count_uninf")  = Mcount_uninf,           
			Named("delta_netBenefit")  = delta_netBenefit,
			Named("delta_winRatio")  = delta_winRatio,
			Named("Delta_netBenefit")  = Delta_netBenefit,
			Named("Delta_winRatio")  = Delta_winRatio,
			Named("n_pairs")  = n_pairs,
			Named("tableScore")  = lsScore
			));
  }
  
}


