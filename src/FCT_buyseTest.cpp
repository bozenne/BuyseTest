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
//' @param endpoint A matrix containing the values of each endpoint (in columns) for each observation (in rows). 
//' @param censoring A matrix containing the values of the censoring variables relative to each endpoint (in columns) for each observation (in rows).
//' @param indexC A list containing, for each strata, which rows of the endpoint and censoring matrices corresponds to the control observations. Not unique when bootstraping.
//' @param posC A list containing, for each strata, the unique identifier of each control observations. 
//' @param indexT A list containing, for each strata, which rows of the endpoint and censoring matrices corresponds to the treatment observations. Not unique when bootstraping.
//' @param posT A list containing, for each strata, the unique identifier of each treatment observations.
//' @param threshold Store the thresholds associated to each endpoint. Must have length D. The threshold is ignored for binary endpoints. 
//' @param indexEndpoint_M1 Store the position of the previous occurence of the endpoint. Must have length D. Only used for time to event endpoints and when method is Peron.
//' @param weight Store the weight associated to each endpoint. Must have length D. 
//' @param method The index of the method used to score the pairs. Must have length D. 1 for continuous, 2 for Gehan, and 3 for Peron.
//' @param D The number of endpoints.
//' @param n_strata The number of strata. 
//' @param n_TTE The number of time-to-event endpoints. 
//' @param n_UTTE The number of unique time-to-event endpoints.
//' @param Wscheme The matrix describing the weighting strategy. For each endpoint (except the first) in column, weights of each pair are initialized at 1 and multiplied by the weight of the endpoints in rows where there is a 1. Must have D lines and D columns.
//' @param index_endpoint The position of the endpoint at each priority in the argument endpoint. Must have length D. 
//' @param index_censoring The position of the censoring at each priority in the argument censoring. Must have length D. 
//' @param index_UTTE The position, among all the unique tte endpoints, of the TTE endpoints. Equals -1 for non tte endpoints. Must have length n_TTE. 
//' @param reanalyzed Will this endpoint be re-analyzed latter with a different threshold.
//' @param list_survTimeC A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows).
//' @param list_survTimeT A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows).
//' @param list_survJumpC A list of matrix containing the survival estimates and survival jumps when the survival for the control arm jumps.
//' @param list_survJumpT A list of matrix containing the survival estimates and survival jumps when the survival for the treatment arm jumps.
//' @param list_lastSurv A list of matrix containing the last survival estimate in each strata (rows) and treatment group (columns).
//' @param p_C Number of nuisance parameter in the survival model for the control group, for each endpoint and strata
//' @param p_T Number of nuisance parameter in the survival model for the treatment group, for each endpoint and strata
//' @param iid_survJumpC A list of matrix containing the iid of the survival estimates in the control group.
//' @param iid_survJumpT A list of matrix containing the iid of the survival estimates in the treatment group.
//' @param correctionUninf Should the uninformative weight be re-distributed to favorable and unfavorable?
//' @param hierarchical Should only the uninformative pairs be analyzed at the lower priority endpoints (hierarchical GPC)? Otherwise all pairs will be compaired for all endpoint (full GPC).
//' @param hprojection Order of the H-projection used to compute the variance.
//' @param neutralAsUninf Should paired classified as neutral be re-analyzed using endpoints of lower priority? 
//' @param keepScore Should the result of each pairwise comparison be kept?
//' @param reserve Should vector storing neutral pairs and uninformative pairs be initialized at their maximum possible length?
//' @param returnIID Should the iid be computed?
//' @keywords function Cpp BuyseTest

// * Function GPC_cpp
//' @name GPC_cpp
//' @export
// [[Rcpp::export]]
List GPC_cpp(arma::mat endpoint,
			 arma::mat censoring,
			 std::vector< arma::uvec > indexC,
			 std::vector< arma::uvec > posC,
			 std::vector< arma::uvec > indexT,
			 std::vector< arma::uvec > posT,
			 std::vector< double > threshold,
			 std::vector< int > indexEndpoint_M1,
			 arma::vec weight,
			 std::vector< int > method,
			 unsigned int D,
			 unsigned int n_strata,
			 unsigned int n_TTE, 
			 int n_UTTE, 
			 arma::mat Wscheme,
			 std::vector<int> index_endpoint, 
			 std::vector<int> index_censoring, 
			 std::vector<int> index_UTTE, 
			 std::vector<bool> reanalyzed, 
			 std::vector< std::vector< arma::mat > > list_survTimeC,
			 std::vector< std::vector< arma::mat > > list_survTimeT,
			 std::vector< std::vector< arma::mat > > list_survJumpC,
			 std::vector< std::vector< arma::mat > > list_survJumpT,
			 std::vector< arma::mat > list_lastSurv,
			 arma::mat p_C,
			 arma::mat p_T,
			 std::vector< std::vector< arma::mat > > iid_survJumpC,
			 std::vector< std::vector< arma::mat > > iid_survJumpT,
			 int correctionUninf,
			 bool hierarchical,
			 int hprojection,
			 bool neutralAsUninf,
			 bool keepScore,
			 bool reserve,
			 int returnIID){

  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.

  /// ** initialization
  // *** object storing the final results
  // esperance
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs by endpoint for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs by endpoint for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs by endpoint for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs by endpoint for each strata
  arma::vec n_pairs(n_strata); // number of pairs sumed over the strats
  arma::vec n_treatment(n_strata); // number of patients in the treatment group over the strats
  arma::vec n_control(n_strata); // number of patients in the control group over the strats
  arma::vec n_cumpairsM1(n_strata); // number of pairs in the previous strata
  
  // variance
  arma::mat iid_favorable; // iid for the favorable scores across endpoints
  arma::mat iid_unfavorable; // iid for the unfavorable scores across endpoints
  arma::colvec iIIDNuisance_favorable; // iid relative to the KM estimates for the favorable scores across endpoints
  arma::colvec iIIDNuisance_unfavorable; // iid relative to the KM estimates for the unfavorable scores across endpoints
  arma::mat Mvar; // total variance-covariance (favorable,unfavorable) scores across endpoints
  if(returnIID>0){
    int n_obs = endpoint.n_rows;

    // iid Pn
    iid_favorable.resize(n_obs,D);
	iid_favorable.fill(0.0);
    iid_unfavorable.resize(n_obs,D);
	iid_unfavorable.fill(0.0);

	// iid nuisance parameters
    if(returnIID>1){  
      iidNuisance_favorable.resize(n_obs,D);
      iidNuisance_favorable.fill(0.0);
      iidNuisance_unfavorable.resize(n_obs,D);
      iidNuisance_unfavorable.fill(0.0);
    }
	
	// variance
    Mvar.resize(D,5);
    Mvar.fill(0.0);

  }

  // *** object necessary for the basic GPC
  // global
  double zeroPlus = 1e-12;

  // at a given strata/endpoint
  bool iMoreEndpoint; // is there any endpoint after this one (if so do not store detailed information about the residual pairs)
  int iMethod; // how the score are computed: Gehan or Peron
  bool iReanalyzed; // for method=Peron: will this endpoint be re-analyzed? If so the score and iid are stored.
  int iIndex_UTTE; // which TTE endpoint 

  std::vector< int > iIndex_control; // index of the pairs in the control arm [current endpoint]
  std::vector< int > iIndex_treatment;  // index of the pairs in the treatment arm [current endpoint]
  std::vector< int > iIndex_control_M1; // index of the pairs in the control arm [previous endpoint]
  std::vector< int > iIndex_treatment_M1;  // index of the pairs in the treatment arm [previous endpoint]


  // residual pairs (RP): pairs with non-0 neutral or uninformative score
  arma::mat iResults_RP; // favorable/unfavorable/(neutral+uninformative) of the residual pairs [current endpoint]
  arma::mat iResults_M1_RP; // favorable/uninformative/weight of the residual pairs. Weights are (neutral+uninformative) cumulated over the previous endpoints.
  arma::uvec iIndex_RP; // index of the new residual pairs among the previous residual pairs.
  int iSize_RP; // number of residual pairs.
  
  // *** object necessary for the Peron correction
  // global

  // at a given strata/endpoint
  double iThreshold_M1;
  arma::mat iSurvTimeC_M1;
  arma::mat iSurvTimeT_M1;
  arma::mat iSurvJumpC_M1;
  arma::mat iSurvJumpT_M1;
  double iLastSurvC_M1;
  double iLastSurvT_M1;

  std::vector< arma::mat > iLs_matIIDC(3);
  std::vector< arma::mat > iLs_matIIDT(3);

  // *** object necessary for the iid
  // global

  // at a given strata/endpoint
  arma::mat iPartialCount_C; // iid in the control group [current endpoint]
  arma::mat iPartialCount_T; // iid in the treatment group [current endpoint]

  // *** object necessary for the iid nuisance
  // global

  // at a given strata/endpoint

  
  // ***  object necessary to store the score of each pair over the endpoints
  // global
  std::vector< arma::mat > lsScore(D);

  // at a given strata/endpoint
  arma::mat iScore;
  arma::mat iMat;
  int iNpairs;
  


  //
  arma::mat iEdSurvC; // temporary derivative regarding the survival parameters 
  arma::mat iEdSurvT; // temporary derivative regarding the survival parameters 
  arma::mat iidNuisance_favorable; 
  arma::mat iidNuisance_unfavorable;

  // favorable/unfavorable/weights of the neutral / uninformative pairs
  std::vector< arma::mat > ls_previousScore(3);  // for all endpoint up to the current endpoint 
  std::vector< arma::mat > ls_previousScore2(3); // (temporary, just to update ls_previousScore)
  
  // for which unique tte endpoint calculation have been stored
  uvec isStored_UTTE = zeros<arma::uvec>(n_UTTE); // weights
  
  // *** others
  arma::uvec iUvec_iter_d(1);
  arma::uvec iUvec_iter_d2(1);
  arma::uvec iUvec_endpoint(1);
  arma::uvec iUvec_censoring(1);
  arma::uvec iUvec;

  
  // ** loop over strata
  // Wscheme.print("Wscheme:");
  // Rcout << "start" << endl;
  for(unsigned int iter_strata=0 ; iter_strata < n_strata ; iter_strata++){

    for(unsigned int iter_d=0 ; iter_d < D; iter_d++){
      // Rcout << endl << "** endpoint " << iter_d << "**" << endl;

      // **** type of endpoint
      iMoreEndpoint = (D>(iter_d+1));
      iMethod = method[iter_d];
      iReanalyzed = reanalyzed[iter_d];
      iIndex_UTTE = index_UTTE[iter_d];
      iUvec_endpoint[0] = index_endpoint[iter_d];
      iUvec_censoring[0] = index_censoring[iter_d];

      // **** compute the current weights of the pairs
      if((iter_d > 0) && hierarchical){
        // Rcout << " compute cumweight: ";
		// matWeight.print("matWeight:");
		// initialize iCumWeight_M1
		iCumWeight_M1.resize(matWeight.n_rows);
		iCumWeight_M1.fill(1.0);
		for(unsigned int iter_endpoint=0 ; iter_endpoint<iter_d ; iter_endpoint++){
		  if(Wscheme(iter_endpoint,iter_d)==1){iCumWeight_M1 %= matWeight.col(iter_endpoint);}
		}
		// Rcout << " total weight M1 = " << sum(iCumWeight_M1) << endl;
      }

      // **** compute scores
      // Rcout << " score" << endl;
      if((iter_d==0) || (hierarchical == false)){
		iScore = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
							  censoring.submat(indexC[iter_strata],iUvec_censoring), censoring.submat(indexT[iter_strata],iUvec_censoring),
							  list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata],
							  list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
							  list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1), 
							  iMethod, correctionUninf,	 p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),
							  Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d),
							  Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d),
							  iIndex_control, iIndex_treatment,
							  iWeight, 
							  iPartialCount_C, iPartialCount_T, iEdSurvC, iEdSurvT,
							  iLs_matIIDC, iLs_matIIDT, returnIID, 
							  neutralAsUninf, keepScore, iMoreEndpoint & hierarchical, iReanalyzed, reserve);
		// add to the total number of pairs the number of pairs found for this endpoint
		if(iter_d==0){
		  n_control[iter_strata] = posC[iter_strata].size();
		  n_treatment[iter_strata] = posT[iter_strata].size();		  
		  n_pairs[iter_strata] = n_control[iter_strata] * n_treatment[iter_strata];
		  //= Mcount_favorable(iter_strata,0) + Mcount_unfavorable(iter_strata,0) + Mcount_neutral(iter_strata,0) + Mcount_uninf(iter_strata,0);
		  if(iter_strata == 0){
			n_cumpairsM1[0] = 0;
		  }else{
			n_cumpairsM1[iter_strata] = n_cumpairsM1[iter_strata-1] + n_control[iter_strata]*n_treatment[iter_strata];
		  }
		}
      }else{

		if(indexEndpoint_M1[iter_d]<0){
		  iThreshold_M1 = -100;
		  iSurvTimeC_M1.resize(0,0);
		  iSurvTimeT_M1.resize(0,0);
		  iSurvJumpC_M1.resize(0,0);
		  iSurvJumpT_M1.resize(0,0);
		  iLastSurvC_M1 = -100;
		  iLastSurvT_M1 = -100;
		}else{
		  iThreshold_M1 = threshold[indexEndpoint_M1[iter_d]];
		  iSurvTimeC_M1 = list_survTimeC[indexEndpoint_M1[iter_d]][iter_strata];
		  iSurvTimeT_M1 = list_survTimeT[indexEndpoint_M1[iter_d]][iter_strata];
		  iSurvJumpC_M1 = list_survJumpC[indexEndpoint_M1[iter_d]][iter_strata];
		  iSurvJumpT_M1 = list_survJumpT[indexEndpoint_M1[iter_d]][iter_strata];
		  iLastSurvC_M1 = list_lastSurv[indexEndpoint_M1[iter_d]](iter_strata,0);
		  iLastSurvT_M1 = list_lastSurv[indexEndpoint_M1[iter_d]](iter_strata,1);
		}
		
		iScore = calcSubsetPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
								 censoring.submat(indexC[iter_strata],iUvec_censoring), censoring.submat(indexT[iter_strata],iUvec_censoring),
								 list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata],
								 list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
								 list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1),
								 iMethod, correctionUninf, p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),	
								 Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d),
								 Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
								 iIndex_control, iIndex_treatment, 
								 iWeight, iIndexWeight_pair, 
								 iPartialCount_C, iPartialCount_T, iEdSurvC, iEdSurvT, returnIID, 
								 neutralAsUninf, keepScore, iMoreEndpoint, iReanalyzed, reserve,
								 iIndex_control_M1, iIndex_treatment_M1, iCumWeight_M1,
								 iThreshold_M1,iSurvTimeC_M1,iSurvTimeT_M1,iSurvJumpC_M1,iSurvJumpT_M1,iLastSurvC_M1,iLastSurvT_M1);
      }
      R_CheckUserInterrupt();
	

      // **** update iid
      // Rcout << " update iid (" << returnIID << ")" << endl;
      if(returnIID>0){
		iUvec_iter_d = {iter_d};
		iid_favorable.submat(posC[iter_strata], iUvec_iter_d) = iPartialCount_C.col(0);
		iid_favorable.submat(posT[iter_strata], iUvec_iter_d) = iPartialCount_T.col(0);
	
		iid_unfavorable.submat(posC[iter_strata], iUvec_iter_d) = iPartialCount_C.col(1);
		iid_unfavorable.submat(posT[iter_strata], iUvec_iter_d) = iPartialCount_T.col(1);

		if(returnIID>1){
		  // new score is w_old * favorable so IF = IF_old * favorable + w_old * IF_favorable
		  // Rcout << "iid nuisance " << endl;
		  if(iMethod == 3){  // w_old * IF_favorable
			iEdSurvC /= n_pairs[iter_strata];
			iEdSurvT /= n_pairs[iter_strata];
			iIIDNuisance_favorable = iid_survJumpC[iIndex_UTTE][iter_strata] * iEdSurvC.col(0) + iid_survJumpT[iIndex_UTTE][iter_strata] * iEdSurvT.col(0);
			iIIDNuisance_unfavorable = iid_survJumpC[iIndex_UTTE][iter_strata] * iEdSurvC.col(1) + iid_survJumpT[iIndex_UTTE][iter_strata] * iEdSurvT.col(1);

			// if(isStored2_UTTE[iIndex_UTTE] && hierarchical){
			//   iidNuisance_favorable.col(iter_d) += (iIIDNuisance_favorable - iidUTTE_favorable.col(iIndex_UTTE));
			//   iidNuisance_unfavorable.col(iter_d) += iIIDNuisance_unfavorable - iidUTTE_unfavorable.col(iIndex_UTTE);
			// }else{ // first time it is analyzed
			iidNuisance_favorable.col(iter_d) += iIIDNuisance_favorable;
			iidNuisance_unfavorable.col(iter_d) += iIIDNuisance_unfavorable;
			// }
		  }
		}
	  }      
	  
      // **** update all Scores
      if(keepScore){
		// Rcout << " update lsScore" << endl;
		iNpairs = iScore.n_rows;
		iMat.resize(iNpairs,4);
		iMat.col(0).fill(iter_strata);

		for(int iPair=0; iPair < iNpairs; iPair++){
		  iMat(iPair,1) = posC[iter_strata](iScore(iPair,0));
		  iMat(iPair,2) = posT[iter_strata](iScore(iPair,1));
		  iMat(iPair,3) = iScore(iPair,0) + iScore(iPair,1)*n_control[iter_strata] + n_cumpairsM1[iter_strata];
		}
		// merge with current table and store
		if(iter_strata==0){
		  lsScore[iter_d] = arma::join_rows(iMat,iScore);
		}else{
		  lsScore[iter_d] = arma::join_cols(lsScore[iter_d], arma::join_rows(iMat,iScore));
		}
      }
    
      // **** check that there remain pairs to be analyzed
      iSize_weight = iWeight.size();
      if(iSize_weight < zeroPlus || (iMoreEndpoint==false) ){
		break;
      }

      if(hierarchical){
		// **** update score/weights associated to the remaing pairs
		// Rcout << " update matWeights" << endl;
		matWeight.resize(iSize_weight, iter_d+1);
		matWeight.col(iter_d) = iWeight; // current endpoint
		if(iter_d > 0){ // previous endpoints
		  iUvec = arma::regspace<uvec>(0, 1, iter_d - 1);
		  matWeight.cols(iUvec) = matWeight_M1.submat(iIndexWeight_pair, iUvec);
		}

		// **** update iid of the score/weigths associated to the remaing pairs
		// iLs_matIIDC, iLs_matIIDT

		  
		// **** update elements for the next step
		matWeight_M1 = matWeight; 
		iIndex_control_M1 = iIndex_control;
		iIndex_treatment_M1 = iIndex_treatment;
	  
		// **** re-initialize all vectors to 0 length
		iIndex_control.resize(0);
		iIndex_treatment.resize(0);
		iWeight.resize(0);
		iIndexWeight_pair.resize(0);
      }
	  
    }
  }
  
  // ** proportion in favor of treatment
  // Rcout << endl << " compute statistics" << endl;
  arma::mat delta_netBenefit(n_strata,D), delta_winRatio(n_strata,D); // matrix containing for each strata and each endpoint the statistic
  arma::vec Delta_netBenefit(D), Delta_winRatio(D); // vector containing for each endpoint the overall statistic

  calcStatistic(delta_netBenefit, delta_winRatio, Delta_netBenefit, Delta_winRatio,
                Mcount_favorable, Mcount_unfavorable,
				iid_favorable, iid_unfavorable, iidNuisance_favorable, iidNuisance_unfavorable,
				Mvar, returnIID,
				posC, posT, 
                D, n_strata, n_pairs, n_control, n_treatment,
				weight, hprojection, lsScore, keepScore);

  // ** export
  return(List::create(Named("count_favorable") = Mcount_favorable,
					  Named("count_unfavorable") = Mcount_unfavorable,
					  Named("count_neutral") = Mcount_neutral,           
					  Named("count_uninf") = Mcount_uninf,
					  Named("delta_netBenefit") = delta_netBenefit,
					  Named("delta_winRatio") = delta_winRatio,
					  Named("Delta_netBenefit") = conv_to< std::vector<double> >::from(Delta_netBenefit),
					  Named("Delta_winRatio") = conv_to< std::vector<double> >::from(Delta_winRatio),
					  Named("n_pairs") = conv_to< std::vector<double> >::from(n_pairs),
					  Named("iid_favorable") = iid_favorable,
					  Named("iid_unfavorable") = iid_unfavorable,
					  Named("iidNuisance_favorable") = iidNuisance_favorable,
					  Named("iidNuisance_unfavorable") = iidNuisance_unfavorable,
					  Named("Mvar") = Mvar,
					  Named("tableScore")  = lsScore
					  ));
}


