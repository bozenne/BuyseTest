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

void updateIID(arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, 
			   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, 
			   const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
			   const arma::mat& iCount_obsC, const arma::mat& iCount_obsT,
			   const vector<int>& activeUTTE, int D_activeUTTE,
			   const arma::mat& iDscore_Dnuisance_C, const arma::mat& iDscore_Dnuisance_T,
			   const std::vector< std::vector< arma::mat > >& iid_survJumpC, const std::vector< std::vector< arma::mat > >& iid_survJumpT,
			   const std::vector<std::vector< arma::mat >> & iPairDweight_Dnuisance_C,
			   const std::vector<std::vector< arma::mat >> & iPairDweight_Dnuisance_T,
			   const arma::vec& vecn_pairs, unsigned int iter_d, int iIndex_UTTE, unsigned int iter_strata, int iMethod, int returnIID);

void updatePairScore(std::vector< arma::mat >& pairScore, arma::mat& iPairScore,
					 unsigned int iter_strata, const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
					 const arma::vec& vecn_control, const arma::vec& vecn_cumpairsM1, unsigned int iter_d);

void updateRP(arma::mat& iRP_score, std::vector< arma::mat >& iRP_Dscore_Dnuisance_C, std::vector< arma::mat >& iRP_Dscore_Dnuisance_T,
			  std::vector<arma::mat>& RP_score, std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_C, std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_T,
			  arma::vec& iPairWeight_nPeron, int iSize_RP, bool neutralAsUninf, int iter_d, int correctionUninf,
			  double zeroPlus, int iIndex_UTTE, int nUTTE_analyzedPeron, int returnIID);

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
//' @param weight Store the weight associated to each endpoint. Must have length D. 
//' @param method The index of the method used to score the pairs. Must have length D. 1 for continuous, 2 for Gehan, and 3 for Peron.
//' @param D The number of endpoints.
//' @param D_UTTE The number of distinct time to event endpoints.
//' @param n_strata The number of strata. 
//' @param nUTTE_analyzedPeron_M1 The number of unique time-to-event endpoints that have been analyzed the Peron scoring rule before the current endpoint. Must have length D.
//' @param index_endpoint The position of the endpoint at each priority in the argument endpoint. Must have length D. 
//' @param index_censoring The position of the censoring at each priority in the argument censoring. Must have length D. 
//' @param index_UTTE The position, among all the unique tte endpoints, of the TTE endpoints. Equals -1 for non tte endpoints. Must have length n_TTE. 
//' @param list_survTimeC A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows).
//' @param list_survTimeT A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows).
//' @param list_survJumpC A list of matrix containing the survival estimates and survival jumps when the survival for the control arm jumps.
//' @param list_survJumpT A list of matrix containing the survival estimates and survival jumps when the survival for the treatment arm jumps.
//' @param list_lastSurv A list of matrix containing the last survival estimate in each strata (rows) and treatment group (columns).
//' @param p_C Number of nuisance parameter in the survival model for the control group, for each endpoint and strata
//' @param p_T Number of nuisance parameter in the survival model for the treatment group, for each endpoint and strata
//' @param iid_survJumpC A list of matrix containing the iid of the survival estimates in the control group.
//' @param iid_survJumpT A list of matrix containing the iid of the survival estimates in the treatment group.
//' @param zeroPlus Value under which doubles are considered 0?
//' @param correctionUninf Should the uninformative weight be re-distributed to favorable and unfavorable?
//' @param hierarchical Should only the uninformative pairs be analyzed at the lower priority endpoints (hierarchical GPC)? Otherwise all pairs will be compaired for all endpoint (full GPC).
//' @param hprojection Order of the H-projection used to compute the variance.
//' @param neutralAsUninf Should paired classified as neutral be re-analyzed using endpoints of lower priority? 
//' @param keepScore Should the result of each pairwise comparison be kept?
//' @param returnIID Should the iid be computed?
//' @param debug Print messages tracing the execution of the function to help debugging. The amount of messages increase with the value of debug (0-5).
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
			 arma::vec weight,
			 arma::vec method,
			 unsigned int D,
			 unsigned int D_UTTE,
			 unsigned int n_strata,
			 arma::vec nUTTE_analyzedPeron_M1,
			 std::vector<unsigned int> index_endpoint, 
			 std::vector<unsigned int> index_censoring, 
			 std::vector<int> index_UTTE, 
			 std::vector< std::vector< arma::mat > > list_survTimeC,
			 std::vector< std::vector< arma::mat > > list_survTimeT,
			 std::vector< std::vector< arma::mat > > list_survJumpC,
			 std::vector< std::vector< arma::mat > > list_survJumpT,
			 std::vector< arma::mat > list_lastSurv,
			 arma::mat p_C,
			 arma::mat p_T,
			 std::vector< std::vector< arma::mat > > iid_survJumpC,
			 std::vector< std::vector< arma::mat > > iid_survJumpT,
			 double zeroPlus,
			 int correctionUninf,
			 bool hierarchical,
			 int hprojection,
			 bool neutralAsUninf,
			 bool keepScore,
			 int returnIID,
			 int debug){

  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following endpoints.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.

  
  /// ** initialization
  int n_obs = endpoint.n_rows;
  
  // *** objects storing the final results
  // score specific to each pair
  std::vector< arma::mat > pairScore(D);
  arma::mat iPairScore; // same but only for the current endpoint
   
  // total score over pairs
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs [all endpoints, strata]
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs [all endpoints, strata]
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs [all endpoints, strata]
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs [all endpoints, strata]

  // number of pairs
  arma::vec vecn_pairs(n_strata); // number of pairs sumed over the strats
  arma::vec vecn_treatment(n_strata); // number of patients in the treatment group over the strats
  arma::vec vecn_control(n_strata); // number of patients in the control group over the strats
  arma::vec vecn_cumpairsM1(n_strata); // number of pairs in the previous strata (used when storing all the pairs in pairScore)
  
  // variance and iid
  arma::mat Mvar; // variance-covariance (favorable,unfavorable) scores  [all endpoints]

  arma::mat iidAverage_favorable; // iid relative to the average over all pairs for the favorable scores [all endpoints]
  arma::mat iidAverage_unfavorable; // iid relative to the average over all pairs for the unfavorable scores [all endpoints]
  arma::mat iCount_obsC; // iidAverage [current endpoint and strata]
  arma::mat iCount_obsT; // iidAverage [current endpoint and strata]

  arma::mat iidNuisance_favorable; // iid relative to the nuisance parameters for the favorable scores [all endpoints]
  arma::mat iidNuisance_unfavorable; // iid relative to the nuisance parameters for the unfavorable scores [all endpoints]
  arma::mat iDscore_Dnuisance_C; // partial derivative regarding nuisance parameters used by iidAverage [current endpoint and strata].
  arma::mat iDscore_Dnuisance_T; // partial derivative regarding nuisance parameters used by iidAverage [current endpoint and strata]. 


  if(returnIID>0){
	Mvar.resize(D,5); // variance(favorable); variance(unfavorable); covariance(favorable,unfavorable); variance(netBenefit); variance(winRatio);
    Mvar.fill(0.0);

	// iid with respect to the averaging over pairs
    iidAverage_favorable.resize(n_obs,D);
	iidAverage_favorable.fill(0.0);
    iidAverage_unfavorable.resize(n_obs,D);
	iidAverage_unfavorable.fill(0.0);

	if(returnIID>1){  
	  // iid with respect to the nuisance parameters
	  iidNuisance_favorable.resize(n_obs,D);
	  iidNuisance_favorable.fill(0.0);
	  iidNuisance_unfavorable.resize(n_obs,D);
	  iidNuisance_unfavorable.fill(0.0);
	}
  }
  
  // ** loop over strata
  for(unsigned int iter_strata=0 ; iter_strata < n_strata ; iter_strata++){

	// *** initialization
	int iSize_RP = 0; // number of residual pairs (initialize to avoid C++warnings). 
	vector<int> activeUTTE(0); // index of the distinct TTE endpoints already analyzed with Peron scoring rule.
	int D_activeUTTE = 0; // number of distinct TTE endpoints already analyzed with Peron scoring rule.
	arma::vec iIndex_control; // position of the controls in the remaining pairs among all controls.
	arma::vec iIndex_treatment; // position of the treated in the remaining pairs among all treated.
	arma::vec iPairWeight_nPeron; // weight associated with each pair for endpoints not analyzed with Peron

	// score of the remaining pairs for the TTE endpoints at the highest previous threshold
	// remaining pairs = pairs with non-0 neutral or uninformative score
	std::vector<arma::mat> RP_score(D_UTTE); // mat(pair;favorable/unfavorable/neutral/uninformative)
	arma::mat iRP_score; // same but only for the current endpoint

	// iid(score) of the remaining pairs for the TTE endpoints at the highest previous threshold
	std::vector< std::vector< arma::mat > > RP_Dscore_Dnuisance_C(D_UTTE); // endpoint favorable/unfavorable/neutral (pair,nuisance parameters)
	std::vector< std::vector< arma::mat > > RP_Dscore_Dnuisance_T(D_UTTE); // endpoint favorable/unfavorable/neutral (pair,nuisance parameters)
	for(unsigned int iter_UTTE=0; iter_UTTE<D_UTTE; iter_UTTE++){
	  RP_Dscore_Dnuisance_C[iter_UTTE].resize(4);
	  RP_Dscore_Dnuisance_T[iter_UTTE].resize(4);
	}
	std::vector< arma::mat > iRP_Dscore_Dnuisance_C(4); //  favorable/unfavorable/neutral/uninformative (pair,nuisance parameters) [current endpoint and strata]
	std::vector< arma::mat > iRP_Dscore_Dnuisance_T(4); //  favorable/unfavorable/neutral/uninformative (pair,nuisance parameters) [current endpoint and strata]

    for(unsigned int iter_d=0 ; iter_d < D; iter_d++){
      if(debug>0){Rcout << endl << "** endpoint " << iter_d << "**" << endl;}

      // *** type of endpoint
      int iMethod = method(iter_d); // how the score are computed: 1: continuous, 2: Gehan or 3: Peron
      int iIndex_UTTE = index_UTTE[iter_d]; // which of the distinct TTE is the current TTE endpoint (-100 if not TTE or if method=Gehan)
	  int iNUTTE_analyzedPeron = nUTTE_analyzedPeron_M1[iter_d]; // number of distinct TTE endpoints already analyzed with Peron scoring rule

	  bool iFirstEndpoint = (iter_d==0) || (hierarchical==false); // should the GPC be performed on all possible pairs or only remaining pairs, accounting for weights.
      bool iMoreEndpoint = (D>(iter_d+1)) && hierarchical; // is there any endpoint after this one (if not, do not store detailed information about the residual pairs)
	  bool iAlreadyAnalyzed = (iIndex_UTTE < iNUTTE_analyzedPeron) && (iIndex_UTTE > (-zeroPlus)); // has the current endpoint already been analyzed using the Peron's scoring rule
	  bool iUpdateIndexNeutral = iMoreEndpoint & neutralAsUninf;
	  bool iUpdateIndexUninf = iMoreEndpoint & (correctionUninf == 0 || (neutralAsUninf && (correctionUninf == 1)));
	  
      // *** compute weights with their iid decomposition
	  arma::vec iPairWeight = iPairWeight_nPeron; // weight associated with each pair over all previous endpoints	
	  std::vector<std::vector< arma::mat >> iPairDweight_Dnuisance_C(4);  
	  std::vector<std::vector< arma::mat >> iPairDweight_Dnuisance_T(4); 
	  std::vector<std::vector< arma::mat >> iPairDweight_Dnuisance_C_M1;  
	  std::vector<std::vector< arma::mat >> iPairDweight_Dnuisance_T_M1; 
	  if(iNUTTE_analyzedPeron > zeroPlus){
   	    if(debug>0){Rcout << " - compute weights("<< iNUTTE_analyzedPeron <<"): " << endl;}
		prepareWeight(iPairWeight, iPairDweight_Dnuisance_C, iPairDweight_Dnuisance_T,
					  activeUTTE, D_activeUTTE,
					  iter_d, iIndex_UTTE, RP_score,
					  RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
					  iNUTTE_analyzedPeron, correctionUninf, zeroPlus, neutralAsUninf, returnIID);
		if(iAlreadyAnalyzed){
		  iPairDweight_Dnuisance_C_M1 = iPairDweight_Dnuisance_C;
		  iPairDweight_Dnuisance_T_M1 = iPairDweight_Dnuisance_T;
		}
	  }
	  // Rcout << iPairWeight << endl;

	  // *** compute scores
	  if(debug>0){Rcout << " - score("<< iFirstEndpoint <<")" << endl;}
      arma::uvec iUvec_endpoint = {index_endpoint[iter_d]};
      arma::uvec iUvec_censoring = {index_censoring[iter_d]};

	  iPairScore = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
								censoring.submat(indexC[iter_strata],iUvec_censoring), censoring.submat(indexT[iter_strata],iUvec_censoring),
								list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata], list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
								list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1), 
								iIndex_control, iIndex_treatment, iPairWeight,
								activeUTTE, D_activeUTTE,
								Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d),
								iRP_score,
								iCount_obsC, iCount_obsT, iDscore_Dnuisance_C, iDscore_Dnuisance_T,
								iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T,
								iPairDweight_Dnuisance_C, iPairDweight_Dnuisance_T,
								zeroPlus, 
								iMethod, returnIID, p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),
								iFirstEndpoint, false, iUpdateIndexNeutral, iUpdateIndexUninf, keepScore, correctionUninf, neutralAsUninf,
								debug);

	  if(iAlreadyAnalyzed){ // substract contribution of the previous analysis
		if(debug>0){Rcout << " - scoreM1("<< iIndex_UTTE <<")" << endl;}
		arma::mat iRP_score_M1= RP_score[iIndex_UTTE]; // scores relative to latest analysis of the same endpoint
		std::vector< arma::mat > iRP_Dscore_Dnuisance_C_M1 = RP_Dscore_Dnuisance_C[iIndex_UTTE]; // iid of the scores relative to latest analysis of the same endpoint
		std::vector< arma::mat > iRP_Dscore_Dnuisance_T_M1 = RP_Dscore_Dnuisance_T[iIndex_UTTE];

		// re-compute using current weights
		double iCount_favorable_M1,iCount_unfavorable_M1,iCount_neutral_M1,iCount_uninf_M1; // initialization necessary
		arma::mat iCount_obsC_M1,iCount_obsT_M1;
		arma::mat iDscore_Dnuisance_C_M1,iDscore_Dnuisance_T_M1;

		// note the values in endpoint, censoring, survTime, survJump, lastSurv are not used (only their dimensions)
		arma::mat iPairScore_M1 = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
									   censoring.submat(indexC[iter_strata],iUvec_censoring), censoring.submat(indexT[iter_strata],iUvec_censoring),
									   list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata], list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
									   list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1), 
									   iIndex_control, iIndex_treatment, iPairWeight,
									   activeUTTE, D_activeUTTE,
									   iCount_favorable_M1, iCount_unfavorable_M1, iCount_neutral_M1, iCount_uninf_M1,
									   iRP_score_M1,
									   iCount_obsC_M1, iCount_obsT_M1, iDscore_Dnuisance_C_M1, iDscore_Dnuisance_T_M1,
									   iRP_Dscore_Dnuisance_C_M1, iRP_Dscore_Dnuisance_T_M1,
									   iPairDweight_Dnuisance_C_M1, iPairDweight_Dnuisance_T_M1,								
									   zeroPlus, 
									   iMethod, returnIID, p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),
									   false, true, false, false, keepScore, correctionUninf, neutralAsUninf,
									   debug);

		Mcount_favorable(iter_strata,iter_d) -= iCount_favorable_M1;
		Mcount_unfavorable(iter_strata,iter_d) -= iCount_unfavorable_M1;
		if(keepScore){
		  iPairScore.col(2) -= iPairScore_M1.col(2); // favorable
		  iPairScore.col(3) -= iPairScore_M1.col(3); // unfavorable corrected
		  iPairScore.col(7) -= iPairScore_M1.col(7); // favorable
		  iPairScore.col(8) -= iPairScore_M1.col(8); // unfavorable corrected
		}
		if(returnIID>0){
		  for(int iCol=0; iCol <2; iCol++){
			iCount_obsC.col(iCol) -= iCount_obsC_M1.col(iCol);
			iCount_obsT.col(iCol) -= iCount_obsT_M1.col(iCol);
			if(returnIID>1){
			  iDscore_Dnuisance_C.col(iCol) -= iDscore_Dnuisance_C_M1.col(iCol);
			  iDscore_Dnuisance_T.col(iCol) -= iDscore_Dnuisance_T_M1.col(iCol);
			  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
				iPairDweight_Dnuisance_C[iCol][activeUTTE[iter_UTTE]] -= iPairDweight_Dnuisance_C_M1[iCol][activeUTTE[iter_UTTE]];
				iPairDweight_Dnuisance_T[iCol][activeUTTE[iter_UTTE]] -= iPairDweight_Dnuisance_T_M1[iCol][activeUTTE[iter_UTTE]];
			  }
			}
		  }
		}
	  }

	  // *** update number of pairs
	  if(iter_d==0){
		if(debug>0){Rcout << " update number of pairs" << endl;}
		vecn_control[iter_strata] = posC[iter_strata].size();
		vecn_treatment[iter_strata] = posT[iter_strata].size();		  
		vecn_pairs[iter_strata] = vecn_control[iter_strata] * vecn_treatment[iter_strata];
		//= Mcount_favorable(iter_strata,0) + Mcount_unfavorable(iter_strata,0) + Mcount_neutral(iter_strata,0) + Mcount_uninf(iter_strata,0);
		if(iter_strata == 0){
		  vecn_cumpairsM1[0] = 0;
		}else{
		  vecn_cumpairsM1[iter_strata] = vecn_cumpairsM1[iter_strata-1] + vecn_control[iter_strata]*vecn_treatment[iter_strata];
		}
	  }

	  // *** update iid
	  if(returnIID>0){
		if(debug>0){Rcout << " update iid (" << returnIID << ")" << endl;}
		updateIID(iidAverage_favorable, iidAverage_unfavorable,
				  iidNuisance_favorable, iidNuisance_unfavorable, 
				  posC, posT,
				  iCount_obsC, iCount_obsT,
				  activeUTTE, D_activeUTTE,
				  iDscore_Dnuisance_C, iDscore_Dnuisance_T,
				  iid_survJumpC, iid_survJumpT,
				  iPairDweight_Dnuisance_C, iPairDweight_Dnuisance_T,
				  vecn_pairs, iter_d, iIndex_UTTE, iter_strata, iMethod, returnIID);
	  } 
	  
	  // *** update pairwise-scores (all pairs)
	  if(keepScore){ // store iPaireScore in pairScore
		if(debug>0){Rcout << " update pairwise scores " << endl;}
		updatePairScore(pairScore, iPairScore,
						iter_strata, posC, posT,
						vecn_control, vecn_cumpairsM1, iter_d);
      }
    
      // *** store scores (and iid) relative to the remaining pairs
      if(iMoreEndpoint){
		// end if no remaining pairs to be analyzed
		iSize_RP = iRP_score.n_rows;
		if(iSize_RP < zeroPlus){break;}
		
		// update position of the remaining pairs among the controls / treated
		iIndex_control = iRP_score.col(1);
		iIndex_treatment = iRP_score.col(2);

		// update iPairWeight_nPeron, RP_score, RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
		// and re-initialize iRP_score, iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T
		if(debug>0){Rcout << " update score/iid for the remaing pairs("<< nUTTE_analyzedPeron_M1[iter_d+1] <<") " << endl;}
		updateRP(iRP_score, iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T,
				 RP_score, RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
				 iPairWeight_nPeron, iSize_RP, neutralAsUninf, iter_d, correctionUninf,
				 zeroPlus, iIndex_UTTE, nUTTE_analyzedPeron_M1[iter_d+1], returnIID);
	  }
  
	  
	} // end endpoint
  } // end strata
  
  // ** proportion in favor of treatment
  // Rcout << endl << " compute statistics" << endl;
  arma::mat delta_netBenefit(n_strata,D), delta_winRatio(n_strata,D); // matrix containing for each strata and each endpoint the statistic
  arma::vec Delta_netBenefit(D), Delta_winRatio(D); // vector containing for each endpoint the overall statistic

  calcStatistic(delta_netBenefit, delta_winRatio, Delta_netBenefit, Delta_winRatio,
                Mcount_favorable, Mcount_unfavorable,
				iidAverage_favorable, iidAverage_unfavorable, iidNuisance_favorable, iidNuisance_unfavorable,
				Mvar, returnIID,
				posC, posT, 
                D, n_strata, vecn_pairs, vecn_control, vecn_treatment,
				weight, hprojection, pairScore, keepScore);

  // ** export
  return(List::create(Named("count_favorable") = Mcount_favorable,
					  Named("count_unfavorable") = Mcount_unfavorable,
					  Named("count_neutral") = Mcount_neutral,           
					  Named("count_uninf") = Mcount_uninf,
					  Named("delta_netBenefit") = delta_netBenefit,
					  Named("delta_winRatio") = delta_winRatio,
					  Named("Delta_netBenefit") = conv_to< std::vector<double> >::from(Delta_netBenefit),
					  Named("Delta_winRatio") = conv_to< std::vector<double> >::from(Delta_winRatio),
					  Named("n_pairs") = conv_to< std::vector<double> >::from(vecn_pairs),
					  Named("iidAverage_favorable") = iidAverage_favorable,
					  Named("iidAverage_unfavorable") = iidAverage_unfavorable,
					  Named("iidNuisance_favorable") = iidNuisance_favorable,
					  Named("iidNuisance_unfavorable") = iidNuisance_unfavorable,
					  Named("Mvar") = Mvar,
					  Named("tableScore")  = pairScore
					  ));
}

// * updateIID
void updateIID(arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, 
			   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, 
			   const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
			   const arma::mat& iCount_obsC, const arma::mat& iCount_obsT,
			   const vector<int>& activeUTTE, int D_activeUTTE,
			   const arma::mat& iDscore_Dnuisance_C, const arma::mat& iDscore_Dnuisance_T,
			   const std::vector< std::vector< arma::mat > >& iid_survJumpC, const std::vector< std::vector< arma::mat > >& iid_survJumpT,
			   const std::vector<std::vector< arma::mat >> & iPairDweight_Dnuisance_C,
			   const std::vector<std::vector< arma::mat >> & iPairDweight_Dnuisance_T,
			   const arma::vec& vecn_pairs, unsigned int iter_d, int iIndex_UTTE, unsigned int iter_strata, int iMethod, int returnIID){

  // ** iid with respect to the average over all pairs
  arma::uvec iUvec_iter_d = {iter_d};
  iidAverage_favorable.submat(posC[iter_strata], iUvec_iter_d) = iCount_obsC.col(0);
  iidAverage_favorable.submat(posT[iter_strata], iUvec_iter_d) = iCount_obsT.col(0);
	
  iidAverage_unfavorable.submat(posC[iter_strata], iUvec_iter_d) = iCount_obsC.col(1);
  iidAverage_unfavorable.submat(posT[iter_strata], iUvec_iter_d) = iCount_obsT.col(1);

  // ** iid with respect to the nuisance parameters
  if(returnIID>1){

	// *** iid of the weights
	arma::rowvec iDweight_Dnuisance_C;
	arma::rowvec iDweight_Dnuisance_T;
	
	for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
	  iDweight_Dnuisance_C = sum(iPairDweight_Dnuisance_C[0][activeUTTE[iter_UTTE]],0);
	  iidNuisance_favorable.col(iter_d) += iid_survJumpC[activeUTTE[iter_UTTE]][iter_strata] * arma::trans(iDweight_Dnuisance_C) / vecn_pairs[iter_strata];
	  iDweight_Dnuisance_T = sum(iPairDweight_Dnuisance_T[0][activeUTTE[iter_UTTE]],0);
	  iidNuisance_favorable.col(iter_d) += iid_survJumpT[activeUTTE[iter_UTTE]][iter_strata] * arma::trans(iDweight_Dnuisance_T) / vecn_pairs[iter_strata];
	  
	  iDweight_Dnuisance_C = sum(iPairDweight_Dnuisance_C[1][activeUTTE[iter_UTTE]],0);
	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[activeUTTE[iter_UTTE]][iter_strata] * arma::trans(iDweight_Dnuisance_C) / vecn_pairs[iter_strata];
	  iDweight_Dnuisance_T = sum(iPairDweight_Dnuisance_T[1][activeUTTE[iter_UTTE]],0);
	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[activeUTTE[iter_UTTE]][iter_strata] * arma::trans(iDweight_Dnuisance_T) / vecn_pairs[iter_strata];
	}

	// *** iid of the proba/score
	if(iMethod == 3){
	  iidNuisance_favorable.col(iter_d) += iid_survJumpC[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_C.col(0)/vecn_pairs[iter_strata];
	  iidNuisance_favorable.col(iter_d) += iid_survJumpT[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_T.col(0)/vecn_pairs[iter_strata];
	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_C.col(1)/vecn_pairs[iter_strata];
	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_T.col(1)/vecn_pairs[iter_strata];
	}

  }

  return;
}


// * updatePairScore
void updatePairScore(std::vector< arma::mat >& pairScore, arma::mat& iPairScore,
					 unsigned int iter_strata, const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
					 const arma::vec& vecn_control, const arma::vec& vecn_cumpairsM1, unsigned int iter_d){

  // ** prepare additional columns for indicating position and strata
  int iNpairs = iPairScore.n_rows;
  arma::mat iMat(iNpairs,4);

  // first column contain strata indicator
  iMat.col(0).fill(iter_strata);

  // following columns contain position relative to control obs, treatment obs, and all pairs
  for(int iPair=0; iPair < iNpairs; iPair++){
	iMat(iPair,1) = posC[iter_strata](iPairScore(iPair,0));
	iMat(iPair,2) = posT[iter_strata](iPairScore(iPair,1));
	iMat(iPair,3) = iPairScore(iPair,0) + iPairScore(iPair,1)*vecn_control[iter_strata] + vecn_cumpairsM1[iter_strata];
  }
  // ** merge additional column with the current table
  if(iter_strata==0){
	pairScore[iter_d] = arma::join_rows(iMat,iPairScore);
  }else{
	pairScore[iter_d] = arma::join_cols(pairScore[iter_d], arma::join_rows(iMat,iPairScore));
  }

  return;
}

// * updateRP
void updateRP(arma::mat& iRP_score, std::vector< arma::mat >& iRP_Dscore_Dnuisance_C, std::vector< arma::mat >& iRP_Dscore_Dnuisance_T,
			  std::vector<arma::mat>& RP_score, std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_C, std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_T,
			  arma::vec& iPairWeight_nPeron, int iSize_RP, bool neutralAsUninf, int iter_d, int correctionUninf,
			  double zeroPlus, int iIndex_UTTE, int nUTTE_analyzedPeron, int returnIID){

  // ** index of the remaining pairs among the analyzed pairs
  arma::uvec iIndex_RP = conv_to<uvec>::from(iRP_score.col(0));

  // ** update/subset score and iid(score) for TTE endpoint analyzed with Peron's scoring rule
  for(int iter_UTTE=0; iter_UTTE<nUTTE_analyzedPeron; iter_UTTE++){
	if(iter_UTTE == iIndex_UTTE){
	  arma::uvec iter_36 = linspace<arma::uvec>(3, 6, 4);
	  RP_score[iter_UTTE] = iRP_score.cols(iter_36);
	  if(returnIID>1){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
		  RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP] = iRP_Dscore_Dnuisance_C[iter_typeRP];
		  RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP] = iRP_Dscore_Dnuisance_T[iter_typeRP];
		}
	  }
	}else{
	  RP_score[iter_UTTE] = RP_score[iter_UTTE].rows(iIndex_RP);
	  if(returnIID>1){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){ 
		  RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP] = RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP].rows(iIndex_RP);
		  RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP] = RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP].rows(iIndex_RP);
		}
	  }
	}
			
  } // end UTTE
	
  // ** update/subset weights of other endpoints
  if(correctionUninf>0){
	if(iter_d>0){
	  iPairWeight_nPeron = iPairWeight_nPeron.rows(iIndex_RP);
	}else{
	  iPairWeight_nPeron.resize(iSize_RP);
	  iPairWeight_nPeron.fill(1.0);
	}
	
	if(iIndex_UTTE < (-zeroPlus)){ // i.e. not analyzed by Peron
	  if(neutralAsUninf){
		iPairWeight_nPeron %= (iRP_score.col(5) + iRP_score.col(6));
	  }else{
		iPairWeight_nPeron %= iRP_score.col(6);
	  }	  
	}
  }else{
	iPairWeight_nPeron.resize(iSize_RP);
	iPairWeight_nPeron.fill(1.0);
  }

  return ;
}
