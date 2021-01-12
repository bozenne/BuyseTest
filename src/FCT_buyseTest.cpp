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

void prepareWeight(arma::vec& iPairWeight, std::vector<std::vector< arma::sp_mat >>& iPairDweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& iPairDweight_Dnuisance_T,
		   std::vector<int>& activeUTTE, int& D_activeUTTE,
		   int iter_d, int iIndex_UTTE, const std::vector<arma::mat>& RP_score,
		   const std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_C, const std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_T,
		   int iNUTTE_analyzedPeron, int correctionUninf, double zeroPlus, std::vector<bool>& neutralAsUninf, int returnIID);

void updateIID(arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, 
	       arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, 
	       const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
	       const arma::mat& iCount_obsC, const arma::mat& iCount_obsT,
	       const std::vector<int>& activeUTTE, int D_activeUTTE,
	       const arma::mat& iDscore_Dnuisance_C, const arma::mat& iDscore_Dnuisance_T,
	       const std::vector< std::vector< arma::mat > >& iid_survJumpC, const std::vector< std::vector< arma::mat > >& iid_survJumpT,
	       const std::vector<std::vector< arma::sp_mat >> & iPairDweight_Dnuisance_C,
	       const std::vector<std::vector< arma::sp_mat >> & iPairDweight_Dnuisance_T,
	       const arma::vec& vecn_pairs, unsigned int iter_d, int iIndex_UTTE, unsigned int iter_strata, int iMethod, int returnIID);

void updatePairScore(std::vector< arma::mat >& pairScore, arma::mat& iPairScore,
		     unsigned int iter_strata, const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
		     const arma::vec& vecn_control, const arma::vec& vecn_cumpairsM1, unsigned int iter_d);

void updateRP(arma::mat& iRP_score, std::vector< arma::sp_mat >& iRP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& iRP_Dscore_Dnuisance_T,
	      std::vector<arma::mat>& RP_score, std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_C, std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_T,
	      arma::vec& iPairWeight_nPeron, int iSize_RP, bool neutralAsUninf, int iter_d, int correctionUninf,
	      double zeroPlus, int iIndex_UTTE, int nUTTE_analyzedPeron, int returnIID);

arma::sp_mat subcol_sp_mat(const arma::sp_mat& X, arma::uvec index);

// * Documentation GPC_cpp
//' @title C++ function performing the pairwise comparison over several endpoints. 
//' @description \code{GPC_cpp} call for each endpoint and each strata the pairwise comparison function suited to the type of endpoint and store the results. 
//' @name GPC_cpp
//' 
//' @param endpoint A matrix containing the values of each endpoint (in columns) for each observation (in rows). 
//' @param status A matrix containing the values of the status variables relative to each endpoint (in columns) for each observation (in rows).
//' @param indexC A list containing, for each strata, which rows of the endpoint and status matrices corresponds to the control observations. Not unique when bootstraping.
//' @param posC A list containing, for each strata, the unique identifier of each control observations. 
//' @param indexT A list containing, for each strata, which rows of the endpoint and status matrices corresponds to the treatment observations. Not unique when bootstraping.
//' @param posT A list containing, for each strata, the unique identifier of each treatment observations.
//' @param threshold Store the thresholds associated to each endpoint. Must have length D. The threshold is ignored for binary endpoints. 
//' @param weight Store the weight associated to each endpoint. Must have length D. 
//' @param method The index of the method used to score the pairs. Must have length D. 1 for continuous, 2 for Gehan, and 3 for Peron.
//' @param op The index of the operator used to score the pairs. Must have length D. 1 for larger is beter, -1 for smaller is better.
//' @param D The number of endpoints.
//' @param D_UTTE The number of distinct time to event endpoints.
//' @param n_strata The number of strata. 
//' @param nUTTE_analyzedPeron_M1 The number of unique time-to-event endpoints that have been analyzed the Peron scoring rule before the current endpoint. Must have length D.
//' @param index_endpoint The position of the endpoint at each priority in the argument endpoint. Must have length D. 
//' @param index_status The position of the status at each priority in the argument status. Must have length D. 
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
//' @param precompute Have the integrals relative to the survival be already computed and stored in list_survTimeC/list_survTimeT and list_survJumpC/list_survJumpT (derivatives)
//' @param returnIID Should the iid be computed?
//' @param debug Print messages tracing the execution of the function to help debugging. The amount of messages increase with the value of debug (0-5).
//'
//' @details GPC_cpp implements GPC looping first over endpoints and then over pairs.
//' To handle multiple endpoints, it stores some of the results which can be memory demanding when considering large sample - especially when computing the iid decomposition.
//' GPC2_cpp implements GPC looping first over pairs and then over endpoints. It has rather minimal memory requirement but does not handle correction for uninformative pairs. 
//'
//' @keywords function Cpp BuyseTest
//'
//' @author Brice Ozenne

// * Function GPC_cpp
//' @name GPC_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List GPC_cpp(arma::mat endpoint,
		   arma::mat status,
		   std::vector< arma::uvec > indexC,
		   std::vector< arma::uvec > posC,
		   std::vector< arma::uvec > indexT,
		   std::vector< arma::uvec > posT,
		   std::vector< double > threshold,
		   arma::vec weight,
		   arma::vec method,
		   std::vector< int > op,
		   unsigned int D,
		   unsigned int D_UTTE,
		   unsigned int n_strata,
		   arma::vec nUTTE_analyzedPeron_M1,
		   std::vector<unsigned int> index_endpoint, 
		   std::vector<unsigned int> index_status, 
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
		   std::vector<bool> neutralAsUninf,
		   bool keepScore,
		   bool precompute,
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
  arma::mat Mcount_favorable(n_strata,D,arma::fill::zeros); // store the total weight of favorable pairs [all endpoints, strata]
  arma::mat Mcount_unfavorable(n_strata,D,arma::fill::zeros); // store the total weight of unfavorable pairs [all endpoints, strata]
  arma::mat Mcount_neutral(n_strata,D,arma::fill::zeros); // store the total weight of neutral pairs [all endpoints, strata]
  arma::mat Mcount_uninf(n_strata,D,arma::fill::zeros); // store the total weight of uninf pairs [all endpoints, strata]

  // number of pairs
  arma::vec vecn_pairs(n_strata); // number of pairs sumed over the strats
  arma::vec vecn_treatment(n_strata); // number of patients in the treatment group over the strats
  arma::vec vecn_control(n_strata); // number of patients in the control group over the strats
  arma::vec vecn_cumpairsM1(n_strata); // number of pairs in the previous strata (used when storing all the pairs in pairScore)
  
  // iid
  arma::mat iidAverage_favorable; // iid relative to the average over all pairs for the favorable scores [all endpoints]
  arma::mat iidAverage_unfavorable; // iid relative to the average over all pairs for the unfavorable scores [all endpoints]
  arma::mat iCount_obsC; // iidAverage [current endpoint and strata]
  arma::mat iCount_obsT; // iidAverage [current endpoint and strata]

  arma::mat iidNuisance_favorable; // iid relative to the nuisance parameters for the favorable scores [all endpoints]
  arma::mat iidNuisance_unfavorable; // iid relative to the nuisance parameters for the unfavorable scores [all endpoints]
  arma::mat iDscore_Dnuisance_C; // partial derivative regarding nuisance parameters used by iidAverage [current endpoint and strata].
  arma::mat iDscore_Dnuisance_T; // partial derivative regarding nuisance parameters used by iidAverage [current endpoint and strata]. 


  if(returnIID>0){
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
    std::vector<int> activeUTTE(0); // index of the distinct TTE endpoints already analyzed with Peron scoring rule.
    std::vector<bool> neutralAsUninfUTTE(D_UTTE); // neutral as uninf option for each of the distinct TTE endpoints already analyzed with Peron scoring rule.
    int D_activeUTTE = 0; // number of distinct TTE endpoints already analyzed with Peron scoring rule.
    arma::vec iIndex_control; // position of the controls in the remaining pairs among all controls.
    arma::vec iIndex_treatment; // position of the treated in the remaining pairs among all treated.
    arma::vec iPairWeight_nPeron; // weight associated with each pair for endpoints not analyzed with Peron

    // score of the remaining pairs for the TTE endpoints at the highest previous threshold
    // remaining pairs = pairs with non-0 neutral or uninformative score
    std::vector<arma::mat> RP_score(D_UTTE); // mat(pair;favorable/unfavorable/neutral/uninformative)
    arma::mat iRP_score; // same but only for the current endpoint

    // iid(score) of the remaining pairs for the TTE endpoints at the smallest previous threshold
    std::vector< std::vector< arma::sp_mat > > RP_Dscore_Dnuisance_C(D_UTTE); // endpoint favorable/unfavorable/neutral (pair,nuisance parameters)
    std::vector< std::vector< arma::sp_mat > > RP_Dscore_Dnuisance_T(D_UTTE); // endpoint favorable/unfavorable/neutral (pair,nuisance parameters)
    std::vector< arma::sp_mat > iRP_Dscore_Dnuisance_C; //  favorable/unfavorable/neutral/uninformative (pair,nuisance parameters) [current endpoint
    std::vector< arma::sp_mat > iRP_Dscore_Dnuisance_T; //  favorable/unfavorable/neutral/uninformative (pair,nuisance parameters) [current endpoint]
    if(returnIID>1){
      for(unsigned int iter_UTTE=0; iter_UTTE<D_UTTE; iter_UTTE++){
	RP_Dscore_Dnuisance_C[iter_UTTE].resize(4);
	RP_Dscore_Dnuisance_T[iter_UTTE].resize(4);
	iRP_Dscore_Dnuisance_C.resize(4);
	iRP_Dscore_Dnuisance_T.resize(4);
      }
    }
	
    for(unsigned int iter_d=0 ; iter_d < D; iter_d++){
      if(debug>0){Rcpp::Rcout << std::endl << "** endpoint " << iter_d << "**" << std::endl;}

      // *** type of endpoint
      int iMethod = method(iter_d); // how the score are computed: 1: continuous, 2: Gehan or 3: Peron
      int iIndex_UTTE = index_UTTE[iter_d]; // which of the distinct TTE is the current TTE endpoint (-100 if not TTE or if method=Gehan)
      int iNUTTE_analyzedPeron = nUTTE_analyzedPeron_M1[iter_d]; // number of distinct TTE endpoints already analyzed with Peron scoring rule

      bool iFirstEndpoint = (iter_d==0) || (hierarchical==false); // should the GPC be performed on all possible pairs or only remaining pairs, accounting for weights.
      bool iMoreEndpoint = (D>(iter_d+1)) && hierarchical; // is there any endpoint after this one (if not, do not store detailed information about the residual pairs)
      bool iAlreadyAnalyzed = (iIndex_UTTE < iNUTTE_analyzedPeron) && (iIndex_UTTE > (-zeroPlus)); // has the current endpoint already been analyzed using the Peron's scoring rule
      bool iUpdateIndexNeutral = iMoreEndpoint & neutralAsUninf[iter_d];
      bool iUpdateIndexUninf = iMoreEndpoint & (correctionUninf == 0 || (neutralAsUninf[iter_d] && (correctionUninf == 1)));
	  
      // *** compute weights with their iid decomposition
      arma::vec iPairWeight = iPairWeight_nPeron; // weight associated with each pair over all previous endpoints	
      std::vector<std::vector< arma::sp_mat >> iPairDweight_Dnuisance_C(4);  
      std::vector<std::vector< arma::sp_mat >> iPairDweight_Dnuisance_T(4); 
      std::vector<std::vector< arma::sp_mat >> iPairDweight_Dnuisance_C_M1;  
      std::vector<std::vector< arma::sp_mat >> iPairDweight_Dnuisance_T_M1; 
      if(iNUTTE_analyzedPeron > zeroPlus){
	if(debug>0){Rcpp::Rcout << " - compute weights("<< iNUTTE_analyzedPeron <<"): " << std::endl;}
	prepareWeight(iPairWeight, iPairDweight_Dnuisance_C, iPairDweight_Dnuisance_T,
		      activeUTTE, D_activeUTTE,
		      iter_d, iIndex_UTTE, RP_score,
		      RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
		      iNUTTE_analyzedPeron, correctionUninf, zeroPlus, neutralAsUninfUTTE, returnIID);
	if(iAlreadyAnalyzed){
	  iPairDweight_Dnuisance_C_M1 = iPairDweight_Dnuisance_C;
	  iPairDweight_Dnuisance_T_M1 = iPairDweight_Dnuisance_T;
	}
      }
      // Rcpp::Rcout << iPairWeight << std::endl;

      // *** compute scores
      if(debug>0){Rcpp::Rcout << " - score("<< iFirstEndpoint <<")" << std::endl;}
      arma::uvec iUvec_endpoint = {index_endpoint[iter_d]};
      arma::uvec iUvec_status = {index_status[iter_d]};

      iPairScore = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
				status.submat(indexC[iter_strata],iUvec_status), status.submat(indexT[iter_strata],iUvec_status),
				list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata], list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
				list_lastSurv[iter_d].row(iter_strata),
				iIndex_control, iIndex_treatment, iPairWeight,
				activeUTTE, D_activeUTTE,
				Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d),
				iRP_score,
				iCount_obsC, iCount_obsT, iDscore_Dnuisance_C, iDscore_Dnuisance_T,
				iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T,
				iPairDweight_Dnuisance_C, iPairDweight_Dnuisance_T,
				zeroPlus, 
				iMethod, op[iter_d], returnIID, p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),
				iFirstEndpoint, false, iUpdateIndexNeutral, iUpdateIndexUninf, keepScore, precompute, correctionUninf, neutralAsUninf[iter_d],
				debug);

      if(iAlreadyAnalyzed){ // substract contribution of the previous analysis
	if(debug>0){Rcpp::Rcout << " - scoreM1("<< iIndex_UTTE <<")" << std::endl;}
	arma::mat iRP_score_M1= RP_score[iIndex_UTTE]; // scores relative to latest analysis of the same endpoint
	std::vector< arma::sp_mat > iRP_Dscore_Dnuisance_C_M1 = RP_Dscore_Dnuisance_C[iIndex_UTTE]; // iid of the scores relative to latest analysis of the same endpoint
	std::vector< arma::sp_mat > iRP_Dscore_Dnuisance_T_M1 = RP_Dscore_Dnuisance_T[iIndex_UTTE];

	// re-compute using current weights
	double iCount_favorable_M1,iCount_unfavorable_M1,iCount_neutral_M1,iCount_uninf_M1; // initialization necessary
	arma::mat iCount_obsC_M1,iCount_obsT_M1;
	arma::mat iDscore_Dnuisance_C_M1,iDscore_Dnuisance_T_M1;

	// note the values in endpoint, status, survTime, survJump, lastSurv are not used (only their dimensions)
	arma::mat iPairScore_M1 = calcAllPairs(endpoint.submat(indexC[iter_strata],iUvec_endpoint), endpoint.submat(indexT[iter_strata],iUvec_endpoint), threshold[iter_d],
					       status.submat(indexC[iter_strata],iUvec_status), status.submat(indexT[iter_strata],iUvec_status),
					       list_survTimeC[iter_d][iter_strata], list_survTimeT[iter_d][iter_strata], list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
					       list_lastSurv[iter_d].row(iter_strata), 
					       iIndex_control, iIndex_treatment, iPairWeight,
					       activeUTTE, D_activeUTTE,
					       iCount_favorable_M1, iCount_unfavorable_M1, iCount_neutral_M1, iCount_uninf_M1,
					       iRP_score_M1,
					       iCount_obsC_M1, iCount_obsT_M1, iDscore_Dnuisance_C_M1, iDscore_Dnuisance_T_M1,
					       iRP_Dscore_Dnuisance_C_M1, iRP_Dscore_Dnuisance_T_M1,
					       iPairDweight_Dnuisance_C_M1, iPairDweight_Dnuisance_T_M1,								
					       zeroPlus, 
					       iMethod, op[iter_d], returnIID, p_C(iter_strata, iter_d), p_T(iter_strata, iter_d),
					       false, true, false, false, keepScore, precompute, correctionUninf, neutralAsUninf[iter_d],
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
	if(debug>0){Rcpp::Rcout << " update number of pairs" << std::endl;}
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
	
	if(debug>0){Rcpp::Rcout << " update iid (" << returnIID << ")" << std::endl;}
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
	if(debug>0){Rcpp::Rcout << " update pairwise scores " << std::endl;}
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

	// update neutral.as.uninf for tte endpoints analysed for the weights (only relevant when using Peron's scoring rule)
	if(iIndex_UTTE>=0){
	  neutralAsUninfUTTE[iIndex_UTTE] = neutralAsUninf[iter_d];
	}
	  
	// update iPairWeight_nPeron, RP_score, RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
	// and re-initialize iRP_score, iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T
	if(debug>0){Rcpp::Rcout << " update score/iid for the remaing pairs("<< nUTTE_analyzedPeron_M1[iter_d+1] <<") " << std::endl;}
	updateRP(iRP_score, iRP_Dscore_Dnuisance_C, iRP_Dscore_Dnuisance_T,
		 RP_score, RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T,
		 iPairWeight_nPeron, iSize_RP, neutralAsUninf[iter_d], iter_d, correctionUninf,
		 zeroPlus, iIndex_UTTE, nUTTE_analyzedPeron_M1[iter_d+1], returnIID);
      }
  
	  
    } // end endpoint
  } // end strata
  
  // ** proportion in favor of treatment
  // Rcpp::Rcout << std::endl << " compute statistics" << std::endl;
  arma::cube delta(n_strata,D,4); // cube containing the endpoint and strata specific statistics: favorable; unfavorable; netBenefit; winRatio;
  arma::mat Delta(D,4); // matrix containing for each endpoint the overall statistics: favorable; unfavorable; netBenefit; winRatio;
  arma::mat Mvar; // variance-covariance of the overall statistics for each endpoint
  if(returnIID > 0){
    Mvar.resize(D,5); // variance(favorable); variance(unfavorable); covariance(favorable,unfavorable); variance(netBenefit); variance(winRatio);
    Mvar.fill(0.0);
  }

  if(debug>0){Rcpp::Rcout << "Compute summary statistics" << std::endl;}
  calcStatistic(delta, Delta, 
                Mcount_favorable, Mcount_unfavorable,
				iidAverage_favorable, iidAverage_unfavorable, iidNuisance_favorable, iidNuisance_unfavorable,
				Mvar, returnIID,
				posC, posT, 
                D, n_strata, vecn_pairs, vecn_control, vecn_treatment,
				weight, hprojection, pairScore, keepScore);

  // ** export
  return(Rcpp::List::create(Rcpp::Named("count_favorable") = Mcount_favorable,
							Rcpp::Named("count_unfavorable") = Mcount_unfavorable,
							Rcpp::Named("count_neutral") = Mcount_neutral,           
							Rcpp::Named("count_uninf") = Mcount_uninf,
							Rcpp::Named("delta") = delta,
							Rcpp::Named("Delta") = Delta,
							Rcpp::Named("n_pairs") = arma::conv_to< std::vector<double> >::from(vecn_pairs),
							Rcpp::Named("iidAverage_favorable") = iidAverage_favorable,
							Rcpp::Named("iidAverage_unfavorable") = iidAverage_unfavorable,
							Rcpp::Named("iidNuisance_favorable") = iidNuisance_favorable,
							Rcpp::Named("iidNuisance_unfavorable") = iidNuisance_unfavorable,
							Rcpp::Named("covariance") = Mvar,
							Rcpp::Named("tableScore")  = pairScore
							));
}

// * Function GPC2_cpp
//' @name GPC_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List GPC2_cpp(arma::mat endpoint,
					arma::mat status,
					std::vector< arma::uvec > indexC,
					std::vector< arma::uvec > posC,
					std::vector< arma::uvec > indexT,
					std::vector< arma::uvec > posT,
					std::vector< double > threshold,
					arma::vec weight,
		    arma::vec method,
		                        std::vector< int > op,
					unsigned int D,
					unsigned int D_UTTE,
					unsigned int n_strata,
					arma::vec nUTTE_analyzedPeron_M1,
					std::vector<unsigned int> index_endpoint, 
					std::vector<unsigned int> index_status, 
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
					int correctionUninf, // not used		   
					bool hierarchical,
					int hprojection,
					std::vector<bool> neutralAsUninf,
					bool keepScore,
					bool precompute,
					int returnIID,
					int debug){
  if(debug>0){Rcpp::Rcout << std::endl;}

  /// ** number of pairs
  int n_obs = endpoint.n_rows;
  bool methodPeron = (method.max()>=4); // if >= 3 then the Peron's scoring rule is used

  // number of pairs
  arma::vec vecn_pairs(n_strata); // number of pairs sumed over the strats
  arma::vec vecn_treatment(n_strata); // number of patients in the treatment group over the strats
  arma::vec vecn_control(n_strata); // number of patients in the control group over the strats
  arma::vec vecn_cumpairsM1(n_strata); // number of pairs in the previous strata (used when storing all the pairs in pairScore)
  for(unsigned int iter_strata=0 ; iter_strata < n_strata ; iter_strata++){
    vecn_control[iter_strata] = posC[iter_strata].size();
    vecn_treatment[iter_strata] = posT[iter_strata].size();		  
    vecn_pairs[iter_strata] = vecn_control[iter_strata] * vecn_treatment[iter_strata];
    if(iter_strata == 0){
      vecn_cumpairsM1[0] = 0;
    }else{
      vecn_cumpairsM1[iter_strata] = vecn_cumpairsM1[iter_strata-1] + vecn_control[iter_strata]*vecn_treatment[iter_strata];
    }
  }
  

  // ** objects storing the final results
  // score specific to each pair
  int n_pairs = arma::sum(vecn_pairs);
  std::vector< std::vector <std::vector <double> > > vecPairScore;
  if(keepScore){
    vecPairScore.resize(D);
    for(unsigned int iter_d=0; iter_d<D; iter_d++){
      vecPairScore[iter_d].resize(15);
      for(int iter_type=0; iter_type<15; iter_type++){
		vecPairScore[iter_d][iter_type].reserve(n_pairs);
      }
    }
  }
  // std::vector< arma::mat > pairScore;
  // if(keepScore){
  //   pairScore.resize(D);
  //   for(unsigned int iter_d=0; iter_d < D; iter_d++){
  //     pairScore[iter_d].resize(n_pairs,15);
  //   }
  // }
	
  // total score over pairs
  arma::mat Mcount_favorable(n_strata,D,arma::fill::zeros); // store the total weight of favorable pairs [all endpoints, strata]
  arma::mat Mcount_unfavorable(n_strata,D,arma::fill::zeros); // store the total weight of unfavorable pairs [all endpoints, strata]
  arma::mat Mcount_neutral(n_strata,D,arma::fill::zeros); // store the total weight of neutral pairs [all endpoints, strata]
  arma::mat Mcount_uninf(n_strata,D,arma::fill::zeros); // store the total weight of uninf pairs [all endpoints, strata]

  // iid
  arma::mat iidAverage_favorable; // iid relative to the average over all pairs for the favorable scores [all endpoints]
  arma::mat iidAverage_unfavorable; // iid relative to the average over all pairs for the unfavorable scores [all endpoints]
  arma::mat iidNuisance_favorable; // iid relative to the nuisance parameters for the favorable scores [all endpoints]
  arma::mat iidNuisance_unfavorable; // iid relative to the nuisance parameters for the unfavorable scores [all endpoints]

  if(returnIID>0){
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

  // ** initialization for the loop
  // over strata
  arma::uvec indexStrataC;
  arma::uvec indexStrataT;
  arma::uvec posStrataC;
  arma::uvec posStrataT;
  unsigned int nStrata_Control;
  unsigned int nStrata_Treatment;

  std::vector< arma::mat > Dfavorable_Dnuisance_strataC; // sum over all pairs of a given strata of the partial derivative regarding nuisance parameters 
  std::vector< arma::mat > Dfavorable_Dnuisance_strataT; // sum over all pairs of a given strata of the partial derivative regarding nuisance parameters 
  std::vector< arma::mat > Dunfavorable_Dnuisance_strataC; // sum over all pairs of a given strata of the partial derivative regarding nuisance parameters 
  std::vector< arma::mat > Dunfavorable_Dnuisance_strataT; // sum over all pairs of a given strata of the partial derivative regarding nuisance parameters 
  if(returnIID>1){
    Dfavorable_Dnuisance_strataC.resize(D_UTTE);
    Dfavorable_Dnuisance_strataT.resize(D_UTTE);
    Dunfavorable_Dnuisance_strataC.resize(D_UTTE);
    Dunfavorable_Dnuisance_strataT.resize(D_UTTE);
  }
  
  // over pair
  int iPair = 0;

  arma::vec iFavorable_UTTE; // store favorable score computed at the last threshold for each UTTE
  arma::vec iUnfavorable_UTTE; // store unfavorable score computed at the last threshold for each UTTE
  arma::vec iWeight_UTTE; // store weight (i.e. complement to 1 of fav+defav) computed at the last threshold for each UTTE
  
  std::vector < arma::vec > iDweight_Dnuisance_C_UTTE; // store derivative of the weight regarding the nuisance parameters computed at the last threshold for each UTTE
  std::vector < arma::vec > iDweight_Dnuisance_T_UTTE; // store derivative of the weight regarding the nuisance parameters computed at the last threshold for each UTTE

  std::vector < arma::mat > iDscore_Dnuisance_C_UTTE; // store derivative of the weight regarding the score computed at the last threshold for each UTTE
  std::vector < arma::mat > iDscore_Dnuisance_T_UTTE; // store derivative of the weight regarding the score computed at the last threshold for each UTTE

  std::vector < arma::mat > iDscore_Dnuisance_C_calcOnePair; // store derivative of the regarding the score computed for each UTTE (input to calcOnePair_SurvPeron)
  std::vector < arma::mat > iDscore_Dnuisance_T_calcOnePair; // store derivative of the regarding the score computed for each UTTE (input to calcOnePair_SurvPeron)
  if(methodPeron){
    iFavorable_UTTE.resize(D_UTTE);
    iUnfavorable_UTTE.resize(D_UTTE);
    iWeight_UTTE.resize(D_UTTE);

    iDweight_Dnuisance_C_UTTE.resize(D_UTTE);
    iDweight_Dnuisance_T_UTTE.resize(D_UTTE);

    iDscore_Dnuisance_C_UTTE.resize(D_UTTE);
    iDscore_Dnuisance_T_UTTE.resize(D_UTTE);

    iDscore_Dnuisance_C_calcOnePair.resize(D);
    iDscore_Dnuisance_T_calcOnePair.resize(D);
  }

  // over endpoint
  int iIndex_UTTE_d; // position of the TTE among the UTTE
  int iMethod; // method used to analyze the endpoint
  double iTempoOperatorDouble; // for switching favorable to unfavorable when op is <0
  arma::colvec iTempoOperatorCol; // for switching favorable to unfavorable when op is <0

  double iCumWeight; // current weight of the pair
  double iNewWeight; // remaining weight to analyze after the current endpoint

  std::vector< double > iPairScore;
  // ** loop over strata
  
  for(unsigned int iter_strata=0 ; iter_strata < n_strata ; iter_strata++){
    if(debug>0){Rcpp::Rcout << "Strata " << iter_strata << std::endl;}

    // subset by strata
    indexStrataC = indexC[iter_strata];
    indexStrataT = indexT[iter_strata];
    posStrataC = posC[iter_strata];
    posStrataT = posT[iter_strata];
    nStrata_Control = indexStrataC.size();
    nStrata_Treatment = indexStrataT.size();

    // prepare d(survival)/d(nuisance)
    if(returnIID > 1){ // 
      for(unsigned int iter_d=0; iter_d < D; iter_d ++){
		if(p_C(iter_strata,iter_d)>0){
		  Dfavorable_Dnuisance_strataC[index_UTTE[iter_d]].resize(p_C(iter_strata,iter_d),D);
		  Dfavorable_Dnuisance_strataC[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata
		  Dunfavorable_Dnuisance_strataC[index_UTTE[iter_d]].resize(p_C(iter_strata,iter_d),D);
		  Dunfavorable_Dnuisance_strataC[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata

		  iDscore_Dnuisance_C_UTTE[index_UTTE[iter_d]].resize(p_C(iter_strata,iter_d),4);
		  iDscore_Dnuisance_C_calcOnePair[iter_d].resize(p_C(iter_strata,iter_d),4);
		}
		if(p_T(iter_strata,iter_d)>0){
		  Dfavorable_Dnuisance_strataT[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),D);
		  Dfavorable_Dnuisance_strataT[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata
		  Dunfavorable_Dnuisance_strataT[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),D);
		  Dunfavorable_Dnuisance_strataT[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata

		  iDscore_Dnuisance_T_UTTE[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),4);
		  iDscore_Dnuisance_T_calcOnePair[iter_d].resize(p_T(iter_strata,iter_d),4);
		}
      }
    }

    // *** loop over pairs
    if(debug>0){Rcpp::Rcout << " - compute scores" << std::endl;}
    for(unsigned int iter_C=0 ; iter_C < nStrata_Control; iter_C++){
      for(unsigned int iter_T=0 ; iter_T < nStrata_Treatment; iter_T++){
		if(debug>1){Rcpp::Rcout << " pair " << iPair << " (" << iter_C << ";" << iter_T << ") ";}

		// **** loop over endpoints
		for(unsigned int iter_d=0 ; iter_d < D; iter_d++){
		  if(debug==3){Rcpp::Rcout << "*" << std::endl;}
		  iIndex_UTTE_d = index_UTTE[iter_d];
		  iMethod = method[iter_d];
  
		  // **** compute weight
		  if(debug>3){Rcpp::Rcout << "w";}
		  iCumWeight = 1;
		  iNewWeight = 0;
		  if(hierarchical && methodPeron){
			for(int iter_UTTE=0 ; iter_UTTE<nUTTE_analyzedPeron_M1[iter_d]; iter_UTTE++){
			  if(iter_UTTE != iIndex_UTTE_d){
				iCumWeight *= iWeight_UTTE[iter_UTTE];
			  }
			}
			if(iCumWeight<zeroPlus){break;}
		  }
	      
		  // **** compute score
		  if(debug>3){Rcpp::Rcout << "s";}
		  if(iMethod == 1){ // continuous or binary endpoint
			iPairScore = calcOnePair_Continuous(endpoint(indexStrataT[iter_T], index_endpoint[iter_d]) - endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
												threshold[iter_d]);
		  }else if(iMethod == 2){ // time to event endpoint with Gehan's scoring rule (right-censored, survival or competing risks)
			iPairScore = calcOnePair_TTEgehan(endpoint(indexStrataT[iter_T], index_endpoint[iter_d]) - endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
											  status(indexStrataC[iter_C], index_status[iter_d]),
											  status(indexStrataT[iter_T], index_status[iter_d]),
											  threshold[iter_d]);
		  }else if(iMethod == 3){ // time to event endpoint with Gehan's scoring rule (left-censored, survival or competing risks)
			iPairScore = calcOnePair_TTEgehan2(endpoint(indexStrataT[iter_T], index_endpoint[iter_d]) - endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
											   status(indexStrataC[iter_C], index_status[iter_d]),
											   status(indexStrataT[iter_T], index_status[iter_d]),
											   threshold[iter_d]);
		  }else if(iMethod == 4){  // time to event endpoint with Peron's scoring rule (right-censored, survival)

			// note: iDscore_Dnuisance_C, iDscore_Dnuisance_T are initalized to 0 in calcOnePair_SurvPeron
			iPairScore = calcOnePair_SurvPeron(endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
											   endpoint(indexStrataT[iter_T], index_endpoint[iter_d]),
											   status(indexStrataC[iter_C], index_status[iter_d]),
											   status(indexStrataT[iter_T], index_status[iter_d]),
											   threshold[iter_d],
											   list_survTimeC[iter_d][iter_strata].row(iter_C), list_survTimeT[iter_d][iter_strata].row(iter_T),
											   list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
											   list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1),
											   iDscore_Dnuisance_C_calcOnePair[iter_d], iDscore_Dnuisance_T_calcOnePair[iter_d],
											   p_C(iter_strata, iter_d), p_T(iter_strata, iter_d), precompute, returnIID);

		  }else if(iMethod == 5){  // time to event endpoint with Peron's scoring rule (right-censored, competing risks)
			iPairScore = calcOnePair_CRPeron(endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
											 endpoint(indexStrataT[iter_T], index_endpoint[iter_d]),
											 status(indexStrataC[iter_C], index_status[iter_d]),
											 status(indexStrataT[iter_T], index_status[iter_d]),
											 threshold[iter_d],
											 list_survTimeC[iter_d][iter_strata].row(iter_C), list_survTimeT[iter_d][iter_strata].row(iter_T),
											 list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],					     
											 list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1), list_lastSurv[iter_d](iter_strata,2), list_lastSurv[iter_d](iter_strata,3),
											 iDscore_Dnuisance_C_calcOnePair[iter_d], iDscore_Dnuisance_T_calcOnePair[iter_d],
											 p_C(iter_strata, iter_d), p_T(iter_strata, iter_d), precompute, returnIID);
		  }
	  
		  // **** operator
		  if(op[iter_d]<0){
		    iTempoOperatorDouble = iPairScore[0];
		    iPairScore[0] = iPairScore[1];
		    iPairScore[1] = iTempoOperatorDouble;

		    if(returnIID > 1 && (iMethod >= 4)){
		      iTempoOperatorCol = iDscore_Dnuisance_C_calcOnePair[iter_d].col(0);
		      iDscore_Dnuisance_C_calcOnePair[iter_d].col(0) = iDscore_Dnuisance_C_calcOnePair[iter_d].col(1);
		      iDscore_Dnuisance_C_calcOnePair[iter_d].col(1) = iTempoOperatorCol;

		      iTempoOperatorCol = iDscore_Dnuisance_T_calcOnePair[iter_d].col(0);
		      iDscore_Dnuisance_T_calcOnePair[iter_d].col(0) = iDscore_Dnuisance_T_calcOnePair[iter_d].col(1);
		      iDscore_Dnuisance_T_calcOnePair[iter_d].col(1) = iTempoOperatorCol;
		    }
		  }

		  // **** remove contribution from previously analyzed threshold of the same endpoint
		  if( (iMethod >= 4) && (nUTTE_analyzedPeron_M1[iter_d]>iIndex_UTTE_d) ){  // endpoint already analyzed 
			iPairScore[0] -= iFavorable_UTTE[iIndex_UTTE_d];
			iPairScore[1] -= iUnfavorable_UTTE[iIndex_UTTE_d];

			if(returnIID>1){
			  iDscore_Dnuisance_C_calcOnePair[iter_d].col(0) -= iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d].col(0);
			  iDscore_Dnuisance_C_calcOnePair[iter_d].col(1) -= iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d].col(1);
			  iDscore_Dnuisance_T_calcOnePair[iter_d].col(0) -= iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d].col(0);
			  iDscore_Dnuisance_T_calcOnePair[iter_d].col(1) -= iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d].col(1);
			}
		  }
		  
		  // **** aggregate favorable score and iid over analyzed pairs
		  // if(iPairScore[0] > zeroPlus){
		  if(debug==4){Rcpp::Rcout << "f";}
		  if(debug>4){Rcpp::Rcout << " favorable=" << iPairScore[0] << " ";}

		  // score
		  Mcount_favorable(iter_strata,iter_d) += iPairScore[0] * iCumWeight;
	  
		  if(returnIID > 0){
			// iid (average)
			iidAverage_favorable(posStrataC[iter_C],iter_d) += iPairScore[0] * iCumWeight;
			iidAverage_favorable(posStrataT[iter_T],iter_d) += iPairScore[0] * iCumWeight;
		  }
		  // iid (nuisance) for the score
		  if( (returnIID > 1) && (iMethod >= 4) ){
			Dfavorable_Dnuisance_strataC[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_C_calcOnePair[iter_d].col(0) * iCumWeight;
			Dfavorable_Dnuisance_strataT[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_T_calcOnePair[iter_d].col(0) * iCumWeight;
		  }
		  // iid (nuisance) for the weight of the pair
		  if( (returnIID > 1) && (nUTTE_analyzedPeron_M1[iter_d] > 0) && hierarchical ){
			for(int iter_UTTE=0 ; iter_UTTE<nUTTE_analyzedPeron_M1[iter_d]; iter_UTTE++){
			  if(iter_UTTE != iIndex_UTTE_d){
				Dfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d) += (iPairScore[0] * iCumWeight / iWeight_UTTE[iter_UTTE]) * iDweight_Dnuisance_C_UTTE[iter_UTTE] ;
				Dfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d) += (iPairScore[0] * iCumWeight / iWeight_UTTE[iter_UTTE]) * iDweight_Dnuisance_T_UTTE[iter_UTTE] ;
			  }
			}
		  }

		  // **** aggregate unfavorable score and iid over analyzed pairs
		  // if(iPairScore[1] > zeroPlus){
		  if(debug==4){Rcpp::Rcout << "d";}
		  if(debug>4){Rcpp::Rcout << " unfavorable=" << iPairScore[1] << " ";}

		  // score
		  Mcount_unfavorable(iter_strata,iter_d) += iPairScore[1] * iCumWeight;
      
		  if(returnIID > 0){
			// iid (average)
			iidAverage_unfavorable(posStrataC[iter_C],iter_d) += iPairScore[1] * iCumWeight;
			iidAverage_unfavorable(posStrataT[iter_T],iter_d) += iPairScore[1] * iCumWeight;
		  }
      
		  // iid (nuisance) for the score
		  if( (returnIID > 1) && (iMethod >= 4) ){
			Dunfavorable_Dnuisance_strataC[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_C_calcOnePair[iter_d].col(1) * iCumWeight;
			Dunfavorable_Dnuisance_strataT[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_T_calcOnePair[iter_d].col(1) * iCumWeight;
		  }

		  // iid (nuisance) for the weight of the pair
		  if( (returnIID > 1) && (nUTTE_analyzedPeron_M1[iter_d] > 0) && hierarchical ){
			for(int iter_UTTE=0 ; iter_UTTE<nUTTE_analyzedPeron_M1[iter_d]; iter_UTTE++){
			  if(iter_UTTE != iIndex_UTTE_d){	      
				Dunfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d) += (iPairScore[1] * iCumWeight / iWeight_UTTE[iter_UTTE]) * iDweight_Dnuisance_C_UTTE[iter_UTTE] ;
				Dunfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d) += (iPairScore[1] * iCumWeight / iWeight_UTTE[iter_UTTE]) * iDweight_Dnuisance_T_UTTE[iter_UTTE] ;
			  }
			}
		  }
    
		  // **** aggregate neutral score and iid over analyzed pairs
		  if(iPairScore[2] > zeroPlus){
			if(debug==4){Rcpp::Rcout << "n";}
			if(debug>4){Rcpp::Rcout << " neutral=" << iPairScore[2] << " ";}
		  
			// score
			Mcount_neutral(iter_strata,iter_d) += iPairScore[2] * iCumWeight;

			// update weight
			if(neutralAsUninf[iter_d]){
			  iNewWeight += iPairScore[2];
			}
		  }
    
		  // **** aggregate uninformative score and iid over analyzed pairs
		  if(iPairScore[3] > zeroPlus){
			if(debug==4){Rcpp::Rcout << "u";}
			if(debug>4){Rcpp::Rcout << " uninformative=" << iPairScore[3] << " " ;}

			// score
			Mcount_uninf(iter_strata,iter_d) += iPairScore[3] * iCumWeight;

			// update weight
			iNewWeight += iPairScore[3];
		  }

		  // **** update pairwise-scores for all pairs
		  if(keepScore){
			if(debug>3){Rcpp::Rcout << " keepScore ";}
			vecPairScore[iter_d][0].push_back(iter_strata);
			vecPairScore[iter_d][1].push_back(posStrataC[iter_C]);
			vecPairScore[iter_d][2].push_back(posStrataT[iter_T]);
			vecPairScore[iter_d][3].push_back(iPair);
			vecPairScore[iter_d][4].push_back(iter_C);
			vecPairScore[iter_d][5].push_back(iter_T);

			vecPairScore[iter_d][6].push_back(iPairScore[0]);
			vecPairScore[iter_d][7].push_back(iPairScore[1]);
			vecPairScore[iter_d][8].push_back(iPairScore[2]);
			vecPairScore[iter_d][9].push_back(iPairScore[3]);
			vecPairScore[iter_d][10].push_back(iCumWeight);
			vecPairScore[iter_d][11].push_back(iPairScore[0] * iCumWeight);
			vecPairScore[iter_d][12].push_back(iPairScore[1] * iCumWeight);
			vecPairScore[iter_d][13].push_back(iPairScore[2] * iCumWeight);
			vecPairScore[iter_d][14].push_back(iPairScore[3] * iCumWeight);
		  } 

		  // **** early stop if nothing left or store weight when TTE endpoint with Peron's scoring rule
		  if(hierarchical){
			if( (iNewWeight < zeroPlus) || (iter_d == (D-1)) ){
			  if(debug>3){Rcpp::Rcout << " exit ";}
			  break;
			}else if(methodPeron && (iIndex_UTTE_d>=0) ){

			  if(debug>3){Rcpp::Rcout << " store ";}
			  if(nUTTE_analyzedPeron_M1[iter_d]>iIndex_UTTE_d){ // endpoint already analyzed (add to previous contributions)
				iFavorable_UTTE[iIndex_UTTE_d] += iPairScore[0];
				iUnfavorable_UTTE[iIndex_UTTE_d] += iPairScore[1];
			  }else{ // restart
				iFavorable_UTTE[iIndex_UTTE_d] = iPairScore[0];
				iUnfavorable_UTTE[iIndex_UTTE_d] = iPairScore[1];
			  }
			  iWeight_UTTE[iIndex_UTTE_d] = iNewWeight;
	      
			  if(returnIID>1){
				if(nUTTE_analyzedPeron_M1[iter_d]>iIndex_UTTE_d){ // endpoint already analyzed (add to previous contributions)
				  iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d] += iDscore_Dnuisance_C_calcOnePair[iter_d];
				  iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d] += iDscore_Dnuisance_T_calcOnePair[iter_d];
				}else{ // restart
				  iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_C_calcOnePair[iter_d];
				  iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_T_calcOnePair[iter_d];
				}
		
				if(neutralAsUninf[iter_d] && (iPairScore[2] > zeroPlus)){
				  iDweight_Dnuisance_C_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_C_calcOnePair[iter_d].col(2) + iDscore_Dnuisance_C_calcOnePair[iter_d].col(3);
				  iDweight_Dnuisance_T_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_T_calcOnePair[iter_d].col(2) + iDscore_Dnuisance_T_calcOnePair[iter_d].col(3);
				}else{
				  iDweight_Dnuisance_C_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_C_calcOnePair[iter_d].col(3);
				  iDweight_Dnuisance_T_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_T_calcOnePair[iter_d].col(3);
				}
			  }
			}
		  }
	  
		}
	
		// *** update pair number
		iPair++;
		if(iPair % 65536 == 0){
		  R_CheckUserInterrupt();
		}
		if(debug>2){Rcpp::Rcout << " done " << std::endl;}
      }
    }
    if(debug>1){Rcpp::Rcout << std::endl;}

    // ** compute iid nuisance    
    if(returnIID>1){
      if(debug>0){Rcpp::Rcout << "compute iid nuisance" << std::endl;}

      for(unsigned int iter_d=0; iter_d < D; iter_d++){
		for(unsigned int iter_UTTE=0 ; iter_UTTE<D_UTTE; iter_UTTE++){	  
		  iidNuisance_favorable.col(iter_d) += iid_survJumpC[iter_UTTE][iter_strata] * Dfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_favorable.col(iter_d) += iid_survJumpT[iter_UTTE][iter_strata] * Dfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[iter_UTTE][iter_strata] * Dunfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[iter_UTTE][iter_strata] * Dunfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
		}
      }
    }

  }
  
  // ** combine pairwise scores into a list of matrices 
  std::vector< arma::mat> pairScore;
  if(keepScore){
    if(debug>0){Rcpp::Rcout << "generate pairScore" << std::endl;}
	int iNpairs;
    pairScore.resize(D);
    for(unsigned int iter_d=0; iter_d<D; iter_d++){
      iNpairs = vecPairScore[iter_d][0].size();
      pairScore[iter_d].resize(iNpairs,15);
      for(int iter_type=0; iter_type<15; iter_type++){
		for(int iter_pair=0; iter_pair<iNpairs; iter_pair++){
	  pairScore[iter_d](iter_pair,iter_type) = vecPairScore[iter_d][iter_type][iter_pair];
	  // pairScore[iter_d].col(iter_type) = arma::conv_to< arma::colvec >::from(vecPairScore[iter_d][iter_type]);
	}
      }
    }
  }
  
  // ** proportion in favor of treatment
  // Rcpp::Rcout << std::endl << " compute statistics" << std::endl;
  arma::cube delta(n_strata,D,4); // cube containing the endpoint and strata specific statistics: favorable; unfavorable; netBenefit; winRatio;
  arma::mat Delta(D,4); // matrix containing for each endpoint the overall statistics: favorable; unfavorable; netBenefit; winRatio;
  arma::mat Mvar; // variance-covariance of the overall statistics for each endpoint
  if(returnIID > 0){
    Mvar.resize(D,5); // variance(favorable); variance(unfavorable); covariance(favorable,unfavorable); variance(netBenefit); variance(winRatio);
    Mvar.fill(0.0);
  }

   if(debug>0){Rcpp::Rcout << "Compute summary statistics" << std::endl;}
   calcStatistic(delta, Delta, 
		Mcount_favorable, Mcount_unfavorable,
		iidAverage_favorable, iidAverage_unfavorable, iidNuisance_favorable, iidNuisance_unfavorable,
		Mvar, returnIID,
		posC, posT, 
		D, n_strata, vecn_pairs, vecn_control, vecn_treatment,
		weight, hprojection, pairScore, keepScore);

  // ** export
  return(Rcpp::List::create(Rcpp::Named("count_favorable") = Mcount_favorable,
			    Rcpp::Named("count_unfavorable") = Mcount_unfavorable,
			    Rcpp::Named("count_neutral") = Mcount_neutral,           
			    Rcpp::Named("count_uninf") = Mcount_uninf,
			    Rcpp::Named("delta") = delta,
			    Rcpp::Named("Delta") = Delta,
			    Rcpp::Named("n_pairs") = arma::conv_to< std::vector<double> >::from(vecn_pairs),
			    Rcpp::Named("iidAverage_favorable") = iidAverage_favorable,
			    Rcpp::Named("iidAverage_unfavorable") = iidAverage_unfavorable,
			    Rcpp::Named("iidNuisance_favorable") = iidNuisance_favorable,
			    Rcpp::Named("iidNuisance_unfavorable") = iidNuisance_unfavorable,
			    Rcpp::Named("covariance") = Mvar,
			    Rcpp::Named("tableScore")  = pairScore
			    ));
}

// * prepareWeight
// author Brice Ozenne
void prepareWeight(arma::vec& iPairWeight, std::vector<std::vector< arma::sp_mat >>& iPairDweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& iPairDweight_Dnuisance_T,
		   std::vector<int>& activeUTTE, int& D_activeUTTE,
		   int iter_d, int iIndex_UTTE, const std::vector<arma::mat>& RP_score,
		   const std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_C, const std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_T,
		   int iNUTTE_analyzedPeron, int correctionUninf, double zeroPlus, std::vector<bool>& neutralAsUninf, int returnIID){


  // ** cumulate weights over previous TTE endpoints
  activeUTTE.resize(0);
  for(int iter_UTTE=0 ; iter_UTTE<iNUTTE_analyzedPeron; iter_UTTE++){
    if(iter_UTTE != iIndex_UTTE){
      activeUTTE.push_back(iter_UTTE);
      if(neutralAsUninf[iter_UTTE]){
	iPairWeight %= (RP_score[iter_UTTE].col(2) + RP_score[iter_UTTE].col(3));
      }else{
	iPairWeight %= RP_score[iter_UTTE].col(3);
      }
    }
  }
  D_activeUTTE = activeUTTE.size();

  // ** compute iid of the weights
  if(returnIID > 1){
    iPairDweight_Dnuisance_C[0].resize(iNUTTE_analyzedPeron);
    iPairDweight_Dnuisance_T[0].resize(iNUTTE_analyzedPeron);
	
    for(int iter_UTTE=0 ; iter_UTTE<iNUTTE_analyzedPeron; iter_UTTE++){
      if(iter_UTTE != iIndex_UTTE){
	if(neutralAsUninf[iter_UTTE]){
	  iPairDweight_Dnuisance_C[0][iter_UTTE] = RP_Dscore_Dnuisance_C[iter_UTTE][2] + RP_Dscore_Dnuisance_C[iter_UTTE][3];
	  for(size_t iCol=0; iCol < iPairDweight_Dnuisance_C[0][iter_UTTE].n_rows; iCol++){
	    iPairDweight_Dnuisance_C[0][iter_UTTE].row(iCol) %= arma::trans(iPairWeight/(RP_score[iter_UTTE].col(2) + RP_score[iter_UTTE].col(3)));
	  }
	  iPairDweight_Dnuisance_T[0][iter_UTTE] = RP_Dscore_Dnuisance_T[iter_UTTE][2] + RP_Dscore_Dnuisance_T[iter_UTTE][3];
	  for(size_t iCol=0; iCol < iPairDweight_Dnuisance_T[0][iter_UTTE].n_rows; iCol++){
	    iPairDweight_Dnuisance_T[0][iter_UTTE].row(iCol) %= arma::trans(iPairWeight/(RP_score[iter_UTTE].col(2) + RP_score[iter_UTTE].col(3)));
	  }
	}else{
	  iPairDweight_Dnuisance_C[0][iter_UTTE] = RP_Dscore_Dnuisance_C[iter_UTTE][3];
	  for(size_t iCol=0; iCol < iPairDweight_Dnuisance_C[0][iter_UTTE].n_rows; iCol++){
	    iPairDweight_Dnuisance_C[0][iter_UTTE].row(iCol) %= arma::trans(iPairWeight/RP_score[iter_UTTE].col(3));
	  }
	  iPairDweight_Dnuisance_T[0][iter_UTTE] = RP_Dscore_Dnuisance_T[iter_UTTE][3];
	  for(size_t iCol=0; iCol < iPairDweight_Dnuisance_T[0][iter_UTTE].n_rows; iCol++){
	    iPairDweight_Dnuisance_T[0][iter_UTTE].row(iCol) %= arma::trans(iPairWeight/RP_score[iter_UTTE].col(3));
	  }
	}
      }else{
	iPairDweight_Dnuisance_C[0][iter_UTTE].resize(0,0);
	iPairDweight_Dnuisance_T[0][iter_UTTE].resize(0,0);
      }
    }
    iPairDweight_Dnuisance_C[1] = iPairDweight_Dnuisance_C[0];
    iPairDweight_Dnuisance_T[1] = iPairDweight_Dnuisance_T[0];
    iPairDweight_Dnuisance_C[2] = iPairDweight_Dnuisance_C[0];
    iPairDweight_Dnuisance_T[2] = iPairDweight_Dnuisance_T[0];
    iPairDweight_Dnuisance_C[3] = iPairDweight_Dnuisance_C[0];
    iPairDweight_Dnuisance_T[3] = iPairDweight_Dnuisance_T[0];
  }  

  return;
}

// * updateIID
// author Brice Ozenne
void updateIID(arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, 
	       arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, 
	       const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT,
	       const arma::mat& iCount_obsC, const arma::mat& iCount_obsT,
	       const std::vector<int>& activeUTTE, int D_activeUTTE,
	       const arma::mat& iDscore_Dnuisance_C, const arma::mat& iDscore_Dnuisance_T,
	       const std::vector< std::vector< arma::mat > >& iid_survJumpC, const std::vector< std::vector< arma::mat > >& iid_survJumpT,
	       const std::vector<std::vector< arma::sp_mat >> & iPairDweight_Dnuisance_C,
	       const std::vector<std::vector< arma::sp_mat >> & iPairDweight_Dnuisance_T,
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
    arma::colvec iDweight_Dnuisance_C;
    arma::colvec iDweight_Dnuisance_T;
    
    for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
      // DOCUMENTATION Armadillo
      // sum: For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1)
      iDweight_Dnuisance_C = sum(iPairDweight_Dnuisance_C[0][activeUTTE[iter_UTTE]],1);
      iidNuisance_favorable.col(iter_d) += iid_survJumpC[activeUTTE[iter_UTTE]][iter_strata] * iDweight_Dnuisance_C / vecn_pairs[iter_strata];
      iDweight_Dnuisance_T = sum(iPairDweight_Dnuisance_T[0][activeUTTE[iter_UTTE]],1);
      iidNuisance_favorable.col(iter_d) += iid_survJumpT[activeUTTE[iter_UTTE]][iter_strata] * iDweight_Dnuisance_T / vecn_pairs[iter_strata];

      iDweight_Dnuisance_C = sum(iPairDweight_Dnuisance_C[1][activeUTTE[iter_UTTE]],1);
      iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[activeUTTE[iter_UTTE]][iter_strata] * iDweight_Dnuisance_C / vecn_pairs[iter_strata];
      iDweight_Dnuisance_T = sum(iPairDweight_Dnuisance_T[1][activeUTTE[iter_UTTE]],1);
      iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[activeUTTE[iter_UTTE]][iter_strata] * iDweight_Dnuisance_T / vecn_pairs[iter_strata];
    }

    // *** iid of the proba/score
    if(iMethod >= 4){
      iidNuisance_favorable.col(iter_d) += iid_survJumpC[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_C.col(0)/vecn_pairs[iter_strata];
      iidNuisance_favorable.col(iter_d) += iid_survJumpT[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_T.col(0)/vecn_pairs[iter_strata];
      iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_C.col(1)/vecn_pairs[iter_strata];
      iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[iIndex_UTTE][iter_strata] * iDscore_Dnuisance_T.col(1)/vecn_pairs[iter_strata];
    }

  }

  return;
}


// * updatePairScore
// author Brice Ozenne
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
// author Brice Ozenne
void updateRP(arma::mat& iRP_score, std::vector< arma::sp_mat >& iRP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& iRP_Dscore_Dnuisance_T,
	      std::vector<arma::mat>& RP_score, std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_C, std::vector< std::vector< arma::sp_mat > >& RP_Dscore_Dnuisance_T,
	      arma::vec& iPairWeight_nPeron, int iSize_RP, bool neutralAsUninf, int iter_d, int correctionUninf,
	      double zeroPlus, int iIndex_UTTE, int nUTTE_analyzedPeron, int returnIID){

  // ** index of the residual pairs among the analyzed pairs
  arma::uvec iIndex_RP = arma::conv_to<arma::uvec>::from(iRP_score.col(0));
  // Rcpp::Rcout << iIndex_RP.size() << std::endl;
  
  // ** update score and iid associated to the residual pairs for TTE endpoint analyzed with Peron's scoring rule
  for(int iter_UTTE=0; iter_UTTE<nUTTE_analyzedPeron; iter_UTTE++){
    if(iter_UTTE == iIndex_UTTE){
      // update with the result of the current endpoint
      arma::uvec iter_36 = arma::linspace<arma::uvec>(3, 6, 4);
      RP_score[iter_UTTE] = iRP_score.cols(iter_36);
      if(returnIID>1){
	for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
	  RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP] = iRP_Dscore_Dnuisance_C[iter_typeRP];
	  RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP] = iRP_Dscore_Dnuisance_T[iter_typeRP];
	}
      }
    }else{	  
      // subset the result of the previous endpoints, only considering the residual pairs
      RP_score[iter_UTTE] = RP_score[iter_UTTE].rows(iIndex_RP);
      if(returnIID>1){
	for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
	  RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP] = subcol_sp_mat(RP_Dscore_Dnuisance_C[iter_UTTE][iter_typeRP], iIndex_RP);
	  RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP] = subcol_sp_mat(RP_Dscore_Dnuisance_T[iter_UTTE][iter_typeRP], iIndex_RP);
	}
      }
    }
			
  }
	
  // ** update weights associated to the residual pairs according to the correction for uninformative pairs
  if(correctionUninf == 0){
    // no correction: all pairs have a weight of 1
    iPairWeight_nPeron.resize(iSize_RP);
    iPairWeight_nPeron.fill(1.0);
  }else{
    // correction
    // 1- gather the weights that were previously used
    if(iter_d>0){
      iPairWeight_nPeron = iPairWeight_nPeron.rows(iIndex_RP);
    }else{
      iPairWeight_nPeron.resize(iSize_RP);
      iPairWeight_nPeron.fill(1.0);
    }

    // 2- update the weights with neutral/uninformative contribution of the current endpoint
    if(iIndex_UTTE < (-zeroPlus)){ // i.e. not analyzed by Peron
      if(neutralAsUninf){
	iPairWeight_nPeron %= (iRP_score.col(5) + iRP_score.col(6));
      }else{
	iPairWeight_nPeron %= iRP_score.col(6);
      }	  
    }
  }

  return ;
}

// * subcol_sp_mat
// inspired from https://stackoverflow.com/questions/40222092/access-and-modify-the-non-zero-elements-of-sparse-matrix-of-class-armasp-mat-u
arma::sp_mat subcol_sp_mat(const arma::sp_mat& X, arma::uvec index){
  int p = X.n_rows;

  // WARNING only works for subsetting over columns since the value in a sp_mat are stored by column
  // also assumes that index is sorted
  int nIndex = index.size();
  bool check = index.is_sorted();
  if(check == false){
    Rcpp::Rcout << "WARNING: index should be sorted when subsetting a sparse matrix by column. Operation may give incorrect results" << std::endl;
  }
  
  // Make const iterator
  arma::sp_mat::const_iterator it = X.begin();
  arma::sp_mat::const_iterator it_end  = X.end();

  // Rcout << X << endl;
  // Calculate number of points
  int n = std::distance(it, it_end);
  if (n <= 0) {
    arma::sp_mat Y(p, nIndex);
    return(Y);
  }
  
  // Collecting locations
  int iIndex = 0;
  unsigned int iCol;
  int iRow;

  std::vector<int> Vrow;
  Vrow.reserve(n);
  std::vector<int> Vcol;
  Vcol.reserve(n);
  std::vector<double> Vvalue;
  Vvalue.reserve(n);

  for(int i = 0; i < n; ++i){
    iCol = it.col();
    iRow = it.row();
    while((nIndex > iIndex) && (iCol > index[iIndex])){
      iIndex++;
    }
    if((iIndex >= nIndex) || (iCol > index[iIndex])){
      break;
    }else if(iCol == index[iIndex]){
      Vrow.push_back(iRow);
      Vcol.push_back(iIndex);
      Vvalue.push_back((*it));
    } 
    ++it;
  }
  
  // Generate matrix
  int nValue = Vvalue.size();
  arma::umat loc(2,nValue);
  for(int iV=0; iV<nValue; iV++){
    loc(0,iV) = Vrow[iV];
    loc(1,iV) = Vcol[iV];
  }
  arma::sp_mat Y(loc, arma::conv_to<arma::colvec>::from(Vvalue), p, nIndex);
  return(Y);
}
