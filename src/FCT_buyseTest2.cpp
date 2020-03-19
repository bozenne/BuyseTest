// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "FCT_calcOnePair.h"
#include "FCT_calcStatistic.h"

// * Function GPC2_cpp
//' @name GPC2_cpp
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
		    bool neutralAsUninf,
		    bool keepScore,
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
  std::vector< arma::mat > pairScore;
  if(keepScore){
    pairScore.resize(D);
    for(unsigned int iter_d=0; iter_d < D; iter_d++){
      pairScore[iter_d].resize(n_pairs,15);
    }
  }
	
  // total score over pairs
  arma::mat Mcount_favorable(n_strata,D,arma::fill::zeros); // store the total weight of favorable pairs [all endpoints, strata]
  arma::mat Mcount_unfavorable(n_strata,D,arma::fill::zeros); // store the total weight of unfavorable pairs [all endpoints, strata]
  arma::mat Mcount_neutral(n_strata,D,arma::fill::zeros); // store the total weight of neutral pairs [all endpoints, strata]
  arma::mat Mcount_uninf(n_strata,D,arma::fill::zeros); // store the total weight of uninf pairs [all endpoints, strata]

  // variance and iid
  arma::mat Mvar; // variance-covariance (favorable,unfavorable) scores  [all endpoints]

  arma::mat iidAverage_favorable; // iid relative to the average over all pairs for the favorable scores [all endpoints]
  arma::mat iidAverage_unfavorable; // iid relative to the average over all pairs for the unfavorable scores [all endpoints]
  arma::mat iidNuisance_favorable; // iid relative to the nuisance parameters for the favorable scores [all endpoints]
  arma::mat iidNuisance_unfavorable; // iid relative to the nuisance parameters for the unfavorable scores [all endpoints]

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

  // ** initialization for the loop
  // over strata
  arma::uvec indexStrataC;
  arma::uvec indexStrataT;
  arma::uvec posStrataC;
  arma::uvec posStrataT;
  int nStrata_Control;
  int nStrata_Treatment;

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
  if(methodPeron){
    iFavorable_UTTE.resize(D_UTTE);
    iUnfavorable_UTTE.resize(D_UTTE);
    iWeight_UTTE.resize(D_UTTE);

    if(returnIID>1){
      iDweight_Dnuisance_C_UTTE.resize(D_UTTE);
      iDweight_Dnuisance_T_UTTE.resize(D_UTTE);

      iDscore_Dnuisance_C_UTTE.resize(D_UTTE);
      iDscore_Dnuisance_T_UTTE.resize(D_UTTE);
    }
  }

  // over endpoint
  int iIndex_UTTE_d; // position of the TTE among the UTTE
  int iMethod; // method used to analyze the endpoint

  double iCumWeight; // current weight of the pair
  double iNewWeight; // remaining weight to analyze after the current endpoint

  std::vector< double > iPairScore;

  arma::mat iDscore_Dnuisance_C;
  arma::mat iDscore_Dnuisance_T;
  
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
	}
	if(p_T(iter_strata,iter_d)>0){
	  Dfavorable_Dnuisance_strataT[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),D);
	  Dfavorable_Dnuisance_strataT[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata
	  Dunfavorable_Dnuisance_strataT[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),D);
	  Dunfavorable_Dnuisance_strataT[index_UTTE[iter_d]].fill(0.0); // will keep the sum over all pairs within strata

	  iDscore_Dnuisance_T_UTTE[index_UTTE[iter_d]].resize(p_T(iter_strata,iter_d),4);
	}
      }
    }

    // *** loop over pairs
    if(debug>0){Rcpp::Rcout << " - compute scores" << std::endl;}
    for(unsigned int iter_C=0 ; iter_C < nStrata_Control; iter_C++){
      for(unsigned int iter_T=0 ; iter_T < nStrata_Treatment; iter_T++){
	if(debug>1){Rcpp::Rcout << " pair " << iPair << " (" << iter_C << ";" << iter_T << ") ";}
    	
	// prepare survival
	if(methodPeron){ // 
	  iFavorable_UTTE.fill(0.0);
	  iUnfavorable_UTTE.fill(0.0);
	  iWeight_UTTE.fill(0.0);
	  if(returnIID>1){
	    iDweight_Dnuisance_C_UTTE.clear();
	    iDweight_Dnuisance_T_UTTE.clear();
	    iDweight_Dnuisance_C_UTTE.resize(D_UTTE);
	    iDweight_Dnuisance_T_UTTE.resize(D_UTTE);

	    for(unsigned int iter_d=0; iter_d < D; iter_d ++){
	      if(p_C(iter_strata,iter_d)>0){
		iDscore_Dnuisance_C_UTTE[index_UTTE[iter_d]].fill(0.0);
	      }
	      if(p_T(iter_strata,iter_d)>0){
		iDscore_Dnuisance_T_UTTE[index_UTTE[iter_d]].fill(0.0);
	      }
	    }
	  }
	}

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

	    if(returnIID>1){
	      iDscore_Dnuisance_C.resize(p_C(iter_strata, iter_d),4);
	      iDscore_Dnuisance_T.resize(p_T(iter_strata, iter_d),4);
	    }
	    
	    // note: iDscore_Dnuisance_C, iDscore_Dnuisance_T are initalized to 0 in calcOneScore_SurvPeron
	    iPairScore = calcOneScore_SurvPeron(endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
						endpoint(indexStrataT[iter_T], index_endpoint[iter_d]),
						status(indexStrataC[iter_C], index_status[iter_d]),
						status(indexStrataT[iter_T], index_status[iter_d]),
						threshold[iter_d],
						list_survTimeC[iter_d][iter_strata].row(iter_C), list_survTimeT[iter_d][iter_strata].row(iter_T), list_survJumpC[iter_d][iter_strata], list_survJumpT[iter_d][iter_strata],
						list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1),
						iDscore_Dnuisance_C, iDscore_Dnuisance_T,
						p_C(iter_strata, iter_d), p_T(iter_strata, iter_d), returnIID);

	  }else if(iMethod == 5){  // time to event endpoint with Peron's scoring rule (right-censored, competing risks)
	    iPairScore = calcOnePair_CRPeron(endpoint(indexStrataC[iter_C], index_endpoint[iter_d]),
					     endpoint(indexStrataT[iter_T], index_endpoint[iter_d]),
					     status(indexStrataC[iter_C], index_status[iter_d]),
					     status(indexStrataT[iter_T], index_status[iter_d]),
					     threshold[iter_d],
					     list_survTimeC[iter_d][iter_strata].row(iter_C), list_survTimeT[iter_d][iter_strata].row(iter_T), list_survJumpC[iter_d][iter_strata],					     
					     list_lastSurv[iter_d](iter_strata,0), list_lastSurv[iter_d](iter_strata,1), list_lastSurv[iter_d](iter_strata,2), list_lastSurv[iter_d](iter_strata,3));
	  }
	  
	  // **** remove contribution from previously analyzed threshold of the same endpoint
	  if( (iMethod >= 4) && (nUTTE_analyzedPeron_M1[iter_d]>0) ){
	    iPairScore[0] -= iFavorable_UTTE[iIndex_UTTE_d];
	    iPairScore[1] -= iUnfavorable_UTTE[iIndex_UTTE_d];

	    if(returnIID>1){
	    iDscore_Dnuisance_C.col(0) -= iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d].col(0);
	    iDscore_Dnuisance_C.col(1) -= iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d].col(1);
	    iDscore_Dnuisance_T.col(0) -= iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d].col(0);
	    iDscore_Dnuisance_T.col(1) -= iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d].col(1);
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
	    if( (returnIID > 1) && (iMethod == 4) ){
	      Dfavorable_Dnuisance_strataC[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_C.col(0) * iCumWeight;
	      Dfavorable_Dnuisance_strataT[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_T.col(0) * iCumWeight;
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
	  // }

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
	    if( (returnIID > 1) && (iMethod == 4) ){
	      Dunfavorable_Dnuisance_strataC[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_C.col(1) * iCumWeight;
	      Dunfavorable_Dnuisance_strataT[iIndex_UTTE_d].col(iter_d) += iDscore_Dnuisance_T.col(1) * iCumWeight;
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
	  // }
    
    
	  // **** aggregate neutral score and iid over analyzed pairs
	  if(iPairScore[2] > zeroPlus){
	    if(debug==4){Rcpp::Rcout << "n";}
	    if(debug>4){Rcpp::Rcout << " neutral=" << iPairScore[2] << " ";}
		  
	    // score
	    Mcount_neutral(iter_strata,iter_d) += iPairScore[2] * iCumWeight;

	    // update weight
	    if(neutralAsUninf){
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
	    pairScore[iter_d](iPair,0) = iter_d;
	    pairScore[iter_d](iPair,1) = posStrataC[iter_C];
	    pairScore[iter_d](iPair,2) = posStrataT[iter_T];
	    pairScore[iter_d](iPair,3) = iPair;
	    pairScore[iter_d](iPair,4) = iter_C;
	    pairScore[iter_d](iPair,5) = iter_T;

	    pairScore[iter_d](iPair,6) = iPairScore[0];
	    pairScore[iter_d](iPair,11) = iPairScore[0] * iCumWeight;
	    pairScore[iter_d](iPair,7) = iPairScore[1];
	    pairScore[iter_d](iPair,12) = iPairScore[1] * iCumWeight;
	    pairScore[iter_d](iPair,8) = iPairScore[2];
	    pairScore[iter_d](iPair,13) = iPairScore[2] * iCumWeight;
	    pairScore[iter_d](iPair,9) = iPairScore[3];
	    pairScore[iter_d](iPair,14) = iPairScore[7] * iCumWeight;
	    pairScore[iter_d](iPair,10) = iCumWeight;
	  } 

	  // **** early stop if nothing left or store weight when TTE endpoint with Peron's scoring rule
	  if(hierarchical){
	    if( (iNewWeight < zeroPlus) || (iter_d == (D-1)) ){
	      if(debug>3){Rcpp::Rcout << " exit ";}
	      break;
	    }else if(methodPeron && (iIndex_UTTE_d>=0) ){
	      if(debug>3){Rcpp::Rcout << " store ";}
	      iFavorable_UTTE[iIndex_UTTE_d] += iPairScore[0];
	      iUnfavorable_UTTE[iIndex_UTTE_d] += iPairScore[1];
	      iWeight_UTTE[iIndex_UTTE_d] = iNewWeight;

	      if(returnIID>1){
		// Rcpp::Rcout  << iDscore_Dnuisance_C_UTTE[0] << std::endl;
		iDscore_Dnuisance_C_UTTE[iIndex_UTTE_d] += iDscore_Dnuisance_C;
		iDscore_Dnuisance_T_UTTE[iIndex_UTTE_d] += iDscore_Dnuisance_T;
		// Rcpp::Rcout  << iDscore_Dnuisance_C_UTTE[0] << std::endl;
		
		if(neutralAsUninf && (iPairScore[2] > zeroPlus)){
		  iDweight_Dnuisance_C_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_C.col(2)+iDscore_Dnuisance_C.col(3);
		  iDweight_Dnuisance_T_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_T.col(2)+iDscore_Dnuisance_T.col(3);
		}else{
		  iDweight_Dnuisance_C_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_C.col(3);
		  iDweight_Dnuisance_T_UTTE[iIndex_UTTE_d] = iDscore_Dnuisance_T.col(3);
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
	for(int iter_UTTE=0 ; iter_UTTE<D_UTTE; iter_UTTE++){
	  iidNuisance_favorable.col(iter_d) += iid_survJumpC[iter_UTTE][iter_strata] * Dfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_favorable.col(iter_d) += iid_survJumpT[iter_UTTE][iter_strata] * Dfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpC[iter_UTTE][iter_strata] * Dunfavorable_Dnuisance_strataC[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
    	  iidNuisance_unfavorable.col(iter_d) += iid_survJumpT[iter_UTTE][iter_strata] * Dunfavorable_Dnuisance_strataT[iter_UTTE].col(iter_d)/vecn_pairs[iter_strata];
	  }
      }
    }

  }

  
  // ** proportion in favor of treatment
  // Rcpp::Rcout << std::endl << " compute statistics" << std::endl;
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
  return(Rcpp::List::create(Rcpp::Named("count_favorable") = Mcount_favorable,
			    Rcpp::Named("count_unfavorable") = Mcount_unfavorable,
			    Rcpp::Named("count_neutral") = Mcount_neutral,           
			    Rcpp::Named("count_uninf") = Mcount_uninf,
			    Rcpp::Named("delta_netBenefit") = delta_netBenefit,
			    Rcpp::Named("delta_winRatio") = delta_winRatio,
			    Rcpp::Named("Delta_netBenefit") = arma::conv_to< std::vector<double> >::from(Delta_netBenefit),
			    Rcpp::Named("Delta_winRatio") = arma::conv_to< std::vector<double> >::from(Delta_winRatio),
			    Rcpp::Named("n_pairs") = arma::conv_to< std::vector<double> >::from(vecn_pairs),
			    Rcpp::Named("iidAverage_favorable") = iidAverage_favorable,
			    Rcpp::Named("iidAverage_unfavorable") = iidAverage_unfavorable,
			    Rcpp::Named("iidNuisance_favorable") = iidNuisance_favorable,
			    Rcpp::Named("iidNuisance_unfavorable") = iidNuisance_unfavorable,
			    Rcpp::Named("Mvar") = Mvar,
			    Rcpp::Named("tableScore")  = pairScore
			    ));
}
