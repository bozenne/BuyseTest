// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:
using namespace Rcpp ;
using namespace std ;
using namespace arma ;

arma::mat calcAllPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
					   arma::colvec deltaC, arma::colvec deltaT, 
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   double lastSurvC, double lastSurvT, 
					   int method, int correctionUninf, int p_C, int p_T,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   std::vector< int >& index_control, std::vector< int >& index_treatment, arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					   bool neutralAsUninf, bool keepScore, bool moreEndpoint, bool reserve);

arma::mat calcSubsetPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
						  arma::colvec deltaC, arma::colvec deltaT, 
						  arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
						  double lastSurvC, double lastSurvT, 
						  int method, int correctionUninf, int p_C, int p_T,
						  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
						  std::vector< int >& index_control, std::vector< int >& index_treatment,  arma::vec& RP_score, arma::uvec& index_RP,
						  arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						  std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
						  bool neutralAsUninf, bool keepScore, bool moreEndpoint, bool reserve,
						  const std::vector< int >& index_control_M1, const vector<int>& index_treatment_M1, const arma::mat& RP_score_M1,
						  std::vector< arma::mat >& RP_Dscore_Dnuisance_C_M1, std::vector< arma::mat >& RP_Dscore_Dnuisance_T_M1);

void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 const std::vector< int >& index_uninfC, const std::vector< int >& index_uninfT, 
					 const std::vector< int >& index_neutralC, const std::vector< int >& index_neutralT,
					 const std::vector< double >& wNeutral, const std::vector< int >& index_wNeutral, std::vector< double >& wUninf, const std::vector< int >& index_wUninf,
					 std::vector< int >& index_control, std::vector< int >& index_treatment, arma::vec& weight, arma::uvec& index_weight,
					 arma::mat& MC_iid, arma::mat& MT_iid, int returnIID,
					 bool firstEndpoint, bool neutralAsUninf, bool moreEndpoint, bool keepScore, arma::mat& matPairScore);

void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   const std::vector< int >& index_neutralC, const std::vector< int >& index_neutralT,
				   const std::vector< double >& wNeutral, const std::vector< int >& index_wNeutral, 
				   std::vector< int >& index_control, std::vector< int >& index_treatment, arma::vec& weight, arma::uvec& index_weight,
				   arma::mat& MC_iid, arma::mat& MT_iid, int returnIID,
				   bool firstEndpoint, bool neutralAsUninf, bool moreEndpoint, bool keepScore, arma::mat& matPairScore);

// * calcAllPairs
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
					   arma::colvec deltaC, arma::colvec deltaT, 
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   double lastSurvC, double lastSurvT, 
					   int method, int correctionUninf, int p_C, int p_T,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					   bool neutralAsUninf, bool keepScore, bool moreEndpoint, bool reserve){

  // ** initialize
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_pair = n_Treatment * n_Control;

  bool updateIndexNeutral = moreEndpoint & neutralAsUninf;
  bool updateIndexUninf = moreEndpoint & (correctionUninf==0);
  bool iUpdateRPNeutral;
  bool iUpdateRPUninf;
  
  // residual pairs (neutral and uninformative)
  int iRP=0; // count the number of RP
  rowvec iRP_score(6);
  if(reserve){
    RP_score.resize(n_pair,6);
  }

  // iid
  if(returnIID > 0){ // 3 because neutral column used when correction
    count_obsC.resize(n_Control, 3); 
    count_obsC.fill(0.0);

    count_obsT.resize(n_Treatment, 3);
    count_obsT.fill(0.0);
  }

  // iid nuisance
  arma::mat iDscore_Dnuisance_C;
  arma::mat iDscore_Dnuisance_T;
  if(returnIID > 1 && method == 3){
    Dscore_Dnuisance_C.resize(p_C, 2);
    Dscore_Dnuisance_C.fill(0.0);

    Dscore_Dnuisance_T.resize(p_T, 2);
    Dscore_Dnuisance_T.fill(0.0);

    iDscore_Dnuisance_C.resize(p_C, 2); // initialized in calcOneScore_TTEperon
	iDscore_Dnuisance_T.resize(p_T, 2); // initialized in calcOneScore_TTEperon

	for(int iter_typeRP=0; iter_typeRP<3; iter_typeRP++){
	  RP_Dscore_Dnuisance_C[iter_typeRP].resize(0,0);
	  RP_Dscore_Dnuisance_T[iter_typeRP].resize(0,0);
	}
  }else{
    Dscore_Dnuisance_C.resize(0, 0);
    Dscore_Dnuisance_T.resize(0, 0);
  }
  
  // pairScore    
  std::vector< double > iPairScore(4); // temporary store results
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }

  // other
  double zeroPlus = pow(10.0,-12.0);
  
  // ** loop over the pairs
  int iter_pair = 0;
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
      // Rcout << iter_pair << endl;
	  iUpdateRPNeutral = false;
	  iUpdateRPUninf = false;
	  
      // score
      if(method == 1){
		iPairScore = calcOnePair_Continuous(Treatment[iter_T] - Control[iter_C], threshold);
      }else if(method == 2){
		iPairScore = calcOnePair_TTEgehan(Treatment[iter_T] - Control[iter_C], deltaC[iter_C], deltaT[iter_T], threshold);
      }else if(method == 3){
		iPairScore = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
										   deltaC[iter_C], deltaT[iter_T], threshold,
										   survTimeC.row(iter_C), survTimeT.row(iter_T),
										   survJumpC, survJumpT, lastSurvC, lastSurvT,
										   iDscore_Dnuisance_C, iDscore_Dnuisance_T, returnIID);
      }	

      // store results
      if(iPairScore[0] > zeroPlus){
		count_favorable += iPairScore[0];
		if(returnIID > 0){
		  count_obsC(iter_C,0) += iPairScore[0];
		  count_obsT(iter_T,0) += iPairScore[0];
		  if(returnIID > 1 && method == 3){
		  	Dscore_Dnuisance_C.col(0) += iDscore_Dnuisance_C.col(0);
		  	Dscore_Dnuisance_T.col(0) += iDscore_Dnuisance_T.col(0);
		  }
		}
      }

      if(iPairScore[1] > zeroPlus){
		count_unfavorable += iPairScore[1];
		if(returnIID > 0){
		  count_obsC(iter_C,1) += iPairScore[1];
		  count_obsT(iter_T,1) += iPairScore[1];
		  if(returnIID > 1 && method == 3){
		  	Dscore_Dnuisance_C.col(1) += iDscore_Dnuisance_C.col(1);
		  	Dscore_Dnuisance_T.col(1) += iDscore_Dnuisance_T.col(1);
		  }
		}
      }

      if(iPairScore[2] > zeroPlus){
		count_neutral += iPairScore[2];
		iUpdateRPNeutral = true;
      }

	  if(iPairScore[3] > zeroPlus){
		count_uninf += iPairScore[3];
		if(returnIID > 0){ // used in the correction at the pair level
		  count_obsC(iter_C,2) += iPairScore[3];
		  count_obsT(iter_T,2) += iPairScore[3];
		}
		iUpdateRPUninf = true;
      }

	  if(iUpdateRPNeutral || iUpdateRPUninf){
		// keep score
		iRP_score[0] = iter_pair;
		iRP_score[1] = iter_C;
		iRP_score[2] = iter_T;
		iRP_score[3] = iPairScore[0];
		iRP_score[4] = iPairScore[1];
		if(iUpdateRPNeutral){
		  iRP_score[5] = iPairScore[2];
		}else{
		  iRP_score[5] = 0;
		}
		if(iUpdateRPUninf){
		  iRP_score[6] = iPairScore[3];
		}else{
		  iRP_score[6] = 0;
		}
		RP_score.insert_rows(iRP,iRP_score);

		// keep iid
		if(returnIID > 1 && method == 3){
		  for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
			RP_Dscore_Dnuisance_C[0][iter_typeRP].insert_rows(iRP,iDscore_Dnuisance_C);
			RP_Dscore_Dnuisance_T[0][iter_typeRP].insert_rows(iRP,iDscore_Dnuisance_T);
		  }
		}
		iRP++;
	  }

      if(keepScore){
		matPairScore.row(iter_pair) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
			  iPairScore[0], // favorable
			  iPairScore[1], // unfavorable
			  iPairScore[2], // neutral
			  iPairScore[3], // uninformative
			  1.0, // weight
			  iPairScore[0], iPairScore[1], iPairScore[2], iPairScore[3] // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
			  });		
      }

      if(iter_pair % 65536 == 0){
		R_CheckUserInterrupt();
      }

      iter_pair++;
    }    
  }
  
  R_CheckUserInterrupt();

  // ** merge neutral and uninformative pairs + correction
  bool firstEndpoint = true;
  arma::uvec index_weight; // empty, just for calling functions
  std::vector< int > index_wNeutral; // empty, just for calling functions
  std::vector< int > index_wUninf; // empty, just for calling functions

  // Rcout << "start merge " << endl;

    // correction possible: if there are uninformative paris and if there are informative pairs
  if((correctionUninf == 1) && (count_uninf > 0)){
      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  count_obsC, count_obsT, returnIID,
					  firstEndpoint, neutralAsUninf, moreEndpoint, keepScore, matPairScore);
      // update by reference
      
    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_neutralC, index_neutralT,
					wNeutral, index_wNeutral,
					index_control, index_treatment, weight, index_weight,
					count_obsC, count_obsT, returnIID,
					firstEndpoint, neutralAsUninf, moreEndpoint, keepScore, matPairScore);
      // updated by reference      
    }
  // Rcout << "end merge " << endl;

  // ** export
  return matPairScore;
  
}


// * calcSubsetPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a given endpoint 
arma::mat calcSubsetPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
						  arma::colvec deltaC, arma::colvec deltaT, 
						  arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
						  double lastSurvC, double lastSurvT, 
						  int method, int correctionUninf, int p_C, int p_T,
						  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
						  std::vector< int >& index_control, std::vector< int >& index_treatment, arma::vec& weight, arma::uvec& index_weight,
 						  arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						  std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& ls_matIIDC, int returnIID,					
						  bool neutralAsUninf, bool keepScore, bool moreEndpoint, bool reserve,
						  const std::vector< int >& index_control_M1, const vector<int>& index_treatment_M1, const arma::vec& cumWeight_M1,
						  double threshold_M1, const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat survJumpC_M1, const arma::mat& survJumpT_M1, double lastSurvC_M1, double lastSurvT_M1
						  ){
  // Rcout << "start calcSubsetPairs " << endl;
  
  // ** initialize
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_pair = index_control_M1.size(); // n_pair may not be n_Control * n_Treatment since some pairs may have been informative regarding the previous endpoint 

  bool updateIndexNeutral = moreEndpoint && neutralAsUninf;
  bool updateIndexUninf = moreEndpoint;

  // neutral
  std::vector< int > index_neutralC(0); // index of the neutral pairs relative to Control
  std::vector< int > index_neutralT(0); // index of the neutral pairs relative to Treatment
  std::vector< int > index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  std::vector< double > wNeutral(0); // weight of the neutral pairs
  if(updateIndexNeutral && reserve){
    index_neutralC.reserve(n_pair);
    index_neutralT.reserve(n_pair);
    index_wNeutral.reserve(n_pair);
    wNeutral.reserve(n_pair);
  }

  // uninf
  std::vector< int > index_uninfC(0); // index of the uninformative pairs relative to Control
  std::vector< int > index_uninfT(0); // index of the uninformative pairs relative to Treatment
  std::vector< int > index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  std::vector< double > wUninf(0); // weight of the uninformative pairs
  if(updateIndexUninf && reserve){
    index_uninfC.reserve(n_pair);
    index_uninfT.reserve(n_pair);
    index_wUninf.reserve(n_pair);
    wUninf.reserve(n_pair);
  }

  // iid
  if(returnIID > 0){
    count_obsC.resize(n_Control, 3);
    count_obsC.fill(0.0);

    count_obsT.resize(n_Treatment, 3);
    count_obsT.fill(0.0);
  }

  // iid nuisance
  arma::mat iDscore_Dnuisance_C;
  arma::mat iDscore_Dnuisance_T;
  arma::mat iDscore_Dnuisance_C_M1;
  arma::mat iDscore_Dnuisance_T_M1;
  if(returnIID > 1){
    Dscore_Dnuisance_C.resize(p_C, 2);
	Dscore_Dnuisance_C.fill(0.0);

	Dscore_Dnuisance_T.resize(p_T, 2);
    Dscore_Dnuisance_T.fill(0.0);
	
	if(method == 3){
	  iDscore_Dnuisance_C.resize(p_C, 2);
	  // initialization done inside calcOneScore_TTEperon
	  iDscore_Dnuisance_C_M1.resize(p_C, 2);

	  iDscore_Dnuisance_T.resize(p_T, 2);
	  // initialization done inside calcOneScore_TTEperon
	  iDscore_Dnuisance_T_M1.resize(p_T, 2);
	}
  }else{
    Dscore_Dnuisance_C.resize(0, 0);
    Dscore_Dnuisance_T.resize(0, 0);
  }

  // store the score of all pairs
  std::vector< double > iPairScore(4); // temporary store results
  std::vector< double > iPairScore_M1(2); // temporary store results
  std::fill(iPairScore_M1.begin(), iPairScore_M1.end(), 0.0);  
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }else{
    matPairScore.resize(0, 0); 
  }
  double weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative;
  double zeroPlus = pow(10.0,-12.0);
  bool alreadyAnalyzed = (threshold_M1>=0);

  // Rcout << "threshold: " << threshold << " | " << "threshold_M1: " << threshold_M1 << " (" << alreadyAnalyzed << ") " <<endl;
  
  // ** loop over the pairs
  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){
    // Rcout << iter_pair << "/" << n_pair << " ";
	
    // *** find index of the pair
    iter_C = index_control_M1[iter_pair];
    iter_T = index_treatment_M1[iter_pair];
	// Rcout << "("<< iter_C << ";" << iter_T << ") ";
	
    // *** score pair
    if(method == 1){
      iPairScore = calcOnePair_Continuous(Treatment[iter_T] - Control[iter_C], threshold);
    }else if(method == 2){
      iPairScore = calcOnePair_TTEgehan(Treatment[iter_T] - Control[iter_C], deltaC[iter_C], deltaT[iter_T], threshold);
    }else if(method == 3){
      iPairScore = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
									 deltaC[iter_C], deltaT[iter_T], threshold,
									 survTimeC.row(iter_C), survTimeT.row(iter_T),
									 survJumpC, survJumpT, lastSurvC, lastSurvT,
									 iDscore_Dnuisance_C, iDscore_Dnuisance_T, returnIID);

      if(alreadyAnalyzed){
		iPairScore_M1 = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
										  deltaC[iter_C], deltaT[iter_T], threshold_M1,
										  survTimeC_M1.row(iter_C), survTimeT_M1.row(iter_T),
										  survJumpC_M1, survJumpT_M1, lastSurvC_M1, lastSurvT_M1,
										  iDscore_Dnuisance_C_M1, iDscore_Dnuisance_T_M1, returnIID);
      }
    }

    // *** compute adjusted for the previous endpoint
    weight_favorable = (iPairScore[0] - iPairScore_M1[0]) * cumWeight_M1(iter_pair);
    weight_unfavorable = (iPairScore[1] - iPairScore_M1[1]) * cumWeight_M1(iter_pair);
    weight_neutral = iPairScore[2] * cumWeight_M1(iter_pair); 
    weight_uninformative = iPairScore[3] * cumWeight_M1(iter_pair);
    // Rcout << cumWeight_M1(iter_pair) << ") " << weight_favorable << " " << weight_unfavorable << " "<< weight_neutral << " " << weight_uninformative << endl;
    // Rcout << "...) " << iPairScore[0] << " " << iPairScore[1] << " "<< iPairScore[2] << " " << iPairScore[3] << endl;

    // *** store results
    if(iPairScore[0] > zeroPlus){
      count_favorable += weight_favorable;
      if(returnIID > 0){
		count_obsC(iter_C,0) += weight_favorable;
		count_obsT(iter_T,0) += weight_favorable;
		if(returnIID > 1 && method == 3){
		  if(iPairScore_M1[1] > zeroPlus){
			Dscore_Dnuisance_C.col(0) += (iDscore_Dnuisance_C.col(0) - iDscore_Dnuisance_C_M1.col(0)) * cumWeight_M1(iter_pair);
			Dscore_Dnuisance_T.col(0) += (iDscore_Dnuisance_T.col(0) - iDscore_Dnuisance_T_M1.col(0)) * cumWeight_M1(iter_pair);
		  }else{
			Dscore_Dnuisance_C.col(0) += iDscore_Dnuisance_C.col(0) * cumWeight_M1(iter_pair);
			Dscore_Dnuisance_T.col(0) += iDscore_Dnuisance_T.col(0) * cumWeight_M1(iter_pair);
		  }
		}
	  }
	}
	
    if(iPairScore[1] > zeroPlus){
      count_unfavorable += weight_unfavorable;
      if(returnIID > 0){
		count_obsC(iter_C,1) += weight_unfavorable;
		count_obsT(iter_T,1) += weight_unfavorable;
		if(returnIID > 1){
		  if(iPairScore_M1[1] > zeroPlus && method == 3){
			Dscore_Dnuisance_C.col(1) += (iDscore_Dnuisance_C.col(1) - iDscore_Dnuisance_C_M1.col(1)) * cumWeight_M1(iter_pair);
			Dscore_Dnuisance_T.col(1) += (iDscore_Dnuisance_T.col(1) - iDscore_Dnuisance_T_M1.col(1)) * cumWeight_M1(iter_pair);
		  }else{
			Dscore_Dnuisance_C.col(1) += iDscore_Dnuisance_C.col(1) * cumWeight_M1(iter_pair);
			Dscore_Dnuisance_T.col(1) += iDscore_Dnuisance_T.col(1) * cumWeight_M1(iter_pair);
		  }
		}
      }
    }

    if(weight_neutral > zeroPlus){
      count_neutral += weight_neutral;
      if(updateIndexNeutral){
		index_neutralC.push_back(iter_C); // index of the pair relative to Control         
		index_neutralT.push_back(iter_T); // index of the pair relative to Treatment
		index_wNeutral.push_back(iter_pair); // index of the pair relative to cumWeight_M1
		wNeutral.push_back(iPairScore[2]); // not weight_neutral since the product is done in BuyseTest.cpp
		if(returnIID > 1 && method == 3){
		}
      }
    }
    if(weight_uninformative > zeroPlus){
      count_uninf += weight_uninformative;
	  if(returnIID > 0){ // used in the correction at the pair level
		count_obsC(iter_C,2) += weight_uninformative;
		count_obsT(iter_T,2) += weight_uninformative;
	  }
      if(updateIndexUninf){
		index_uninfC.push_back(iter_C); // index of the pair relative to Control    
		index_uninfT.push_back(iter_T); // index of the pair relative to Treatment
		index_wUninf.push_back(iter_pair); // index of the pair relative to cumWeight_M1
		wUninf.push_back(iPairScore[3]); // not weight_uninformative since the product is done in BuyseTest.cpp
		if(returnIID > 1 && method == 3){
		}
      }
    }

    if(keepScore){
      matPairScore.row(iter_pair) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
					    iPairScore[0] - iPairScore_M1[0], // favorable
					    iPairScore[1] - iPairScore_M1[1], // unfavorable
					    iPairScore[2], // neutral
					    iPairScore[3], // uninformative
					    cumWeight_M1(iter_pair), // weight
					    weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
	});
    }

    if(iter_pair % 65536 == 0){
      R_CheckUserInterrupt();
    }

  }

  R_CheckUserInterrupt();

  // ** merge neutral and uninformative pairs + correction
  bool firstEndpoint = false;
  // Rcout << "start merge " << endl;
      
  if(correctionUninf > 0 && count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
    // correction possible: if there are uninformative paris and if there are informative pairs

    if(correctionUninf == 1){
      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfC, index_uninfT, 
		      index_neutralC, index_neutralT, 
		      wNeutral, index_wNeutral, wUninf, index_wUninf,
		      index_control, index_treatment, weight, index_weight,
		      count_obsC, count_obsT, returnIID,
		      firstEndpoint, neutralAsUninf, moreEndpoint, keepScore, matPairScore);
      // updated by reference
    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		    index_neutralC, index_neutralT,
		    wNeutral, index_wNeutral,
		    index_control, index_treatment, weight, index_weight,
		    count_obsC, count_obsT, returnIID,
		    firstEndpoint, neutralAsUninf, moreEndpoint, keepScore, matPairScore);
      // updated by reference      
    }
  }else{

    if(moreEndpoint){
      noCorrection(index_uninfC, index_uninfT, 
		   index_neutralC, index_neutralT,
		   wNeutral, index_wNeutral, wUninf, index_wUninf,
		   index_control, index_treatment, weight, index_weight,
		   method, firstEndpoint, neutralAsUninf);
    }
    // Rcout << " | " << index_control.size() << " " << index_treatment.size() << " " << weight.size() << " " << index_weight.size() << endl;
  }
  // Rcout << "end merge " << endl;

  // ** export 
  return matPairScore;
  
}

// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, bool moreEndpoint,
					 arma::mat& count_obsC, arma::mat& count_obsT, int returnIID,
					 bool keepScore, arma::mat& matPairScore){

  // compute factor
  double factorFavorable;
  double factorUnfavorable;
  double factorNeutral;
  if(factorFavorable + factorUnfavorable + factorNeutral>0){
	factorFavorable = (count_favorable)/(count_favorable + count_unfavorable + count_neutral); 
	factorUnfavorable = (count_unfavorable)/(count_favorable + count_unfavorable + count_neutral); 
	factorNeutral  = (count_neutral)/(count_favorable + count_unfavorable + count_neutral); 
  }else{
	factorFavorable = 1/3;
	factorUnfavorable = 1/3;
	factorNeutral = 1/3;
  }
  
  // update global score
  count_favorable += factorFavorable * count_uninf;    
  count_unfavorable += factorUnfavorable * count_uninf;    
  count_neutral += factorNeutral * count_uninf;   
  count_uninf = 0;

  if(returnIID > 0){
    count_obsC.col(0) += factorFavorable * count_obsC.col(2);
    count_obsC.col(1) += factorUnfavorable * count_obsC.col(2);
    count_obsT.col(0) += factorFavorable * count_obsT.col(2);
    count_obsT.col(1) += factorUnfavorable * count_obsT.col(2);
  }
  
  // update pair score
  if(moreEndpoint){
	RP_score.col(3) = factorFavorable * RP_score.col(6);
	RP_score.col(4) = factorUnfavorable * RP_score.col(6);
	RP_score.col(5) = factorNeutral * RP_score.col(6);
	RP_score.col(6) = 0;
  }
	
  // update keep scores
  if(keepScore){
    matPairScore.col(7) = matPairScore.col(7) + factorFavorable * matPairScore.col(10);
    matPairScore.col(8) = matPairScore.col(8) + factorUnfavorable * matPairScore.col(10);
    matPairScore.col(9) = matPairScore.col(9) + factorNeutral * matPairScore.col(10);
    (matPairScore.col(10)).fill(0.0);
  }
}

// * correctionIPW
void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		   const std::vector< int >& index_neutralC, const std::vector< int >& index_neutralT,
		   const std::vector< double >& wNeutral, const std::vector< int >& index_wNeutral, 
		   std::vector< int >& index_control, std::vector< int >& index_treatment, arma::vec& weight, arma::uvec& index_weight,
		   arma::mat& count_obsC, arma::mat& count_obsT, int returnIID,
		   bool firstEndpoint, bool neutralAsUninf, bool moreEndpoint, bool keepScore, arma::mat& matPairScore){

  // compute factor
  double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
  // update global score
  count_favorable *= factor;    
  count_unfavorable *= factor;    
  count_neutral *= factor;   
  count_uninf = 0;

  if(returnIID > 0){
    count_obsC *= factor;
    count_obsT *= factor;
  }
  
  // new index/weights
  if(moreEndpoint && neutralAsUninf){
    int size = wNeutral.size();
    
    weight.resize(size);
    for(int iIndex=0; iIndex < size; iIndex++){
      weight(iIndex) = wNeutral[iIndex] * factor;
    }

    index_control = index_neutralC;
    index_treatment = index_neutralT;
    if(firstEndpoint==false){
      index_weight = arma::conv_to<uvec>::from(index_wNeutral);
    }
  }

  // update keep scores
  if(keepScore){
    matPairScore.col(7) = matPairScore.col(7) * factor;
    matPairScore.col(8) = matPairScore.col(8) * factor;
    matPairScore.col(9) = matPairScore.col(9) * factor;
    (matPairScore.col(10)).fill(0.0);
  }
}


