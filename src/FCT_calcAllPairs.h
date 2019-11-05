// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:
using namespace Rcpp;
using namespace std;
using namespace arma;

arma::mat calcAllPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
					   arma::colvec deltaC, arma::colvec deltaT, 
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   double lastSurvC, double lastSurvT, 
					   int method, double zeroPlus, int correctionUninf, int p_C, int p_T,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					   bool neutralAsUninf, bool keepScore, bool moreEndpoint, int reserve);

arma::mat calcSubsetPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
						  arma::colvec deltaC, arma::colvec deltaT, 
						  arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
						  double lastSurvC, double lastSurvT, 
						  int method, double zeroPlus, int correctionUninf, int p_C, int p_T,
						  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
						  arma::mat& RP_score, 
						  arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						  std::vector<arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
						  bool neutralAsUninf, bool keepScore, bool moreEndpoint, int reserve,
						  const arma::vec& index_control, const arma::vec& index_treatment, const arma::vec& weight,
						  const arma::mat& RP_score_M1, const std::vector<int>& activeUTTE, 
						  std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
						  const std::vector<arma::mat >& RP_Dscore_Dnuisance_C_M1, const std::vector< arma::mat >& RP_Dscore_Dnuisance_T_M1);

void correctionPairs(int method, double zeroPlus,
					 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, arma::mat& matPairScore,
					 arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					 const vector<int>& activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					 bool neutralAsUninf, bool keepScore, bool moreEndpoint);

void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   const vector<int>& activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool moreEndpoint);

// * calcAllPairs
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
					   arma::colvec deltaC, arma::colvec deltaT, 
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   double lastSurvC, double lastSurvT, 
					   int method, double zeroPlus, int correctionUninf, int p_C, int p_T,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					   bool neutralAsUninf, bool keepScore, bool moreEndpoint, int reserve){

  // ** initialize
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_pair = n_Treatment * n_Control;

  bool updateIndexNeutral = moreEndpoint & neutralAsUninf;
  bool updateIndexUninf = moreEndpoint & (correctionUninf < 2); // if correctionUninf == 2, only keep pair when non-0 neutral score
  bool iUpdateRPNeutral;
  bool iUpdateRPUninf;
  
  // residual pairs (neutral and uninformative)
  int n_RP=0;
  vector<int> vec_indexPair(0);
  vector<int> vec_indexC(0);
  vector<int> vec_indexT(0);
  vector<double> vec_favorable(0);
  vector<double> vec_unfavorable(0);
  vector<double> vec_neutral(0);
  vector<double> vec_uninformative(0);
  
  // iid
  if(returnIID > 0){ 
    count_obsC.resize(n_Control, 4); 
    count_obsC.fill(0.0);

    count_obsT.resize(n_Treatment, 4);
    count_obsT.fill(0.0);
  }

  // iid nuisance
  arma::mat iDscore_Dnuisance_C;
  arma::mat iDscore_Dnuisance_T;
  int space = 0;
  if(returnIID > 1 && method == 3){
    Dscore_Dnuisance_C.resize(p_C, 4);
    Dscore_Dnuisance_C.fill(0.0);

    Dscore_Dnuisance_T.resize(p_T, 4);
    Dscore_Dnuisance_T.fill(0.0);

    iDscore_Dnuisance_C.resize(p_C, 4); // initialized in calcOneScore_TTEperon
	iDscore_Dnuisance_T.resize(p_T, 4); // initialized in calcOneScore_TTEperon
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

      // *** store results
      if(iPairScore[0] > zeroPlus){
		// favorable score
		count_favorable += iPairScore[0];

		// iid(favorable score)
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
		// unfavorable score
		count_unfavorable += iPairScore[1];

		// iid(unfavorable score)
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
		// neutral score
		count_neutral += iPairScore[2];
		iUpdateRPNeutral = updateIndexNeutral;

		// iid(neutral score)
		if(returnIID > 0){ // used in the correction at the pair level
		  count_obsC(iter_C,2) += iPairScore[2];
		  count_obsT(iter_T,2) += iPairScore[2];
		  if(returnIID > 1 && method == 3){
		  	Dscore_Dnuisance_C.col(2) += iDscore_Dnuisance_C.col(2);
		  	Dscore_Dnuisance_T.col(2) += iDscore_Dnuisance_T.col(2);
		  }
		}
      }

	  if(iPairScore[3] > zeroPlus){
		// uninformative score
		count_uninf += iPairScore[3];
		iUpdateRPUninf = updateIndexUninf;

		// iid(uninformative score)
		if(returnIID > 0){ // used in the correction at the pair level
		  count_obsC(iter_C,3) += iPairScore[3];
		  count_obsT(iter_T,3) += iPairScore[3];
		  if(returnIID > 1 && method == 3){
		  	Dscore_Dnuisance_C.col(3) += iDscore_Dnuisance_C.col(3);
		  	Dscore_Dnuisance_T.col(3) += iDscore_Dnuisance_T.col(3);
		  }
		}
      }

	  // keep information relative to the residual pairs
	  if(iUpdateRPNeutral || iUpdateRPUninf){
		vec_indexPair.push_back(iter_pair);
		vec_indexC.push_back(iter_C);
		vec_indexT.push_back(iter_T);
		vec_favorable.push_back(iPairScore[0]);
		vec_unfavorable.push_back(iPairScore[1]);
		vec_neutral.push_back(iPairScore[2]);
		vec_uninformative.push_back(iPairScore[3]);

		if(returnIID > 1 && method == 3){
		  
		  for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
			if(space == 0){
			  RP_Dscore_Dnuisance_C[iter_typeRP].resize(n_RP+reserve,p_C);
			  RP_Dscore_Dnuisance_T[iter_typeRP].resize(n_RP+reserve,p_T);
			  space = reserve;
			}
		
			RP_Dscore_Dnuisance_C[iter_typeRP].row(iter_pair) = trans(iDscore_Dnuisance_C.col(iter_typeRP));
			RP_Dscore_Dnuisance_T[iter_typeRP].row(iter_pair) = trans(iDscore_Dnuisance_T.col(iter_typeRP));
			space --;
		  }
		}
		n_RP++;
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
  if(returnIID > 1 && method == 3 && space > 0){
	for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
	  RP_Dscore_Dnuisance_C[iter_typeRP].resize(n_RP,p_C);
	  RP_Dscore_Dnuisance_T[iter_typeRP].resize(n_RP,p_T);
	}
  }
  
  // ** merge information about the residual pairs
  // Rcout << "start merge " << endl;
  RP_score.resize(n_RP,7);
  if(n_RP>0){
	RP_score.col(0) = conv_to<colvec>::from(vec_indexPair);
	RP_score.col(1) = conv_to<colvec>::from(vec_indexC);
	RP_score.col(2) = conv_to<colvec>::from(vec_indexT);
	RP_score.col(3) = conv_to<colvec>::from(vec_favorable);
	RP_score.col(4) = conv_to<colvec>::from(vec_unfavorable);
	RP_score.col(5) = conv_to<colvec>::from(vec_neutral);
	RP_score.col(6) = conv_to<colvec>::from(vec_uninformative);
  }
  // Rcout << "end merge " << endl;
  
  // ** correction for uninformative pairs
  // (updated by reference)
  // Rcout << "start correction " << endl;
  
  if(count_uninf > 0){
	vector<int> activeUTTE;
	std::vector<std::vector< arma::mat >> Dweight_Dnuisance_C;
	std::vector<std::vector< arma::mat >> Dweight_Dnuisance_T;
	  
	if(correctionUninf == 1){
      correctionPairs(method, zeroPlus,
					  count_favorable, count_unfavorable, count_neutral, count_uninf,
					  RP_score, matPairScore,
					  count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					  activeUTTE, 0, Dweight_Dnuisance_C, Dweight_Dnuisance_T,
					  RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					  neutralAsUninf, keepScore, moreEndpoint);
      
    }else if(correctionUninf == 2){
	  correctionIPW(method, zeroPlus,
					count_favorable, count_unfavorable, count_neutral, count_uninf,
					RP_score, matPairScore,
					count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					activeUTTE, 0, Dweight_Dnuisance_C, Dweight_Dnuisance_T,
					RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					neutralAsUninf, keepScore, moreEndpoint);
    }
  }
  // Rcout << "end correction " << endl;

  // ** export
  return matPairScore;
  
}


// * calcSubsetPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a given endpoint 
arma::mat calcSubsetPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
						  arma::colvec deltaC, arma::colvec deltaT, 
						  arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
						  double lastSurvC, double lastSurvT, 
						  int method, double zeroPlus, int correctionUninf, int p_C, int p_T,
						  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
						  arma::mat& RP_score, 
						  arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						  std::vector<arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
						  bool neutralAsUninf, bool keepScore, bool moreEndpoint, int reserve,
						  const arma::vec& index_control, const arma::vec& index_treatment, const arma::vec& weight,
						  const arma::mat& RP_score_M1, const std::vector<int>& activeUTTE, int D_activeUTTE,
						  std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
						  const std::vector<arma::mat >& RP_Dscore_Dnuisance_C_M1, const std::vector< arma::mat >& RP_Dscore_Dnuisance_T_M1){
  // Rcout << "start calcSubsetPairs " << endl;
  
  // ** initialize
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_pair = weight.size(); // n_pair may not be n_Control * n_Treatment since some pairs may have been informative regarding the previous endpoint 
  bool alreadyAnalyzed = (RP_score_M1.n_rows > 0);
  
  bool updateIndexNeutral = moreEndpoint & neutralAsUninf;
  bool updateIndexUninf = moreEndpoint & (correctionUninf < 2); // if correctionUninf == 2, only keep pair when non-0 neutral score
  bool iUpdateRPNeutral = false; // (initialize to avoid C++warnings)
  bool iUpdateRPUninf = false; // (initialize to avoid C++warnings)
  
  // residual pairs (neutral and uninformative)
  int n_RP=0;
  vector<int> vec_indexPair(0);
  vector<int> vec_indexC(0);
  vector<int> vec_indexT(0);
  vector<double> vec_favorable(0);
  vector<double> vec_unfavorable(0);
  vector<double> vec_neutral(0);
  vector<double> vec_uninformative(0);

  // iid
  if(returnIID > 0){ 
    count_obsC.resize(n_Control, 4); 
    count_obsC.fill(0.0);

    count_obsT.resize(n_Treatment, 4);
    count_obsT.fill(0.0);
  }

  // iid nuisance
  arma::mat iDscore_Dnuisance_C;
  arma::mat iDscore_Dnuisance_T;
  arma::mat iDscore_Dnuisance_C_M1;
  arma::mat iDscore_Dnuisance_T_M1;
  int space = 0;
  if(returnIID > 1 && method == 3){
    Dscore_Dnuisance_C.resize(p_C, 4);
    Dscore_Dnuisance_C.fill(0.0);

    Dscore_Dnuisance_T.resize(p_T, 4);
    Dscore_Dnuisance_T.fill(0.0);

    iDscore_Dnuisance_C.resize(p_C, 4); // initialized in calcOneScore_TTEperon
	iDscore_Dnuisance_T.resize(p_T, 4); // initialized in calcOneScore_TTEperon

	iDscore_Dnuisance_C_M1.resize(p_C, 4);
	iDscore_Dnuisance_C_M1.fill(0.0);
	iDscore_Dnuisance_T_M1.resize(p_T, 4);
	iDscore_Dnuisance_T_M1.fill(0.0);
  }else{
    Dscore_Dnuisance_C.resize(0, 0);
    Dscore_Dnuisance_T.resize(0, 0);
  }
  
  // pairScore
  double weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative;
  std::vector< double > iPairScore(4); // temporary store results
  rowvec iPairScore_M1(4); // temporary store results
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }

  // ** loop over the pairs
  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){
    // Rcout << iter_pair << "/" << n_pair << " ";
	
    // *** find index of the pair
    iter_C = index_control[iter_pair];
    iter_T = index_treatment[iter_pair];
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
		iPairScore_M1 = RP_score_M1.row(iter_pair);
		if(returnIID > 1){
		  for(int iType=0; iType<4; iType++){
			iDscore_Dnuisance_C_M1.col(iType) = RP_Dscore_Dnuisance_C_M1[iType].row(iter_pair);
			iDscore_Dnuisance_T_M1.col(iType) = RP_Dscore_Dnuisance_T_M1[iType].row(iter_pair);
		  }
		}
      }
    }

    // *** compute adjusted for the previous endpoint
    weight_favorable = (iPairScore[0] - iPairScore_M1(0)) * weight(iter_pair);
    weight_unfavorable = (iPairScore[1] - iPairScore_M1(1)) * weight(iter_pair);
    weight_neutral = iPairScore[2] * weight(iter_pair); 
    weight_uninformative = iPairScore[3] * weight(iter_pair);
    // Rcout << weight(iter_pair) << ") " << weight_favorable << " " << weight_unfavorable << " "<< weight_neutral << " " << weight_uninformative << endl;
    // Rcout << "...) " << iPairScore[0] << " " << iPairScore[1] << " "<< iPairScore[2] << " " << iPairScore[3] << endl;

    // *** store results
    if(iPairScore[0] > zeroPlus){
	  // favorable score
      count_favorable += weight_favorable;

	  // iid(favorable score)
      if(returnIID > 0){
		count_obsC(iter_C,0) += weight_favorable;
		count_obsT(iter_T,0) += weight_favorable;
		
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
			  Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]].row(iter_pair) *= weight_favorable;
		  }

		  if(method == 3){
			if(iPairScore_M1(0) > zeroPlus){
			  Dscore_Dnuisance_C.col(0) += (iDscore_Dnuisance_C.col(0) - iDscore_Dnuisance_C_M1.col(0)) * weight(iter_pair);
			  Dscore_Dnuisance_T.col(0) += (iDscore_Dnuisance_T.col(0) - iDscore_Dnuisance_T_M1.col(0)) * weight(iter_pair);
			}else{
			  Dscore_Dnuisance_C.col(0) += iDscore_Dnuisance_C.col(0) * weight(iter_pair);
			  Dscore_Dnuisance_T.col(0) += iDscore_Dnuisance_T.col(0) * weight(iter_pair);
			}
		  }
		} // end returnIID > 1
		
	  }
	}
	
    if(iPairScore[1] > zeroPlus){
	  // unfavorable score
      count_unfavorable += weight_unfavorable;

	  // iid(unfavorable score)
      if(returnIID > 0){
		count_obsC(iter_C,1) += weight_unfavorable;
		count_obsT(iter_T,1) += weight_unfavorable;
		
		if(returnIID > 1){
		  
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
			Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]].row(iter_pair) *= weight_unfavorable;
		  }

		  if(method == 3){
			if(iPairScore_M1(1) > zeroPlus){
			  Dscore_Dnuisance_C.col(1) += (iDscore_Dnuisance_C.col(1) - iDscore_Dnuisance_C_M1.col(1)) * weight(iter_pair);
			  Dscore_Dnuisance_T.col(1) += (iDscore_Dnuisance_T.col(1) - iDscore_Dnuisance_T_M1.col(1)) * weight(iter_pair);
			}else{
			  Dscore_Dnuisance_C.col(1) += iDscore_Dnuisance_C.col(1) * weight(iter_pair);
			  Dscore_Dnuisance_T.col(1) += iDscore_Dnuisance_T.col(1) * weight(iter_pair);
			}
		  }
		} // end returnIID > 1
		
	  }
    }

    if(weight_neutral > zeroPlus){
	  // neutral score
      count_neutral += weight_neutral;
	  iUpdateRPNeutral = updateIndexNeutral;

	  // iid(neutral score)
      if(returnIID > 0){
		count_obsC(iter_C,2) += weight_neutral;
		count_obsT(iter_T,2) += weight_neutral;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
			Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]].row(iter_pair) *= weight_neutral;
		  }

		  if(method==3){
			Dscore_Dnuisance_C.col(2) += iDscore_Dnuisance_C.col(2) * weight(iter_pair);
			Dscore_Dnuisance_T.col(2) += iDscore_Dnuisance_T.col(2) * weight(iter_pair);
		  }
		} // end returnIID
		
	  }
	}
	
    if(weight_uninformative > zeroPlus){
	  // uninformative score
      count_uninf += weight_uninformative;
	  iUpdateRPUninf = updateIndexUninf;

	  // iid(uninformative score)
	  if(returnIID > 0){ // used in the correction at the pair level
		count_obsC(iter_C,3) += weight_uninformative;
		count_obsT(iter_T,3) += weight_uninformative;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
			Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].row(iter_pair) *= weight_uninformative;
		  }

		  if(method == 3){
			Dscore_Dnuisance_C.col(3) += iDscore_Dnuisance_C.col(3) * weight(iter_pair);
			Dscore_Dnuisance_T.col(3) += iDscore_Dnuisance_T.col(3) * weight(iter_pair);
		  }
		} // end returnIID
		
	  }
	}
  
	// keep information relative to the residual pairs
	if(iUpdateRPNeutral || iUpdateRPUninf){
	  vec_indexPair.push_back(iter_pair);
	  vec_indexC.push_back(iter_C);
	  vec_indexT.push_back(iter_T);
	  vec_favorable.push_back(iPairScore[0]);
	  vec_unfavorable.push_back(iPairScore[1]);
	  vec_neutral.push_back(iPairScore[2]);
	  vec_uninformative.push_back(iPairScore[3]);

	  if(returnIID > 1 && method == 3){
		  
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
		  if(space == 0){
			RP_Dscore_Dnuisance_C[iter_typeRP].resize(n_RP+reserve,p_C);
			RP_Dscore_Dnuisance_T[iter_typeRP].resize(n_RP+reserve,p_T);
			space = reserve;
		  }
		
		  RP_Dscore_Dnuisance_C[iter_typeRP].row(iter_pair) = trans(iDscore_Dnuisance_C.col(iter_typeRP));
		  RP_Dscore_Dnuisance_T[iter_typeRP].row(iter_pair) = trans(iDscore_Dnuisance_T.col(iter_typeRP));
		  space --;
		}
	  }
	  n_RP++;
	}

    if(keepScore){
      matPairScore.row(iter_pair) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
			iPairScore[0] - iPairScore_M1(0), // favorable
			iPairScore[1] - iPairScore_M1(1), // unfavorable
			iPairScore[2], // neutral
			iPairScore[3], // uninformative
			weight(iter_pair), // weight
			weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
			});
    }

    if(iter_pair % 65536 == 0){
      R_CheckUserInterrupt();
    }

  }

  R_CheckUserInterrupt();

  // ** merge information about the residual pairs
  // Rcout << "start merge " << endl;
  RP_score.resize(n_RP,7);
  if(n_RP>0){
	RP_score.col(0) = conv_to<colvec>::from(vec_indexPair);
	RP_score.col(1) = conv_to<colvec>::from(vec_indexC);
	RP_score.col(2) = conv_to<colvec>::from(vec_indexT);
	RP_score.col(3) = conv_to<colvec>::from(vec_favorable);
	RP_score.col(4) = conv_to<colvec>::from(vec_unfavorable);
	RP_score.col(5) = conv_to<colvec>::from(vec_neutral);
	RP_score.col(6) = conv_to<colvec>::from(vec_uninformative);
  }
  // Rcout << "end merge " << endl;

  // ** correction for uninformative pairs
  // (updated by reference)
  // Rcout << "start correction " << endl;
  
  if(count_uninf > 0){
	if(correctionUninf == 1){
      correctionPairs(method, zeroPlus,
					  count_favorable, count_unfavorable, count_neutral, count_uninf,
					  RP_score, matPairScore,
					  count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					  activeUTTE, D_activeUTTE, Dweight_Dnuisance_C, Dweight_Dnuisance_T, 
					  RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					  neutralAsUninf, keepScore, moreEndpoint);
      
    }else if(correctionUninf == 2){
	  correctionIPW(method, zeroPlus,
					count_favorable, count_unfavorable, count_neutral, count_uninf,
					RP_score, matPairScore,
					count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					activeUTTE, D_activeUTTE, Dweight_Dnuisance_C, Dweight_Dnuisance_T, 
					RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					neutralAsUninf, keepScore, moreEndpoint);
    }
  }
  // Rcout << "end correction " << endl;

  // ** export 
  return matPairScore;
  
}

// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
void correctionPairs(int method, double zeroPlus,
					 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, arma::mat& matPairScore,
					 arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					 const vector<int>& activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					 bool neutralAsUninf, bool keepScore, bool moreEndpoint){

  // compute factor
  double factorFavorable;
  double factorUnfavorable;
  double factorNeutral;
  if(count_favorable + count_unfavorable + count_neutral > zeroPlus){
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
    count_obsC.col(0) += factorFavorable * count_obsC.col(3);
    count_obsC.col(1) += factorUnfavorable * count_obsC.col(3);
    count_obsC.col(2) += factorNeutral * count_obsC.col(3);
	(count_obsC.col(3)).fill(0.0);
	
    count_obsT.col(0) += factorFavorable * count_obsT.col(3);
    count_obsT.col(1) += factorUnfavorable * count_obsT.col(3);
    count_obsT.col(2) += factorNeutral * count_obsT.col(3);
	(count_obsT.col(3)).fill(0.0);

	if(returnIID>1){

	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]] += factorFavorable * Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]] += factorUnfavorable * Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]] += factorNeutral * Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].fill(0.0);

		Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]] += factorFavorable * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]] += factorUnfavorable * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]] += factorNeutral * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].fill(0.0);
	  }
	  
	  if(method==3){
		Dscore_Dnuisance_C.col(0) += factorFavorable * Dscore_Dnuisance_C.col(3);
		Dscore_Dnuisance_C.col(1) += factorUnfavorable * Dscore_Dnuisance_C.col(3);
		Dscore_Dnuisance_C.col(2) += factorNeutral * Dscore_Dnuisance_C.col(3);
		(Dscore_Dnuisance_C.col(3)).fill(0.0);

		Dscore_Dnuisance_T.col(0) += factorFavorable * Dscore_Dnuisance_T.col(3);
		Dscore_Dnuisance_T.col(1) += factorUnfavorable * Dscore_Dnuisance_T.col(3);
		Dscore_Dnuisance_T.col(2) += factorNeutral * Dscore_Dnuisance_T.col(3);
		(Dscore_Dnuisance_T.col(3)).fill(0.0);
	  }
	}
  }
  
  // update pair score
  if(moreEndpoint){
	if(factorNeutral > zeroPlus){ 
	  RP_score.col(3) += factorFavorable * RP_score.col(6);
	  RP_score.col(4) += factorUnfavorable * RP_score.col(6);
	  RP_score.col(5) += factorNeutral * RP_score.col(6);
	  RP_score.col(6) = 0;
	  if(returnIID>1 && method==3){
		RP_Dscore_Dnuisance_C[0] += factorFavorable * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[1] += factorUnfavorable * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[2] += factorNeutral * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[3].fill(0.0); 

		RP_Dscore_Dnuisance_T[0] += factorFavorable * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[1] += factorUnfavorable * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[2] += factorNeutral * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[3].fill(0.0); 
	  }
	}else{
	  RP_score.resize(0,0);
	  if(returnIID>1 && method==3){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
		  RP_Dscore_Dnuisance_C[iter_typeRP].resize(0,0);
		  RP_Dscore_Dnuisance_T[iter_typeRP].resize(0,0);
		}
	  }
	}
	
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
void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   const vector<int>& activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool moreEndpoint){

  // compute factor
  double factor;
  if(count_favorable + count_unfavorable + count_neutral > zeroPlus){
	factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
  }else{
	factor = 0;
  }

  // update global score
  count_favorable *= factor;    
  count_unfavorable *= factor;    
  count_neutral *= factor;   
  count_uninf = 0;

  if(returnIID > 0){
    count_obsC.col(0) *= factor;
    count_obsC.col(1) *= factor;
    count_obsC.col(2) *= factor;
	(count_obsC.col(3)).fill(0.0);
	
    count_obsT.col(0) *= factor;
    count_obsT.col(1) *= factor;
    count_obsT.col(2) *= factor;
	(count_obsT.col(3)).fill(0.0);

	if(returnIID>1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].fill(0.0);

		Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].fill(0.0);
	  }

	  if(method==3){
		Dscore_Dnuisance_C.col(0) *= factor;
		Dscore_Dnuisance_C.col(1) *= factor;
		Dscore_Dnuisance_C.col(2) *= factor;
		(Dscore_Dnuisance_C.col(3)).fill(0.0);

		Dscore_Dnuisance_T.col(0) *= factor;
		Dscore_Dnuisance_T.col(1) *= factor;
		Dscore_Dnuisance_T.col(2) *= factor;
		(Dscore_Dnuisance_T.col(3)).fill(0.0);
	  }
	}
  }
  
  // update pair score
  if(moreEndpoint){

	if(neutralAsUninf && factor > zeroPlus){
	  RP_score.col(3) *= factor;
	  RP_score.col(4) *= factor;
	  RP_score.col(5) *= factor;
	  RP_score.col(6) = 0;
	  if(returnIID>1 && method==3){
		RP_Dscore_Dnuisance_C[0] *= factor; 
		RP_Dscore_Dnuisance_C[1] *= factor; 
		RP_Dscore_Dnuisance_C[2] *= factor; 
		RP_Dscore_Dnuisance_C[3].fill(0.0); 

		RP_Dscore_Dnuisance_T[0] *= factor; 
		RP_Dscore_Dnuisance_T[1] *= factor; 
		RP_Dscore_Dnuisance_T[2] *= factor; 
		RP_Dscore_Dnuisance_T[3].fill(0.0); 
	  }
	}else{
	  RP_score.resize(0,0);
	  if(returnIID>1 && method==3){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
		  RP_Dscore_Dnuisance_C[iter_typeRP].resize(0,0);
		  RP_Dscore_Dnuisance_T[iter_typeRP].resize(0,0);
		}
	  }
	}
  }
  
  // update keep scores
  if(keepScore){
    matPairScore.col(7) *= factor;
    matPairScore.col(8) *= factor;
    matPairScore.col(9) *= factor;
	(matPairScore.col(10)).fill(0.0);
  }
}


