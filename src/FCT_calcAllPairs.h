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
					   arma::colvec statusC, arma::colvec statusT,					   
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   arma::rowvec lastSurv, 					   
					   arma::vec index_control, arma::vec index_treatment, arma::vec weight,					   
					   vector<int> activeUTTE, int D_activeUTTE,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,					   
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T,
					   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					   double zeroPlus, 
					   int method, int returnIID, int p_C, int p_T, 
					   bool firstEndpoint, bool evalM1, bool updateIndexNeutral, bool updateIndexUninf, bool keepScore, int correctionUninf, bool neutralAsUninf,
					   int debug);

void prepareWeight(arma::vec& iPairWeight, std::vector<std::vector< arma::mat >>& iPairDweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& iPairDweight_Dnuisance_T,
				   std::vector<int>& activeUTTE, int& D_activeUTTE,
				   int iter_d, int iIndex_UTTE, const std::vector<arma::mat>& RP_score,
				   const std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_C, const std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_T,
				   int iNUTTE_analyzedPeron, int correctionUninf, double zeroPlus, bool neutralAsUninf, int returnIID);

void correctionPairs(int method, double zeroPlus,
					 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, arma::mat& matPairScore,
					 arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					 vector<int> activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					 bool neutralAsUninf, bool keepScore, bool updateRP);

void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   vector<int> activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool updateRP);

// * calcAllPairs
// perform pairwise comparisons over all possible pairs for a continuous endpoints
// WARNING: count_favorable, count_unfavorable, count_neutral, count_uninf are not initialized to 0 (to be able to sum over strata)
// WARNING: when evalM1 is true then moreEndpoint should be false
//
// OUTPUT: count_favorable, count_unfavorable, count_neutral, count_uninf,
//         RP_score
//         count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T
//         RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T
//         matPairScore
// author Brice Ozenne
arma::mat calcAllPairs(arma::colvec Control, arma::colvec Treatment, double threshold,
					   arma::colvec statusC, arma::colvec statusT,					   
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   arma::rowvec lastSurv,
					   arma::vec index_control, arma::vec index_treatment, arma::vec weight,					   
					   vector<int> activeUTTE, int D_activeUTTE,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,					   
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T,
					   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					   double zeroPlus, 
					   int method, int returnIID, int p_C, int p_T, 
					   bool firstEndpoint, bool evalM1, bool updateIndexNeutral, bool updateIndexUninf, bool keepScore, int correctionUninf, bool neutralAsUninf,
					   int debug){

  // ** initialize
  int iter_C; // index of the treated patient of the pair in the treatment arm
  int iter_T; // index of the control patient of the pair in the control arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_pair; // number of pairs
  if(firstEndpoint){	
	n_pair = n_Treatment * n_Control;
	iter_C = -1;
	iter_T = 0;
  }else{
	n_pair = weight.size();
  }
  
  double iWeight=1;
  bool iUpdateRPNeutral;
  bool iUpdateRPUninf;

  // counts
  count_favorable = 0;
  count_unfavorable = 0;
  count_neutral = 0;
  count_uninf = 0;

  // pairScore  
  std::vector< double > iPairScore(4); // temporary store results
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }
  vector<int> vec_indexPair(0);
  vector<int> vec_indexC(0);
  vector<int> vec_indexT(0);
  vector<double> vec_favorable(0);
  vector<double> vec_unfavorable(0);
  vector<double> vec_neutral(0);
  vector<double> vec_uninformative(0);

  // residual pairs (neutral and uninformative)
  int n_RP=0;
  if(evalM1==false){
	RP_score.resize(0,0);
  }
  
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
  int reserve = 0; 
  if(returnIID > 1 && method == 4){
    Dscore_Dnuisance_C.resize(p_C, 4);
    Dscore_Dnuisance_C.fill(0.0);

    Dscore_Dnuisance_T.resize(p_T, 4);
    Dscore_Dnuisance_T.fill(0.0);

    iDscore_Dnuisance_C.resize(p_C, 4); // initialized in calcOneScore_SurvPeron
	iDscore_Dnuisance_T.resize(p_T, 4); // initialized in calcOneScore_SurvPeron

	if(evalM1){	
	  for(int iType=0; iType<4; iType++){
		RP_Dscore_Dnuisance_C[iType] = arma::trans(RP_Dscore_Dnuisance_C[iType]);
		RP_Dscore_Dnuisance_T[iType] = arma::trans(RP_Dscore_Dnuisance_T[iType]);
	  }
	}else{
	  // find the number of RP by running simplified GPC
	  int iter_CC = -1;
	  int iter_TT = 0;
	  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){
		if(firstEndpoint){
		  if(iter_CC < (n_Control-1)){
			iter_CC ++;
		  }else if(iter_T < (n_Treatment-1)){
			iter_CC = 0;
			iter_TT ++;
		  }else{
			break;
		  }
		}else{
		  iter_CC = index_control[iter_pair];
		  iter_TT = index_treatment[iter_pair];
		  iWeight = weight(iter_pair);
		}
		iPairScore = calcOnePair_TTEgehan(Treatment[iter_TT] - Control[iter_CC], statusC[iter_CC], statusT[iter_TT], threshold);
		if(((iPairScore[2] > zeroPlus) && updateIndexNeutral) || ((iPairScore[3] > zeroPlus) && updateIndexUninf)){
		  reserve++;
		}
	  }
	  for(int iType=0; iType<4; iType++){
		RP_Dscore_Dnuisance_C[iType].resize(reserve,p_C);
		RP_Dscore_Dnuisance_T[iType].resize(reserve,p_T);
	  }
	}
  }else{
    Dscore_Dnuisance_C.resize(0, 0);
    Dscore_Dnuisance_T.resize(0, 0);
  }
  
  // ** loop over the pairs
  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){
	if(debug>1){Rcout << iter_pair << "/" << n_pair << " ";}
	iUpdateRPNeutral = false;
	iUpdateRPUninf = false;
	
	// *** find index of the pair
	if(firstEndpoint){
	  if(iter_C < (n_Control-1)){
		iter_C ++;
	  }else if(iter_T < (n_Treatment-1)){
		iter_C = 0;
		iter_T ++;
	  }else{
		break;
	  }
	}else{
	  iter_C = index_control[iter_pair];
	  iter_T = index_treatment[iter_pair];
	  iWeight = weight(iter_pair);
	}
	if(debug>2){Rcout << "("<< iter_C << "/" << n_Control <<";" << iter_T << "/" << n_Treatment << ") ";}
	if(debug>3){Rcout << "("<< Treatment[iter_T] << ";" << Control[iter_C] << ";" << threshold <<") ";}
    
	// *** compute score
	if(method == 1){
	  iPairScore = calcOnePair_Continuous(Treatment[iter_T] - Control[iter_C], threshold);
	}else if(method == 2){
	  iPairScore = calcOnePair_TTEgehan(Treatment[iter_T] - Control[iter_C], statusC[iter_C], statusT[iter_T], threshold);
	}else if(method == 3){
	  iPairScore = calcOnePair_TTEgehan2(Treatment[iter_T] - Control[iter_C], statusC[iter_C], statusT[iter_T], threshold);
	}else if(method > 3){
	  if(evalM1){
		iPairScore[0] = RP_score(iter_pair,0);
		iPairScore[1] = RP_score(iter_pair,1);
		iPairScore[2] = RP_score(iter_pair,2);
		iPairScore[3] = RP_score(iter_pair,3);
		if(returnIID > 1){
		  for(int iType=0; iType<4; iType++){
			iDscore_Dnuisance_C.col(iType) = RP_Dscore_Dnuisance_C[iType].col(iter_pair);
			iDscore_Dnuisance_T.col(iType) = RP_Dscore_Dnuisance_T[iType].col(iter_pair);
		  }
		}
	  }else if(method == 4){
		iPairScore = calcOneScore_SurvPeron(Control[iter_C], Treatment[iter_T], 
											statusC[iter_C], statusT[iter_T], threshold,
											survTimeC.row(iter_C), survTimeT.row(iter_T),
											survJumpC, survJumpT, lastSurv(0), lastSurv(1),
											iDscore_Dnuisance_C, iDscore_Dnuisance_T, returnIID);
	  }else if(method == 5){
		iPairScore = calcOnePair_CRPeron(Control[iter_C], Treatment[iter_T],
										 statusC[iter_C], statusT[iter_T], threshold,
										 survTimeC.row(iter_C), survTimeT.row(iter_T), survJumpC,
										 lastSurv(0), lastSurv(1), lastSurv(2), lastSurv(3));
	  }
	}	

	// *** store results
	if(iPairScore[0] > zeroPlus){
	  if(debug>4){Rcout << " favorable=" << iPairScore[0] << " ";}

	  // favorable score
	  count_favorable += iPairScore[0] * iWeight;

	  // iid(favorable score)
	  if(returnIID > 0){
		count_obsC(iter_C,0) += iPairScore[0] * iWeight;
		count_obsT(iter_T,0) += iPairScore[0] * iWeight;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
			Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[0] * iWeight;
			Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[0] * iWeight;
		  }
		  if(method == 4){
			Dscore_Dnuisance_C.col(0) += iDscore_Dnuisance_C.col(0) * iWeight;
			Dscore_Dnuisance_T.col(0) += iDscore_Dnuisance_T.col(0) * iWeight;
		  }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
		(Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
	  }
	}

	if(iPairScore[1] > zeroPlus){
	  if(debug>4){Rcout << " unfavorable=" << iPairScore[1] << " ";}

	  // unfavorable score
	  count_unfavorable += iPairScore[1] * iWeight;

	  // iid(unfavorable score)
	  if(returnIID > 0){
		count_obsC(iter_C,1) += iPairScore[1] * iWeight;
		count_obsT(iter_T,1) += iPairScore[1] * iWeight;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
            Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[1] * iWeight;
            Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[1] * iWeight;
          }
		  if(method == 4){
			Dscore_Dnuisance_C.col(1) += iDscore_Dnuisance_C.col(1) * iWeight;
			Dscore_Dnuisance_T.col(1) += iDscore_Dnuisance_T.col(1) * iWeight;
		  }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
		(Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
	  }
	}

	if(iPairScore[2] > zeroPlus){
  	  if(debug>4){Rcout << " neutral=" << iPairScore[2] << " ";}
		  
	  // neutral score
	  count_neutral += iPairScore[2] * iWeight;
	  iUpdateRPNeutral = updateIndexNeutral;

	  // iid(neutral score)
	  if(returnIID > 0){ // used in the correction at the pair level
		count_obsC(iter_C,2) += iPairScore[2] * iWeight;
		count_obsT(iter_T,2) += iPairScore[2] * iWeight;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
            Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[2] * iWeight;
            Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[2] * iWeight;
          }
		  if(method == 4){
			Dscore_Dnuisance_C.col(2) += iDscore_Dnuisance_C.col(2) * iWeight;
			Dscore_Dnuisance_T.col(2) += iDscore_Dnuisance_T.col(2) * iWeight;
		  }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
		(Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
	  }
	}

	if(iPairScore[3] > zeroPlus){
  	  if(debug>4){Rcout << " uninformative=" << iPairScore[3] << " " ;}

	  // uninformative score
	  count_uninf += iPairScore[3] * iWeight;
	  iUpdateRPUninf = updateIndexUninf;

	  // iid(uninformative score)
	  if(returnIID > 0){ // used in the correction at the pair level
		count_obsC(iter_C,3) += iPairScore[3] * iWeight;
		count_obsT(iter_T,3) += iPairScore[3] * iWeight;
		if(returnIID > 1){
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
            Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[3] * iWeight;
            Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].row(iter_pair) *= iPairScore[3] * iWeight;
          }
		  if(method == 4){
			Dscore_Dnuisance_C.col(3) += iDscore_Dnuisance_C.col(3) * iWeight;
			Dscore_Dnuisance_T.col(3) += iDscore_Dnuisance_T.col(3) * iWeight;
		  }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
		(Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].row(iter_pair)).fill(0.0);
	  }
	}

	// keep information relative to the residual pairs
	if(iUpdateRPNeutral || iUpdateRPUninf){
	  if(debug>4){Rcout << " updateRP " ;}

	  vec_indexPair.push_back(iter_pair);
	  vec_indexC.push_back(iter_C);
	  vec_indexT.push_back(iter_T);
	  vec_favorable.push_back(iPairScore[0]);
	  vec_unfavorable.push_back(iPairScore[1]);
	  vec_neutral.push_back(iPairScore[2]);
	  vec_uninformative.push_back(iPairScore[3]);

	  if((returnIID > 1) && (method == 4) && (evalM1 == false)){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){		
		  RP_Dscore_Dnuisance_C[iter_typeRP].row(n_RP) = arma::trans(iDscore_Dnuisance_C.col(iter_typeRP));
		  RP_Dscore_Dnuisance_T[iter_typeRP].row(n_RP) = arma::trans(iDscore_Dnuisance_T.col(iter_typeRP));
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
			iWeight, // weight
			iPairScore[0] * iWeight, iPairScore[1] * iWeight, iPairScore[2] * iWeight, iPairScore[3] * iWeight // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
			});		
	}

	if(iter_pair % 65536 == 0){
	  R_CheckUserInterrupt();
	}

  } // end  iter_pair
  
  R_CheckUserInterrupt();

  // ** resize RP_Dscore_Dnuisance for the case where Peron could decidedly classify the pair
  if((returnIID > 1) && (method == 4) && (evalM1==false) && (n_RP != reserve)){	
	for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
	  RP_Dscore_Dnuisance_C[iter_typeRP].resize(n_RP,p_C);
	  RP_Dscore_Dnuisance_T[iter_typeRP].resize(n_RP,p_T);
	}
  }
 
  // ** merge information about the residual pairs
  if(debug>0){Rcout << "collect RP from vectors to matrix " << endl;}
  if((n_RP>0) && (evalM1==false)){
	RP_score.resize(n_RP,7);
	RP_score.col(0) = conv_to<colvec>::from(vec_indexPair);
	RP_score.col(1) = conv_to<colvec>::from(vec_indexC);
	RP_score.col(2) = conv_to<colvec>::from(vec_indexT);
	RP_score.col(3) = conv_to<colvec>::from(vec_favorable);
	RP_score.col(4) = conv_to<colvec>::from(vec_unfavorable);
	RP_score.col(5) = conv_to<colvec>::from(vec_neutral);
	RP_score.col(6) = conv_to<colvec>::from(vec_uninformative);
  }  
  
  // ** correction for uninformative pairs
  if(count_uninf > 0){
	if(correctionUninf == 1){
	  if(debug>0){Rcout << "correction(pair level)" << endl;}
	  correctionPairs(method, zeroPlus,
					  count_favorable, count_unfavorable, count_neutral, count_uninf,
					  RP_score, matPairScore,
					  count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					  activeUTTE, D_activeUTTE, Dweight_Dnuisance_C, Dweight_Dnuisance_T, 
					  RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					  neutralAsUninf, keepScore, updateIndexNeutral && (evalM1==false));
      
	}else if(correctionUninf == 2){
	  if(debug>0){Rcout << "correction(trial level)" << endl;}
	  correctionIPW(method, zeroPlus,
					count_favorable, count_unfavorable, count_neutral, count_uninf,
					RP_score, matPairScore,
					count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					activeUTTE, D_activeUTTE, Dweight_Dnuisance_C, Dweight_Dnuisance_T, 
					RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					neutralAsUninf, keepScore, updateIndexNeutral && (evalM1==false));
	}
  }

  // ** export
  return matPairScore;
  
}

// * prepareWeight
// author Brice Ozenne
void prepareWeight(arma::vec& iPairWeight, std::vector<std::vector< arma::mat >>& iPairDweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& iPairDweight_Dnuisance_T,
				   std::vector<int>& activeUTTE, int& D_activeUTTE,
				   int iter_d, int iIndex_UTTE, const std::vector<arma::mat>& RP_score,
				   const std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_C, const std::vector< std::vector< arma::mat > >& RP_Dscore_Dnuisance_T,
				   int iNUTTE_analyzedPeron, int correctionUninf, double zeroPlus, bool neutralAsUninf, int returnIID){


  // ** cumulate weights over previous TTE endpoints
  activeUTTE.resize(0);
  for(int iter_UTTE=0 ; iter_UTTE<iNUTTE_analyzedPeron; iter_UTTE++){
	if(iter_UTTE != iIndex_UTTE){
	  activeUTTE.push_back(iter_UTTE);
	  if(neutralAsUninf){
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
		if(neutralAsUninf){
		  iPairDweight_Dnuisance_C[0][iter_UTTE] = RP_Dscore_Dnuisance_C[iter_UTTE][2] + RP_Dscore_Dnuisance_C[iter_UTTE][3];
		  iPairDweight_Dnuisance_C[0][iter_UTTE].each_col() %= iPairWeight/(RP_score[iter_UTTE].col(2) + RP_score[iter_UTTE].col(3));

		  iPairDweight_Dnuisance_T[0][iter_UTTE] = RP_Dscore_Dnuisance_T[iter_UTTE][2] + RP_Dscore_Dnuisance_T[iter_UTTE][3];
		  iPairDweight_Dnuisance_T[0][iter_UTTE].each_col() %= iPairWeight/(RP_score[iter_UTTE].col(2) + RP_score[iter_UTTE].col(3));
		}else{
		  iPairDweight_Dnuisance_C[0][iter_UTTE] = RP_Dscore_Dnuisance_C[iter_UTTE][3];
		  iPairDweight_Dnuisance_C[0][iter_UTTE].each_col() %= iPairWeight/RP_score[iter_UTTE].col(3);

		  iPairDweight_Dnuisance_T[0][iter_UTTE] = RP_Dscore_Dnuisance_T[iter_UTTE][3];
		  iPairDweight_Dnuisance_T[0][iter_UTTE].each_col() %= iPairWeight/RP_score[iter_UTTE].col(3);
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

// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
//
// OUTPUT: count_favorable, count_unfavorable, count_neutral, count_uninf,
//         RP_score, matPairScore,
//         count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T, Dweight_Dnuisance_C, Dweight_Dnuisance_T
//         RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T
// author Brice Ozenne
void correctionPairs(int method, double zeroPlus,
					 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, arma::mat& matPairScore,
					 arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					 vector<int> activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
					 bool neutralAsUninf, bool keepScore, bool updateRP){

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
	  
	  if(method==4){
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
  if(updateRP){
	if(factorNeutral > zeroPlus){
	  // used to substract contribution of the threhold at smaller thresholds
	  RP_score.col(3) += factorFavorable * RP_score.col(6);
	  RP_score.col(4) += factorUnfavorable * RP_score.col(6);

	  // used for the weights and to substract contribution of the threhold at smaller thresholds
	  RP_score.col(5) += factorNeutral * RP_score.col(6); 
	  (RP_score.col(6)).fill(0.0);  
		
	  if(returnIID>1 && method==4){
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
	  if(returnIID>1 && method==4){
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
// author Brice Ozenne
void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   vector<int> activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::mat >& RP_Dscore_Dnuisance_C, std::vector< arma::mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool updateRP){

  // compute factor
  double factor;
  if(count_favorable + count_unfavorable + count_neutral > zeroPlus){
	factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
  }else{
	factor = 0;
  }
  // Rcout << factor << endl;

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

	  if(method==4){
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
  if(updateRP){

	if(neutralAsUninf && factor > zeroPlus){
	  // used to substract contribution of the current threshold at smaller thresholds
	  RP_score.col(3) *= factor; 
	  RP_score.col(4) *= factor; 

	  // used for the weights and to substract contribution of the threhold at smaller thresholds
	  RP_score.col(5) *= factor; 
	  (RP_score.col(6)).fill(0.0);
	  
	  if(returnIID>1 && method==4){
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
	  if(returnIID>1 && method==4){
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

  return;
}


