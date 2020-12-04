// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:

arma::mat calcAllPairs(arma::colvec endpointC, arma::colvec endpointT, double threshold,
					   arma::colvec statusC, arma::colvec statusT,					   
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   arma::rowvec lastSurv, 					   
					   arma::vec index_control, arma::vec index_treatment, arma::vec weight,					   
					   std::vector<int> activeUTTE, int D_activeUTTE,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,					   
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T,
					   std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
					   double zeroPlus, 
					   int method, int returnIID, int p_C, int p_T, 
					   bool firstEndpoint, bool evalM1, bool updateIndexNeutral, bool updateIndexUninf, bool keepScore, bool precompute, int correctionUninf, bool neutralAsUninf,
					   int debug);

void correctionPairs(int method, double zeroPlus,
					 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 arma::mat& RP_score, arma::mat& matPairScore,
					 arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					 std::vector<int> activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T, int returnIID,
					 bool neutralAsUninf, bool keepScore, bool updateRP);

void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   std::vector<int> activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool updateRP);

void add4vec(std::vector<int>& Vrow,
			 std::vector<int>& Vcol,
			 std::vector<double>& Vvalue0,
			 std::vector<double>& Vvalue1,
			 std::vector<double>& Vvalue2,
			 std::vector<double>& Vvalue3,
			 int n,
			 const arma::mat M);

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
arma::mat calcAllPairs(arma::colvec endpointC, arma::colvec endpointT, double threshold,
					   arma::colvec statusC, arma::colvec statusT,					   
					   arma::mat survTimeC, arma::mat survTimeT, arma::mat survJumpC, arma::mat survJumpT,
					   arma::rowvec lastSurv,
					   arma::vec index_control, arma::vec index_treatment, arma::vec weight,					   
					   std::vector<int> activeUTTE, int D_activeUTTE,
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,					   
					   arma::mat& RP_score,
					   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
					   std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T,
					   std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
					   double zeroPlus, 
					   int method, int returnIID, int p_C, int p_T, 
					   bool firstEndpoint, bool evalM1, bool updateIndexNeutral, bool updateIndexUninf, bool keepScore, bool precompute, int correctionUninf, bool neutralAsUninf,
					   int debug){

  // ** initialize
  int iter_C; // index of the treated patient of the pair in the treatment arm
  int iter_T; // index of the control patient of the pair in the control arm
  int n_Control = endpointC.size(); // number of patients from the control arm
  int n_Treatment = endpointT.size(); // number of patients from the treatment arm
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
  std::vector<int> vec_indexPair(0);
  std::vector<int> vec_indexC(0);
  std::vector<int> vec_indexT(0);
  std::vector<double> vec_favorable(0);
  std::vector<double> vec_unfavorable(0);
  std::vector<double> vec_neutral(0);
  std::vector<double> vec_uninformative(0);

  // residual pairs (neutral and uninformative)
  int n_RP=0;
  if(evalM1==false){
	RP_score.resize(0,0);
  }
  
  // iid (average) aggregated over all pairs
  if(returnIID > 0){ 
    count_obsC.resize(n_Control, 4); 
    count_obsC.fill(0.0);

    count_obsT.resize(n_Treatment, 4);
    count_obsT.fill(0.0);
  }

  // iid (nuisance)
  arma::mat iDscore_Dnuisance_C, iDscore_Dnuisance_T;
  std::vector<int> vecRow_Dnuisance_C, vecRow_Dnuisance_T;
  std::vector<int> vecCol_Dnuisance_C, vecCol_Dnuisance_T;
  std::vector<std::vector<double>> vecValue_Dnuisance_C(4), vecValue_Dnuisance_T(4);
			 
  if(returnIID > 1 && method >= 4){
    // aggregated over all pairs
    Dscore_Dnuisance_C.resize(p_C, 4);
    Dscore_Dnuisance_C.fill(0.0);

    Dscore_Dnuisance_T.resize(p_T, 4);
    Dscore_Dnuisance_T.fill(0.0);

    // pair specific
    iDscore_Dnuisance_C.resize(p_C, 4); // initialized to 0 in calcOnePair_SurvPeron
    iDscore_Dnuisance_T.resize(p_T, 4); // initialized to 0 in calcOnePair_SurvPeron
  }
  
  // ** loop over the pairs
  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){
    if(debug>1){Rcpp::Rcout << iter_pair << "/" << n_pair << " ";}
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
	if(debug>2){Rcpp::Rcout << "("<< iter_C << "/" << n_Control <<";" << iter_T << "/" << n_Treatment << ") ";}
	if(debug>3){Rcpp::Rcout << "("<< endpointT[iter_T] << ";" << endpointC[iter_C] << ";" << threshold <<") ";}
    
	// *** compute score
	if(evalM1 && method > 3){ // Previously analyzed endpoint with Peron's scoring rule. Extract previous score for the residual pairs.
	  iPairScore[0] = RP_score(iter_pair,0);
	  iPairScore[1] = RP_score(iter_pair,1);
	  iPairScore[2] = RP_score(iter_pair,2);
	  iPairScore[3] = RP_score(iter_pair,3);
	  if(returnIID > 1){
		for(int iType=0; iType<4; iType++){ // convertion from sparse to non-sparse
		  arma::colvec tempo_C(RP_Dscore_Dnuisance_C[iType].col(iter_pair));
		  iDscore_Dnuisance_C.col(iType) = tempo_C;
		  arma::colvec tempo_T(RP_Dscore_Dnuisance_T[iType].col(iter_pair));
		  iDscore_Dnuisance_T.col(iType) = tempo_T;
		}
	  }
	}else if(method == 1){ // continuous or binary endpoint
	  iPairScore = calcOnePair_Continuous(endpointT[iter_T] - endpointC[iter_C], threshold);
	}else if(method == 2){ // time to event endpoint with Gehan's scoring rule (right-censored, survival or competing risks)
	  iPairScore = calcOnePair_TTEgehan(endpointT[iter_T] - endpointC[iter_C], statusC[iter_C], statusT[iter_T], threshold);
	}else if(method == 3){ // time to event endpoint with Gehan's scoring rule (left-censored, survival or competing risks)
	  iPairScore = calcOnePair_TTEgehan2(endpointT[iter_T] - endpointC[iter_C], statusC[iter_C], statusT[iter_T], threshold);
	}else if(method == 4){  // time to event endpoint with Peron's scoring rule (right-censored, survival)
	  iPairScore = calcOnePair_SurvPeron(endpointC[iter_C], endpointT[iter_T], statusC[iter_C], statusT[iter_T], threshold,
										 survTimeC.row(iter_C), survTimeT.row(iter_T),
										 survJumpC, survJumpT, lastSurv(0), lastSurv(1),
										 iDscore_Dnuisance_C, iDscore_Dnuisance_T,
										 p_C, p_T, precompute, returnIID);

	}else if(method == 5){  // time to event endpoint with Peron's scoring rule (right-censored, competing risks)
	  iPairScore = calcOnePair_CRPeron(endpointC[iter_C], endpointT[iter_T], statusC[iter_C], statusT[iter_T], threshold,
									   survTimeC.row(iter_C), survTimeT.row(iter_T),
									   survJumpC, survJumpT, lastSurv(0), lastSurv(1), lastSurv(2), lastSurv(3),
									   iDscore_Dnuisance_C, iDscore_Dnuisance_T,
									   p_C, p_T, precompute, returnIID);
	}

	// *** aggregate favorable score and iid over analyzed pairs
	if(iPairScore[0] > zeroPlus){
	  if(debug>4){Rcpp::Rcout << " favorable=" << iPairScore[0] << " ";}

	  // score
	  count_favorable += iPairScore[0] * iWeight;

	  if(returnIID > 0){
	    // iid (average)
	    count_obsC(iter_C,0) += iPairScore[0] * iWeight;
	    count_obsT(iter_T,0) += iPairScore[0] * iWeight;
	    if(returnIID > 1){
	      // iid (nuisance) for the score
	      if(method >= 4){
		    Dscore_Dnuisance_C.col(0) += iDscore_Dnuisance_C.col(0) * iWeight;
		    Dscore_Dnuisance_T.col(0) += iDscore_Dnuisance_T.col(0) * iWeight;
		  }
		  // iid (nuisance) for the weight of the pair
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		    Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[0] * iWeight;
		    Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[0] * iWeight;
		  }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
	    (Dweight_Dnuisance_C[0][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
	    (Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
	  }
	}

	// *** aggregate unfavorable score and iid over analyzed pairs
	if(iPairScore[1] > zeroPlus){
	  if(debug>4){Rcpp::Rcout << " unfavorable=" << iPairScore[1] << " ";}

	  // score
	  count_unfavorable += iPairScore[1] * iWeight;

	  if(returnIID > 0){
	    // iid (average
	    count_obsC(iter_C,1) += iPairScore[1] * iWeight;
	    count_obsT(iter_T,1) += iPairScore[1] * iWeight;
	    if(returnIID > 1){
	      // iid (nuisance) for the score
	      if(method >= 4){
		Dscore_Dnuisance_C.col(1) += iDscore_Dnuisance_C.col(1) * iWeight;
		Dscore_Dnuisance_T.col(1) += iDscore_Dnuisance_T.col(1) * iWeight;
	      }
	      // iid (nuisance) for the weight of the pair
	      for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[1] * iWeight;
		Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[1] * iWeight;
	      }
	    }
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[1][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
		(Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
	  }
	}

	// *** aggregate neutral score and iid over analyzed pairs
	if(iPairScore[2] > zeroPlus){
  	  if(debug>4){Rcpp::Rcout << " neutral=" << iPairScore[2] << " ";}
		  
	  // score
	  count_neutral += iPairScore[2] * iWeight;
	  iUpdateRPNeutral = updateIndexNeutral;

	  if(returnIID > 0){ 
		// iid (average)
		count_obsC(iter_C,2) += iPairScore[2] * iWeight;
		count_obsT(iter_T,2) += iPairScore[2] * iWeight;
		if(returnIID > 1){
		  // iid (nuisance) for the score
		  if(method >= 4){
			Dscore_Dnuisance_C.col(2) += iDscore_Dnuisance_C.col(2) * iWeight;
			Dscore_Dnuisance_T.col(2) += iDscore_Dnuisance_T.col(2) * iWeight;
		  }
		  // iid (nuisance) for the weight of the pair
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
            Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[2] * iWeight;
            Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[2] * iWeight;
          }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[2][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
		(Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
	  }
	}

	// *** aggregate uninformative score and iid over analyzed pairs
	if(iPairScore[3] > zeroPlus){
  	  if(debug>4){Rcpp::Rcout << " uninformative=" << iPairScore[3] << " " ;}

	  // score
	  count_uninf += iPairScore[3] * iWeight;
	  iUpdateRPUninf = updateIndexUninf;

	  // iid(uninformative score)
	  if(returnIID > 0){
		// iid (average)
		count_obsC(iter_C,3) += iPairScore[3] * iWeight;
		count_obsT(iter_T,3) += iPairScore[3] * iWeight;
		if(returnIID > 1){
		  // iid (nuisance) for the score
		  if(method >= 4){
			Dscore_Dnuisance_C.col(3) += iDscore_Dnuisance_C.col(3) * iWeight;
			Dscore_Dnuisance_T.col(3) += iDscore_Dnuisance_T.col(3) * iWeight;
		  }
		  // iid (nuisance) for the weight of the pair
		  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
            Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[3] * iWeight;
            Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].col(iter_pair) *= iPairScore[3] * iWeight;
          }
		}
	  }
	}else if(returnIID > 1){
	  for(int iter_UTTE=0; iter_UTTE<D_activeUTTE; iter_UTTE++){
		(Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
		(Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]].col(iter_pair)) *= 0;
	  }
	}

	// *** keep score + iid relative to the residual pairs
	if(iUpdateRPNeutral || iUpdateRPUninf){
	  if(debug>4){Rcpp::Rcout << " updateRP " ;}

	  vec_indexPair.push_back(iter_pair);
	  vec_indexC.push_back(iter_C);
	  vec_indexT.push_back(iter_T);
	  vec_favorable.push_back(iPairScore[0]);
	  vec_unfavorable.push_back(iPairScore[1]);
	  vec_neutral.push_back(iPairScore[2]);
	  vec_uninformative.push_back(iPairScore[3]);

	  if((returnIID > 1) && (method >= 4) && (evalM1 == false)){
		add4vec(vecRow_Dnuisance_C, vecCol_Dnuisance_C,
				vecValue_Dnuisance_C[0], vecValue_Dnuisance_C[1], vecValue_Dnuisance_C[2], vecValue_Dnuisance_C[3],
				n_RP, iDscore_Dnuisance_C);
		add4vec(vecRow_Dnuisance_T, vecCol_Dnuisance_T,
				vecValue_Dnuisance_T[0], vecValue_Dnuisance_T[1], vecValue_Dnuisance_T[2], vecValue_Dnuisance_T[3],
				n_RP, iDscore_Dnuisance_T);
	  }
	  n_RP++;
	}

	// *** keep score relative to all pairs
	if(keepScore){
	  matPairScore.row(iter_pair) = arma::rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
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

  // ** merge information about the residual pairs
  if(debug>0){Rcpp::Rcout << "collect RP from vectors to matrix " << std::endl;}
  if((n_RP>0) && (evalM1==false)){
	// *** score
	RP_score.resize(n_RP,7);
	RP_score.col(0) = arma::conv_to<arma::colvec>::from(vec_indexPair);
	RP_score.col(1) = arma::conv_to<arma::colvec>::from(vec_indexC);
	RP_score.col(2) = arma::conv_to<arma::colvec>::from(vec_indexT);
	RP_score.col(3) = arma::conv_to<arma::colvec>::from(vec_favorable);
	RP_score.col(4) = arma::conv_to<arma::colvec>::from(vec_unfavorable);
	RP_score.col(5) = arma::conv_to<arma::colvec>::from(vec_neutral);
	RP_score.col(6) = arma::conv_to<arma::colvec>::from(vec_uninformative);

	// *** iid
	// assemble RP_Dscore_Dnuisance from the vectors
	if((returnIID > 1) && (method >= 4)){
  	  // DOCUMENTATION armadillo
	  // (locations, values, n_rows, n_cols)
	  // locations is a dense matrix of type umat, with a size of 2 x N, where N is the number of values to be inserted;
	  // the location of the i-th element is specified by the contents of the i-th column of the locations matrix, where the row is in locations(0,i), and the column is in locations(1,i)

	  // Rcpp::Rcout << "start C: (" << p_C << ";" << n_RP << ";" << Vvalue0_C.size() << ";" << Vvalue1_C.size() << ";" << Vvalue2_C.size() << ";" << Vvalue3_C.size()<<")" << std::endl;
	  int nval_C = vecRow_Dnuisance_C.size();
	  arma::umat iLoc_C(2,nval_C);
	  for(int iV=0; iV<nval_C; iV++){
		iLoc_C(0,iV) = vecRow_Dnuisance_C[iV];
		iLoc_C(1,iV) = vecCol_Dnuisance_C[iV];
	  }
	  arma::colvec iColvecValue_Dnuisance_C;
	  arma::uvec iTestN0_C;
	  for(int iType=0; iType<4; iType++){
		iColvecValue_Dnuisance_C = arma::conv_to<arma::colvec>::from(vecValue_Dnuisance_C[iType]);
		iTestN0_C = find(iColvecValue_Dnuisance_C); // find locations with non-0 score
		RP_Dscore_Dnuisance_C[iType] = arma::sp_mat(iLoc_C.cols(iTestN0_C), iColvecValue_Dnuisance_C.rows(iTestN0_C), p_C, n_RP);
	  }
	  // Rcpp::Rcout << "end C " << std::endl;

	  // Rcpp::Rcout << "start T: (" << p_T << ";" << n_RP << ";" << Vvalue0_T.size() << ";" << Vvalue1_T.size() << ";" << Vvalue2_T.size() << ";" << Vvalue3_T.size()<<")" << std::endl;
  	  int nval_T = vecRow_Dnuisance_T.size();
	  arma::umat iLoc_T(2,nval_T);
	  for(int iV=0; iV<nval_T; iV++){
		iLoc_T(0,iV) = vecRow_Dnuisance_T[iV];
		iLoc_T(1,iV) = vecCol_Dnuisance_T[iV];
	  }
	  arma::colvec iTolvecValue_Dnuisance_T;
	  arma::uvec iTestN0_T;
	  for(int iType=0; iType<4; iType++){
		iTolvecValue_Dnuisance_T = arma::conv_to<arma::colvec>::from(vecValue_Dnuisance_T[iType]);
		iTestN0_T = find(iTolvecValue_Dnuisance_T); // find locations with non-0 score
		RP_Dscore_Dnuisance_T[iType] = arma::sp_mat(iLoc_T.cols(iTestN0_T), iTolvecValue_Dnuisance_T.rows(iTestN0_T), p_T, n_RP);
	  }
	  // Rcpp::Rcout << "end T " << std::endl;
	}
  }
 
  // ** correction for uninformative pairs
  if(count_uninf > 0){
	if(correctionUninf == 1){
	  if(debug>0){Rcpp::Rcout << "correction(pair level)" << std::endl;}
	  correctionPairs(method, zeroPlus,
					  count_favorable, count_unfavorable, count_neutral, count_uninf,
					  RP_score, matPairScore,
					  count_obsC, count_obsT, Dscore_Dnuisance_C, Dscore_Dnuisance_T,
					  activeUTTE, D_activeUTTE, Dweight_Dnuisance_C, Dweight_Dnuisance_T, 
					  RP_Dscore_Dnuisance_C, RP_Dscore_Dnuisance_T, returnIID,
					  neutralAsUninf, keepScore, updateIndexNeutral && (evalM1==false));
      
	}else if(correctionUninf == 2){
	  if(debug>0){Rcpp::Rcout << "correction(trial level)" << std::endl;}
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
					 std::vector<int> activeUTTE, int D_activeUTTE,
					 std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
					 std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T, int returnIID,
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
		Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]] *= 0;

		Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]] += factorFavorable * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]] += factorUnfavorable * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]] += factorNeutral * Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]];
		Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]] *= 0;
	  }
	  
	  if(method>=4){
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
		
	  if(returnIID>1 && method>=4){
		RP_Dscore_Dnuisance_C[0] += factorFavorable * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[1] += factorUnfavorable * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[2] += factorNeutral * RP_Dscore_Dnuisance_C[3]; 
		RP_Dscore_Dnuisance_C[3] *= 0; 

		RP_Dscore_Dnuisance_T[0] += factorFavorable * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[1] += factorUnfavorable * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[2] += factorNeutral * RP_Dscore_Dnuisance_T[3]; 
		RP_Dscore_Dnuisance_T[3] *= 0; 
	  }
	}else{
	  RP_score.resize(0,0);
	  if(returnIID>1 && method>=4){
		for(int iter_typeRP=0; iter_typeRP<4; iter_typeRP++){
		  RP_Dscore_Dnuisance_C[iter_typeRP].resize(0,0);
		  RP_Dscore_Dnuisance_T[iter_typeRP].resize(0,0);
		}
	  }
	}
  }
	
  // update keep scores
  if(keepScore){
    matPairScore.col(7) += factorFavorable * matPairScore.col(10);
    matPairScore.col(8) += factorUnfavorable * matPairScore.col(10);
    matPairScore.col(9) += factorNeutral * matPairScore.col(10);
    (matPairScore.col(10)).fill(0.0);
  }
}

// * correctionIPW
// author Brice Ozenne
void correctionIPW(int method, double zeroPlus,
				   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   arma::mat& RP_score, arma::mat& matPairScore,
				   arma::mat& count_obsC, arma::mat& count_obsT, arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
				   std::vector<int> activeUTTE, int D_activeUTTE,
				   std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_C, std::vector<std::vector< arma::sp_mat >>& Dweight_Dnuisance_T,
				   std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_C, std::vector< arma::sp_mat >& RP_Dscore_Dnuisance_T, int returnIID,
				   bool neutralAsUninf, bool keepScore, bool updateRP){

  // compute factor
  double factor;
  if(count_favorable + count_unfavorable + count_neutral > zeroPlus){
	factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
  }else{
	factor = 0;
  }
  // Rcpp::Rcout << factor << std::endl;

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
		Dweight_Dnuisance_C[3][activeUTTE[iter_UTTE]] *= 0;

		Dweight_Dnuisance_T[0][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[1][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[2][activeUTTE[iter_UTTE]] *= factor;
		Dweight_Dnuisance_T[3][activeUTTE[iter_UTTE]] *= 0;
	  }

	  if(method>=4){
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
	  
	  if(returnIID>1 && method>=4){
		RP_Dscore_Dnuisance_C[0] *= factor; 
		RP_Dscore_Dnuisance_C[1] *= factor; 
		RP_Dscore_Dnuisance_C[2] *= factor; 
		RP_Dscore_Dnuisance_C[3] *= 0; 

		RP_Dscore_Dnuisance_T[0] *= factor; 
		RP_Dscore_Dnuisance_T[1] *= factor; 
		RP_Dscore_Dnuisance_T[2] *= factor; 
		RP_Dscore_Dnuisance_T[3] *= 0; 
	  }
	}else{
	  RP_score.resize(0,0);
	  if(returnIID>1 && method>=4){
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

// * add4vec
void add4vec(std::vector<int>& Vrow,
			 std::vector<int>& Vcol,
			 std::vector<double>& Vvalue0,
			 std::vector<double>& Vvalue1,
			 std::vector<double>& Vvalue2,
			 std::vector<double>& Vvalue3,
			 int n,
			 const arma::mat M){

  // ** find parameters for which the derivative is not 0
  // DOCUMENTATION Armadillo
  // For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) 
  arma::umat N = abs(M)>0;
  arma::uvec indexN0 = find(sum(N,1)); // sum by
  int nN0 = indexN0.size();

  // ** extend vectors
  int l = Vrow.size();
  Vrow.resize(l+nN0);
  Vcol.resize(l+nN0);
  Vvalue0.resize(l+nN0);
  Vvalue1.resize(l+nN0);
  Vvalue2.resize(l+nN0);
  Vvalue3.resize(l+nN0);

  // ** update
  int iRow;
  for(int i=0; i<nN0; i++){
	iRow = indexN0[i];
	Vrow[l+i] = iRow;
	Vcol[l+i] = n;
	Vvalue0[l+i] = M(iRow, 0);
	Vvalue1[l+i] = M(iRow, 1);
	Vvalue2[l+i] = M(iRow, 2);
	Vvalue3[l+i] = M(iRow, 3);
  }

}


// * sp_mat4vec [no used]
// adapted from https://stackoverflow.com/questions/40222092/access-and-modify-the-non-zero-elements-of-sparse-matrix-of-class-armasp-mat-u
// convert sparse to vector
void sp_mat4vec(const arma::sp_mat& X, int RP,
				int& n,
				std::vector<int>& Vrow,
				std::vector<int>& Vcol,
				std::vector<double>& Vvalue0,
				std::vector<double>& Vvalue1,
				std::vector<double>& Vvalue2,
				std::vector<double>& Vvalue3){
    
  // Make const iterator
  arma::sp_mat::const_iterator it = X.begin();
  arma::sp_mat::const_iterator it_end = X.end();

  // Rcout << X << endl;
  // Calculate number of points
  n = std::distance(it, it_end);

  // Build a location storage matrix
  Vrow.resize(n);
  Vcol.resize(n);
  std::fill(Vcol.begin(), Vcol.end(), RP);
  Vvalue0.resize(n);
  std::fill(Vvalue0.begin(), Vvalue0.end(), 0.0);
  Vvalue1.resize(n);
  std::fill(Vvalue1.begin(), Vvalue1.end(), 0.0);
  Vvalue2.resize(n);
  std::fill(Vvalue2.begin(), Vvalue2.end(), 0.0);
  Vvalue3.resize(n);
  std::fill(Vvalue3.begin(), Vvalue3.end(), 0.0);
  
  // Collecting locations
  int iCol;
  for(int i = 0; i < n; ++i){
  	Vrow[i] = it.row();
  	iCol = it.col();
	if(iCol == 0){
	  Vvalue0[i] = (*it);
	}else if(iCol == 1){
	  Vvalue1[i] = (*it);
	}else if(iCol == 2){
	  Vvalue2[i] = (*it);
	}else if(iCol == 3){
	  Vvalue3[i] = (*it);
	}
  	++it; // increment
	// Rcout << Vrow[i] << endl;  
	// Rcout << Vcol[i] << endl;  
  }
  return;
}

