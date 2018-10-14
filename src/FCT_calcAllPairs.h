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

arma::mat calcAllPairs(const arma::colvec& Control, const arma::colvec& Treatment, double threshold,
		       const arma::colvec& deltaC, const arma::colvec& deltaT, 
		       const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
		       double lastSurvC, double lastSurvT, 
		       int method, int correctionUninf,
		       double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		       vector<int>& index_neutralC, vector<int>& index_neutralT,
		       vector<int>& index_uninfC, vector<int>& index_uninfT, 
		       arma::vec& weight,
		       bool neutralAsUninf, bool keepScore, bool moreEndpoint);

arma::mat calcSubsetPairs(const arma::colvec& Control, const arma::colvec& Treatment, double threshold, 
			  const arma::colvec& deltaC, const arma::colvec& deltaT, 
			  const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
			  double lastSurvC, double lastSurvT, 
			  const vector<int>& index_neutralC_M1, const vector<int>& index_neutralT_M1, 
			  const vector<int>& index_uninfC_M1, const vector<int>& index_uninfT_M1,
			  const arma::vec& cumWeight_M1, const double threshold_M1,
			  const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
			  int method, int correctionUninf,
			  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
			  vector<int>& index_neutralC, vector<int>& index_neutralT, 
			  vector<int>& index_uninfC, vector<int>& index_uninfT,
			  arma::vec& weight, arma::uvec& index_weight,
			  bool neutralAsUninf, bool keepScore, bool moreEndpoint);

void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		     vector<int>& index_uninfC, vector<int>& index_uninfT, 
		     vector<int>& index_neutralC, vector<int>& index_neutralT,
		     const vector<double>& wNeutral, const vector<int>& index_wNeutral, const vector<double>& wUninf, const vector<int>& index_wUninf,
		     arma::vec& weight, arma::uvec& index_weight,
		     bool updateWeight, bool updateIndex, bool keepScore, arma::mat& matPairScore);

void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		   vector<int>& index_uninfC, vector<int>& index_uninfT, 
		   const vector<double>& wNeutral, 
		   arma::vec& weight,  
		   bool updateWeight, bool updateIndex, bool keepScore, arma::mat& matPairScore);


// * calcAllPairs
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs(const arma::colvec& Control, const arma::colvec& Treatment, double threshold,
		       const arma::colvec& deltaC, const arma::colvec& deltaT, 
		       const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
		       double lastSurvC, double lastSurvT, 
		       int method, int correctionUninf,
		       double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		       vector<int>& index_neutralC, vector<int>& index_neutralT,
		       vector<int>& index_uninfC, vector<int>& index_uninfT, 
		       arma::vec& weight,
		       bool neutralAsUninf, bool keepScore, bool moreEndpoint){

  // ** initialize
  int n_Treatment = Treatment.size(); // number of patients from the treatment arm
  int n_Control = Control.size(); // number of patients from the control arm
  int n_pair = n_Treatment * n_Control;

  vector<double> wNeutral(0); // weight of the neutral pairs
  vector<double> wUninf(0); // weight of the uninformative pairs
  
  if(wNeutral.max_size() < (n_pair/10.0)){
    index_neutralC.reserve(n_pair);
    index_neutralT.reserve(n_pair);
    index_uninfC.reserve(n_pair);
    index_uninfT.reserve(n_pair);
    wNeutral.reserve(n_pair);
    wUninf.reserve(n_pair);
  }

  bool updateIndexNeutral = moreEndpoint && neutralAsUninf;
  bool updateIndexUninf = moreEndpoint;
    
  vector<double> iScore(4); // temporary store results
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }else{
    matPairScore.resize(0, 0); 
  }

  double zeroPlus = pow(10.0,-12.0);
  
  // ** loop over the pairs
  int iter_pair = 0;
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients

      // score
      if(method == 1){
	iScore = calcOnePair_Continuous(Control[iter_C], Treatment[iter_T], threshold, 1.0);
      }else if(method == 2){
	iScore = calcOnePair_TTEgehan(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, 1.0);
      }else if(method == 3){
	iScore = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
				       deltaC[iter_C], deltaT[iter_T], threshold,
				       survTimeC.row(iter_C), survTimeT.row(iter_T),
				       survJumpC, survJumpT, lastSurvC, lastSurvT);      
      }	

      // store results
      count_favorable += iScore[0];
      count_unfavorable += iScore[1];
      count_neutral += iScore[2];
      count_uninf += iScore[3];
	
      if( updateIndexNeutral && (iScore[2] > zeroPlus) ){
	index_neutralC.push_back(iter_C);     
	index_neutralT.push_back(iter_T);
        wNeutral.push_back(iScore[2]);
      }
      if( updateIndexUninf && (iScore[3] > zeroPlus) ){
	index_uninfC.push_back(iter_C);     
	index_uninfT.push_back(iter_T);
        wUninf.push_back(iScore[3]); 
      }
      if(keepScore){
        matPairScore.row(iter_pair) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
					      iScore[0], // favorable
					      iScore[1], // unfavorable
					      iScore[2], // neutral
					      iScore[3], // uninformative
					      1.0, // weight
					      0.0, 0.0, 0.0, 0.0 // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
	  });
	iter_pair++;
      }
      
    }
  }

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(correctionUninf > 0 && count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
      bool updateWeight = neutralAsUninf && moreEndpoint;
      bool updateIndex = false;

      if(correctionUninf == 1){
      vector<int> index_wNeutral;
      vector<int> index_wUninf;
      arma::uvec index_weight;

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfC, index_uninfT, 
		      index_neutralC, index_neutralT, 
		      wNeutral, index_wNeutral, wUninf, index_wUninf,
		      weight, index_weight,
		      updateWeight, updateIndex, keepScore, matPairScore);
      // update by reference
      
    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		    index_uninfC, index_uninfT, 
		    wNeutral,
		    weight,
		    updateWeight, updateIndex, keepScore, matPairScore);
      // update by reference
      
    }
  }else{
    wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
    weight = arma::conv_to<vec>::from(wNeutral);
    // NOTE: no index_weight to update since it is the first outcome
  }
	
  // ** export
  return matPairScore;
  
}


// * calcSubsetPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a given endpoint 
arma::mat calcSubsetPairs(const arma::colvec& Control, const arma::colvec& Treatment, double threshold, 
			  const arma::colvec& deltaC, const arma::colvec& deltaT, 
			  const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
			  double lastSurvC, double lastSurvT, 
			  const vector<int>& index_neutralC_M1, const vector<int>& index_neutralT_M1,
			  const vector<int>& index_uninfC_M1, const vector<int>& index_uninfT_M1,
			  const arma::vec& cumWeight_M1, const double threshold_M1,
			  const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
			  int method, int correctionUninf,
			  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
			  vector<int>& index_neutralC, vector<int>& index_neutralT, 
			  vector<int>& index_uninfC, vector<int>& index_uninfT,
			  arma::vec& weight, arma::uvec& index_weight,
			  bool neutralAsUninf, bool keepScore, bool moreEndpoint){

  // ** initialize
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  int nNeutral_pairs_M1 = index_neutralC_M1.size();
  int nUninf_pairs_M1 = index_uninfC_M1.size();
  int n_pair = nNeutral_pairs_M1 + nUninf_pairs_M1;
  
  index_neutralT.reserve(n_pair);
  index_neutralC.reserve(n_pair);
  index_uninfC.reserve(n_pair);
  index_uninfT.reserve(n_pair);

  vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  index_wNeutral.reserve(n_pair);
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  index_wUninf.reserve(n_pair);
  vector<double> wNeutral(0); // weight of the neutral pairs
  wNeutral.reserve(n_pair);
  vector<double> wUninf(0); // weight of the uninformative pairs
  wUninf.reserve(n_pair);
  
  vector<double> iScore(4); // temporary store results
  vector<double> iScoreM1(4); // temporary store results
  arma::mat matPairScore; // score of all pairs
  if(keepScore){
    matPairScore.resize(n_pair, 11); // store results from all scores
  }else{
    matPairScore.resize(0, 0); 
  }

  double weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative;
  
  bool updateIndexNeutral = moreEndpoint && neutralAsUninf;
  bool updateIndexUninf = moreEndpoint;
  bool test_tauM1 = survTimeT_M1.n_cols>1;
  double zeroPlus = pow(10.0,-12.0);
  
  // ** loop over the neutral pairs
  for(int iter_pair=0; iter_pair<n_pair ; iter_pair++){

    // find index of the pair
    if(iter_pair<nNeutral_pairs_M1){
      iter_T = index_neutralT_M1[iter_pair];
      iter_C = index_neutralC_M1[iter_pair];
    }else{
      iter_T = index_uninfT_M1[iter_pair];
      iter_C = index_uninfC_M1[iter_pair];
    }

    // score
    if(method == 1){
      iScore = calcOnePair_Continuous(Control[iter_C], Treatment[iter_T], threshold, cumWeight_M1(iter_pair));
      std::fill(iScoreM1.begin(), iScoreM1.end(), 0.0);
    }else if(method == 2){
      iScore = calcOnePair_TTEgehan(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, cumWeight_M1(iter_pair));
      std::fill(iScoreM1.begin(), iScoreM1.end(), 0.0);
    }else if(method == 3){
      iScore = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
				     deltaC[iter_C], deltaT[iter_T], threshold,
				     survTimeC.row(iter_C), survTimeT.row(iter_T),
				     survJumpC, survJumpT, lastSurvC, lastSurvT);

      if(test_tauM1){ // useless if pairs from a different outcome
	iScoreM1 = calcOneScore_TTEperon(Control[iter_C], Treatment[iter_T], 
					 deltaC[iter_C], deltaT[iter_T], threshold_M1,
					 survTimeC_M1.row(iter_C), survTimeT_M1.row(iter_T),
					 survJumpC_M1, survJumpT_M1, lastSurvC, lastSurvT);
      }else{
      std::fill(iScoreM1.begin(), iScoreM1.end(), 0.0);
      }
    }	

    // store
    weight_favorable = (iScore[0] - iScoreM1[0]) * cumWeight_M1(iter_pair);
    weight_unfavorable = (iScore[1] - iScoreM1[1]) * cumWeight_M1(iter_pair);
    weight_neutral = iScore[2] * cumWeight_M1(iter_pair); 
    weight_uninformative = iScore[3] * cumWeight_M1(iter_pair);

    count_favorable += weight_favorable;
    count_unfavorable += weight_unfavorable;
    count_neutral += weight_neutral;
    count_uninf += weight_uninformative;
      
    if(updateIndexNeutral && iScore[2] > zeroPlus){
      index_neutralC.push_back(iter_C); // index of the pair relative to Control         
      index_neutralT.push_back(iter_T); // index of the pair relative to Treatment
      index_wNeutral.push_back(iter_pair); // index of the pair relative to cumWeight_M1
      wNeutral.push_back(iScore[2]); // not weight_neutral since the product is done in BuyseTest.cpp
    }
    if(updateIndexUninf && iScore[3] > zeroPlus){
      index_uninfC.push_back(iter_C); // index of the pair relative to Control    
      index_uninfT.push_back(iter_T); // index of the pair relative to Treatment
      index_wUninf.push_back(iter_pair); // index of the pair relative to cumWeight_M1
      wUninf.push_back(iScore[3]); // not weight_uninformative since the product is done in BuyseTest.cpp
    }
    if(keepScore){
      matPairScore.row(iter_pair) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
					    iScore[0] - iScoreM1[0], // favorable
					    iScore[1] - iScoreM1[1], // unfavorable
					    iScore[2], // neutral
					    iScore[3], // uninformative
					    cumWeight_M1(iter_pair), // weight
					    weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative // favorable corrected, unfavorable corrected, neutral corrected, uninformative corrected
	});
    }      
      
  }

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(correctionUninf > 0 && count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
      bool updateWeight = neutralAsUninf && moreEndpoint;
      bool updateIndex = moreEndpoint;

    if(correctionUninf == 1){

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfC, index_uninfT, 
		      index_neutralC, index_neutralT, 
		      wNeutral, index_wNeutral, wUninf, index_wUninf,
		      weight, index_weight,
		      updateWeight, updateIndex, keepScore, matPairScore);
      // updated by reference

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		    index_uninfC, index_uninfT, 
		    wNeutral,
		    weight,
		    updateWeight, updateIndex, keepScore, matPairScore);
      index_weight = arma::conv_to<uvec>::from(index_wNeutral);
      // updated by reference      
    }
  }else{
    wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
    weight = arma::conv_to<vec>::from(wNeutral);
    index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
    index_weight = arma::conv_to<uvec>::from(index_wNeutral);
  }

  // ** export 
  return matPairScore;
  
}



// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		     vector<int>& index_uninfC, vector<int>& index_uninfT, 
		     vector<int>& index_neutralC, vector<int>& index_neutralT,
		     const vector<double>& wNeutral, const vector<int>& index_wNeutral, const vector<double>& wUninf, const vector<int>& index_wUninf,
		     arma::vec& weight, arma::uvec& index_weight,
		     bool updateWeight, bool updateIndex, bool keepScore, arma::mat& matPairScore){

  // compute factor
  double factorFavorable = (count_favorable)/(count_favorable + count_unfavorable + count_neutral); 
  double factorUnfavorable = (count_unfavorable)/(count_favorable + count_unfavorable + count_neutral); 
  double factorNeutral  = (count_neutral)/(count_favorable + count_unfavorable + count_neutral); 

  // update global score
  count_favorable += factorFavorable * count_uninf;    
  count_unfavorable += factorUnfavorable * count_uninf;    
  count_neutral += factorNeutral * count_uninf;   
  count_uninf = 0;

  // update neutral pairs
  // note it is enough to match the pairs on the treatment arm
  // if neutralAsUninf is false then the neutral pairs are not used at the following endpoints so no need to update
  if(updateWeight){

    int n_neutralT = index_neutralT.size();
    int n_uninfT = index_uninfT.size();
    int newSize = n_uninfT + n_neutralT; // too long
    vector<int> indexNew_neutralC(newSize);
    vector<int> indexNew_neutralT(newSize);
    weight.resize(newSize);
    if(updateIndex){
      index_weight.resize(newSize);
    }

    int iPair = 0;
    int iNeutral = 0;
    int iUninf = 0;

    while(iUninf < n_uninfT || iNeutral < n_neutralT){
	  
      if(iUninf >= n_uninfT || (iNeutral < n_neutralT && index_neutralT[iNeutral]<index_uninfT[iUninf]) || (iNeutral < n_neutralT && index_neutralT[iNeutral]==index_uninfT[iUninf] && index_neutralC[iNeutral]<index_uninfC[iUninf])){
	// current pair is neutral
	indexNew_neutralC[iPair] = index_neutralC[iNeutral];
	indexNew_neutralT[iPair] = index_neutralT[iNeutral];
	weight(iPair) = wNeutral[iNeutral];
	if(updateIndex){
    	  index_weight(iPair) = index_wNeutral[iNeutral];
	}
	iNeutral ++;
		
      }else if(iUninf < n_uninfT && iNeutral < n_neutralT && index_neutralT[iNeutral]==index_uninfT[iUninf] && index_neutralC[iNeutral]==index_uninfC[iUninf]){
	// current pair is both neutral and uninformative: increase the weight
	indexNew_neutralC[iPair] = index_neutralC[iNeutral];
	indexNew_neutralT[iPair] = index_neutralT[iNeutral];
	weight(iPair) = wNeutral[iNeutral] + factorNeutral * wUninf[iUninf];
	if(updateIndex){
    	  index_weight(iPair) = index_wNeutral[iNeutral];
	}
	iNeutral ++;
	iUninf ++;

      }else{
	// current pair is uninformative
	indexNew_neutralC[iPair] = index_uninfC[iUninf];
	indexNew_neutralT[iPair] = index_uninfT[iUninf];
	weight(iPair) = factorNeutral * wUninf[iUninf];
	if(updateIndex){
    	  index_weight(iPair) = index_wUninf[iUninf];
	}
	iUninf ++;
      }
	  
      iPair ++;
    }

    // remove empty values in the vector
    // they were all initialized assuming no pair was both neutral and uninf
    if(iPair < newSize){
      indexNew_neutralC.erase(indexNew_neutralC.begin() + iPair, indexNew_neutralC.end());
      indexNew_neutralT.erase(indexNew_neutralT.begin() + iPair, indexNew_neutralT.end());
      weight.resize(iPair);
      if(updateIndex){
	index_weight.resize(iPair);
      }
    }
	
    // update index
    index_neutralC = indexNew_neutralC;
    index_neutralT = indexNew_neutralT;
  }
	
  // remove uninformative pairs and associated weights
  index_uninfT.resize(0);
  index_uninfC.resize(0);
  
  // update keep scores
  if(keepScore){
    matPairScore.col(7) = matPairScore.col(4) + factorNeutral * matPairScore.col(5);
    matPairScore.col(8) = matPairScore.col(2) + factorFavorable * matPairScore.col(5);
    matPairScore.col(9) = matPairScore.col(3) + factorUnfavorable * matPairScore.col(5);
    matPairScore.col(10) = matPairScore.col(4) + factorNeutral * matPairScore.col(5);
  }
}

// * correctionIPW
void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		   vector<int>& index_uninfC, vector<int>& index_uninfT, 
		   const vector<double>& wNeutral, 
		   arma::vec& weight,  
		   bool updateWeight, bool updateIndex, bool keepScore, arma::mat& matPairScore){


  // compute factor
  double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    
  // update global score
  count_favorable *= factor;    
  count_unfavorable *= factor;    
  count_neutral *= factor;   
  count_uninf = 0;
    
  // new weights
  if(updateWeight){
    int size = wNeutral.size();
    
    weight.resize(size);
    for(int iIndex=0; iIndex < size; iIndex++){
      weight(iIndex) = wNeutral[iIndex] * factor;
    }
  }
  
  // remove uninformative pairs and associated weights
  if(updateIndex){
    index_uninfT.resize(0);
    index_uninfC.resize(0);
  }
  
  // update keep scores
  if(keepScore){
    matPairScore.col(7) = matPairScore.col(2) * factor;
    matPairScore.col(8) = matPairScore.col(3) * factor;
    matPairScore.col(9) = matPairScore.col(4) * factor;
    matPairScore.col(10) = 0.0;
  }
}
