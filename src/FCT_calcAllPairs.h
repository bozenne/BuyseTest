// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

arma::mat calcAllPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                   bool keepScore);

arma::mat calcSubsetPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                      double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                      vector<int>& index_neutralT,  vector<int>& index_neutralC, const int nNeutral_pairs, 
                                      vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                      const arma::vec& Wpairs, vector<int>& index_wNeutral,
                                      bool keepScore);

arma::mat calcAllPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                const arma::colvec& deltaT, const arma::colvec& deltaC,
                                const int correctionTTE,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                vector<double>& wNeutral, 
                                bool keepScore);
arma::mat calcAllPairs_TTEperon( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                 const arma::colvec& deltaT, const arma::colvec& deltaC,
								 const arma::mat& survTimeT, const arma::mat& survTimeC, const arma::mat& survJumpC, const arma::mat& survJumpT,
                                 const int correctionTTE,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                 vector<double>& wUninf,
                                 bool keepScore);

arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
                                   const int correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepScore);
arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
								   const arma::mat& survTimeT, const arma::mat& survTimeC, const arma::mat& survJumpC, const arma::mat& survJumpT,
                                   const int correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1,
								   const arma::mat& survTimeT_M1, const arma::mat& survTimeC_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepScore);

void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 vector<int>& index_uninfT, vector<int>& index_uninfC,
					 vector<int>& index_neutralT, vector<int>& index_neutralC,
					 vector<double>& wNeutral, vector<double>& wUninf,
					 bool keepScore, arma::mat& score);

void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   vector<int>& index_uninfT, vector<int>& index_uninfC,
				   vector<int>& index_neutralT, vector<int>& index_neutralC,
				   vector<double>& wNeutral, vector<double>& wUninf,
				   bool keepScore, arma::mat& score);


// * calcAllPairs_Continuous
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                   bool keepScore){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  arma::rowvec iRow; // temporary store results
  arma::mat score;
  if(keepScore){
    score.set_size(n_Treatment * n_Control, 10); // store results from all scores
  }else{
    score.set_size(0, 0); 
  }
  
  
  // ** loop over the pairs
  int iter_pairs=0;
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
      
      iRow = calcOnePair_Continuous(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  1, -1, 
                                    count_favorable, count_unfavorable, count_neutral,count_uninf,
                                    index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                                    NULL1_vector, NULL2_vector,
                                    keepScore);
      
      if(keepScore){
        score.row(iter_pairs) = iRow;
        iter_pairs++;
      }
      
    }
  }
  
  // ** export
  return score;
  
}


// * calcSubsetPairs_Continuous
// perform pairwise comparisons over the neutral and uniformative pairs for a continuous endpoint 
arma::mat calcSubsetPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                      double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                      vector<int>& index_neutralT,  vector<int>& index_neutralC, const int nNeutral_pairs, 
                                      vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                      const arma::vec& Wpairs, vector<int>& index_wNeutral,
                                      bool keepScore){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  // vector<int> index_wNeutral(0); // index of the neutral and uninformative pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the neutral and uninformative pairs relative to Wpairs
  
  arma::rowvec iRow; // temporary store results
  arma::mat score;
  if(keepScore){
    score.set_size(nNeutral_pairs+nUninf_pairs, 10); // store results from all scores
  }else{
    score.set_size(0, 0); 
  }
  
  // ** loop over the neutral pairs
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      iRow = calcOnePair_Continuous(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  Wpairs(iter_pairs), iter_pairs, 
                                    count_favorable, count_unfavorable, count_neutral, count_uninf,
                                    indexNew_neutralT, indexNew_neutralC, indexNew_uninfT, indexNew_uninfC, 
                                    index_wNeutral, index_wUninf,
                                    keepScore);
      
      if(keepScore){
        score.row(iter_pairs) = iRow;
      }
      
    }
    
  }
  
  // ** loop over the uninformative pairs
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      iRow = calcOnePair_Continuous(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,
                                    Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs, 
                                    count_favorable, count_unfavorable, count_neutral,count_uninf,
                                    indexNew_neutralT, indexNew_neutralC, indexNew_uninfT, indexNew_uninfC, 
                                    index_wNeutral, index_wUninf,
                                    keepScore);      
      if(keepScore){
        score.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
    
  }
  
  // ** export 
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;
  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // index_wNeutral is returned by reference
  
  return score;
  
}



// * calcAllPairs_TTEgehan
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                const arma::colvec& deltaT, const arma::colvec& deltaC,
                                const int correctionTTE,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                vector<double>& wNeutral, 
                                bool keepScore){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  
  arma::mat score;
  int iter_pairs = 0;
  if(keepScore){
    score.set_size(n_Treatment * n_Control, 10); // store results from all scores
  }else{
    score.set_size(0, 0); 
  }
  //vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs
  arma::rowvec iRow; // temporary store results
  
  // ** loop over the pairs
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments  
  
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients

      iRow = calcOnePair_TTEgehan(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C, 1, -1,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                                  NULL1_vector, NULL2_vector,
                                  keepScore);
      
      if(keepScore){
        score.row(iter_pairs) = iRow;
        iter_pairs++;
      }
    }
  }
  
  // ** update weights
  wNeutral.resize(index_neutralT.size());
  std::fill(wNeutral.begin(),wNeutral.end(),1.0);
  wUninf.resize(index_uninfT.size());
  std::fill(wUninf.begin(),wUninf.end(),1.0);
  
  // ** correction
  if(correctionTTE == 1){
	correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfT, index_uninfC,
					index_neutralT, index_neutralC,
					wNeutral, wUninf,
					keepScore, score);

  }else if(correctionTTE == 3){
	correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
				  index_uninfT, index_uninfC,
				  index_neutralT, index_neutralC,
				  wNeutral, wUninf,
				  keepScore, score);
  }
  
  
  // ** export
  // NOTE: no index_wNeutral to update since it is the first outcome
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return score;
}


// * calcAllPairs_TTEperon
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEperon( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                 const arma::colvec& deltaT, const arma::colvec& deltaC,
								 const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
                                 const int correctionTTE,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                 vector<double>& wNeutral, 
                                 bool keepScore){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  arma::mat score;
  if(keepScore){
    score.set_size(n_Treatment * n_Control, 10); // store results from all scores
  }else{
    score.set_size(0, 0);
  }
  //vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs
  arma::rowvec iRow; // temporary store results
  int iter_pairs = 0;
  
  // ** loop over the pairs
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  double zeroPlus = pow(10.0,-12.0);

  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
      /* Rcout << iter_T << ";" << iter_C ; */
	  
	  proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C],
											  deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
											  survTimeC, survTimeT, survJumpC, survJumpT);
      
      /* Rcout << endl; */
	  count_favorable += proba_threshold[0];    
	  count_unfavorable += proba_threshold[1];

	  if(proba_threshold[2] > zeroPlus){ // i.e. test neutral == 1
        index_neutralT.push_back(iter_T);
        index_neutralC.push_back(iter_C);
        
        wNeutral.push_back(proba_threshold[2]);
        count_neutral += proba_threshold[2];        
	  }
        
	  if(proba_threshold[3] > zeroPlus){
		index_uninfT.push_back(iter_T);
		index_uninfC.push_back(iter_C);
          
		wUninf.push_back(proba_threshold[3]);
		count_uninf += proba_threshold[3];
	  }
        
        
	  if(keepScore){
		score.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C, // indexT, indexC
			  proba_threshold[0], // favorable
			  proba_threshold[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3], // uninformative
			  1, // weight
              NA_REAL, NA_REAL, NA_REAL // unfavorable corrected, neutral corrected, uninformative corrected
	    });
		}

		iter_pairs++;
        }
  }
  
 // ** correction
  if(correctionTTE == 1){
     correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		     keepScore, score);

  }else if(correctionTTE == 3){
     correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		   keepScore, score);
   }
  
  // ** export
  // NOTE: no index_wNeutral to update since it is the first outcome
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return score;
}




// * calcSubsetPairs_TTEgehan
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
                                   const int correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepScore){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  //vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs 
  // vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  arma::rowvec iRow; // temporary store results
  arma::mat score;
  if(keepScore){
    score.set_size(nNeutral_pairs+nUninf_pairs, 10); // store results from all scores
  }else{
    score.set_size(0, 0);
  }
  
  
  // ** Neutral pairs
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
      
      iRow = calcOnePair_TTEgehan(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C, Wpairs(iter_pairs), iter_pairs,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
                                  index_wNeutral, index_wUninf,
                                  keepScore);
      
      if(keepScore){
        score.row(iter_pairs) = iRow;
      }
    }
  }
  // ** Uninformative pairs
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
      
      iRow = calcOnePair_TTEgehan(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                  Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
                                  index_wNeutral, index_wUninf,
                                  keepScore);
      
      if(keepScore){
        score.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
  }
  
  
  // ** update weights
  wNeutral.resize(index_neutralT.size());
  std::fill(wNeutral.begin(),wNeutral.end(),1.0);
  wUninf.resize(index_uninfT.size());
  std::fill(wUninf.begin(),wUninf.end(),1.0);

  // ** update index
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;

  // ** correction
    if(correctionTTE == 1){
     correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		     keepScore, score);

  }else if(correctionTTE == 3){
     correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		   keepScore, score);
   }

  // ** export  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  
  return score;
}



// * calcSubsetPairs_TTEperon
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
                                   const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,								   
                                   const int correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1,
						           const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepScore){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  //vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs 
  // vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  arma::rowvec iRow; // temporary store results
  arma::mat score;
  if(keepScore){
    score.set_size(nNeutral_pairs+nUninf_pairs, 10); // store results from all scores
  }else{
    score.set_size(0, 0);
  }
  
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  vector<double> proba_thresholdM1(4); // probaF, probaUF, test.neutral and test.uninformative for the previous threhold
  double weight_favorable, weight_unfavorable, weight_neutral, weight_uninformative;  
  bool test_tauM1 = survTimeT_M1.n_cols>1; // test whether it is the first time that the endpoint is used
  double zeroPlus = pow(10.0,-12.0);

  
  // ** Neutral pairs
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix

	  // *** Compute probas
       proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                               survTimeC, survTimeT, survJumpC, survJumpT);
      
      if(test_tauM1){ // useless if pairs from a different outcome
          proba_thresholdM1 = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    survTimeC_M1, survTimeT_M1, survJumpC_M1, survJumpT_M1);
      }else{
        proba_thresholdM1[0] = 0;
        proba_thresholdM1[1] = 0;      
      }
	  
	  // *** update weights
	  // note: proba_thresholdM1[2] is 0 except in the case of no censoring.
	  //       In this case where all probs are 0 or 1 and there is no need to remove proba_thresholdM1[3].
        weight_favorable = (proba_threshold[0] - proba_thresholdM1[0]) * Wpairs(iter_pairs);
        weight_unfavorable = (proba_threshold[1] - proba_thresholdM1[1]) * Wpairs(iter_pairs);
        weight_neutral = proba_threshold[2] * Wpairs(iter_pairs); 
	    weight_uninformative = proba_threshold[3] * Wpairs(iter_pairs);

        count_favorable += weight_favorable;
        count_unfavorable += weight_unfavorable;

      if(weight_neutral > zeroPlus){
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
        index_wNeutral.push_back(iter_pairs);
        
        wNeutral.push_back(proba_threshold[2]); // not weight_neutral since the product is done in BuyseTest.cpp
        count_neutral += weight_neutral;
      }

	  if(weight_uninformative > zeroPlus){
		indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);
        index_wUninf.push_back(iter_pairs);        
         
        wUninf.push_back(proba_threshold[3]);  // not weight_uninformative since the product is done in BuyseTest.cpp
        count_uninf += weight_uninformative;
      }
        
	  if(keepScore){
		score.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C, // indexT, indexC
			  proba_threshold[0] - proba_thresholdM1[0], // favorable
			  proba_threshold[1] - proba_thresholdM1[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3],  // uninformative,
			  Wpairs(iter_pairs), // weight
			  NA_REAL, NA_REAL, NA_REAL // unfavorable corrected, neutral corrected, uninformative corrected					       
			  });        
      }
    }
  }
  
  // *** loop over the uninformative pairs
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix

 	  // *** Compute probas
        proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                survTimeC, survTimeT, survJumpC, survJumpT);            
      
      if(test_tauM1){
          proba_thresholdM1 = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    survTimeC_M1, survTimeT_M1, survJumpC_M1, survJumpT_M1);  
      }else{
        proba_thresholdM1[0] = 0;
        proba_thresholdM1[1] = 0;      
      }

  	  // *** update weights
	  // note: proba_thresholdM1[2] is 0 except in the case of no censoring.
	  //       In this case where all probs are 0 or 1 and there is no need to remove proba_thresholdM1[3].

        weight_favorable = (proba_threshold[0] - proba_thresholdM1[0]) * Wpairs(nNeutral_pairs+iter_pairs);    
        weight_unfavorable = (proba_threshold[1] - proba_thresholdM1[1]) * Wpairs(nNeutral_pairs+iter_pairs);
        weight_neutral = proba_threshold[2] * Wpairs(nNeutral_pairs+iter_pairs); // note: it is not a mistake that proba_thresholdM1 does not appear here; 
        weight_uninformative = proba_threshold[3] * Wpairs(nNeutral_pairs+iter_pairs); // note: it is not a mistake that proba_thresholdM1 does not appear here; 

        count_favorable += weight_favorable;    
        count_unfavorable += weight_unfavorable;

      if(weight_neutral > zeroPlus){
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
        index_wNeutral.push_back(nNeutral_pairs+iter_pairs);
        
        wNeutral.push_back(proba_threshold[2]); // not weight_neutral since the product is done in BuyseTest.cpp 
        count_neutral += weight_neutral;
      }

	  if(weight_uninformative > zeroPlus){
          indexNew_uninfT.push_back(iter_T);
          indexNew_uninfC.push_back(iter_C);
          index_wUninf.push_back(nNeutral_pairs+iter_pairs);
          
          wUninf.push_back(proba_threshold[3]);  // not weight_uninformative since the product is done in BuyseTest.cpp
          count_uninf += weight_uninformative;
        }
        
	  if(keepScore){
		score.row(nNeutral_pairs+iter_pairs) = rowvec({(double)iter_T, (double)iter_C,  // indexT, indexC
			  proba_threshold[0] - proba_thresholdM1[0], // favorable
			  proba_threshold[1] - proba_thresholdM1[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3],  // uninformative,
			  Wpairs(nNeutral_pairs + iter_pairs), // weight
			  NA_REAL, NA_REAL, NA_REAL // unfavorable corrected, neutral corrected, uninformative corrected
			  });        
      }
    }
  }    

  // Rcout << endl;
  
  // ** update index
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;

 // ** correction
  if(correctionTTE == 1){
     correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		     keepScore, score);

  }else if(correctionTTE == 3){
     correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
		      index_uninfT, index_uninfC,
		      index_neutralT, index_neutralC,
		      wNeutral, wUninf,
		   keepScore, score);
   }
  
  // ** export
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());

  return score;
}

// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		      vector<int>& index_uninfT, vector<int>& index_uninfC,
		      vector<int>& index_neutralT, vector<int>& index_neutralC,
		      vector<double>& wNeutral, vector<double>& wUninf,
		      bool keepScore, arma::mat& score){

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
    int n_neutralT = index_neutralT.size();
    int n_uninfT = index_uninfT.size();
    int iNeutral = 0;
    int iUninf = 0;
    
    while(iNeutral < n_neutralT && iUninf < n_uninfT){
      if(index_neutralT[iNeutral]<index_uninfT[iUninf]){
	iNeutral ++;
      }else if(index_neutralT[iNeutral]>index_uninfT[iUninf]){
	iUninf ++;	
      }else if(index_neutralT[iNeutral]==index_uninfT[iUninf]){
        wNeutral[iNeutral] += factorNeutral * wUninf[iUninf];
	iNeutral ++;
	iUninf ++;		
      }
	
    }

    // remove uninformative pairs and associated weights
    index_uninfT.resize(0);
    index_uninfC.resize(0);
    wUninf.resize(0);

    // update keep scores
    if(keepScore){
       score.col(6) = score.col(2) + factorFavorable * score.col(5);
       score.col(7) = score.col(3) + factorUnfavorable * score.col(5);
       score.col(8) = score.col(4) + factorNeutral * score.col(5);
    }
    
}
 
// * correctionIPW
void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
		   vector<int>& index_uninfT, vector<int>& index_uninfC,
		   vector<int>& index_neutralT, vector<int>& index_neutralC,
		   vector<double>& wNeutral, vector<double>& wUninf,
		   bool keepScore, arma::mat& score){


    // compute factor
    double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    
    // update global score
    count_favorable *= factor;    
    count_unfavorable *= factor;    
    count_neutral *= factor;   
    count_uninf = 0;
    
    // update neutral pairs
    int n_wNeutral = wNeutral.size();
    for(int iNeutral = 0; iNeutral< n_wNeutral; iNeutral++){      
        wNeutral[iNeutral] *= factor;
    }
    
    // remove uninformative pairs and associated weights
    index_uninfT.resize(0);
    index_uninfC.resize(0);
    wUninf.resize(0);

    // update keep scores
    if(keepScore){
       score.col(6) = score.col(2) * factor;
       score.col(7) = score.col(3) * factor;
       score.col(8) = score.col(4) * factor;
    }
}
