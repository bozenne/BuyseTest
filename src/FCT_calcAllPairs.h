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

arma::mat calcAllPairs_Continuous(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
								  const int correctionUninf,
								  double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
								  vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT,
								  bool neutralAsUninf, bool keepScore);

arma::mat calcSubsetPairs_Continuous(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
									 const int correctionUninf,
									 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
									 vector<int>& index_neutralC, vector<int>& index_neutralT, const int nNeutral_pairs, 
									 vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
									 const arma::vec& Wpairs, vector<int>& index_wNeutral,
									 bool neutralAsUninf, bool keepScore);

arma::mat calcAllPairs_TTEgehan(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
                                const arma::colvec& deltaC, const arma::colvec& deltaT, 
                                const int correctionUninf,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT, 
                                vector<double>& wNeutral, 
                                bool neutralAsUninf, bool keepScore);
arma::mat calcAllPairs_TTEperon( const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
                                 const arma::colvec& deltaC, const arma::colvec& deltaT, 
								 const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
								 double lastSurvC, double lastSurvT, 
                                 const int correctionUninf,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT, 
                                 vector<double>& wUninf,
                                 bool neutralAsUninf, bool keepScore);

arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold, 
                                   const arma::colvec& deltaC, const arma::colvec& deltaT, 
                                   const int correctionUninf,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralC, vector<int>& index_neutralT, const int nNeutral_pairs, 
                                   vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool neutralAsUninf, bool keepScore);

arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold, 
                                   const arma::colvec& deltaC, const arma::colvec& deltaT, 
								   const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
                                   double lastSurvC, double lastSurvT,
								   const int correctionUninf,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralC, vector<int>& index_neutralT, const int nNeutral_pairs, 
                                   vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1,
								   const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool neutralAsUninf, bool keepScore);

void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 vector<int>& index_uninfC, vector<int>& index_uninfT,
					 vector<int>& index_neutralC, vector<int>& index_neutralT, 
					 vector<double>& wNeutral, vector<int>& index_wNeutral, vector<double>& wUninf, vector<int>& index_wUninf,
					 bool neutralAsUninf, bool updateW, bool updateIndex, bool keepScore, arma::mat& score);

void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   vector<int>& index_uninfC, vector<int>& index_uninfT, 
				   vector<int>& index_neutralC, vector<int>& index_neutralT, 
				   vector<double>& wNeutral, vector<double>& wUninf,
				   bool updateW, bool keepScore, arma::mat& score);


// * calcAllPairs_Continuous
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs_Continuous( const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
								   const int correctionUninf,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT, 
                                   bool neutralAsUninf, bool keepScore){
  
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
      
      iRow = calcOnePair_Continuous(Control[iter_C], Treatment[iter_T], threshold, iter_C, iter_T, 1, -1, 
                                    count_favorable, count_unfavorable, count_neutral,count_uninf,
                                    index_neutralC, index_neutralT, index_uninfC, index_uninfT, 
                                    NULL1_vector, NULL2_vector,
                                    neutralAsUninf, keepScore);
      
      if(keepScore){
        score.row(iter_pairs) = iRow;
        iter_pairs++;
      }
      
    }
  }

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
	vector<double> wNeutral(0);
	vector<double> wUninf(0);

	  if(correctionUninf == 1){
      vector<int> index_wNeutral(0);
      vector<int> index_wUninf(0);

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT, 
					  index_neutralC, index_neutralT, 
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, false, false, keepScore, score);

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					false, keepScore, score);
    }
  }
  
  // ** export
  return score;
  
}


// * calcSubsetPairs_Continuous
// perform pairwise comparisons over the neutral and uniformative pairs for a continuous endpoint 
arma::mat calcSubsetPairs_Continuous( const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
									  const int correctionUninf,
                                      double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                      vector<int>& index_neutralC, vector<int>& index_neutralT,  const int nNeutral_pairs, 
                                      vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
                                      const arma::vec& Wpairs, vector<int>& index_wNeutral,
                                      bool neutralAsUninf, bool keepScore){
  
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
      
      iRow = calcOnePair_Continuous(Control[iter_C], Treatment[iter_T], threshold, iter_C, iter_T,
									Wpairs(iter_pairs), iter_pairs, 
                                    count_favorable, count_unfavorable, count_neutral, count_uninf,
                                    indexNew_neutralC, indexNew_neutralT, indexNew_uninfC, indexNew_uninfT,
                                    index_wNeutral, index_wUninf,
                                    neutralAsUninf, keepScore);
      
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
      
      iRow = calcOnePair_Continuous(Control[iter_C], Treatment[iter_T], threshold, iter_C, iter_T, 
                                    Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs, 
                                    count_favorable, count_unfavorable, count_neutral,count_uninf,
                                    indexNew_neutralC, indexNew_neutralT, indexNew_uninfC, indexNew_uninfT, 
                                    index_wNeutral, index_wUninf,
                                    neutralAsUninf, keepScore);      
      if(keepScore){
        score.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
    
  }

  // ** update index
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
	vector<double> wNeutral(0);
	vector<double> wUninf(0);

	if(correctionUninf == 1){

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT, 
					  index_neutralC, index_neutralT, 
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, false, true, keepScore, score);

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					false, keepScore, score);
    }
  }

  // ** export 
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // index_wNeutral is returned by reference

  return score;
  
}



// * calcAllPairs_TTEgehan
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEgehan(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
                                const arma::colvec& deltaC, const arma::colvec& deltaT, 
                                const int correctionUninf,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
								vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT,
                                vector<double>& wNeutral, 
                                bool neutralAsUninf, bool keepScore){
  
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

      iRow = calcOnePair_TTEgehan(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T, 1, -1,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  index_neutralC, index_neutralT, index_uninfC, index_uninfT, 
                                  NULL1_vector, NULL2_vector,
                                  neutralAsUninf, keepScore);
      
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
  
  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
    if(correctionUninf == 1){
      vector<int> index_wNeutral(0);
      vector<int> index_wUninf(0);

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT, 
					  index_neutralC, index_neutralT, 
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, true, false, keepScore, score);

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					true, keepScore, score);
    }
  }  
  
  // ** export
  // NOTE: no index_wNeutral to update since it is the first outcome
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return score;
}


// * calcAllPairs_TTEperon
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEperon( const arma::colvec& Control, const arma::colvec& Treatment, const double threshold,
                                 const arma::colvec& deltaC, const arma::colvec& deltaT, 
								 const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
								 double lastSurvC, double lastSurvT, 
                                 const int correctionUninf,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralC, vector<int>& index_neutralT, vector<int>& index_uninfC, vector<int>& index_uninfT,
                                 vector<double>& wNeutral, 
                                 bool neutralAsUninf, bool keepScore){
  
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
	  
      proba_threshold = calcOneProba_TTEperon(Control[iter_C], Treatment[iter_T], 
											  deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T,
											  survTimeC, survTimeT, survJumpC, survJumpT,
											  lastSurvC, lastSurvT);
      
      /* Rcout << endl; */
      count_favorable += proba_threshold[0];    
      count_unfavorable += proba_threshold[1];

      if(proba_threshold[2] > zeroPlus){ // i.e. test neutral == 1
        count_neutral += proba_threshold[2];        

		if(neutralAsUninf){
		  index_neutralC.push_back(iter_C);
		  index_neutralT.push_back(iter_T);
		  wNeutral.push_back(proba_threshold[2]);
		}
      }
        
      if(proba_threshold[3] > zeroPlus){
		wUninf.push_back(proba_threshold[3]);

		count_uninf += proba_threshold[3];
		index_uninfC.push_back(iter_C);
		index_uninfT.push_back(iter_T);
      }
        
        
      if(keepScore){
		score.row(iter_pairs) = rowvec({(double)iter_C, (double)iter_T, // indexC, indexT
			  proba_threshold[0], // favorable
			  proba_threshold[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3], // uninformative
			  1, // weight
			  proba_threshold[0], proba_threshold[1], proba_threshold[2] // unfavorable corrected, neutral corrected, uninformative corrected
			  });
      }

      iter_pairs++;
    }
  }
  
  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  // Rcout << "correction: " << correctionUninf << " | uninf=" << count_uninf << endl;
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
	
    if(correctionUninf == 1){
	  vector<int> index_wNeutral(0);
      vector<int> index_wUninf(0);

      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT,
					  index_neutralC, index_neutralT, 
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, true, false, keepScore, score);

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					true, keepScore, score);
	}
  }
  
  // ** export
  // NOTE: no index_wNeutral to update since it is the first outcome
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return score;
}




// * calcSubsetPairs_TTEgehan
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold, 
                                   const arma::colvec& deltaC, const arma::colvec& deltaT,
                                   const int correctionUninf,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralC, vector<int>& index_neutralT, const int nNeutral_pairs, 
                                   vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool neutralAsUninf, bool keepScore){

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
      
      iRow = calcOnePair_TTEgehan(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T,
								  Wpairs(iter_pairs), iter_pairs,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  indexNew_uninfC, indexNew_uninfT, indexNew_neutralC, indexNew_neutralT, 
                                  index_wNeutral, index_wUninf,
                                  neutralAsUninf, keepScore);
      
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
      
      iRow = calcOnePair_TTEgehan(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T,
                                  Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs,  
                                  count_favorable, count_unfavorable, count_neutral, count_uninf,
                                  indexNew_uninfC, indexNew_uninfT, indexNew_neutralC, indexNew_neutralT, 
                                  index_wNeutral, index_wUninf,
                                  neutralAsUninf, keepScore);
      
      if(keepScore){
        score.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
  }
  
  // ** update index
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;
  
  // ** update weights
  wNeutral.resize(index_neutralT.size());
  std::fill(wNeutral.begin(),wNeutral.end(),1.0);
  wUninf.resize(index_uninfT.size());
  std::fill(wUninf.begin(),wUninf.end(),1.0);

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
    if(correctionUninf == 1){
      correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT, 
					  index_neutralC, index_neutralT, 
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, true, true, keepScore, score);

    }else if(correctionUninf == 2){
      correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					true, keepScore, score);
    }
  }

  // ** export  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  
  return score;
}



// * calcSubsetPairs_TTEperon
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Control, const arma::colvec& Treatment, const double threshold, 
                                   const arma::colvec& deltaC, const arma::colvec& deltaT, 
                                   const arma::mat& survTimeC, const arma::mat& survTimeT, const arma::mat& survJumpC, const arma::mat& survJumpT,
								   double lastSurvC, double lastSurvT, 
                                   const int correctionUninf,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralC, vector<int>& index_neutralT, const int nNeutral_pairs, 
                                   vector<int>& index_uninfC, vector<int>& index_uninfT, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1,
								   const arma::mat& survTimeC_M1, const arma::mat& survTimeT_M1, const arma::mat& survJumpC_M1, const arma::mat& survJumpT_M1,
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool neutralAsUninf, bool keepScore){
  
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
      proba_threshold = calcOneProba_TTEperon(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T,
											  survTimeC, survTimeT, survJumpC, survJumpT,
											  lastSurvC, lastSurvT);
      
      if(test_tauM1){ // useless if pairs from a different outcome
		proba_thresholdM1 = calcOneProba_TTEperon(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold_M1, iter_C, iter_T, 
												  survTimeC_M1, survTimeT_M1, survJumpC_M1, survJumpT_M1,
												  lastSurvC, lastSurvT);
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
        count_neutral += weight_neutral;

		if(neutralAsUninf){
		  indexNew_neutralC.push_back(iter_C);
		  indexNew_neutralT.push_back(iter_T);
		  index_wNeutral.push_back(iter_pairs);
		  wNeutral.push_back(proba_threshold[2]); // not weight_neutral since the product is done in BuyseTest.cpp
		}
      }

      if(weight_uninformative > zeroPlus){
		count_uninf += weight_uninformative;

		indexNew_uninfC.push_back(iter_C);
		indexNew_uninfT.push_back(iter_T);
		index_wUninf.push_back(iter_pairs);        
        wUninf.push_back(proba_threshold[3]);  // not weight_uninformative since the product is done in BuyseTest.cpp
      }
        
      if(keepScore){
		score.row(iter_pairs) = rowvec({(double)iter_C, (double)iter_T, // indexT, indexC
			  proba_threshold[0] - proba_thresholdM1[0], // favorable
			  proba_threshold[1] - proba_thresholdM1[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3],  // uninformative,
			  Wpairs(iter_pairs), // weight
			  weight_favorable, weight_unfavorable, weight_neutral // unfavorable corrected, neutral corrected, uninformative corrected					       
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
      proba_threshold = calcOneProba_TTEperon(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold, iter_C, iter_T,
											  survTimeC, survTimeT, survJumpC, survJumpT,
											  lastSurvC, lastSurvT);            
      
      if(test_tauM1){
		proba_thresholdM1 = calcOneProba_TTEperon(Control[iter_C], Treatment[iter_T], deltaC[iter_C], deltaT[iter_T], threshold_M1, iter_C, iter_T, 
												  survTimeC_M1, survTimeT_M1, survJumpC_M1, survJumpT_M1,
												  lastSurvC, lastSurvT);  
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
        count_neutral += weight_neutral;

		if(neutralAsUninf){
		  indexNew_neutralC.push_back(iter_C);
		  indexNew_neutralT.push_back(iter_T);
		  index_wNeutral.push_back(nNeutral_pairs+iter_pairs);
		  wNeutral.push_back(proba_threshold[2]); // not weight_neutral since the product is done in BuyseTest.cpp
		}
      }

      if(weight_uninformative > zeroPlus){
		count_uninf += weight_uninformative;

		indexNew_uninfC.push_back(iter_C);
		indexNew_uninfT.push_back(iter_T);
		index_wUninf.push_back(nNeutral_pairs+iter_pairs);
		wUninf.push_back(proba_threshold[3]);  // not weight_uninformative since the product is done in BuyseTest.cpp
      }
        
      if(keepScore){
		score.row(nNeutral_pairs+iter_pairs) = rowvec({(double)iter_C, (double)iter_T,  // indexT, indexC
			  proba_threshold[0] - proba_thresholdM1[0], // favorable
			  proba_threshold[1] - proba_thresholdM1[1], // unfavorable
			  proba_threshold[2], // neutral
			  proba_threshold[3],  // uninformative,
			  Wpairs(nNeutral_pairs + iter_pairs), // weight
			  weight_favorable, weight_unfavorable, weight_neutral // unfavorable corrected, neutral corrected, uninformative corrected
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

  // ** correction for uninformative pairs
  // correction possible: if there are uninformative paris
  //                      if there are informative pairs
  if(count_uninf > 0 && (count_favorable + count_unfavorable + count_neutral) > 0){
	if(correctionUninf == 1){
	  correctionPairs(count_favorable, count_unfavorable, count_neutral, count_uninf,
					  index_uninfC, index_uninfT,
					  index_neutralC, index_neutralT,
					  wNeutral, index_wNeutral, wUninf, index_wUninf,
					  neutralAsUninf, true, true, keepScore, score);

	}else if(correctionUninf == 2){
	  correctionIPW(count_favorable, count_unfavorable, count_neutral, count_uninf,
					index_uninfC, index_uninfT, 
					index_neutralC, index_neutralT, 
					wNeutral, wUninf,
					true, keepScore, score);
	}
  }
  
  // ** export
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());

  return score;
}

// * correctionPairs
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
void correctionPairs(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					 vector<int>& index_uninfC, vector<int>& index_uninfT, 
					 vector<int>& index_neutralC, vector<int>& index_neutralT,
					 vector<double>& wNeutral, vector<int>& index_wNeutral, vector<double>& wUninf, vector<int>& index_wUninf,
					 bool neutralAsUninf, bool updateW, bool updateIndex, bool keepScore, arma::mat& score){

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
  // if neutralAsUninf is false then the neutral pairs are not used so no need to update
  if(neutralAsUninf && (updateW || updateIndex)){

	int n_neutralT = index_neutralT.size();
	int n_uninfT = index_uninfT.size();
	int newSize = n_uninfT + n_neutralT; // too long
	vector<int> indexNew_neutralC(newSize);
	vector<int> indexNew_neutralT(newSize);
	vector<double> wNeutralNew;
	if(updateW){
	  wNeutralNew.resize(newSize);
	}
	vector<int> indexNew_wNeutral;
	if(updateIndex){
	  indexNew_wNeutral.resize(newSize);
	}

	int iPair = 0;
	int iNeutral = 0;
	int iUninf = 0;

	while(iUninf < n_uninfT || iNeutral < n_neutralT){
	  
	  if(iUninf >= n_uninfT || (iNeutral < n_neutralT && index_neutralT[iNeutral]<index_uninfT[iUninf]) || (iNeutral < n_neutralT && index_neutralT[iNeutral]==index_uninfT[iUninf] && index_neutralC[iNeutral]<index_uninfC[iUninf])){
		// current pair is neutral
		indexNew_neutralC[iPair] = index_neutralC[iNeutral];
		indexNew_neutralT[iPair] = index_neutralT[iNeutral];
		if(updateW){
		  wNeutralNew[iPair] = wNeutral[iNeutral];
		}
		if(updateIndex){
    	  indexNew_wNeutral[iPair] = index_wNeutral[iNeutral];
		}
		iNeutral ++;
		
	  }else if(iUninf < n_uninfT && iNeutral < n_neutralT && index_neutralT[iNeutral]==index_uninfT[iUninf] && index_neutralC[iNeutral]==index_uninfC[iUninf]){
		// current pair is both neutral and uninformative: increase the weight
		indexNew_neutralC[iPair] = index_neutralC[iNeutral];
		indexNew_neutralT[iPair] = index_neutralT[iNeutral];
		if(updateW){
		  wNeutralNew[iPair] = wNeutral[iNeutral] + factorNeutral * wUninf[iUninf];
		}
		if(updateIndex){
    	  indexNew_wNeutral[iPair] = index_wNeutral[iNeutral];
		}
		iNeutral ++;
		iUninf ++;

      }else{
		// current pair is uninformative
		indexNew_neutralC[iPair] = index_uninfC[iUninf];
		indexNew_neutralT[iPair] = index_uninfT[iUninf];
		if(updateW){
		  wNeutralNew[iPair] = factorNeutral * wUninf[iUninf];
		}
		if(updateIndex){
    	  indexNew_wNeutral[iPair] = index_wUninf[iUninf];
		}
		iUninf ++;
      }
	  
	  iPair ++;
    }

	if(iPair < newSize){
	  indexNew_neutralC.erase(indexNew_neutralC.begin() + iPair, indexNew_neutralC.end());
	  indexNew_neutralT.erase(indexNew_neutralT.begin() + iPair, indexNew_neutralT.end());
	  wNeutralNew.erase(wNeutralNew.begin() + iPair, wNeutralNew.end());
	  indexNew_wNeutral.erase(indexNew_wNeutral.begin() + iPair, indexNew_wNeutral.end());
	}
	
	// update index
	index_neutralC = indexNew_neutralC;
	index_neutralT = indexNew_neutralT;
	if(updateW){
	  wNeutral = wNeutralNew;
	}
	if(updateIndex){
	  index_wNeutral = indexNew_wNeutral;
	}
  }
	
  // remove uninformative pairs and associated weights
  index_uninfT.resize(0);
  index_uninfC.resize(0);
  if(updateW){
	wUninf.resize(0);
  }
  if(updateIndex){
	index_wUninf.resize(0);
  }
  
  // update keep scores
  if(keepScore){
	score.col(7) = score.col(2) + factorFavorable * score.col(5);
	score.col(8) = score.col(3) + factorUnfavorable * score.col(5);
	score.col(9) = score.col(4) + factorNeutral * score.col(5);
  }
}

// * correctionIPW
void correctionIPW(double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				   vector<int>& index_uninfC, vector<int>& index_uninfT, 
				   vector<int>& index_neutralC, vector<int>& index_neutralT, 
				   vector<double>& wNeutral, vector<double>& wUninf,
				   bool updateW, bool keepScore, arma::mat& score){


  // compute factor
  double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    
  // update global score
  count_favorable *= factor;    
  count_unfavorable *= factor;    
  count_neutral *= factor;   
  count_uninf = 0;
    
  // update neutral pairs
  if(updateW){
	int n_wNeutral = wNeutral.size();
	for(int iNeutral = 0; iNeutral< n_wNeutral; iNeutral++){      
	  wNeutral[iNeutral] *= factor;
	}
  }
    
  // remove uninformative pairs and associated weights
  index_uninfT.resize(0);
  index_uninfC.resize(0);
  if(updateW){
	wUninf.resize(0);
  }
  
  // update keep scores
  if(keepScore){
	score.col(7) = score.col(2) * factor;
	score.col(8) = score.col(3) * factor;
	score.col(9) = score.col(4) * factor;
  }
}
