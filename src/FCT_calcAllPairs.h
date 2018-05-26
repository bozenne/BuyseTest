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
                                   bool keepComparison);

arma::mat calcSubsetPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                      double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                      vector<int>& index_neutralT,  vector<int>& index_neutralC, const int nNeutral_pairs, 
                                      vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                      const arma::vec& Wpairs, vector<int>& index_wNeutral,
                                      bool keepComparison);

arma::mat calcAllPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                const arma::colvec& deltaT, const arma::colvec& deltaC,
                                const bool correctionTTE,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                vector<double>& wNeutral, 
                                bool keepComparison);
arma::mat calcAllPairs_TTEperon( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                 const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
                                 const int methodTTE, const bool correctionTTE,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                 vector<double>& wUninf,
                                 bool keepComparison);

arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
                                   const bool correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepComparison);
arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
                                   const int methodTTE, const bool correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1, const arma::mat& matKMT_M1, const arma::mat& matKMC_M1,
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepComparison);

// * calcAllPairs_Continuous
// perform pairwise comparisons over all possible pairs for a continuous endpoints
arma::mat calcAllPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                   bool keepComparison){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  arma::rowvec iRow; // temporary store results
  arma::mat comparison;
  if(keepComparison){
    comparison.set_size(n_Treatment * n_Control, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0); 
  }
  
  
  // ** loop over the pairs
  int iter_pairs=0;
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
      
      iRow = calcOnePair_Continuous(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  1, -1, 
                                    count_favorable, count_unfavorable, count_neutral,count_uninf,
                                    index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
                                    NULL1_vector, NULL2_vector,
                                    keepComparison);
      
      if(keepComparison){
        comparison.row(iter_pairs) = iRow;
        iter_pairs++;
      }
      
    }
  }
  
  // ** export
  return comparison;
  
}


// * calcSubsetPairs_Continuous
// perform pairwise comparisons over the neutral and uniformative pairs for a continuous endpoint 
arma::mat calcSubsetPairs_Continuous( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                      double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                      vector<int>& index_neutralT,  vector<int>& index_neutralC, const int nNeutral_pairs, 
                                      vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                      const arma::vec& Wpairs, vector<int>& index_wNeutral,
                                      bool keepComparison){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  // vector<int> index_wNeutral(0); // index of the neutral and uninformative pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the neutral and uninformative pairs relative to Wpairs
  
  arma::rowvec iRow; // temporary store results
  arma::mat comparison;
  if(keepComparison){
    comparison.set_size(nNeutral_pairs+nUninf_pairs, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0); 
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
                                    keepComparison);
      
      if(keepComparison){
        comparison.row(iter_pairs) = iRow;
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
                                    keepComparison);      
      if(keepComparison){
        comparison.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
    
  }
  
  // ** export 
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;
  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // index_wNeutral is returned by reference
  
  return comparison;
  
}



// * calcAllPairs_TTEgehan
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                const arma::colvec& deltaT, const arma::colvec& deltaC,
                                const bool correctionTTE,
                                double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                vector<double>& wNeutral, 
                                bool keepComparison){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  
  arma::mat comparison;
  int iter_pairs = 0;
  if(keepComparison){
    comparison.set_size(n_Treatment * n_Control, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0); 
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
                                  keepComparison);
      
      if(keepComparison){
        comparison.row(iter_pairs) = iRow;
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
  if(correctionTTE){
    double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    
    count_favorable *= factor;    
    count_unfavorable *= factor;    
    count_neutral *= factor;   
    count_uninf = 0;
    
    index_uninfT.resize(0);
    index_uninfC.resize(0);

    std::fill(wNeutral.begin(),wNeutral.end(),factor);
    wUninf.resize(0);

    // if(keepComparison){
	    // comparison.col(2) *= factor;
	    // comparison.col(3) *= factor;
	    // comparison.col(4) *= factor;
	    // comparison.col(5) *= 0;
	  // }
  }
  
  // ** export
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return comparison;
}


// * calcAllPairs_TTEperon
// perform pairwise comparisons over all possible pairs for a TTE endpoint
arma::mat calcAllPairs_TTEperon( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
                                 const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
                                 const int methodTTE, const bool correctionTTE,
                                 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC,
                                 vector<double>& wNeutral, 
                                 bool keepComparison){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  arma::mat comparison;
  if(keepComparison){
    comparison.set_size(n_Treatment * n_Control, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0);
  }
  //vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs
  arma::rowvec iRow; // temporary store results
  int iter_pairs = 0;
  
  // ** loop over the pairs
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  double weight_residual;   
  
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
    for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
      
      if(methodTTE == 1){
        proba_threshold = calcOneProba_TTEpeto(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                               matKMT, matKMC);  
      } else if(methodTTE == 2){
        proba_threshold = calcOneProba_TTEefron(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);
      } else{ // methodTTE == 3
        proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);
      }
      
      if(proba_threshold[2]>0.5){ // i.e. test neutral == 1
        index_neutralT.push_back(iter_T);
        index_neutralC.push_back(iter_C);
        
        wNeutral.push_back(1);
        count_neutral++;
        
        if(keepComparison){
          comparison.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C, 0, 0, 1, 0});
          iter_pairs++;
        }
        
      }else{
        
        weight_residual = 1 - (proba_threshold[0] + proba_threshold[1]);
        
        if(proba_threshold[3]>0.5 && weight_residual>pow(10.0,-12.0)){ // i.e. test uninformative == 1
          index_uninfT.push_back(iter_T);
          index_uninfC.push_back(iter_C);
          
          wUninf.push_back(weight_residual); 
          count_uninf += weight_residual;
        }
        
        count_favorable += proba_threshold[0];    
        count_unfavorable += proba_threshold[1];
        
        if(keepComparison){
          comparison.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C,
                         proba_threshold[0], proba_threshold[1], 0, weight_residual});
          iter_pairs++;
        }
      }
      
    }
  }
  
  
  // ** correction
  if(correctionTTE){
    double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);

    count_favorable *= factor;
    count_unfavorable *= factor;
    count_neutral *= factor;
    count_uninf = 0;
    
    index_uninfT.resize(0);
    index_uninfC.resize(0);
    
    for(unsigned int iNeutral=0; iNeutral<wNeutral.size() ; iNeutral++){
      wNeutral[iNeutral] *= factor;
    }
    wUninf.resize(0);

    // if(keepComparison){
	    // comparison.col(2) *= factor;
	    // comparison.col(3) *= factor;
	    // comparison.col(4) *= factor;
	    // comparison.col(5) *= 0;
	  // }
  }
  
  // ** export
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return comparison;
}




// * calcSubsetPairs_TTEgehan
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEgehan(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC,
                                   const bool correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepComparison){
  
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
  arma::mat comparison;
  if(keepComparison){
    comparison.set_size(nNeutral_pairs+nUninf_pairs, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0);
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
                                  keepComparison);
      
      if(keepComparison){
        comparison.row(iter_pairs) = iRow;
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
                                  keepComparison);
      
      if(keepComparison){
        comparison.row(nNeutral_pairs + iter_pairs) = iRow;
      }
    }
  }
  
  
  // ** update weights
  wNeutral.resize(index_neutralT.size());
  std::fill(wNeutral.begin(),wNeutral.end(),1.0);
  wUninf.resize(index_uninfT.size());
  std::fill(wUninf.begin(),wUninf.end(),1.0);
  
  // ** correction
  if(correctionTTE){
    double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    
    count_favorable *= factor;
    count_unfavorable *= factor;
    count_neutral *= factor;
    count_uninf = 0;
	
    index_uninfT.resize(0);
    index_uninfC.resize(0);

    std::fill(wNeutral.begin(),wNeutral.end(),factor);     
    wUninf.resize(0);

    // if(keepComparison){
	    // comparison.col(2) *= factor;
	    // comparison.col(3) *= factor;
	    // comparison.col(4) *= factor;
	    // comparison.col(5) *= 0;
	  // }
  }
  
  // ** export
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;
  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  
  return comparison;
}



// * calcSubsetPairs_TTEperon
// perform pairwise comparisons over the neutral and uniformative pairs for a TTE endpoint
arma::mat calcSubsetPairs_TTEperon(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
                                   const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
                                   const int methodTTE, const bool correctionTTE,
                                   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
                                   vector<int>& index_neutralT, vector<int>& index_neutralC, const int nNeutral_pairs, 
                                   vector<int>& index_uninfT, vector<int>& index_uninfC, const int nUninf_pairs,
                                   const arma::vec& Wpairs, const double threshold_M1, const arma::mat& matKMT_M1, const arma::mat& matKMC_M1, 
                                   vector<double>& wNeutral, vector<int>& index_wNeutral,
                                   bool keepComparison){
  
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
  arma::mat comparison;
  if(keepComparison){
    comparison.set_size(nNeutral_pairs+nUninf_pairs, 6); // store results from all comparisons
  }else{
    comparison.set_size(0, 0);
  }
  
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  vector<double> proba_thresholdM1(4); // probaF, probaUF, test.neutral and test.uninformative for the previous threhold
  double weight_residual, weight_favorable, weight_unfavorable;  
  bool test_tauM1 = matKMT_M1.n_cols>1; // test whether it is the first time that the endpoint is used
  
  
  // ** Neutral pairs
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(methodTTE == 1){
        proba_threshold = calcOneProba_TTEpeto(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                               matKMT, matKMC);        
      }else if(methodTTE == 2){
        proba_threshold = calcOneProba_TTEefron(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);
      }else { // methodTTE == 3
        proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);
      }
      
      if(test_tauM1){ // useless if pairs from a different outcome
        if(methodTTE == 1){
          proba_thresholdM1 = calcOneProba_TTEpeto(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                   matKMT_M1, matKMC_M1);
        }else if(methodTTE == 2){
          proba_thresholdM1 = calcOneProba_TTEefron(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    matKMT_M1, matKMC_M1);
        }else { // methodTTE == 3
          proba_thresholdM1 = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    matKMT_M1, matKMC_M1);
        }
      }else{
        proba_thresholdM1[0] = 0;
        proba_thresholdM1[1] = 0;      
      }
      
      if(proba_threshold[2]>0.5){ // i.e. test neutral == 1
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
        index_wNeutral.push_back(iter_pairs);
        
        wNeutral.push_back(1); 
        count_neutral+=Wpairs(iter_pairs);
        if(keepComparison){
          comparison.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C, 0, 0, 1, 0});
        }
      }else{
        
        weight_residual = 1-(proba_threshold[0] + proba_threshold[1]); // note: it is not a mistake that proba_thresholdM1 does not appear here; 
        weight_favorable = (proba_threshold[0] - proba_thresholdM1[0])*Wpairs(iter_pairs);
        weight_unfavorable = (proba_threshold[1] - proba_thresholdM1[1])*Wpairs(iter_pairs);
        
        if(proba_threshold[3]>0.5 && weight_residual>pow(10.0,-12.0)){ // i.e. test uninformative == 1
          indexNew_uninfT.push_back(iter_T);
          indexNew_uninfC.push_back(iter_C);
          index_wUninf.push_back(iter_pairs);        
          
          wUninf.push_back(weight_residual); 
          count_uninf += Wpairs(iter_pairs)*weight_residual;
        }
        
        count_favorable += weight_favorable;
        count_unfavorable += weight_unfavorable;
        if(keepComparison){
          comparison.row(iter_pairs) = rowvec({(double)iter_T, (double)iter_C,
                         weight_favorable,
                         weight_unfavorable,
                         0,
                         Wpairs(iter_pairs)*weight_residual});
        }
      }
      
    }
  }
  
  // *** loop over the uninformative pairs
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(methodTTE == 1){
        proba_threshold = calcOneProba_TTEpeto(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                               matKMT, matKMC);  
      }else if(methodTTE == 2){
        proba_threshold = calcOneProba_TTEefron(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);  
      }else { // methodTTE == 3
        proba_threshold = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
                                                matKMT, matKMC);            
      }
      
      if(test_tauM1){
        if(methodTTE == 1){
          proba_thresholdM1 = calcOneProba_TTEpeto(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                   matKMT_M1, matKMC_M1);        
        }else if(methodTTE == 2){
          proba_thresholdM1 = calcOneProba_TTEefron(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    matKMT_M1, matKMC_M1);
        }else { // methodTTE == 3
          proba_thresholdM1 = calcOneProba_TTEperon(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
                                                    matKMT_M1, matKMC_M1);  
        }
      }else{
        proba_thresholdM1[0] = 0;
        proba_thresholdM1[1] = 0;      
      }
      
      if(proba_threshold[2]>0.5){ // i.e. test neutral == 1
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
        index_wNeutral.push_back(nNeutral_pairs+iter_pairs);
        
        wNeutral.push_back(1); 
        count_neutral += Wpairs(nNeutral_pairs+iter_pairs);
        if(keepComparison){
          comparison.row(nNeutral_pairs+iter_pairs) = rowvec({(double)iter_T, (double)iter_C, 0, 0, 1, 0});
        }
      }else{
        
        weight_residual = 1-(proba_threshold[0] + proba_threshold[1]); // note: it is not a mistake that proba_thresholdM1 does not appear here; 
        weight_favorable = (proba_threshold[0] - proba_thresholdM1[0])*Wpairs(nNeutral_pairs+iter_pairs);    
        weight_unfavorable = (proba_threshold[1] - proba_thresholdM1[1])*Wpairs(nNeutral_pairs+iter_pairs);
        
        if(proba_threshold[3]>0.5 && weight_residual>pow(10.0,-12.0)){ // i.e. test uninformative == 1
          indexNew_uninfT.push_back(iter_T);
          indexNew_uninfC.push_back(iter_C);
          index_wUninf.push_back(nNeutral_pairs+iter_pairs);
          
          wUninf.push_back(weight_residual); 
          count_uninf += Wpairs(nNeutral_pairs+iter_pairs)*weight_residual;
        }
        
        count_favorable += weight_favorable;    
        count_unfavorable += weight_unfavorable;
        if(keepComparison){
          comparison.row(nNeutral_pairs+iter_pairs) = rowvec({(double)iter_T, (double)iter_C,
                         weight_favorable,
                         weight_unfavorable,
                         0,
                         Wpairs(nNeutral_pairs+iter_pairs)*weight_residual});
        }
      }
    }
  }    
  
  // ** correction
  if(correctionTTE){
    double factor = (count_favorable + count_unfavorable + count_neutral + count_uninf)/(count_favorable + count_unfavorable + count_neutral);
    count_favorable *= factor;
    count_unfavorable *= factor;
    count_neutral *= factor;
    count_uninf = 0;

    index_uninfT.resize(0);
    index_uninfC.resize(0);
    
    for(unsigned int iNeutral=0; iNeutral<wNeutral.size() ; iNeutral++){
      wNeutral[iNeutral] *= factor;
    }
    wUninf.resize(0);

    // if(keepComparison){
	    // comparison.col(2) *= factor;
	    // comparison.col(3) *= factor;
	    // comparison.col(4) *= factor;
	    // comparison.col(5) *= 0;
	  // }
    
  }
  
  // ** export
  index_neutralT = indexNew_neutralT;
  index_neutralC = indexNew_neutralC;
  index_uninfT = indexNew_uninfT;
  index_uninfC = indexNew_uninfC;
  
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
  return comparison;
}
