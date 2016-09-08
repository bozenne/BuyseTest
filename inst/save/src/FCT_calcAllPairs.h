// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

struct structAllPairs {
  double count_favorable;
  double count_unfavorable;
  double count_neutral;
  double count_uninf;
  vector<int> index_neutralT;
  vector<int> index_neutralC;
  vector<int> index_uninfT;
  vector<int> index_uninfC;
  vector<double> w;
  vector<int> index_w;
};

//// SUMMARY : calcAll functions ////
structAllPairs calcAllPairs_BinaryOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control);

structAllPairs calcSubsetPairs_BinaryOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs);

structAllPairs calcSubsetPairs_WeightedBinaryOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs);

structAllPairs calcAllPairs_ContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold);

structAllPairs calcSubsetPairs_ContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs);

structAllPairs calcSubsetPairs_WeightedContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs);

structAllPairs calcAllPairs_TTEOutcome_Gehan_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC);

structAllPairs calcSubsetPairs_TTEOutcome_Gehan_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, const arma::colvec& deltaT, const arma::colvec& deltaC, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs);

structAllPairs calcAllPairs_TTEOutcome_PetoEfronPeron_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold,
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
const int PEP);

structAllPairs calcSubsetPairs_TTEOutcome_PetoEfronPeron_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
const vector<int>& index_neutralT, const vector<int>& index_neutralC, int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, double threshold_M1,  const arma::mat& matKMT_M1, const arma::mat& matKMC_M1, const int PEP);

//  fct1 : perform pairwise comparisons over all possible pairs for a binary endpoint //////////////////////////////////////////
inline structAllPairs calcAllPairs_BinaryOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  calcOnePair_BinaryOutcome_cpp(Treatment[iter_T],Control[iter_C],iter_T,iter_C,1,-1,
  count_favorable,count_unfavorable,count_neutral,count_uninf,
  index_uninfT,index_uninfC,index_neutralT,index_neutralC,
  NULL1_vector, NULL2_vector);
  
  }
  }
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = index_neutralT;
  res.index_neutralC = index_neutralC;
  res.index_uninfT = index_uninfT;
  res.index_uninfC = index_uninfC;
    
  return(res);
    
}

//  fct1bis : perform pairwise comparisons over a prespecified subset of pairs for a binary endpoint ////////////////////////////////////////////////////////////
inline structAllPairs calcSubsetPairs_BinaryOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_BinaryOutcome_cpp(Treatment[iter_T],Control[iter_C],iter_T,iter_C,1,-1,
      count_favorable,count_unfavorable,count_neutral,count_uninf,
      indexNew_uninfT,indexNew_uninfC,indexNew_neutralT,indexNew_neutralC,
      NULL1_vector, NULL2_vector);        
      
    }
    
  }
  
  //// loop over the uninformative pairs ////
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_BinaryOutcome_cpp(Treatment[iter_T],Control[iter_C],iter_T,iter_C,1,-1,
      count_favorable,count_unfavorable,count_neutral,count_uninf,
      indexNew_uninfT,indexNew_uninfC,indexNew_neutralT,indexNew_neutralC,
      NULL1_vector, NULL2_vector);
      
    }
    
  }
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = indexNew_neutralT;
  res.index_neutralC = indexNew_neutralC;
  res.index_uninfT = indexNew_uninfT;
  res.index_uninfC = indexNew_uninfC;
  
  return(res);
    
}

//  fct1ter : perform a weighted pairwise comparisons over a prespecified subset of pairs for a binary endpoint ///////////////////////////
inline structAllPairs calcSubsetPairs_WeightedBinaryOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs,
const arma::vec& Wpairs){
  
  int iter_T,iter_C;  // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_BinaryOutcome_cpp(Treatment[iter_T], Control[iter_C], iter_T, iter_C, Wpairs(iter_pairs), iter_pairs,
      count_favorable, count_unfavorable, count_neutral, count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      index_wUninf, index_wNeutral);
    }
    
  }
  
  //// loop over the uninformative pairs ////
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs];  // index of the control patient of the pair in the Control matrix
      
      calcOnePair_BinaryOutcome_cpp(Treatment[iter_T], Control[iter_C], iter_T, iter_C, Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs,
      count_favorable, count_unfavorable, count_neutral, count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      index_wUninf, index_wNeutral);                  
      
    }
    
  }
  
  //// export ////
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // merging
  
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = indexNew_neutralT;
  res.index_neutralC = indexNew_neutralC;
  res.index_uninfT = indexNew_uninfT;
  res.index_uninfC = indexNew_uninfC;
  //res.w = -1.0;
  res.index_w = index_wNeutral;
  
  return(res);
    
}


//  fct2 : perform pairwise comparisons over all possible pairs for a continuous endpoint ///////////////////////////////
inline structAllPairs calcAllPairs_ContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  calcOnePair_ContinuousOutcome_cpp(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  1, -1, 
  count_favorable, count_unfavorable, count_neutral,count_uninf,
  index_uninfT, index_uninfC, index_neutralT, index_neutralC,
  NULL1_vector, NULL2_vector);
  
  }
  }
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = index_neutralT;
  res.index_neutralC = index_neutralC;
  res.index_uninfT = index_uninfT;
  res.index_uninfC = index_uninfC;
  
  return(res);
    
}

//  fct2bis : perform pairwise comparisons over a prespecified subset of pairs for a continuous endpoint//////////////////////
inline structAllPairs calcSubsetPairs_ContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_ContinuousOutcome_cpp(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  1, -1, 
      count_favorable, count_unfavorable, count_neutral,count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      NULL1_vector, NULL2_vector);
      
    }
    
  }
  
  //// loop over the uninformative pairs ////
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix 
      
      calcOnePair_ContinuousOutcome_cpp(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  1, -1, 
      count_favorable, count_unfavorable, count_neutral,count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      NULL1_vector, NULL2_vector);
      
    }
    
  }
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = indexNew_neutralT;
  res.index_neutralC = indexNew_neutralC;
  res.index_uninfT = indexNew_uninfT;
  res.index_uninfC = indexNew_uninfC;
  
  return(res);
    
}


//  fct2ter : perform a weighted pairwise comparisons over a prespecified subset of pairs for a continuous endpoint //////////////////
inline structAllPairs calcSubsetPairs_WeightedContinuousOutcome_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs,
const arma::vec& Wpairs){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  vector<int> index_wNeutral(0); // index of the neutral and uninformative pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the neutral and uninformative pairs relative to Wpairs
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_ContinuousOutcome_cpp(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  Wpairs(iter_pairs), iter_pairs, 
      count_favorable, count_unfavorable, count_neutral,count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      index_wNeutral, index_wUninf);
      
    }
    
  }
  
  //// loop over the uninformative pairs ////
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      calcOnePair_ContinuousOutcome_cpp(Treatment[iter_T], Control[iter_C], threshold, iter_T, iter_C,  Wpairs(nNeutral_pairs+iter_pairs), nNeutral_pairs+iter_pairs, 
      count_favorable, count_unfavorable, count_neutral,count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      index_wNeutral, index_wUninf);
      
    }
    
  }
  
  //// export ////
  index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
  
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = indexNew_neutralT;
  res.index_neutralC = indexNew_neutralC;
  res.index_uninfT = indexNew_uninfT;
  res.index_uninfC = indexNew_uninfC;
  //res.w = -1.0;
  res.index_w = index_wNeutral;
  
  return(res);
  
}


//  fct3 : perform pairwise comparisons over all possible pairs for a TTE endpoint ///////////////////////////////////////////
inline structAllPairs calcAllPairs_TTEOutcome_Gehan_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  calcOnePair_TTEOutcome_Gehan_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C, 1, -1,  
  count_favorable, count_unfavorable, count_neutral, count_uninf,
  index_uninfT, index_uninfC, index_neutralT, index_neutralC,
  NULL1_vector, NULL2_vector);
  
  }}
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = index_neutralT;
  res.index_neutralC = index_neutralC;
  res.index_uninfT = index_uninfT;
  res.index_uninfC = index_uninfC;
   
  return(res);
    
}

//  fct3bis : perform pairwise comparisons over a prespecified subset of pairs for a TTE endpoint ///////////////////////////
inline structAllPairs calcSubsetPairs_TTEOutcome_Gehan_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, 
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs,
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0;  // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninformative pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<int> NULL1_vector(0); // only to match function arguments
  vector<int> NULL2_vector(0); // only to match function arguments
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
      
      calcOnePair_TTEOutcome_Gehan_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C, 1, -1,  
      count_favorable, count_unfavorable, count_neutral, count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      NULL1_vector, NULL2_vector);
    }
  }
  
  //// loop over the uninformative pairs ////
  if(nUninf_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
      iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
      iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
      
      calcOnePair_TTEOutcome_Gehan_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C, 1, -1,  
      count_favorable, count_unfavorable,  count_neutral, count_uninf,
      indexNew_uninfT, indexNew_uninfC, indexNew_neutralT, indexNew_neutralC,
      NULL1_vector, NULL2_vector);
    }
  }
  
  //// export ////    
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = indexNew_neutralT;
  res.index_neutralC = indexNew_neutralC;
  res.index_uninfT = indexNew_uninfT;
  res.index_uninfC = indexNew_uninfC;
 
  return(res);
    
}

//  fct3ter : perform a weighted pairwise comparisons over all possible pairs for a TTE endpoint //////////////////////
inline structAllPairs calcAllPairs_TTEOutcome_PetoEfronPeron_cpp( const arma::colvec& Treatment, const arma::colvec& Control, const double threshold,
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
const int PEP){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  vector<double> wUninf(0); // weights of the uninformative pairs
  
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  double weight_residual;  
 
 //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  if(PEP == 1){
  proba_threshold = calcOneProba_TTEOutcome_Peto_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
  matKMT, matKMC);  
  } else if(PEP == 2){
  proba_threshold = calcOneProba_TTEOutcome_Efron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
  matKMT, matKMC);
  } else{ // PEP == 3
  proba_threshold = calcOneProba_TTEOutcome_Peron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
  matKMT, matKMC);
  }
  
  if(proba_threshold[2]>0.5){ // i.e. test neutral == 1
  index_neutralT.push_back(iter_T);
  index_neutralC.push_back(iter_C);
  count_neutral=count_neutral+1; 
  }else{
    
    weight_residual=1-(proba_threshold[0]+proba_threshold[1]);
    
    if(proba_threshold[3]>0.5 && weight_residual>0){ // i.e. test uninformative == 1
    index_uninfT.push_back(iter_T);
    index_uninfC.push_back(iter_C);
        
    wUninf.push_back(weight_residual); 
    count_uninf=count_uninf+weight_residual;
    }
    
    count_favorable = count_favorable + proba_threshold[0];    
    count_unfavorable = count_unfavorable + proba_threshold[1];  
    
  }
  
  }
  }
  
  
  //// export ////
  structAllPairs res;
  res.count_favorable = count_favorable;
  res.count_unfavorable = count_unfavorable;
  res.count_neutral = count_neutral;
  res.count_uninf = count_uninf;
  res.index_neutralT = index_neutralT;
  res.index_neutralC = index_neutralC;
  res.index_uninfT = index_uninfT;
  res.index_uninfC = index_uninfC;
  res.w = wUninf;
  
  return(res);
}

//  fct3quater : perform a weighted pairwise comparisons over a prespecified subset of pairs for a TTE endpoint //////////////////
inline structAllPairs calcSubsetPairs_TTEOutcome_PetoEfronPeron_cpp(const arma::colvec& Treatment, const arma::colvec& Control, const double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
const vector<int>& index_neutralT, const vector<int>& index_neutralC, const int nNeutral_pairs, 
const vector<int>& index_uninfT, const vector<int>& index_uninfC, const int nUninf_pairs,
const arma::vec& Wpairs, const double threshold_M1, const arma::mat& matKMT_M1, const arma::mat& matKMC_M1, const int PEP){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  double count_favorable=0; // number of favorable pairs
  double count_unfavorable=0; // number of unfavorable pairs
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs 
  vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  vector<double> proba_threshold(4); // probaF, probaUF, test.neutral and test.uninformative for the current threhold
  vector<double> proba_thresholdM1(4); // probaF, probaUF, test.neutral and test.uninformative for the previous threhold
  double weight_residual;  
  
  bool test_tauM1 = matKMT_M1.n_cols>1; // test whether it is the first time that the endpoint is used
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(PEP == 1){
      proba_threshold = calcOneProba_TTEOutcome_Peto_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
      matKMT, matKMC);        
      }else if(PEP == 2){
      proba_threshold = calcOneProba_TTEOutcome_Efron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
      matKMT, matKMC);
      }else { // PEP == 3
      proba_threshold = calcOneProba_TTEOutcome_Peron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
      matKMT, matKMC);
      }
      
      if(test_tauM1){ // no useless if pairs from a different outcome
        if(PEP == 1){
            proba_thresholdM1 = calcOneProba_TTEOutcome_Peto_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
            matKMT_M1, matKMC_M1);
        }else if(PEP == 2){
            proba_thresholdM1 = calcOneProba_TTEOutcome_Efron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
            matKMT_M1, matKMC_M1);
        }else { // PEP == 3
            proba_thresholdM1 = calcOneProba_TTEOutcome_Peron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
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
      count_neutral=count_neutral+Wpairs(iter_pairs);    
      }else{

        weight_residual=1-(proba_threshold[0]+proba_threshold[1]);
        
         if(proba_threshold[3]>0.5 && weight_residual>0){ // i.e. test uninformative == 1
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);
        index_wUninf.push_back(iter_pairs);        
        
        wUninf.push_back(weight_residual); 
        count_uninf=count_uninf+Wpairs(iter_pairs)*weight_residual;
       }
       
        count_favorable = count_favorable + (proba_threshold[0] - proba_thresholdM1[0])*Wpairs(iter_pairs);    
        count_unfavorable = count_unfavorable + (proba_threshold[1] - proba_thresholdM1[1])*Wpairs(iter_pairs);  
   
//        if(proba_threshold[0]+proba_threshold[1]>1.00001 || proba_threshold[0]<-0.00001 || proba_threshold[1]<-0.00001){Rcout <<iter_pairs << "(Neutral) : Pf=" << proba_threshold[0] - proba_thresholdM1[0] << " Puf" << proba_threshold[1] - proba_thresholdM1[1] << endl;}
 //       if(weight_residual>1 || weight_residual < -0.00001){Rcout <<iter_pairs << "(Neutral) : w=" << weight_residual << endl;}

      }
      
    }
  }
  
     
    //// loop over the uninformative pairs ////
    if(nUninf_pairs>0){
      
      for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
        iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
        iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
        
        if(PEP == 1){
        proba_threshold = calcOneProba_TTEOutcome_Peto_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
        matKMT, matKMC);  
        }else if(PEP == 2){
        proba_threshold = calcOneProba_TTEOutcome_Efron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
        matKMT, matKMC);  
        }else { // PEP == 3
        proba_threshold = calcOneProba_TTEOutcome_Peron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold, iter_T, iter_C,
        matKMT, matKMC);            
        }
        
      if(test_tauM1){
        if(PEP == 1){
        proba_thresholdM1 = calcOneProba_TTEOutcome_Peto_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
        matKMT_M1, matKMC_M1);        
        }else if(PEP == 2){
        proba_thresholdM1 = calcOneProba_TTEOutcome_Efron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
        matKMT_M1, matKMC_M1);
        }else { // PEP == 3
        proba_thresholdM1 = calcOneProba_TTEOutcome_Peron_cpp(Treatment[iter_T], Control[iter_C], deltaT[iter_T], deltaC[iter_C], threshold_M1, iter_T, iter_C,
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
      count_neutral=count_neutral+Wpairs(nNeutral_pairs+iter_pairs);    
      }else{
        
        weight_residual=1-(proba_threshold[0]+proba_threshold[1]);
        if(proba_threshold[3]>0.5 && weight_residual>0){ // i.e. test uninformative == 1
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);
        index_wUninf.push_back(nNeutral_pairs+iter_pairs);
        
        wUninf.push_back(weight_residual); 
        count_uninf=count_uninf+Wpairs(nNeutral_pairs+iter_pairs)*weight_residual;
        }
        
        count_favorable = count_favorable + (proba_threshold[0] - proba_thresholdM1[0])*Wpairs(nNeutral_pairs+iter_pairs);    
        count_unfavorable = count_unfavorable + (proba_threshold[1] - proba_thresholdM1[1])*Wpairs(nNeutral_pairs+iter_pairs);  
        
     //    if(proba_threshold[0]+proba_threshold[1]>1.00001 || proba_threshold[0]<-0.00001 || proba_threshold[1]<-0.00001){Rcout <<iter_pairs << "(Uninf) : Pf=" << proba_threshold[0] - proba_thresholdM1[0] << " Puf" << proba_threshold[1] - proba_thresholdM1[1] << endl;}
  //      if(weight_residual>1 || weight_residual < -0.00001){Rcout <<iter_pairs << "(Uninf) : w=" << weight_residual <<" " << proba_threshold[0] << " " << proba_threshold[1] << endl;}

      }
        
    }
  }
      
      //// export ////
      wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end());
      index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());
      
      structAllPairs res;
      res.count_favorable = count_favorable;
      res.count_unfavorable = count_unfavorable;
      res.count_neutral = count_neutral;
      res.count_uninf = count_uninf;
      res.index_neutralT = indexNew_neutralT;
      res.index_neutralC = indexNew_neutralC;
      res.index_uninfT = indexNew_uninfT;
      res.index_uninfC = indexNew_uninfC;
      res.w = wNeutral;
      res.index_w = index_wNeutral;
      
      return(res);
}
