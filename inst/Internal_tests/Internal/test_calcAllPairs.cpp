#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "../src/FCT_buyseTest.cpp"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

//// SUMMARY : test functions ////
void Test_Display_cpp(double count_favorable, double count_unfavorable, double count_neutral, double count_uninf,
vector<int> index_uninfT, vector<int> index_uninfC, vector<int> index_neutralT, vector<int> index_neutralC,
vector<double> weights, vector<int>  index_weights, int trace);


/////// fct1 - binary endpoint ///////////////////////////////////////////////////

// [[Rcpp::export]]
List Test_calcAllPairs_BinaryOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, int trace=2){
  
  List res;
  
  res = calcAllPairs_BinaryOutcome_cpp(Treatment, Control);
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
  }
  
  // export
  return(res);
}


// [[Rcpp::export]]
List Test_calcSubsetPairs_BinaryOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs, int trace=2){
  
  List res;
  
  res = calcSubsetPairs_BinaryOutcome_cpp(Treatment, Control,
  index_neutralT, index_neutralC, nNeutral_pairs,
  index_uninfT, index_uninfC,nUninf_pairs);
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
  }
  
  // export
  return(res);
}


// [[Rcpp::export]]
List Test_calcSubsetPairs_WeightedBinaryOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, int trace=2){
  
  List res;
  
  res = calcSubsetPairs_WeightedBinaryOutcome_cpp(Treatment, Control,
  index_neutralT, index_neutralC, nNeutral_pairs,
  index_uninfT, index_uninfC,nUninf_pairs,
  Wpairs);
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  res[9],vector<int>(0),trace);
  }
  
  // export
  return(res);
}


/////// fct2 - continuous endpoint ///////////////////////////////////////////////////

// [[Rcpp::export]]
List Test_calcAllPairs_ContinuousOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold, int trace=2){
  
  List res;
  
  res = calcAllPairs_ContinuousOutcome_cpp(Treatment, Control, threshold);
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
  }
  
  // export
  return(res);
}


// [[Rcpp::export]]
List Test_calcSubsetPairs_ContinuousOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs, int trace=2){
  
  List res;
  
  res = calcSubsetPairs_ContinuousOutcome_cpp(Treatment, Control, threshold,
  index_neutralT, index_neutralC, nNeutral_pairs,
  index_uninfT, index_uninfC,nUninf_pairs);
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
  }
  
  // export
  return(res);
}


// [[Rcpp::export]]
List Test_calcSubsetPairs_WeightedContinuousOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, int trace=2){
  
  List res;
  
    res = calcSubsetPairs_WeightedContinuousOutcome_cpp(Treatment, Control, threshold,
    index_neutralT, index_neutralC, nNeutral_pairs,
    index_uninfT, index_uninfC,nUninf_pairs,
    Wpairs);    
  
  // display
  if(trace>0){   
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  res[9],vector<int>(0),trace);
  }
  
  // export
  return(res);
}

/////// fct3 - TTE endpoint ///////////////////////////////////////////////////

// [[Rcpp::export]]
List Test_calcAllPairs_TTEOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold,
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC, int type, int trace=2){
  
  List res;
  
  if(type==0){
  res = calcAllPairs_TTEOutcome_Gehan_cpp(Treatment, Control, threshold,
  deltaT, deltaC);
  }else{
  res = calcAllPairs_TTEOutcome_PetoEfronPeron_cpp(Treatment, Control, threshold,
  deltaT, deltaC, matKMT, matKMC,
  type);
  }
  
  // display
  if(trace>0){   
    if(type==0){
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
    }else{
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  res[8],vector<int>(0),trace);
    }
  }
  
  // export
  return(res);
}


// [[Rcpp::export]]
List Test_calcSubsetPairs_TTEOutcome_cpp(const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC,
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, double threshold_M1, const arma::mat& matKMT_M1, const arma::mat& matKMC_M1, int type, int trace=2){
  
  List res;
  
  if(type==0){
  res = calcSubsetPairs_TTEOutcome_Gehan_cpp(Treatment, Control, threshold, 
  deltaT, deltaC,
  index_neutralT, index_neutralC, nNeutral_pairs,
  index_uninfT, index_uninfC,nUninf_pairs);
  }else{
    res = calcSubsetPairs_TTEOutcome_PetoEfron_cpp(Treatment, Control, threshold, 
    deltaT, deltaC, matKMT, matKMC,
    index_neutralT, index_neutralC, nNeutral_pairs, 
    index_uninfT, index_uninfC, nUninf_pairs,
    Wpairs, threshold_M1, matKMT_M1, matKMC_M1, type);
  }

  // display
  if(trace>0){   
    if(type==0){
  Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  vector<double>(0),vector<int>(0),trace);
    }else{
     Test_Display_cpp(as<double>(res[0]), as<double>(res[1]), as<double>(res[2]), as<double>(res[3]),
  res[4], res[5], res[6], res[7],
  res[8],res[9],trace);
    }
  }
  
  // export
  return(res);
}

/////// fct4  - Display ///////////////////////////////////////////////////

void Test_Display_cpp(double count_favorable, double count_unfavorable, double count_neutral, double count_uninf,
vector<int> index_neutralT, vector<int> index_neutralC, vector<int> index_uninfT, vector<int> index_uninfC,
vector<double> weights, vector<int>  index_weights, int trace){
  
  // display count
  Rcout << "count  (favorable / unfavorable / neutral / uninf) : " ;
  Rcout << count_favorable <<" " << count_unfavorable << " " << count_neutral << " " << count_uninf <<  endl;
  
  // display index
  Rcout << "index size (neutralT / neutralC / uninfT/ uninfC) : " ;
  Rcout << index_neutralT.size() << " " << index_neutralC.size() << " " <<  index_uninfT.size() <<" " << index_uninfC.size() << endl;
  
  if(index_neutralT.size()>0 && trace>1){   
    Rcout << "     # index_neutralT : " ;
    for(unsigned int iter=0 ; iter<index_neutralT.size() ; iter++){
      Rcout << index_neutralT[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_neutralC.size()>0 && trace>1){
    Rcout << "     # index_neutralC : " ;
    for(unsigned int iter=0 ; iter<index_neutralC.size() ; iter++){
      Rcout << index_neutralC[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_uninfT.size()>0 && trace>1){    
    Rcout << "     # index_uninfT : " ;
    for(unsigned int iter=0 ; iter<index_uninfT.size() ; iter++){
      Rcout << index_uninfT[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_uninfC.size()>0 && trace>1){     
    Rcout << "     # index_uninfC : " ;
    for(unsigned int iter=0 ; iter<index_uninfC.size() ; iter++){
      Rcout << index_uninfC[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  // display wIndex
  Rcout << "size (weights/ index_weights) : " ;
  Rcout << weights.size() <<" " << index_weights.size() <<  endl;
  
  if(weights.size()>0 && trace>1){  
    Rcout << "     # weights : " ;
    for(unsigned int iter=0 ; iter<weights.size() ; iter++){
      Rcout << weights[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_weights.size()>0 && trace>1){   
    Rcout << "     # index_weights : " ;
    for(unsigned int iter=0 ; iter<index_weights.size() ; iter++){
      Rcout << index_weights[iter]  <<" " ;
    }    
    Rcout << endl;
  }
  
  
  
  return;
}