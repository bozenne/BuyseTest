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
vector<int> index_wUninf, vector<int>  index_wNeutral, int trace);


/////// fct1 - binary endpoint ///////////////////////////////////////////////////
// [[Rcpp::export]]
List Test_calcOnePair_BinaryOutcome_cpp(vector<int> endpoint_T, vector<int> endpoint_C, vector<double> Wpair, vector<int> iter_pair, int trace=2){
  
  // initialization
  int index_T=0, index_C=0 ;
  
  double count_favorable=0, count_unfavorable=0, count_neutral=0, count_uninf=0;
  vector<int> index_uninfT(0), index_uninfC(0), index_neutralT(0), index_neutralC(0),index_wUninf(0), index_wNeutral(0);
  
  
  // calc
  for(unsigned int iter_ex=0; iter_ex<endpoint_T.size() ; iter_ex++){
    calcOnePair_BinaryOutcome_cpp( endpoint_T[iter_ex],  endpoint_C[iter_ex],  index_T,  index_C,  Wpair[iter_ex],  iter_pair[iter_ex],
    count_favorable, count_unfavorable, count_neutral,  count_uninf,
    index_uninfT, index_uninfC, index_neutralT, index_neutralC,
    index_wUninf,  index_wNeutral);
    
    index_T++;
    index_C++;
  }
  
  // display
  if(trace>0){
  Test_Display_cpp(count_favorable, count_unfavorable, count_neutral, count_uninf,
  index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
  index_wUninf, index_wNeutral,trace);
  } 
 
  // export
  return(List::create(
    Named("endpoint_T")  = endpoint_T,
    Named("endpoint_C")  = endpoint_C,
    //      Named("index_T")  = index_T,
    //      Named("index_C")  = index_C,           
    Named("Wpair")  = Wpair,           
    Named("iter_pair")  = iter_pair,
    Named("count_favorable")  = count_favorable,
    Named("count_unfavorable")  = count_unfavorable,
    Named("count_neutral")  = count_neutral,
    Named("count_uninf")  = count_uninf,
    Named("index_uninfT")  = index_uninfT,
    Named("index_uninfC")  = index_uninfC,
    Named("index_neutralT")  = index_neutralT,
    Named("index_neutralC")  = index_neutralC,
    Named("index_wUninf")  = index_wUninf,
    Named("index_wNeutral")  = index_wNeutral
    ));
}


/////// fct2 - continuous endpoint ///////////////////////////////////////////////////
// [[Rcpp::export]]
List Test_calcOnePair_ContinuousOutcome_cpp(vector<double> endpoint_T, vector<double> endpoint_C, double threshold, vector<double> Wpair, vector<int> iter_pair, bool trace=2){
  
  // initialization
  int index_T=0, index_C=0 ;
  
  double count_favorable=0, count_unfavorable=0, count_neutral=0, count_uninf=0;
  vector<int> index_uninfT(0), index_uninfC(0), index_neutralT(0), index_neutralC(0),index_wUninf(0), index_wNeutral(0);
  
  
  // calc
  for(unsigned int iter_ex=0; iter_ex<endpoint_T.size() ; iter_ex++){
    calcOnePair_ContinuousOutcome_cpp( endpoint_T[iter_ex],  endpoint_C[iter_ex], threshold,  index_T,  index_C,  Wpair[iter_ex],  iter_pair[iter_ex],
    count_favorable, count_unfavorable, count_neutral,  count_uninf,
    index_uninfT, index_uninfC, index_neutralT, index_neutralC,
    index_wUninf,  index_wNeutral);
    
    index_T++;
    index_C++;
  }
  
  // display
  if(trace>0){
  Test_Display_cpp(count_favorable, count_unfavorable, count_neutral, count_uninf,
  index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
  index_wUninf, index_wNeutral,trace);
  }
  
  // export
  return(List::create(
    Named("endpoint_T")  = endpoint_T,
    Named("endpoint_C")  = endpoint_C,
    //      Named("index_T")  = index_T,
    //      Named("index_C")  = index_C,           
    Named("Wpair")  = Wpair,           
    Named("iter_pair")  = iter_pair,
    Named("count_favorable")  = count_favorable,
    Named("count_unfavorable")  = count_unfavorable,
    Named("count_neutral")  = count_neutral,
    Named("count_uninf")  = count_uninf,
    Named("index_uninfT")  = index_uninfT,
    Named("index_uninfC")  = index_uninfC,
    Named("index_neutralT")  = index_neutralT,
    Named("index_neutralC")  = index_neutralC,
    Named("index_wUninf")  = index_wUninf,
    Named("index_wNeutral")  = index_wNeutral
    ));
}


/////// fct3 - TTE endpoint ///////////////////////////////////////////////////
// [[Rcpp::export]]
List Test_calcOnePair_TTEOutcome_cpp(vector<double> endpoint_T, vector<double> endpoint_C, vector<double> delta_T, vector<double> delta_C,
arma::mat survival_T, arma::mat survival_C,
double threshold, vector<double> Wpair, vector<int> iter_pair,int type, bool trace=2){
  
  // initialization
  int index_T=0, index_C=0 ;
  
  double count_favorable=0, count_unfavorable=0, count_neutral=0, count_uninf=0;
  vector<int> index_uninfT(0), index_uninfC(0), index_neutralT(0), index_neutralC(0),index_wUninf(0), index_wNeutral(0);
  
  vector<double> res(4);
  
  // calc
  if(type==0){
 
    for(unsigned int iter_ex=0; iter_ex<endpoint_T.size() ; iter_ex++){
      
      calcOnePair_TTEOutcome_Gehan_cpp( endpoint_T[iter_ex],  endpoint_C[iter_ex], delta_T[iter_ex] , delta_C[iter_ex], threshold,  index_T,  index_C,  Wpair[iter_ex],  iter_pair[iter_ex],
      count_favorable, count_unfavorable, count_neutral,  count_uninf,
      index_neutralT, index_neutralC, index_uninfT, index_uninfC, 
      index_wUninf,  index_wNeutral);    
      
      index_T++;
      index_C++;
    }
  }else if(type==1){ // iter_pair stands for iter C and iter T   
    res  = calcOneProba_TTEOutcome_Peto_cpp(endpoint_T[0],  endpoint_C[0], delta_T[0] , delta_C[0], threshold, iter_pair[0], iter_pair[1],
    survival_T, survival_C);
  }else if(type==2){ // iter_pair stands for iter C and iter T
     res  = calcOneProba_TTEOutcome_Efron_cpp(endpoint_T[0],  endpoint_C[0], delta_T[0] , delta_C[0], threshold, iter_pair[0], iter_pair[1],
    survival_T, survival_C);
  }else if(type==3){
      res  = calcOneProba_TTEOutcome_Peron_cpp(endpoint_T[0],  endpoint_C[0], delta_T[0] , delta_C[0], threshold, iter_pair[0], iter_pair[1],
    survival_T, survival_C);
  }

  // display
   if(trace>0){
   if(type==0){
  Test_Display_cpp(count_favorable, count_unfavorable, count_neutral, count_uninf,
  index_uninfT, index_uninfC, index_neutralT, index_neutralC,
  index_wUninf, index_wNeutral,trace);
  }else{
    Rcout << "output (probaF / probaUF / test.neutral / test.uninf) : " ;
    Rcout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << endl;
  }
  }
  
  // export
  if(type==0){
    return(List::create(
      Named("endpoint_T")  = endpoint_T,
      Named("endpoint_C")  = endpoint_C,
      //      Named("index_T")  = index_T,
      //      Named("index_C")  = index_C,           
      Named("Wpair")  = Wpair,           
      Named("iter_pair")  = iter_pair,
      Named("count_favorable")  = count_favorable,
      Named("count_unfavorable")  = count_unfavorable,
      Named("count_neutral")  = count_neutral,
      Named("count_uninf")  = count_uninf,
      Named("index_uninfT")  = index_uninfT,
      Named("index_uninfC")  = index_uninfC,
      Named("index_neutralT")  = index_neutralT,
      Named("index_neutralC")  = index_neutralC,
      Named("index_wUninf")  = index_wUninf,
      Named("index_wNeutral")  = index_wNeutral
      ));
  }else{
    return(List::create(
      Named("res") = res
      ));
  }
}



/////// fct4  - Display ///////////////////////////////////////////////////

void Test_Display_cpp(double count_favorable, double count_unfavorable, double count_neutral, double count_uninf,
vector<int> index_neutralT, vector<int> index_neutralC, vector<int> index_uninfT, vector<int> index_uninfC, 
vector<int> index_wUninf, vector<int>  index_wNeutral, int trace){
  
  // display count
  Rcout << "count  (favorable / unfavorable / neutral / uninf) : " ;
  Rcout << count_favorable <<" " << count_unfavorable << " " << count_neutral << " " << count_uninf <<  endl;
  
  // display index
  Rcout << "index size (neutralT / neutralC / uninfT/ uninfC) : " ;
  Rcout << index_neutralT.size() << " " << index_neutralC.size() << " " <<  index_uninfT.size() <<" " << index_uninfC.size() << endl;
  
  if(index_neutralT.size()>0 && trace>0){   
    Rcout << "     # index_neutralT : " ;
    for(unsigned int iter=0 ; iter<index_neutralT.size() ; iter++){
      Rcout << index_neutralT[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_neutralC.size()>0 && trace>0){
    Rcout << "     # index_neutralC : " ;
    for(unsigned int iter=0 ; iter<index_neutralC.size() ; iter++){
      Rcout << index_neutralC[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_uninfT.size()>0 && trace>0){    
    Rcout << "     # index_uninfT : " ;
    for(unsigned int iter=0 ; iter<index_uninfT.size() ; iter++){
      Rcout << index_uninfT[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_uninfC.size()>0 && trace>0){     
    Rcout << "     # index_uninfC : " ;
    for(unsigned int iter=0 ; iter<index_uninfC.size() ; iter++){
      Rcout << index_uninfC[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  // display wIndex
  Rcout << "index size (wNeutral/ wUninf) : " ;
  Rcout << index_wNeutral.size() <<" " << index_wUninf.size() <<  endl;
  
  if(index_wNeutral.size()>0 && trace>0){  
    Rcout << "     # index_wUNeutral : " ;
    for(unsigned int iter=0 ; iter<index_wNeutral.size() ; iter++){
      Rcout << index_wNeutral[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  if(index_wUninf.size()>0 && trace>0){   
    Rcout << "     # index_wUninf : " ;
    for(unsigned int iter=0 ; iter<index_wUninf.size() ; iter++){
      Rcout << index_wUninf[iter] <<" " ;
    }    
    Rcout << endl;
  }
  
  
  
  return;
}