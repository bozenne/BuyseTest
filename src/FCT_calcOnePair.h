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

inline arma::rowvec calcOnePair_Continuous(const double endpoint_T, const double endpoint_C, const double threshold,
					   const int index_T, const int index_C, const double Wpair, const int iter_pair, 
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
					   vector<int>& index_wUninf, vector<int>&  index_wNeutral,
					   bool keepComparison);
 
inline arma::rowvec calcOnePair_TTEgehan(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
		       		 const int index_T, const int index_C, const double Wpair, const int iter_pair, 
				 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
				 vector<int>& index_wUninf, vector<int>&  index_wNeutral,
				 bool keepComparison);
 
inline vector<double> calcOneProba_TTEpeto(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					   const arma::mat& survival_T, const arma::mat& survival_C);
 
inline vector<double> calcOneProba_TTEefron(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					    const arma::mat& survival_T, const arma::mat& survival_C);
 
inline vector<double> calcOneProba_TTEperon(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					    const arma::mat& survival_T, const arma::mat& survival_C);

inline double calcIntegralProba(const double time_begin, const arma::mat& survival, const int index_col);




// * calcOnePair_Continuous
inline arma::rowvec calcOnePair_Continuous(const double endpoint_T, const double endpoint_C, const double threshold,
					   const int index_T, const int index_C, const double Wpair, const int iter_pair, 
					   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
					   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
					   vector<int>& index_wUninf, vector<int>&  index_wNeutral,
					   bool keepComparison){

  arma::rowvec iRow;
  
  if(R_IsNA(endpoint_T) || R_IsNA(endpoint_C)){
    index_uninfT.push_back(index_T);
    index_uninfC.push_back(index_C);     
    count_uninf+=Wpair;
    
    if(iter_pair>=0){
      index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs 
    }
    if(keepComparison){
      iRow = {(double)index_T, (double)index_C, 0, 0, 0, Wpair};
    }
  }else{    
    double diff = endpoint_T-endpoint_C; // difference between the endpoints from the treatment and control patients of the pair
    
    if(diff>=threshold && diff>pow(10.0,-12.0)){ // diff>0 for threshold == 0
      count_favorable+=Wpair;
      if(keepComparison){
	iRow = {(double)index_T, (double)index_C, Wpair, 0, 0, 0};
      }
    }else if(diff<= -threshold && diff<pow(10.0,-12.0)){ // diff<0 for threshold == 0
      count_unfavorable+=Wpair;
      if(keepComparison){
	iRow = {(double)index_T, (double)index_C, 0, Wpair, 0, 0};
      }

    }else{
      index_neutralT.push_back(index_T);
      index_neutralC.push_back(index_C);
      count_neutral+=Wpair;
      
      if(iter_pair>=0){
        index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs      
      }
      if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, 0, 0, Wpair, 0};
      }

    }
    
  }
  
  return iRow;
  
}

// * calcOnePair_TTEgehan
inline arma::rowvec calcOnePair_TTEgehan(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
				 const int index_T, const int index_C, const double Wpair, const int iter_pair, 
				 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
				 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
				 vector<int>& index_wUninf, vector<int>&  index_wNeutral,
				 bool keepComparison){
  
  double diff = endpoint_T-endpoint_C; // difference between the endpoints from the treatment and control patients of the pair
  arma::rowvec iRow;

  if(delta_T==1){
    if(delta_C==1){
      
      if(diff>=threshold && diff>pow(10.0,-12.0)){  // (1,1) >= tau    : favorable
        count_favorable+=Wpair;
	if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, Wpair, 0, 0, 0};
        }
      }else if(diff<= -threshold && diff<pow(10.0,-12.0)){              // (1,1) <= -tau   : unfavorable
	count_unfavorable+=Wpair;
	if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, 0, Wpair, 0, 0};
        }
      }else{                                  //  (1,1)  ]-tau;tau[ : uninformative
	index_neutralT.push_back(index_T);
	index_neutralC.push_back(index_C);
	count_neutral+=Wpair;
      
	if(iter_pair>=pow(10.0,-12.0)){
	  index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs    
	}
	if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, 0, 0, Wpair, 0};
        }

      }      
      
    }else{ // delta_C==0
    
      if(diff<= -threshold && diff<pow(10.0,-12.0)){              // (1,0) <= -tau   : unfavorable
	count_unfavorable+=Wpair;           
	if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, 0, Wpair, 0, 0};
        }
      }else{                             //  (1,0)  ]-tau;+Inf[ : uninformative
	index_uninfT.push_back(index_T);
	index_uninfC.push_back(index_C);
	count_uninf+=Wpair;
    
	if(iter_pair>=0){
	  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs 
	}
	if(keepComparison){
	  iRow = {(double)index_T, (double)index_C, 0, 0, 0, Wpair};
        }

      } 
    }
    
  }else{ // delta_T==0
    if(delta_C){
    
      if(diff>=threshold && diff>pow(10.0,-12.0)){   // (0,1) > tau    : favorable
	count_favorable+=Wpair;      
	if(keepComparison){
	 iRow = {(double)index_T, (double)index_C, Wpair, 0, 0, 0};
        }
      }else{                 //  (1,0)  ]-Inf;+tau[ : uninformative
	index_uninfT.push_back(index_T);
	index_uninfC.push_back(index_C);
	count_uninf+=Wpair;
    
	if(iter_pair>=0){
	  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs
	}
	if(keepComparison){
	  iRow = {(double)index_T, (double)index_C, 0, 0, 0, Wpair};
        }
      }
    
    }else{ // delta_C==0                (0,0) ]-Inf;+Inf[ : uninformative
      index_uninfT.push_back(index_T);
      index_uninfC.push_back(index_C); 
      count_uninf+=Wpair;
  
      if(iter_pair>=0){
	index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs 
      }
	if(keepComparison){
	  iRow = {(double)index_T, (double)index_C, 0, 0, 0, Wpair};
        }
    }
  }
  /* Rcout << iRow[0] << " " << iRow[1] << " " <<  iRow[2] << " " << iRow[3] << " " << iRow[4] << " " << iRow[5]  << endl; */
  return iRow;
  
}

// * calcOneProba_TTEpeto
inline vector<double> calcOneProba_TTEpeto(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					   const arma::mat& survival_T, const arma::mat& survival_C){
  
  // survival_T [S(x_i-tau,x_i,x_i+tau)  (survival at treatment event times)
  // survival_C [S(y_i-tau,y_i,y_i+tau) (survival at control event times)
  double diff = endpoint_T-endpoint_C;
  vector<double> proba(4,0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative
   
  if(delta_T==1){
    if(delta_C==1){
      
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (1,1) >= tau
      
	proba[0] = 1.0; // favorable
	//      proba[1] = 0.0; // unfavorable
      
      }else if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,1) <= -tau
      
	//      proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable 
      
      }else{  // (1,1) [-tau,tau]
      
	//      proba[0] = 0.0; // favorable   
	//      proba[1] = 0.0; // unfavorable 
	proba[2] = 1.0; // neutral
      
      }      
      
    }else{ // deltaC[iter_C]==0
    
    
      if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,0) <= -tau
    
	//    proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable  
    
      }else{
	proba[3] = 1.0; // uninformative
       
	if(diff>threshold){ // (1,0) >= tau
      
	  if(R_IsNA(survival_T(index_T,2)) == false){
	    proba[0]=1.0-survival_T(index_T,0)/survival_C(index_C,1); // favorable
	    proba[1]=survival_T(index_T,2)/survival_C(index_C,1); // unfavorable             
	  }// else {
	  //      proba[0]=0.0; // favorable
	  //      proba[1]=0.0; // unfavorable 
	  // }
      
	}else{ // (1,0) [-tau ; tau]
      
	  if(R_IsNA(survival_T(index_T,2)) == false){      
	    //      proba[0]=0.0; // favorable
	    proba[1]=survival_T(index_T,2)/survival_C(index_C,1); // unfavorable  
	  }// else {
	  //      proba[0]=0.0; // favorable
	  //      proba[1]=0.0; // unfavorable 
	  // }
      
	}    
      
      } 
    }
    
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 
  
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (0,1) >= taux
    
	proba[0] = 1.0; // favorable     
	//    proba[1] = 0.0; // unfavorable    
    
      }else{
	proba[3] = 1.0; // uninformative
      
	if(diff< -threshold){  // (0,1) < tau
      
	  if(R_IsNA(survival_C(index_C,2)) == false){ 
	    proba[0]=survival_C(index_C,2)/survival_T(index_T,1); // favorable
	    proba[1]=1-survival_C(index_C,0)/survival_T(index_T,1); // unfavorable
	  }// else {
	  //      proba[0]=0.0; // favorable
	  //      proba[1]=0.0; // unfavorable 
	  // }
      
	}else{  // (0,1) [-tau,tau]
	  //      Rcout << " proba : "  << index_C << " " << index_T << " | " << survival_C(index_C,2) << " " << survival_T(index_T,1) << endl;
	  if(R_IsNA(survival_C(index_C,2)) == false){
	    proba[0]=survival_C(index_C,2)/survival_T(index_T,1); // favorable
	    //      proba[1]=0.0; // unfavorable             
	  }// else {
	  //      proba[0]=0.0; // favorable
	  //      proba[1]=0.0; // unfavorable 
	  // }
	}
      } 
    
    }else{ // delta_C==0
      proba[3] = 1.0; // uninformative
  
      if(diff>threshold){ // (0,0) > tau
    
	if(R_IsNA(survival_T(index_T,2)) == false){
	  proba[0]=1.0-0.5*survival_T(index_T,0)/survival_C(index_C,1); // favorable
	  proba[1]=0.5*survival_T(index_T,2)/survival_C(index_C,1); // unfavorable  
	}// else {
	//      proba[0]=0.0; // favorable
	//      proba[1]=0.0; // unfavorable 
	// }
    
      }else if(diff< -threshold){  // (0,0) < -tau
    
	if(R_IsNA(survival_C(index_C,2)) == false){
	  proba[0]=0.5*survival_C(index_C,2)/survival_T(index_T,1); // favorable
	  proba[1]=1.0-0.5*survival_C(index_C,0)/survival_T(index_T,1); // unfavorable        
	}// else {
	//      proba[0]=0.0; // favorable
	//      proba[1]=0.0; // unfavorable 
	// }
  
      }else{  // (0,0) [-tau;tau]
  
	if( (R_IsNA(survival_C(index_C,2))==false) && (R_IsNA(survival_T(index_T,2))==false)){
	  proba[0]=0.5*survival_C(index_C,2)/survival_T(index_T,1); // favorable
	  proba[1]=0.5*survival_T(index_T,2)/survival_C(index_C,1); // unfavorable  
	}// else {
	//      proba[0]=0.0; // favorable
	//      proba[1]=0.0; // unfavorable 
	// }
  
      }
  
    }}
  
  return(proba);  
}

// * calcOneProba_TTEefron
inline vector<double> calcOneProba_TTEefron(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					    const arma::mat& survival_T, const arma::mat& survival_C){
  
  // survival_T [S_T(x_i-tau,x_i,x_i+tau) S_C(x_i-tau,x_i,x_i+tau)] (survival at treatment event times)
  // survival_C [S_C(y_i-tau,y_i,y_i+tau) S_T(y_i-tau,y_i,y_i+tau)] (survival at control event times)
  
  double diff = endpoint_T-endpoint_C;
  vector<double> proba(4,0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative

  if(delta_T==1){
    if(delta_C==1){
      
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (1,1) >= tau
	proba[0] = 1.0; // favorable
	//      proba[1] = 0.0; // unfavorable
      
      }else if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,1) <= -tau
	//      proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable 
      
      }else{  // (1,1) [-tau,tau]
	proba[2] = 1.0; // neutral
      
      }      
      
    }else{ // deltaC[iter_C]==0
    
    
      if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,0) <= -tau
	//    proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable  
    
      }else{
      
	proba[3] = 1.0; // uninformative
            
	if(diff>threshold){ // (1,0) > tau

	  proba[0]=1.0-survival_T(index_T,3)/survival_C(index_C,1); // favorable 1-[Sc(x_i-taux)/Sc(y_j)]
	  proba[1]=survival_T(index_T,5)/survival_C(index_C,1); // unfavorable  [Sc(x_i+taux)/Sc(y_j)]          
        
	}else{ // (1,0) [-tau ; tau]
 
	  //      proba[0]=0.0; // favorable
	  proba[1]=survival_T(index_T,5)/survival_C(index_C,1); // unfavorable   [Sc(x_i+taux)/Sc(y_j)] 
      
	} 
      }
    }
    
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 
    
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (0,1) >= taux
	proba[0] = 1.0; // favorable     
	//    proba[1] = 0.0; // unfavorable    
    
      }else{
	proba[3] = 1.0; // uninformative
      
	if(diff< -threshold){  // (0,1) < tau
 
	  proba[0]=survival_C(index_C,5)/survival_T(index_T,1); // favorable [St(y_j+taux)/Sc(x_i)] 
	  proba[1]=1-survival_C(index_C,3)/survival_T(index_T,1); // unfavorable [St(y_j-taux)/St(x_i)] 
          
	}else{  // (1,0) [-tau,tau]

	  proba[0]=survival_C(index_C,5)/survival_T(index_T,1); // favorable [St(y_j+taux)/Sc(x_i)] 
	  //      proba[1]=0.0; // unfavorable
      
	}
      } 
    
    }else{ // delta_C==0
      proba[3] = 1.0; // uninformative

      if(diff>threshold){ // (0,0) > tau

	proba[0]=1.0-survival_T(index_T,3)/survival_C(index_C,1)-calcIntegralProba(endpoint_T-threshold, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable 1-[St(x_i-taux)/Sc(y_j)]-I/(St(x_i)*Sc(y_j)) 
	proba[1]=survival_T(index_T,5)/survival_C(index_C,1)+calcIntegralProba(endpoint_T+threshold, survival_C, 7)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable  
  

      }else if(diff< -threshold){  // (0,0) < -tau
	proba[0]=-calcIntegralProba(endpoint_C, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable
	proba[1]=1.0+calcIntegralProba(endpoint_C, survival_C, 7)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable        
  

      }else{  // (0,0) [-tau;tau]

	proba[0]=-calcIntegralProba(endpoint_C, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable
	proba[1]=survival_T(index_T,5)/survival_C(index_C,1)+calcIntegralProba(endpoint_T+threshold, survival_C, 7)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable  
  

      }
  
    }}
  
  return(proba);  
}

// * calcOneProba_TTEperon
inline vector<double> calcOneProba_TTEperon(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold, const int index_T, const int index_C,
					    const arma::mat& survival_T, const arma::mat& survival_C){
  
  // survival_T : survival at treatment event times
  //          [0-2] unordered St (treatment time - tau, treatment time, treatment time + tau)
  //          [3-5] unordered Sc (treatment time - tau, treatment time, treatment time + tau)
  //          [6] ordered St (treatment time)
  //          [7-8] ordered Sc (treatment time - tau, treatment time + tau)
  //          [9] orderd treatment survival times
  //          [10] orderd event type (1 : event, 0 : censoring)
  
  // survival_C : survival at control event times
  //          [0-2] unordered Sc (control time - tau, control time, control time + tau)
  //          [3-5] unordered St (control time - tau, control time, control time + tau)
  //          [6] ordered Sc (control time)
  //          [7-8] ordered St (control time - tau, contol time + tau)
  //          [9] orderd control survival times
  //          [10] orderd event type (1 : event, 0 : censoring)
    
  
  double diff = endpoint_T-endpoint_C;
  vector<double> proba(4,0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative

  if(delta_T==1){
    if(delta_C==1){
      
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (1,1) >= tau
	proba[0] = 1.0; // favorable
	//      proba[1] = 0.0; // unfavorable
      
      }else if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,1) <= -tau
	//      proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable 
      
      }else{  // (1,1) [-tau,tau]
	proba[2] = 1.0; // neutral
      
      }      
      
    }else{ // deltaC[iter_C]==0
        
      if(diff<= -threshold && diff<pow(10.0,-12.0)){ // (1,0) <= -tau
	//    proba[0] = 0.0; // favorable   
	proba[1] = 1.0; // unfavorable  
    
      }else{
      
	proba[3] = 1.0; // uninformative
            
	if(diff>threshold){ // (1,0) > tau

	  if(R_IsNA(survival_T(index_T,3))==false){
	    proba[0]=1.0-survival_T(index_T,3)/survival_C(index_C,1); // favorable 1-[Sc(x_i-taux)/Sc(y_j)]
	  }// else {
	  //      proba[0]=0.0; // favorable
	  // }
      
	  if(R_IsNA(survival_T(index_T,5))==false){ 
	    proba[1]=survival_T(index_T,5)/survival_C(index_C,1); // unfavorable  [Sc(x_i+taux)/Sc(y_j)]          
	  }// else {
	  //      proba[1]=0.0; // unfavorable 
	  // }
      
	}else{ // (1,0) [-tau ; tau]

	  //      proba[0]=0.0; // favorable

	  if(R_IsNA(survival_T(index_T,5))==false){       
	    proba[1]=survival_T(index_T,5)/survival_C(index_C,1); // unfavorable   [Sc(x_i+taux)/Sc(y_j)] 
	  }// else {
	  //      proba[1]=0.0; // unfavorable 
	  // }    
      
	} 
      }
    }
    
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 
    
      if(diff>=threshold && diff>pow(10.0,-12.0)){ // (0,1) >= taux
	proba[0] = 1.0; // favorable     
	//    proba[1] = 0.0; // unfavorable    
    
      }else{
	proba[3] = 1.0; // uninformative
      
	if(diff< -threshold){  // (0,1) < tau

	  if(R_IsNA(survival_C(index_C,5))==false){
	    proba[0]=survival_C(index_C,5)/survival_T(index_T,1); // favorable [St(y_j+taux)/Sc(x_i)]   
	  } // else {
	  //      proba[0]=0.0; // favorable
	  // }  
      
	  if(R_IsNA(survival_C(index_C,3))==false){       
	    proba[1]=1-survival_C(index_C,3)/survival_T(index_T,1); // unfavorable [St(y_j-taux)/St(x_i)] 
	  } // else {
	  //      proba[1]=0.0; // unfavorable 
	  // }  
      
	}else{  // (1,0) [-tau,tau]
	  if(R_IsNA(survival_C(index_C,5))==false){       
	    proba[0]=survival_C(index_C,5)/survival_T(index_T,1); // favorable [St(y_j+taux)/Sc(x_i)] 
	  }// else {
	  //      proba[0]=0.0; // favorable
	  // }

	  //      proba[1]=0.0; // unfavorable 

	}
      } 
    
    }else{ // delta_C==0
      proba[3] = 1.0; // uninformative

      if(diff>threshold){ // (0,0) > tau

	if(R_IsNA(survival_T(index_T,3))==false){
	  proba[0]=1.0-survival_T(index_T,3)/survival_C(index_C,1) - calcIntegralProba(endpoint_T-threshold, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable 1-[St(x_i-taux)/Sc(y_j)]-I/(St(x_i)*Sc(y_j)) 
  
	}// else {
	//      proba[0]=0.0; // favorable
	// }
  
	proba[1]=-calcIntegralProba(endpoint_T, survival_T, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable  
   
      }else if(diff< -threshold){  // (0,0) < -tau
	proba[0]=-calcIntegralProba(endpoint_C, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable

	if(R_IsNA(survival_C(index_C,3))==false){
	  proba[1]=1.0-survival_C(index_C,3)/survival_T(index_T,1) - calcIntegralProba(endpoint_C-threshold, survival_T, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable        
  
	}// else {
	//      proba[1]=0.0; // unfavorable 
	// }
  
      }else{  // (0,0) [-tau;tau]

	proba[0]=-calcIntegralProba(endpoint_C, survival_C, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // favorable
	proba[1]=-calcIntegralProba(endpoint_T, survival_T, 8)/(survival_T(index_T,1)*survival_C(index_C,1)); // unfavorable  
  

      }
  
    }}
  
  return(proba);  
}

// * calcIntegralProba
inline double calcIntegralProba(const double time_begin, const arma::mat& survival, const int index_col){
  
  // index col must be 7 or 8
  
  // survival [0-2] unordered Sc (control time - tau, control time, control time + tau)
  //          [3-5] unordered St (control time - tau, control time, control time + tau)
  //          [6] ordered Sc (control time)
  //          [7-8] ordered St (control time - tau, contol time + tau)
  //          [9] orderd survival times
  //          [10] orderd event type (1 : event, 0 : censoring)
  
  double integral = 0;  
  int  n_survival=survival.n_rows;
  
  if(survival(0,10)>0.5 && survival(0,9)>time_begin && R_IsNA(survival(0,index_col))==false){
    integral += survival(0,index_col)*(survival(0,6)-1);
  }
      
    
  if(n_survival>0){
    
    for(int iter_time=1 ; iter_time<n_survival ; iter_time++){
      
      if(survival(iter_time,10)>0.5 && survival(iter_time,9)>time_begin && R_IsNA(survival(iter_time,index_col))==false){
        integral += survival(iter_time,index_col)*(survival(iter_time,6)-survival(iter_time-1,6));
      }
      
    } 
    
  }
  
  return(integral);
  
}
