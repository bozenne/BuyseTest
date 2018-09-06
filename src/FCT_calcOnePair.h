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
 
inline vector<double> calcOneProba_TTEperon(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
											const int index_T, const int index_C,
                    					    const arma::mat& survTimeC, const arma::mat& survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT);

inline double calcIntegralProba(const arma::mat& survival, double start);

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
      iRow = {(double)index_T, (double)index_C, // indexT, indexC
			  0, // favorable
			  0, // unfavorable
			  0, // neutral
			  Wpair, // uninformative
			  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
      };
    }
  }else{    
    double diff = endpoint_T-endpoint_C; // difference between the endpoints from the treatment and control patients of the pair
    
    if(diff >= threshold){ // diff>0 for threshold == 0
      count_favorable+=Wpair;
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				Wpair, // favorable
				0, // unfavorable
				0, // neutral
				0, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		};
      }
    }else if(diff <= -threshold){ // diff<0 for threshold == 0
      count_unfavorable+=Wpair;
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				Wpair, // unfavorable
				0, // neutral
				0, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		}; 
      }

    }else{
      index_neutralT.push_back(index_T);
      index_neutralC.push_back(index_C);
      count_neutral+=Wpair;
      
      if(iter_pair>=0){
        index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs      
      }
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				Wpair, // neutral
				0, // uninformative
				NA_REAL, NA_REAL, NA_REAL  // favorable corrected, unfavorable corrected, neutral corrected
		};
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
      
      if(diff >= threshold){  // (1,1) >= tau    : favorable
        count_favorable+=Wpair;
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				  Wpair, // favorable
				  0, // unfavorable
				  0, // neutral
				  0, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else if(diff <= -threshold){              // (1,1) <= -tau   : unfavorable
		count_unfavorable+=Wpair;
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  Wpair, // unfavorable
				  0, // neutral 
				  0, // uninformative
    	          NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                                  //  (1,1)  ]-tau;tau[ : uninformative
		index_neutralT.push_back(index_T);
		index_neutralC.push_back(index_C);
		count_neutral+=Wpair;
      
		if(iter_pair > 0){
		  index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs, only stored at the second outcome
		}
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  Wpair, // neutral
				  0, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }

      }      
      
    }else if(delta_C==0){
	
      if(diff <= -threshold){              // (1,0) <= -tau   : unfavorable
		count_unfavorable+=Wpair;           
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  Wpair, // unfavorable
				  0, // neutral 
				  0, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                             //  (1,0)  ]-tau;+Inf[ : uninformative
		index_uninfT.push_back(index_T);
		index_uninfC.push_back(index_C);
		count_uninf+=Wpair;
    
		if(iter_pair > 0){
		  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs ; only stored for the second outcome
		}
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  0, // neutral 
				  Wpair, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }

      } 
    }else if(delta_C==2){
	  count_unfavorable+=Wpair;           
	  if(keepComparison){
		iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				0, // favorable
				Wpair, // unfavorable
				0, // neutral 
				0, // uninformative	
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected	  
		};
	  }      
    }
    
  }else if(delta_T==0){
    
    if(delta_C==1){
    
      if(diff >= threshold){   // (0,1) > tau    : favorable
		count_favorable+=Wpair;      
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  Wpair, // favorable
				  0, // unfavorable
				  0, // neutral 
				  0, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                 //  (1,0)  ]-Inf;+tau[ : uninformative
		index_uninfT.push_back(index_T);
		index_uninfC.push_back(index_C);
		count_uninf+=Wpair;
    
		if(iter_pair > 0){
		  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs ; only stored for the second outcome
		}
		if(keepComparison){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  0, // neutral 
				  Wpair, // uninformative
				  NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }
    
    }else{ // delta_C==0 or 2              // (0,0) ]-Inf;+Inf[ : uninformative
      
      index_uninfT.push_back(index_T);
      index_uninfC.push_back(index_C); 
      count_uninf+=Wpair;
  
      if(iter_pair > 0){
		index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome
      }
	  if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				0, // neutral 
				Wpair, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		};
	  }
    }
    
  }else if(delta_T==2){

    if(delta_C==1){
      count_favorable+=Wpair;      
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				Wpair, // favorable
				0, // unfavorable
				0, // neutral 
				0, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		};
	  }
    }else if(delta_C==2){

      index_neutralT.push_back(index_T);
      index_neutralC.push_back(index_C);
      count_neutral+=Wpair;
      
      if(iter_pair > 0){
		index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome   
      }
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				Wpair, // neutral 
				0, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		};
      }

    }else if(delta_C==0){
      index_uninfT.push_back(index_T);
      index_uninfC.push_back(index_C); 
      count_uninf+=Wpair;
  
      if(iter_pair > 0){
		index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome
      }
      if(keepComparison){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				0, // neutral 
				Wpair, // uninformative
				NA_REAL, NA_REAL, NA_REAL // favorable corrected, unfavorable corrected, neutral corrected
		};
       	// (indexT, indexC, favorable, unfavorable, neutral, uninformative, weight for IPWC)
      }
    }
    
  }
  
  /* Rcout << iRow[0] << " " << iRow[1] << " " <<  iRow[2] << " " << iRow[3] << " " << iRow[4] << " " << iRow[5]  << endl; */
  return iRow;
  
}

// * calcOneProba_TTEperon
inline vector<double> calcOneProba_TTEperon(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
											const int index_T, const int index_C,
											const arma::mat& survTimeC, const arma::mat& survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT){
  
  // survTimeC and survTimeT: survival at control/treatment observation times
  //        [0]    time 
  //        [1-3]  survival estimated in the control arm: time - tau, time, time + tau
  //        [4-6]  survival estimated in the treatment arm: time - tau, time, time + tau
  
  // survcCJumpC and survTJumpT: survival at jump times
  //        [0]  time
  //        [1]  survival estimated at time + tau (treatment arm and control arm)
  //        [2]  d(survival) estimated at time (control arm and treatment arm)
    
  
  double diff = endpoint_T-endpoint_C;
  vector<double> proba(4,0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative
  /* Rcout << " (" << delta_T << ";" << delta_C << ")"; */
  if(delta_T==1){
    if(delta_C==1){
      
      if(diff >= threshold){ // (1,1) >= tau
		proba[0] = 1.0; // favorable
		//      proba[1] = 0.0; // unfavorable
      
      }else if(diff <= -threshold){ // (1,1) <= -tau
		//      proba[0] = 0.0; // favorable   
		proba[1] = 1.0; // unfavorable 
      
      }else{  // (1,1) [-tau,tau]
		proba[2] = 1.0; // neutral
      
      }      
      
    }else{ // deltaC[iter_C]==0
        
      if(diff <= -threshold){ // (1,0) <= -tau
		//    proba[0] = 0.0; // favorable
		proba[1] = 1.0; // unfavorable  
    
      }else{

		// favorable
		if(diff >= threshold && R_IsNA(survTimeT(index_T,1))==false){ // (1,0) > tau
		  proba[0]=1.0-survTimeT(index_T,1)/survTimeC(index_C,2); // 1-[Sc(x_i-taux)/Sc(y_j)]
  	    }// else {
		//      proba[0]=0.0; 
		// }

		// unfavorable
		// for both diff > threshold and |diff|<threshold
		if(R_IsNA(survTimeT(index_T,3))==false){
		  proba[1]=survTimeT(index_T,3)/survTimeC(index_C,2); //  [Sc(x_i+taux)/Sc(y_j)]          
		}// else {
		//      proba[1]=0.0; 
		// }

		// neutral
		// proba[2] = 0.0;
	
		// uninformative
		proba[3] = 1.0 - (proba[1] + proba[0]); 
            
      }
    }
    
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 
    
      if(diff >= threshold){ // (0,1) >= taux
		proba[0] = 1.0; // favorable     
		//    proba[1] = 0.0; // unfavorable    
    
      }else{

		// favorable
		// for both diff > threshold and |diff|<threshold
		if(R_IsNA(survTimeC(index_C,6))==false){
		  proba[0]=survTimeC(index_C,6)/survTimeT(index_T,5); // [St(y_j+taux)/St(x_i)]   
		} // else {
		//      proba[0]=0.0;
		// }
	
		if(diff <= -threshold && R_IsNA(survTimeC(index_C,4))==false){  // (0,1) < tau
		  proba[1]=1-survTimeC(index_C,4)/survTimeT(index_T,5); // unfavorable [St(y_j-taux)/St(x_i)] 
		}// else {
		//      proba[1]=0.0; // unfavorable 
		// }
		
		// neutral
		// proba[2] = 0.0;
	
		// uninformative
		proba[3] = 1.0 - (proba[1] + proba[0]); 

      } 
    
    }else{ // delta_C==0

	  double denom = survTimeT(index_T,5)*survTimeC(index_C,2);
	  double intFavorable; 
	  double intUnfavorable;
	  if(diff >= threshold){
		intFavorable = calcIntegralProba(survJumpC, endpoint_T-threshold) / denom;
	  }else{
		intFavorable = calcIntegralProba(survJumpC, endpoint_C) / denom;
	  }
	  if(diff <= -threshold){
		intUnfavorable = calcIntegralProba(survJumpT, endpoint_C-threshold) / denom;
	  }else{
		intUnfavorable = calcIntegralProba(survJumpT, endpoint_T) / denom;
	  }
	  // Rcout << denom << " " << intFavorable << " " << intUnfavorable << endl;
	  
	  // favorable
      if(diff>threshold){ // (0,0) > tau
		// favorable 1-[Sc(x_i-taux)/Sc(y_j)]-I/(St(x_i)*Sc(y_j)) 
	    proba[0] = 1.0 - survTimeT(index_T,1)/survTimeC(index_C,2) - intFavorable;
	  }else{
		proba[0]=  -intFavorable;
	  }

	  // unfavorable        
      if(diff < -threshold){  // (0,0) < -tau
		// unfavorable 1-[St(y_j-taux)/St(x_i)]-I/(St(x_i)*Sc(y_j))		
		proba[1]= 1.0 - survTimeC(index_C,4)/survTimeT(index_T,5) - intUnfavorable; 
	  }else{
		proba[1]= -intUnfavorable;
	  }

	  // neutral
	  // proba[2] = 0.0;
	
	  // uninformative
	  proba[3] = 1.0 - (proba[1] + proba[0]); 

  
    }}
  
  return(proba);  
}

// * calcIntegralProba
inline double calcIntegralProba(const arma::mat& survival, double start){
  // computes \int_t>tau S dS

  // survival contains:
  // [0] times
  // [1] S
  // [2] dS
  
  double integral = 0;  
  int nJump = survival.n_rows;
  
  if(nJump>0){    
    for(int iter_time=0 ; iter_time<nJump ; iter_time++){

	  if(R_IsNA(survival(iter_time,1))){break;}
	  
      if(survival(iter_time,0) >= start){
        integral += survival(iter_time,1)*survival(iter_time,2);
      }

	}     
  }
  
  return(integral);
}





























