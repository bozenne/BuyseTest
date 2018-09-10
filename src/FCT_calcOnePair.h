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
										   bool neutralAsUninf, bool keepScore);
 
inline arma::rowvec calcOnePair_TTEgehan(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
										 const int index_T, const int index_C, const double Wpair, const int iter_pair, 
										 double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
										 vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
										 vector<int>& index_wUninf, vector<int>&  index_wNeutral,
										 bool neutralAsUninf, bool keepScore);
 
inline vector<double> calcOneProba_TTEperon(const double endpoint_T, const double endpoint_C, const double delta_T, const double delta_C, const double threshold,
											const int index_T, const int index_C,
                    					    const arma::mat& survTimeC, const arma::mat& survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT);

double calcIntegralProba(const arma::mat& survival, double start);

// * calcOnePair_Continuous
inline arma::rowvec calcOnePair_Continuous(const double endpoint_T, const double endpoint_C, const double threshold,
										   const int index_T, const int index_C, const double Wpair, const int iter_pair, 
										   double& count_favorable, double& count_unfavorable, double& count_neutral, double& count_uninf,
										   vector<int>& index_neutralT, vector<int>& index_neutralC, vector<int>& index_uninfT, vector<int>& index_uninfC, 
										   vector<int>& index_wUninf, vector<int>&  index_wNeutral,
										   bool neutralAsUninf, bool keepScore){

  arma::rowvec iRow;
  
  if(R_IsNA(endpoint_T) || R_IsNA(endpoint_C)){
    count_uninf+=Wpair;
    
    index_uninfT.push_back(index_T);
    index_uninfC.push_back(index_C);     
    if(iter_pair>=0){
      index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs 
    }
    if(keepScore){
      iRow = {(double)index_T, (double)index_C, // indexT, indexC
			  0, // favorable
			  0, // unfavorable
			  0, // neutral
			  1, // uninformative
			  Wpair, // weight
			  0, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
      };
    }
  }else{    
    double diff = endpoint_T-endpoint_C; // difference between the endpoints from the treatment and control patients of the pair
    
    if(diff >= threshold){ // diff>0 for threshold == 0
      count_favorable+=Wpair;
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				1, // favorable
				0, // unfavorable
				0, // neutral
				0, // uninformative
				Wpair, // weight
				Wpair, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		};
      }
    }else if(diff <= -threshold){ // diff<0 for threshold == 0
      count_unfavorable+=Wpair;
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				1, // unfavorable
				0, // neutral
				0, // uninformative
				Wpair, // weight
				0, Wpair, 0 // favorable corrected, unfavorable corrected, neutral corrected
		}; 
      }

    }else{
      count_neutral+=Wpair;

	  if(neutralAsUninf){
		index_neutralT.push_back(index_T);
		index_neutralC.push_back(index_C);
		if(iter_pair>=0){
		  index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs      
		}
	  }
	  
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				1, // neutral
				0, // uninformative
				Wpair, // weight
				0, 0, Wpair  // favorable corrected, unfavorable corrected, neutral corrected
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
										 bool neutralAsUninf, bool keepScore){
  
  double diff = endpoint_T-endpoint_C; // difference between the endpoints from the treatment and control patients of the pair
  arma::rowvec iRow;
  
  if(delta_T==1){
    
    if(delta_C==1){
      
      if(diff >= threshold){  // (1,1) >= tau    : favorable
        count_favorable+=Wpair;
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				  1, // favorable
				  0, // unfavorable
				  0, // neutral
				  0, // uninformative
				  Wpair, // weight
				  Wpair, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else if(diff <= -threshold){              // (1,1) <= -tau   : unfavorable
		count_unfavorable+=Wpair;
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  1, // unfavorable
				  0, // neutral 
				  0, // uninformative
				  Wpair, // weight
    	          0, Wpair, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                                  //  (1,1)  ]-tau;tau[ : uninformative
		count_neutral+=Wpair;

		if(neutralAsUninf){
		  index_neutralT.push_back(index_T);
		  index_neutralC.push_back(index_C);
		  if(iter_pair > 0){
			index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs, only stored at the second outcome
		  }
		}
		
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  1, // neutral
				  0, // uninformative
				  Wpair, // weight
				  0, 0, Wpair // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }

      }      
      
    }else if(delta_C==0){
	
      if(diff <= -threshold){              // (1,0) <= -tau   : unfavorable
		count_unfavorable+=Wpair;           
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  1, // unfavorable
				  0, // neutral 
				  0, // uninformative
				  Wpair, // weight
				  0, Wpair, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                             //  (1,0)  ]-tau;+Inf[ : uninformative
		count_uninf+=Wpair;
    
		index_uninfT.push_back(index_T);
		index_uninfC.push_back(index_C);
		if(iter_pair > 0){
		  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs ; only stored for the second outcome
		}
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  0, // neutral 
				  1, // uninformative
				  Wpair, // weight
				  0, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }

      } 
    }else if(delta_C==2){
	  count_unfavorable+=Wpair;           
	  if(keepScore){
		iRow = {(double)index_T, (double)index_C,  // indexT, indexC
				0, // favorable
				1, // unfavorable
				0, // neutral 
				0, // uninformative	
				Wpair, // weight
				0, Wpair, 0 // favorable corrected, unfavorable corrected, neutral corrected	  
		};
	  }      
    }
    
  }else if(delta_T==0){
    
    if(delta_C==1){
    
      if(diff >= threshold){   // (0,1) > tau    : favorable
		count_favorable+=Wpair;      
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  1, // favorable
				  0, // unfavorable
				  0, // neutral 
				  0, // uninformative
				  Wpair, // weight
				  Wpair, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }else{                 //  (1,0)  ]-Inf;+tau[ : uninformative
		count_uninf+=Wpair;
    
		index_uninfT.push_back(index_T);
		index_uninfC.push_back(index_C);
		if(iter_pair > 0){
		  index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs ; only stored for the second outcome
		}
		if(keepScore){
		  iRow = {(double)index_T, (double)index_C, // indexT, indexC
				  0, // favorable
				  0, // unfavorable
				  0, // neutral 
				  1, // uninformative
				  Wpair, // weight
				  0, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		  };
        }
      }
    
    }else{ // delta_C==0 or 2              // (0,0) ]-Inf;+Inf[ : uninformative
      
      count_uninf+=Wpair;
  
      index_uninfT.push_back(index_T);
      index_uninfC.push_back(index_C); 
      if(iter_pair > 0){
		index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome
      }
	  if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				0, // neutral 
				1, // uninformative
				Wpair, // weight
				0, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
		};
	  }
    }
    
  }else if(delta_T==2){

    if(delta_C==1){
      count_favorable+=Wpair;      
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				1, // favorable
				0, // unfavorable
				0, // neutral 
				0, // uninformative
				Wpair, // weight
				Wpair, 0, 0  // favorable corrected, unfavorable corrected, neutral corrected
		};
	  }
    }else if(delta_C==2){

      count_neutral+=Wpair;

	  if(neutralAsUninf){
		index_neutralT.push_back(index_T);
		index_neutralC.push_back(index_C);
		if(iter_pair > 0){
		  index_wNeutral.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome   
		}
	  }
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				1, // neutral 
				0, // uninformative
				Wpair, // weight
				0, 0, Wpair // favorable corrected, unfavorable corrected, neutral corrected
		};
      }

    }else if(delta_C==0){
      count_uninf+=Wpair;
  
      index_uninfT.push_back(index_T);
      index_uninfC.push_back(index_C); 
      if(iter_pair > 0){
		index_wUninf.push_back(iter_pair); // index of the pair relative to Wpairs  ; only stored for the second outcome
      }
      if(keepScore){
		iRow = {(double)index_T, (double)index_C, // indexT, indexC
				0, // favorable
				0, // unfavorable
				0, // neutral 
				1, // uninformative
				Wpair, // weight
				0, 0, 0 // favorable corrected, unfavorable corrected, neutral corrected
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
      
      if(diff >= threshold){ 
		proba[0] = 1.0; // favorable
		// proba[1] = 0.0; // unfavorable  
		// proba[2] = 0.0; // neutral
		// proba[3] = 0.0; // uniformative
      
      }else if(diff <= -threshold){ 
		// proba[0] = 0.0; // favorable
		proba[1] = 1.0; // unfavorable 
		// proba[2] = 0.0; // neutral      
		// proba[3] = 0.0; // uniformative
      }else{ 
		// proba[0] = 0.0; // favorable
		// proba[1] = 0.0; // unfavorable  
		proba[2] = 1.0; // neutral
   		// proba[3] = 0.0; // uniformative
      }      
      
    }else{ // deltaC[iter_C]==0
        
      if(diff >= threshold){ //
		if(R_IsNA(survTimeT(index_T,1))==false && R_IsNA(survTimeT(index_T,3))==false){ // (1,0) > tau
		  // favorable
		  proba[0] = 1.0 - survTimeT(index_T,1)/survTimeC(index_C,2); // 1-[Sc(x_i-taux)/Sc(y_j)]
		  // unfavorable
		  proba[1] = survTimeT(index_T,3)/survTimeC(index_C,2); //  [Sc(x_i+taux)/Sc(y_j)]
		  // neutral
		  proba[2] = 1 - (proba[0] + proba[1]); //  [Sc(x_i+taux)/Sc(y_j)]
		  // uniformative
		  // proba[3] = 0.0;
		}else{
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
		   proba[3] = 1.0; // uniformative
		}
		
      }else if(diff <= -threshold){
		// proba[0] = 0.0; // favorable
		proba[1] = 1.0; // unfavorable  
		// proba[2] = 0.0; // neutral
   		// proba[3] = 0.0; // uniformative
	  }else{ // |diff| < threshold

		if(R_IsNA(survTimeT(index_T,3))==false){
		  // favorable
		  // proba[0]=0.0; 
		  // unfavorable		  
		  proba[1] = survTimeT(index_T,3)/survTimeC(index_C,2); //  [Sc(x_i+taux)/Sc(y_j)]
		  // neutral
		  proba[2] = 1 - (proba[0] + proba[1]); //  1-[Sc(x_i+taux)/Sc(y_j)]
		  // uniformative
		  // proba[3] = 0.0;
		}else {
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
   	      proba[3] = 1.0; // uniformative
		}
      }
	  
    }
    
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 
    
      if(diff >= threshold){ // 
		proba[0] = 1.0; // favorable     
		// proba[1] = 0.0; // unfavorable  
		// proba[2] = 0.0; // neutral
   		// proba[3] = 0.0; // uniformative
      }else if(diff <= -threshold){
		if(R_IsNA(survTimeC(index_C,6))==false && R_IsNA(survTimeC(index_C,4))==false){
		  // favorable
		  proba[0] = survTimeC(index_C,6)/survTimeT(index_T,5); // [St(y_j+taux)/St(x_i)]
		  // unfavorable
		  proba[1] = 1 - survTimeC(index_C,4)/survTimeT(index_T,5); // 1-[St(y_j-taux)/St(x_i)]
		  // neutral
		  proba[2] = 1 - (proba[0] + proba[1]); // [St(y_j-taux)-St(y_j+taux)]/St(x_i)
		  // uniformative
		  // proba[3] = 0.0;
		} else {
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
   	      proba[3] = 1.0; // uniformative
	    }
	  }else{
		if(R_IsNA(survTimeC(index_C,6))==false){
		  // favorable
		  proba[0] = survTimeC(index_C,6)/survTimeT(index_T,5); // [St(y_j+taux)/St(x_i)]   
		  // unfavorable
		  // proba[1]=0.0;
		  // neutral
		  proba[2] = 1 - (proba[0] + proba[1]); // 1 - St(y_j+taux)/St(x_i)
		  // uniformative
		  // proba[3] = 0.0;
		} else {
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
		  proba[3] = 1.0; // uninformative
		 }


	  } 
    
	}else{ // delta_C==0

	  double denom = survTimeT(index_T,5)*survTimeC(index_C,2);
	  double intFavorable; 
	  double intUnfavorable;
	  
	  if(diff >= threshold){
		if(R_IsNA(survTimeT(index_T,1))==false){
		  intFavorable = calcIntegralProba(survJumpC, endpoint_T-threshold) / denom;
		  intUnfavorable = calcIntegralProba(survJumpT, endpoint_T) / denom;
		  // favorable
		  proba[0] = 1.0 - survTimeT(index_T,1)/survTimeC(index_C,2) - intFavorable;
		  // unfavorable  
		  proba[1]= -intUnfavorable;
		  // neutral  
		  proba[2]= 1 - (proba[0] + proba[1]);
		}else{
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
		  proba[3] = 1.0; // uninformative
		}
	  }else if(diff <= -threshold){
		if(R_IsNA(survTimeC(index_C,4))==false){
		  intFavorable = calcIntegralProba(survJumpC, endpoint_C) / denom;
		  intUnfavorable = calcIntegralProba(survJumpT, endpoint_C-threshold) / denom;

		  // favorable
		  proba[0]= -intFavorable;
		  // unfavorable
		  proba[1]= 1.0 - survTimeC(index_C,4)/survTimeT(index_T,5) - intUnfavorable;
		  // neutral
		  proba[2]= 1 - (proba[0] + proba[1]);
		  // uninformative
		  // proba[3] = 0.0;
		}else{
		  // proba[0] = 0.0; // favorable
		  // proba[1] = 0.0; // unfavorable  
		  // proba[2] = 0.0; // neutral
		  proba[3] = 1.0; // uninformative
		}
	  }else{
		intFavorable = calcIntegralProba(survJumpC, endpoint_C) / denom;
		intUnfavorable = calcIntegralProba(survJumpT, endpoint_T) / denom;
		// favorable
		proba[0]= -intFavorable;
		// unfavorable
		proba[1]= -intUnfavorable;
		// neutral
		proba[2]= 1 - (proba[0] + proba[1]);
	    // uninformative
		// proba[3] = 0.0;
	  }
	  // Rcout << "{" << denom << " " << intFavorable << " " << intUnfavorable << "}" << endl;
	  
    }}

  return(proba);  
}

// * calcIntegralProba
//' @title C++ Function Computing the Integral Terms for the Peron Method. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored.
//' @name calcIntegralProba
//' 
//' @param survival [matrix] Contains the jump times in the first column,
//' the survival in the other arm at times plus threshold in the second column,
//' and the jump in survival in the third column.
//' @param start [numeric] Time at which to start the integral.
//'
//' @keywords function Cpp internal
//' @export
// [[Rcpp::export]]
double calcIntegralProba(const arma::mat& survival, double start){
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





























