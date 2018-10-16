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

inline vector<double> calcOnePair_Continuous(double diff, double threshold);
 
inline vector<double> calcOnePair_TTEgehan(double diff, double delta_C, double delta_T, double threshold);
 
inline vector<double> calcOneScore_TTEperon(double endpoint_C, double endpoint_T, double delta_C, double delta_T, double threshold,
											arma::rowvec survTimeC, arma::rowvec survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT,
											double lastSurvC, double lastSurvT);

double calcIntegralScore_cpp(const arma::mat& survival, double start);

// * calcOnePair_Continuous
inline vector<double> calcOnePair_Continuous(double diff, double threshold){

  // ** initialize
  std::vector<double> score(4,0.0);

  // ** score
  if(R_IsNA(diff)){ // missing data: uninformative
    score[3] = 1.0;
  }else{    
    
    if(diff >= threshold){ // favorable
      score[0] = 1.0;
    }else if(diff <= -threshold){ // unfavorable
      score[1] = 1.0;
    }else{ // neutral
      score[2] = 1.0;
    }
    
  }

  // ** export
  // Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << endl;
  return(score);
  
}

// * calcOnePair_TTEgehan
inline vector<double> calcOnePair_TTEgehan(double diff, double delta_C, double delta_T, double threshold){
  
  // ** initialize
  std::vector<double> score(4,0.0);
  // Rcout << diff << " " << delta_T << " " << delta_C << " " << threshold << endl;

  // ** score
  if(delta_T==1){
    
    if(delta_C==1){ // (treatment event, control event)
      
      if(diff >= threshold){         // >= tau    : favorable
        score[0] = 1.0;
      }else if(diff <= -threshold){ // <= -tau    : unfavorable
		score[1] = 1.0;
      }else{                        // ]-tau;tau[ : neutral
		score[2] = 1.0;
      }
      
    }else if(delta_C==0){ // (treatment event, control censored)
	
      if(diff <= -threshold){ // <= -tau   : unfavorable
		score[1] = 1.0;
      }else{                  // otherwise : uninformative
		score[3] = 1.0;
      }
      
    }else if(delta_C==2){ // (treatment event, control competing risk)
      score[1] = 1.0; //  unfavorable
    }
    
  }else if(delta_T==0){
    
    if(delta_C==1){ // (treatment censored, control event)
    
      if(diff >= threshold){ // > tau    : favorable
		score[0] = 1.0;
      }else{                 // otherwise: uninformative
		score[3] = 1.0;
      }
    
    }else{ // (treatment censored, control censored/competing risk): uninformative
	  score[3] = 1.0;
    }
    
  }else if(delta_T==2){ 

    if(delta_C==1){ // (treatment competing risk, control event): favorable
	  score[0] = 1.0;
    }else if(delta_C==2){ // (treatment competing risk, control competing risk): neutral
	  score[2] = 1.0;
    }else if(delta_C==0){ // (treatment competing risk, control censored): uninformative
	  score[3] = 1.0;
    }
    
  }

  // ** export
  // Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << endl;
  return(score);
  
}

// * calcOneScore_TTEperon
inline vector<double> calcOneScore_TTEperon(double endpoint_C, double endpoint_T, double delta_C, double delta_T, double threshold,
											arma::rowvec survTimeC, arma::rowvec survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT,
											double lastSurvC, double lastSurvT){
  
  // survTimeC and survTimeT: survival at control/treatment observation times
  //        [0]    time 
  //        [1-3]  survival estimated in the control arm: time - tau, time, time + tau
  //        [4-6]  survival estimated in the treatment arm: time - tau, time, time + tau
  
  // survcCJumpC and survTJumpT: survival at jump times
  //        [0]  time
  //        [1]  survival estimated at time + tau (treatment arm and control arm)
  //        [2]  d(survival) estimated at time (control arm and treatment arm)
    
  
  // ** initialize
  double diff = endpoint_T-endpoint_C;
  vector<double> score(4,0.0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative
  double upperFavorable;
  double upperUnfavorable;
  
  /* Rcout << " (" << delta_T << ";" << delta_C << ")"; */
  // ** compute favorable and unfavorable
  if(delta_T==1){
    if(delta_C==1){
      
      if(diff >= threshold){ 
		score[0] = 1.0; // favorable
		// score[1] = 0.0; // unfavorable  
		// score[2] = 0.0; // neutral
		// score[3] = 0.0; // uniformative
      
      }else if(diff <= -threshold){ 
		// score[0] = 0.0; // favorable
		score[1] = 1.0; // unfavorable 
		// score[2] = 0.0; // neutral      
		// score[3] = 0.0; // uniformative
      }else{ 
		// score[0] = 0.0; // favorable
		// score[1] = 0.0; // unfavorable  
		score[2] = 1.0; // neutral
		// score[3] = 0.0; // uniformative
      }      

      upperFavorable = score[0];
      upperUnfavorable = score[1];
	
    }else{ // deltaC[iter_C]==0

      // favorable
      if(diff >= threshold){
		if(survTimeC(2)==0){ // since 0 = lastSurvC >= survTimeT(1)/lastSurv >= 0 we know that survTimeT(1)/lastSurv = 0
		  score[0] = 0.0;
		  upperFavorable = score[0];
		}else if(R_IsNA(survTimeT(1))==false){
		  score[0] = 1.0 - survTimeT(1)/survTimeC(2); // 1-[Sc(x_i-tau)/Sc(y_j)]
		  upperFavorable = score[0];
		}else{
		  score[0] = 1.0 - lastSurvC/survTimeC(2); // 1-[Sc(max)/Sc(y_j)] (lower bound)
		  upperFavorable = 1.0;  // (upper bound)
		}	
      }else {
		// score[0] = 0.0;
		upperFavorable = score[0];
      }

      // unfavorable
      if(diff <= -threshold || survTimeC(2)==0){ //
		score[1] = 1.0;
		upperUnfavorable = score[1];
      }else if(R_IsNA(survTimeT(3))==false){
		score[1] = survTimeT(3)/survTimeC(2); //  [Sc(x_i+tau)/Sc(y_j)]
		upperUnfavorable = score[1];
      }else {
		// score[1] = 0.0 // (lower bound)
		upperUnfavorable = lastSurvC/survTimeC(2); // (upper bound)
      }

    }
  }else{ // deltaT[iter_T]==0
    if(delta_C==1){ 

      // favorable
      if(diff >= threshold || survTimeT(5) == 0){ // 
		score[0] = 1.0;
		upperFavorable = score[0];
      }else if(R_IsNA(survTimeC(6))==false){
		score[0] = survTimeC(6)/survTimeT(5); // [St(y_j+tau)/St(x_i)]
		upperFavorable = score[0];
      }else{
		// score[0] = 0.0 // lower bound
		upperFavorable = lastSurvT/survTimeT(5); // upper bound
      }

      // unfavorable
      if(diff <= -threshold){
		if(survTimeT(5) == 0){
		  score[1] = 0.0;
		  upperUnfavorable = score[1];
		}else if(R_IsNA(survTimeC(4))==false){
		  score[1] = 1.0 - survTimeC(4)/survTimeT(5); // 1-[St(y_j-tau)/St(x_i)]
		  upperUnfavorable = score[1];
		}else{
		  score[1] = 1.0 - lastSurvT/survTimeT(5); // 1-[St(max)/St(x_i)] (lower bound)
		  upperUnfavorable = 1.0; // (upper bound)
		}
      }else{
		// score[1] = 0.0;
		upperUnfavorable = score[1];
      }
      
    }else{ // delta_C==0

      double denom = survTimeT(5)*survTimeC(2);
      double intFavorable; 
      double intUnfavorable;

      // favorable
      if(diff >= threshold){
		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_T-threshold) / denom;  // -intFavorable is already the lower bound

		if(R_IsNA(survTimeT(1))==false){
		  score[0] = 1.0 - survTimeT(1)/survTimeC(2) - intFavorable; // (lower bound)
		  upperFavorable = score[0] + lastSurvC*lastSurvT/denom; // (upper bound)
		}else{
		  score[0] = 1.0 - lastSurvC/survTimeC(2) - intFavorable; // (lower bound)
		  upperFavorable = 1.0 - intFavorable + lastSurvC*lastSurvT/denom; // (upper bound)
		}

      }else{
		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_C) / denom; // -intFavorable is already the lower bound

		score[0] = -intFavorable; // (lower bound)
		upperFavorable = score[0] + lastSurvC*lastSurvT/denom; // (upper bound)
      }
      
      // unfavorable
      if(diff <= -threshold){	
		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_C-threshold) / denom; // -intUnfavorable is already the lower bound
	
		if(R_IsNA(survTimeC(4))==false){
		  score[1] = 1.0 - survTimeC(4)/survTimeT(5) - intUnfavorable; // (lower bound)
		  upperUnfavorable = score[1] + lastSurvC*lastSurvT/denom;  // (upper bound)
		}else{
		  score[1] = 1.0 - lastSurvT/survTimeT(5) - intUnfavorable; // (lower bound)
		  upperUnfavorable = 1.0 - intUnfavorable  + lastSurvC*lastSurvT/denom;  // (upper bound)
		}

      }else{
		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_T) / denom; // -intUnfavorable is already the lower bound
		score[1]= -intUnfavorable; // (lower bound)
		upperUnfavorable = score[1] + lastSurvC*lastSurvT/denom;  // (upper bound)
      }
      
    }}

  // ** compute neutral and uninformative
  // neutral
  score[2] = std::max(1 - upperFavorable - upperUnfavorable, 0.0);

  // uninformative
  score[3] = std::max(1 - (score[0] + score[1] + score[2]), 0.0);

  // ** export
  // Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << " (upper) " << upperFavorable << " " << upperUnfavorable << endl;
  return(score);  
}

// * calcIntegralScore_cpp
//' @title C++ Function Computing the Integral Terms for the Peron Method. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored.
//' @name calcIntegralScore_cpp
//' 
//' @param survival [matrix] Contains the jump times in the first column,
//' the survival in the other arm at times plus threshold in the second column,
//' and the jump in survival in the third column.
//' @param start [numeric] Time at which to start the integral.
//'
//' @keywords function Cpp internal
//' @export
// [[Rcpp::export]]
double calcIntegralScore_cpp(const arma::mat& survival, double start){
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
	  
      if(survival(iter_time,0) > start){ // strict
        integral += survival(iter_time,1)*survival(iter_time,2);
      }

    }     
  }

  return(integral);
}
