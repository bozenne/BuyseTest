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

inline std::vector< double > calcOnePair_Continuous(double diff, double threshold);
 
inline std::vector< double > calcOnePair_TTEgehan(double diff, double delta_C, double delta_T, double threshold);
 
inline std::vector< double > calcOneScore_TTEperon(double endpoint_C, double endpoint_T, double delta_C, double delta_T, double threshold,
											arma::rowvec survTimeC, arma::rowvec survTimeT,
											const arma::mat& survJumpC, const arma::mat& survJumpT,
											double lastSurvC, double lastSurvT);

inline std::vector< double > CalcOnePair_TTEperon_CR(double endpoint_T, double endpoint_C, double delta_T, double delta_C, 
                                                  double tau, int index_T, int index_C, const arma::mat& cifTimeT, 
                                                  const arma::mat& cifTimeC, const arma::mat& cifJumpC, 
                                                  double lastCif1C, double lastCif2C, double lastCif1T, 
                                                  double lastCif2T)

std::vector<double> calcIntegralScore_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv);

double CalcIntegral_Peron_CR(const arma::mat& cif, double start_val, double stop_val, double CIF_t, double lastCIF, int type)

// * calcOnePair_Continuous
inline std::vector< double > calcOnePair_Continuous(double diff, double threshold){

  // ** initialize
  std::vector< double > score(4,0.0);

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
inline std::vector< double > calcOnePair_TTEgehan(double diff, double delta_C, double delta_T, double threshold){
  
  // ** initialize
  std::vector< double > score(4,0.0);
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
inline std::vector< double > calcOneScore_TTEperon(double endpoint_C, double endpoint_T, double delta_C, double delta_T, double threshold,
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
  std::vector< double > score(4,0.0); // [0] favorable, [1] unfavorable, [2] test neutral [3] test uniformative
  double upperFavorable;
  double upperUnfavorable;

  // ** deal with null survival
  // according to the survival the observation will die immediatly after the observation time.
  // so we treat it as if was an event
  if(delta_C==0 && (survTimeC(2) == 0) ){
	delta_C = 1;
  }
  if(delta_T==0 && (survTimeT(5) == 0) ){
	delta_T = 1;
  }
  	
  /* Rcout << " (" << delta_T << ";" << delta_C << ")"; */
  // ** compute favorable and unfavorable
  if(delta_T==1){
    if(delta_C==1){
      
      if(diff >= threshold){ 
		score[0] = 1.0; // favorable
		// score[1] = 0.0; // unfavorable        
      }else if(diff <= -threshold){ 
		// score[0] = 0.0; // favorable
		score[1] = 1.0; // unfavorable 
      }else{ 
		// score[0] = 0.0; // favorable
		// score[1] = 0.0; // unfavorable  
      }      

      upperFavorable = score[0];
      upperUnfavorable = score[1];
	
    }else{ // deltaC[iter_C]==0

      // favorable
      if(diff >= threshold){
		if(R_IsNA(survTimeT(1))==false){
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
      if(diff <= -threshold){ 
		score[1] = 1.0;
		upperUnfavorable = score[1];
	  }else {
		if(R_IsNA(survTimeT(3))==false){
		  score[1] = survTimeT(3)/survTimeC(2); //  [Sc(x_i+tau)/Sc(y_j)]
		  upperUnfavorable = score[1];
		}else {
		  // score[1] = 0.0 // (lower bound)
		  upperUnfavorable = lastSurvC/survTimeC(2); // (upper bound)
		}
	  }

    }
	
  }else{ // deltaT[iter_T]==0

	if(delta_C==1){ 

      // favorable
      if(diff >= threshold){ // 
		score[0] = 1.0;
		upperFavorable = score[0];
      }else {
		if(R_IsNA(survTimeC(6))==false){
		  score[0] = survTimeC(6)/survTimeT(5); // [St(y_j+tau)/St(x_i)]
		  upperFavorable = score[0];
		}else{
		  // score[0] = 0.0 // lower bound
		  upperFavorable = lastSurvT/survTimeT(5); // upper bound
		}
	  }

      // unfavorable
      if(diff <= -threshold){
		if(R_IsNA(survTimeC(4))==false){
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
	  std::vector< double > intFavorable; 
      std::vector< double > intUnfavorable;

      // favorable
      if(diff >= threshold){
		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_T-threshold, lastSurvT, lastSurvC);  // -intFavorable is already the lower bound

		if(R_IsNA(survTimeT(1))==false){
		  score[0] = 1.0 - survTimeT(1)/survTimeC(2) - intFavorable[0] / denom; // (lower bound)
		  upperFavorable = 1.0 - survTimeT(1)/survTimeC(2) - intFavorable[1] / denom; // (upper bound)
		}else{
		  score[0] = 1.0 - lastSurvC/survTimeC(2) - intFavorable[0] / denom; // (lower bound)
		  upperFavorable = 1.0 - intFavorable[1] / denom; // (upper bound)
		}

      }else{
		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_C, lastSurvT, lastSurvC); // -intFavorable is already the lower bound
		score[0] = -intFavorable[0] / denom; // (lower bound)
		upperFavorable = -intFavorable[1] / denom; // (upper bound)
		// Rcout << intFavorable[0] << " " << intFavorable[1] << " " << lastSurvC*lastSurvT << "/" << denom << endl;
      }
      
      // unfavorable
      if(diff <= -threshold){	
		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_C-threshold, lastSurvC, lastSurvT); // -intUnfavorable is already the lower bound
	
		if(R_IsNA(survTimeC(4))==false){
		  score[1] = 1.0 - survTimeC(4)/survTimeT(5) - intUnfavorable[0] / denom; // (lower bound)
		  upperUnfavorable = 1.0 - survTimeC(4)/survTimeT(5) - intUnfavorable[1] / denom; // (upper bound)
		  // Rcout << intUnfavorable[0] << " " << intUnfavorable[1] << " " << survTimeC(4)/survTimeT(5) << "/" << denom << endl;
		}else{
		  score[1] = 1.0 - lastSurvT/survTimeT(5) - intUnfavorable[0] / denom; // (lower bound)
		  upperUnfavorable = 1.0 - intUnfavorable[1] / denom; // (upper bound)
		}
		

      }else{
		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_T, lastSurvC, lastSurvT); // -intUnfavorable is already the lower bound
		score[1]= -intUnfavorable[0] / denom; // (lower bound)
		upperUnfavorable = -intUnfavorable[1] / denom;  // (upper bound)
		// Rcout << intUnfavorable[0] << " " << intUnfavorable[1] << " " << lastSurvC*lastSurvT << "/" << denom << endl;
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

// * calcOnePair_Peron
// [[Rcpp::export]]
inline std::vector< double > CalcOnePair_Peron_CR(double endpoint_T, double endpoint_C, double delta_T, double delta_C, 
                                               double tau, int index_T, int index_C, const arma::mat& cifTimeT, 
                                               const arma::mat& cifTimeC, const arma::mat& cifJumpC, 
                                               double lastCif1C, double lastCif2C, double lastCif1T, 
                                               double lastCif2T) {
  
  // cifTimeC and cifTimeT: cumulative incidence at control/treatment observation times
  //        [1]    times 
  //        [2-4]  cif of event of interest estimated in the control arm: time - tau, time, time + tau
  //        [5-7]  cif of event of interest estimated in the treatment arm: time - tau, time, time + tau
  //        [8]  cif of competing event estimated in the control arm: time
  //        [9]  cif of competing event estimated in the treatment arm: time
  
  // cifJumpC: cumulative incidence of event of interest in control group at jump times
  //        [1]  jump times in control group (unique values in ascending order)
  //        [2-3]  cif of the treatment group at times-tau and times+tau
  //        [4]  d(cif) of control group estimated at time 
  
  double diff = endpoint_T - endpoint_C;
  double Cif1T_t = cifTimeT(index_T,5);
  std::vector< double > proba(5, 0.0); // [0] favorable, [1] unfavorable, [2] neutral competing, [3] neutral event, 
  // [4] uninformative
  double denomC = 1 - cifTimeC(index_C, 2) - cifTimeC(index_C, 7);
  double denomT = 1 - cifTimeT(index_T, 5) - cifTimeT(index_T, 8);
  
  if(delta_T == 2) {
    if(delta_C == 2) { // (2,2)
      proba[2] = 1.0; // systematically neutral competing
    } else if(delta_C == 1){ // (2,1)
      proba[0] = 1.0; // systematically favorable
    } else if(delta_C == 0) { // (2,0)
      proba[0] = (lastCif1C - cifTimeC(index_C, 2))/denomC;
      // proba[1] = 0
      proba[2] = (lastCif2C - cifTimeC(index_C, 7))/denomC;
      // proba[3] = 0 // since the treated patient had the competing event
      proba[4] = 1 - (proba[0] + proba[2]);
    }
  } else if(delta_T == 1){
    if(delta_C == 2){ // (1,2)
      proba[1] = 1.0; // systematically defavorable
    } else if(delta_C == 1){ // (1,1)
      if(diff >= tau) {
        proba[0] = 1.0;
      } else if(diff <= -tau) {
        proba[1] = 1.0;
      } else { // |diff| < tau
        proba[3] = 1.0;
      }
    } else if(delta_C == 0) { // (1,0)
      if(diff >= tau) {
        if(R_IsNA(cifTimeT(index_T,1)) == false) {
          proba[0] = (cifTimeT(index_T,1) - cifTimeC(index_C,2))/denomC;
        } else {
          proba[0] = (lastCif1C - cifTimeC(index_C,2))/denomC;
        }
        if(R_IsNA(cifTimeT(index_T,4)) == false) {
          proba[1] = (lastCif1C - cifTimeT(index_T,3) + lastCif2C - cifTimeC(index_C,7))/denomC;
        } else {
          proba[1] = (lastCif2C - cifTimeC(index_C,7))/denomC;
        }
        if((R_IsNA(cifTimeT(index_T,3)) == false) & (R_IsNA(cifTimeT(index_T,1)) == false)) {
          proba[3] = (cifTimeT(index_T,3) - cifTimeT(index_T,1))/denomC;
        } else if ((R_IsNA(cifTimeT(index_T,4)) == true) & (R_IsNA(cifTimeT(index_T,2)) == false)) {
          proba[3] = (lastCif1C - cifTimeT(index_T,1))/denomC;
        } else {
          proba[3] = 0.0;
        }
      } else if(diff <= -tau) {
        proba[1] = 1.0;
      } else { // |diff| < tau
        if(R_IsNA(cifTimeT(index_T,3)) == false) {
          proba[1] = (lastCif1C - cifTimeT(index_T,3) + lastCif2C - cifTimeC(index_C,7))/denomC;
          proba[3] = (cifTimeT(index_T,3) - cifTimeC(index_C,2))/denomC;
        } else {
          proba[1] = (lastCif2C - cifTimeC(index_C,7))/denomC;
          proba[3] = (lastCif1C - cifTimeC(index_C,2))/denomC;
        }
      }
      proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    }
  } else { // delta_T == 0
    if(delta_C == 2) { // (0,2)
      proba[1] = (lastCif1T - cifTimeT(index_T,5))/denomT;
      proba[2] = (lastCif2T - cifTimeT(index_T,8))/denomT;
      proba[4] = 1 - (proba[1] + proba[2]);
    } else if(delta_C == 1) { // (0,1)
      if(diff >= tau) {
        proba[0] = 1.0;
      } else if(diff <= -tau) {
        if(R_IsNA(cifTimeC(index_C,6)) == false) {
          proba[0] = (lastCif1T - cifTimeC(index_C,6) + lastCif2T - cifTimeT(index_T,8))/denomT;
        } else {
          proba[0] = (lastCif2T - cifTimeT(index_T,8))/denomT;
        }
        if(R_IsNA(cifTimeC(index_C,4)) == false) {
          proba[1] = (cifTimeC(index_C,4) - cifTimeT(index_T,5))/denomT;
        } else {
          proba[1] = (lastCif1T - cifTimeT(index_T,5))/denomT;
        }
        if((R_IsNA(cifTimeC(index_C,6)) == false) & (R_IsNA(cifTimeC(index_C,4)) == false)) {
          proba[3] = (cifTimeC(index_C,6) - cifTimeC(index_C,4))/denomT;
        } else if ((R_IsNA(cifTimeC(index_C,6)) == true) & (R_IsNA(cifTimeC(index_C,4)) == false)) {
          proba[3] = (lastCif1T - cifTimeC(index_C,4))/denomT;
        } else {
          proba[3] = 0.0;
        }
      } else { // |diff| < tau
        if(R_IsNA(cifTimeC(index_C,6)) == false) {
          proba[0] = (lastCif1T - cifTimeC(index_C,6) + lastCif2T - cifTimeT(index_T,8))/denomT;
          proba[3] = (cifTimeC(index_C,6) - cifTimeT(index_T,5))/denomT;
        } else {
          proba[0] = (lastCif2T - cifTimeT(index_T,8))/denomT;
          proba[3] = (lastCif1T - cifTimeT(index_T,5))/denomT;
        }
      }
      proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    } else if (delta_C == 0) { // (0,0)
      double prob21 = (lastCif2T - cifTimeT(index_T,8))*(lastCif1C - cifTimeC(index_C,2))/(denomT*denomC);
      double prob12 = (lastCif2C - cifTimeC(index_C,7))*(lastCif1T - cifTimeT(index_T,5))/(denomT*denomC);
      double prob22 = (lastCif2C - cifTimeC(index_C,7))*(lastCif2T - cifTimeT(index_T,8))/(denomT*denomC);
      if(diff >= tau) {
        // lower bound of each integral
        double intFav = CalcIntegral_Peron_CR(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = CalcIntegral_Peron_CR(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent1 = CalcIntegral_Peron_CR(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T_t, lastCif1T, 3)/(denomT*denomC);
        double intNeutralEvent2 = CalcIntegral_Peron_CR(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        if (R_IsNA(cifTimeT(index_T,1)) == false) {
          proba[0] = ((cifTimeT(index_T,1) - cifTimeC(index_C,2))/denomC)*((lastCif1T - cifTimeT(index_T,5))/denomT) + 
            intFav + prob21;          
        } else {
          proba[0] = ((lastCif1C - cifTimeC(index_C,2))/denomC)*((lastCif1T - cifTimeT(index_T,5))/denomT) + 
            intFav + prob21;
        }
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent1 + intNeutralEvent2;
      } else if(diff <= -tau) {
        double intFav = CalcIntegral_Peron_CR(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = CalcIntegral_Peron_CR(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent = CalcIntegral_Peron_CR(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        proba[0] = intFav + prob21;
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent;
      } else { // |diff| < tau
        double intFav = CalcIntegral_Peron_CR(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = CalcIntegral_Peron_CR(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent1 = CalcIntegral_Peron_CR(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 3)/(denomT*denomC);
        double intNeutralEvent2 = CalcIntegral_Peron_CR(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        proba[0] = intFav + prob21;
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent1 + intNeutralEvent2;
      }
      proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    }
  }
  return proba;
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
//' @param lastSurv [numeric,>0] last survival value for the survival function in the second column.
//' @param lastdSurv [numeric,>0] last survival value for the survival function in the third column.
//'
//' @keywords function Cpp internal
//' @export
// [[Rcpp::export]]
std::vector< double > calcIntegralScore_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv){
  // computes \int_t>tau S dS

  // survival contains:
  // [0] times
  // [1] S
  // [2] dS
  
  std::vector< double > integral(2,0.0); // lower and upper bound
  bool stopdSurv = true;
  int nJump = survival.n_rows;
  
  if(nJump>0){    
    for(int iter_time=0 ; iter_time<nJump ; iter_time++){

      if(R_IsNA(survival(iter_time,1))){
		integral[1] = integral[0] + lastSurv*survival(iter_time,2); // upper bound
		stopdSurv = false;
		break;
	  }
	  
      if(survival(iter_time,0) > start){ // strict
        integral[0] += survival(iter_time,1)*survival(iter_time,2); // increment lower bound
      }

    }

	if(stopdSurv){
	  integral[1] = integral[0] - survival(nJump-1,1)*lastdSurv; // upper bound (minus because dSurv = (0 - surv(tmax))
	}
  }

  return(integral);
}


// * CalcIntegral_Peron_CR
//' @title C++ Function Computing the Integral Terms for the Peron Method in the presence of competing risks (CR). 
//' @description Compute the integral with respect to the jump in CIF for pairs where both outcomes are censored.
//' @name CalcIntegral_Peron_CR
//' 
//' @param cif [matrix] cif[1] = jump times in control group (event of interest), cif[2-3] = CIF of event of interest in group 
//' T at times - tau and times + tau, cif[4] : jump in cif of control group at times (event of interest).
//' @param start_val [numeric] Time at which to start the integral.
//' @param stop_val [numeric] Time at which to stop the integral.
//' @param CIF_t [numeric] CIF of event of interest in group T evaluated at observed time of treatment patient.
//' @param lastCIF [numeric, >0] last value of CIF of event type 1 in group T.
//' @param type [numeric] Indicates the type of integral to compute (1 for wins, 2 for losses, 3 for neutral pairs with two 
//' events of interest - integral with t+tau and xi - and 4 for neutral pairs with two events of interest - integral with 
//' t+tau and t-tau).
//'
//' @keywords function Cpp internal
//' @export
// [[Rcpp::export]]
double CalcIntegral_Peron_CR(const arma::mat& cif, double start_val, double stop_val, double CIF_t, 
                    double lastCIF, int type){
  
  double integral = 0.0;
  int nJump = cif.n_rows;
  
  if (nJump > 0) {
    if(type == 1) {
      for(int i = 0; i<nJump; i++){
        if(R_IsNA(cif(i,2))) {break;}
        if(cif(i,0) > start_val) {
          integral = integral + (lastCIF - cif(i, 2))*cif(i, 3);
        }
      }
    } else if(type == 2) {
      for(int i = 0; i<nJump; i++){
        double lb = cif(i,1);
        if(R_IsNA(cif(i,1))) {lb = lastCIF;}
        if(cif(i,0) > start_val) {
          integral = integral + (lb - CIF_t)*cif(i, 3);
        }
      }
    } else if(type == 3) {
      for(int i = 0; i<nJump; i++){
        double lb = cif(i,2);
        if(R_IsNA(cif(i,2))) {lb = lastCIF;}
        if((cif(i,0) > start_val) & (cif(i,0) <= stop_val)) {
          integral = integral + (lb - CIF_t)*cif(i, 3);
        }
      }
    } else if(type == 4) {
      for(int i = 0; i<nJump; i++){
        double lb = cif(i, 2);
        if(R_IsNA(cif(i,2))) {lb = lastCIF;}
        if(R_IsNA(cif(i,1))) {break;}
        if(cif(i,0) > start_val) {
          integral = integral + (lb - cif(i, 1))*cif(i, 3);
        }
      }
    }
  }
  
  return integral; 
  
}
