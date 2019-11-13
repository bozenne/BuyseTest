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
												   double lastSurvC, double lastSurvT,
												   arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T, int returnIID);

std::vector<double> calcIntegralScore_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
										  bool returnDeriv, int column, arma::mat& derivSurv, arma::mat& derivSurvD);

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
												   double lastSurvC, double lastSurvT,
												   arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T, int returnIID){
  
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

  if(returnIID > 1){
  	Dscore_Dnuisance_C.fill(0.0);
  	Dscore_Dnuisance_T.fill(0.0);
  }
  
  // ** deal with null survival
  // according to the survival the observation will die immediatly after the observation time.
  // so we treat it as if was an event
  if(delta_C==0 && (survTimeC(2) == 0) ){
    delta_C = 1;
  }
  if(delta_T==0 && (survTimeT(5) == 0) ){
    delta_T = 1;
  }
  	
  // Rcout << " (" << delta_T << ";" << delta_C << ")";
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
		// Rcout << "(2+a) ";
		if(R_IsNA(survTimeT(1))==false){
		  score[0] = 1.0 - survTimeT(1)/survTimeC(2); // 1-[Sc(x_i-tau)/Sc(y_j)]
		  upperFavorable = score[0];
		  if(returnIID>1){ // if((returnIID>1) && (survTimeT(1)>0)){
			Dscore_Dnuisance_C(survTimeT(7),0) -= 1/survTimeC(2); // derivative regarding Sc(x_i-tau)
			Dscore_Dnuisance_C(survTimeC(8),0) += survTimeT(1)/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  }
		}else{
		  score[0] = 1.0 - lastSurvC/survTimeC(2); // 1-[Sc(max)/Sc(y_j)] (lower bound)
		  upperFavorable = 1.0;  // (upper bound)
		  if(returnIID>1){ //if((returnIID>1) && (lastSurvC>0)){
			Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,0) -= 1/survTimeC(2); // derivative regarding Sc(max)
			Dscore_Dnuisance_C(survTimeC(8),0) += lastSurvC/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  }
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
		// Rcout << "(2-b) " << "";
		if((R_IsNA(survTimeT(3))==false) & (survTimeT(3) > 0)){
		  score[1] = survTimeT(3)/survTimeC(2); //  [Sc(x_i+tau)/Sc(y_j)]
		  upperUnfavorable = score[1];
		  if(returnIID>1){ // if((returnIID>1) && (survTimeT(3)>0)){
			Dscore_Dnuisance_C(survTimeT(9),1) += 1/survTimeC(2); // derivative regarding Sc(x_i+tau)
			Dscore_Dnuisance_C(survTimeC(8),1) -= survTimeT(3)/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  }
		}else {
		  // score[1] = 0.0 // (lower bound)
		  upperUnfavorable = lastSurvC/survTimeC(2); // (upper bound)
		  // if(returnIID > 1){
		  // Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,2) += 1/survTimeC(2); // derivative regarding Sc(x_i+tau)
		  // Dscore_Dnuisance_C(survTimeC(8),2) -= lastSurvC/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  // }
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
		// Rcout << "(3+b) ";
		if((R_IsNA(survTimeC(6))==false) && (survTimeC(6) > 0)){
		  score[0] = survTimeC(6)/survTimeT(5); // [St(y_j+tau)/St(x_i)]
		  upperFavorable = score[0];
		  if(returnIID>1){ // if((returnIID>1) && (survTimeC(6)>0)){
			Dscore_Dnuisance_T(survTimeC(12),0) += 1/survTimeT(5); // derivative regarding St(y_j+tau)
			Dscore_Dnuisance_T(survTimeT(11),0) -= survTimeC(6)/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  }
		}else{
		  // score[0] = 0.0 // lower bound
		  upperFavorable = lastSurvT/survTimeT(5); // upper bound
		  // if(returnIID > 1){
		  // Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,2) += 1/survTimeT(5); // derivative regarding St(y_j+tau)
		  // Dscore_Dnuisance_T(survTimeT(11),2) -= lastSurvT/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  // }
		}
      }

      // unfavorable
      if(diff <= -threshold){
		// Rcout << "(3-a) ";
		if(R_IsNA(survTimeC(4))==false){
		  score[1] = 1.0 - survTimeC(4)/survTimeT(5); // 1-[St(y_j-tau)/St(x_i)]
		  upperUnfavorable = score[1];
		  if(returnIID>1){ // if((returnIID>1) && (survTimeC(4)>0)){
			Dscore_Dnuisance_T(survTimeC(10),1) -= 1/survTimeT(5); // derivative regarding St(y_j-tau)
			Dscore_Dnuisance_T(survTimeT(11),1) += survTimeC(4)/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  }
		}else{
		  score[1] = 1.0 - lastSurvT/survTimeT(5); // 1-[St(max)/St(x_i)] (lower bound)
		  upperUnfavorable = 1.0; // (upper bound)
		  if(returnIID>1){ // if((returnIID>1) && (lastSurvT>0)){
			Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1, 1) -= 1/survTimeT(5); // derivative regarding St(max)
			Dscore_Dnuisance_T(survTimeT(11),1) += lastSurvT/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  }
		}
      }else{
		// score[1] = 0.0;
		upperUnfavorable = score[1];
      }
      
    }else{ // delta_C==0

      double denom = survTimeT(5)*survTimeC(2);
      std::vector< double > intFavorable; 
      std::vector< double > intUnfavorable;
   	  arma::mat intDscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows,4); // initialized in calcIntegralScore_cpp
	  arma::mat intDscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows,4); // initialized in calcIntegralScore_cpp

	  // favorable
	  if(diff >= threshold){
		// Rcout << "(4+a) ";
		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_T-threshold, lastSurvT, lastSurvC,
											 (returnIID > 1), 0, intDscore_Dnuisance_T, intDscore_Dnuisance_C);  // -intFavorable is already the lower bound
		if(R_IsNA(survTimeT(1))==false){		  
		  score[0] = 1.0 - survTimeT(1)/survTimeC(2) - intFavorable[0] / denom; // (lower bound)
		  upperFavorable = 1.0 - survTimeT(1)/survTimeC(2) - intFavorable[1] / denom; // (upper bound)
		  if(returnIID>1){ // if((returnIID>1) && (survTimeT(1)>0)){
			Dscore_Dnuisance_C(survTimeT(7),0) -= 1/survTimeC(2); // derivative regarding Sc(x_i-tau)
			Dscore_Dnuisance_C(survTimeC(8),0) += survTimeT(1)/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  }
		}else{
		  score[0] = 1.0 - lastSurvC/survTimeC(2) - intFavorable[0] / denom; // (lower bound)
		  upperFavorable = 1.0 - intFavorable[1] / denom; // (upper bound)
		  if(returnIID>1){ // if((returnIID>1) && (lastSurvC>0)){
			Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,0) -= 1/survTimeC(2); // derivative regarding Sc(max)
			Dscore_Dnuisance_C(survTimeC(8),0) += lastSurvC/pow(survTimeC(2),2); // derivative regarding Sc(y_j)
		  }
		}

		if(returnIID > 1){ //		  if(returnIID > 1 && intFavorable[0]>0){
		  Dscore_Dnuisance_C(survTimeC(8),0) += intFavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
		  Dscore_Dnuisance_T(survTimeT(11),0) += intFavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
		  Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
		  Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
		}
		// Rcout << "end " << endl;

	  }else{
		// Rcout << "(4+b) ";

		intFavorable = calcIntegralScore_cpp(survJumpC, endpoint_C, lastSurvT, lastSurvC,
											 (returnIID > 1), 0, intDscore_Dnuisance_T, intDscore_Dnuisance_C); // -intFavorable is already the lower bound

		score[0] = -intFavorable[0] / denom; // (lower bound)
		upperFavorable = -intFavorable[1] / denom; // (upper bound)
		if(returnIID>1){ //		if((returnIID>1) && (intFavorable[0]>0)){
		  Dscore_Dnuisance_C(survTimeC(8),0) += intFavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
		  Dscore_Dnuisance_T(survTimeT(11),0) += intFavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
		  Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
		  Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
		}
		// Rcout << "end" << endl;
	  }
      
	  // unfavorable
	  if(diff <= -threshold){	
		// Rcout << "(4-a) ";

		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_C-threshold, lastSurvC, lastSurvT,
											   (returnIID > 1), 1, intDscore_Dnuisance_C, intDscore_Dnuisance_T); // -intUnfavorable is already the lower bound

		if(R_IsNA(survTimeC(4))==false){
		  score[1] = 1.0 - survTimeC(4)/survTimeT(5) - intUnfavorable[0] / denom; // (lower bound)
		  upperUnfavorable = 1.0 - survTimeC(4)/survTimeT(5) - intUnfavorable[1] / denom; // (upper bound)
		  if(returnIID>1){ // if((returnIID>1) && (survTimeC(4)>0)){
			Dscore_Dnuisance_T(survTimeC(10),1) -= 1/survTimeT(5); // derivative regarding St(y_j-tau)
			Dscore_Dnuisance_T(survTimeT(11),1) += survTimeC(4)/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  }
		}else{
		  score[1] = 1.0 - lastSurvT/survTimeT(5) - intUnfavorable[0] / denom; // (lower bound)
		  upperUnfavorable = 1.0 - intUnfavorable[1] / denom; // (upper bound)
		  if(returnIID>1){ //		  if((returnIID>1) && (lastSurvT>0)){
			Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,1) -= 1/survTimeT(5); // derivative regarding St(y_j-tau)
			Dscore_Dnuisance_T(survTimeT(11),1) += lastSurvT/pow(survTimeT(5),2); // derivative regarding St(x_i)
		  }
		}

		if(returnIID>1){ //		if((returnIID>1) && (intUnfavorable[0]>0)){
		  Dscore_Dnuisance_C(survTimeC(8),1) += intUnfavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
		  Dscore_Dnuisance_T(survTimeT(11),1) += intUnfavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
		  Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
		  Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
		}
		
		// Rcout << "end ";

	  }else{
		// Rcout << "(4-b) ";
		intUnfavorable = calcIntegralScore_cpp(survJumpT, endpoint_T, lastSurvC, lastSurvT,
											   (returnIID > 1), 1, intDscore_Dnuisance_C, intDscore_Dnuisance_T); // -intUnfavorable is already the lower bound
		
		score[1]= -intUnfavorable[0] / denom; // (lower bound)
		upperUnfavorable = -intUnfavorable[1] / denom;  // (upper bound)
		if(returnIID>1){ //		if((returnIID>1) && (intUnfavorable[0]>0)){
		  Dscore_Dnuisance_C(survTimeC(8),1) += intUnfavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
		  Dscore_Dnuisance_T(survTimeT(11),1) += intUnfavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
		  Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
		  Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
		}
		// Rcout << "end ";
	  }
	  // Rcout << endl;
    }}

  // ** compute neutral and uninformative
  // neutral
  double lowerNeutral = 1 - upperFavorable - upperUnfavorable;
  if(lowerNeutral >= 0.0){ // otherwise 0
	score[2] = lowerNeutral;
	if(returnIID>1){
	  Dscore_Dnuisance_C.col(2) = - Dscore_Dnuisance_C.col(0) - Dscore_Dnuisance_C.col(1);
	  Dscore_Dnuisance_T.col(2) = - Dscore_Dnuisance_T.col(0) - Dscore_Dnuisance_T.col(1);
	}
  }//else{ // initialized at 0
  //	score[2] = 0.0;
  //}

 
  // uninformative
  double upperUninformative = 1 - (score[0] + score[1] + score[2]);
  if(upperUninformative >= 0.0){ // otherwise 0
	score[3] = upperUninformative;
	if(returnIID>1){
	  Dscore_Dnuisance_C.col(3) = - Dscore_Dnuisance_C.col(0) - Dscore_Dnuisance_C.col(1) - Dscore_Dnuisance_C.col(2);
	  Dscore_Dnuisance_T.col(3) = - Dscore_Dnuisance_T.col(0) - Dscore_Dnuisance_T.col(1) - Dscore_Dnuisance_T.col(2);
	}
  }//else{ // initialized at 0
  //	score[3] = 0.0;
  //}

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
//' @param start [integer] time at which to start the integral.
//' @param lastSurv [numeric,>0] last survival value for the survival function in the second column.
//' @param lastdSurv [numeric,>0] last survival value for the survival function in the third column.
//' @param returnDeriv [logical] should the derivative regarding the survival parameters be return. 
//' @param column [integer] column of \code{derivSurv} and \code{derivSurvD} to be filled.
//' @param derivSurv [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival.
//' @param derivSurvD [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival used to compute the jumps.
//'
//' @keywords function Cpp internal
//' @export
// [[Rcpp::export]]
std::vector< double > calcIntegralScore_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
											bool returnDeriv, int column, arma::mat& derivSurv, arma::mat& derivSurvD){

  // WARNING here lastdSurv is actually abs(lastdSurv)
  // computes \int_t>tau S dS

  // survival contains:
  // [0] times
  // [1] S
  // [2] dS
  
  std::vector< double > integral(2,0.0); // lower and upper bound
  bool stopdSurv = true;
  int nJump = survival.n_rows;
  derivSurv.fill(0.0);
  derivSurvD.fill(0.0);

  if(nJump>0){    
	for(int iter_time=0 ; iter_time<nJump ; iter_time++){
	  // Rcout << "Jump: " << iter_time << "/" << nJump << endl;

	  // check whether we are at or beyond the last point of the survivla curve
	  if(R_IsNA(survival(iter_time,1))){
		integral[1] = integral[0] + lastSurv*survival(iter_time,2); // upper bound
		if(returnDeriv){
		  derivSurv(derivSurv.n_rows-1, column) += survival(iter_time,2); // derivative regarding S(t+tau)
		  derivSurvD(survival(iter_time,4), column) -= lastSurv; // derivative regarding dS(t-)
		  derivSurvD(survival(iter_time,5), column) += lastSurv; // derivative regarding dS(t+)
		}
		stopdSurv = false;
		break;
	  }else if(survival(iter_time,1) == 0){
		integral[1] = integral[0];
		stopdSurv = false;
		break;
	  }

	  // compute contribution at the current timepoint
	  if(survival(iter_time,0) > start){ // strict
		integral[0] += survival(iter_time,1)*survival(iter_time,2); // increment lower bound
		if(returnDeriv){
		  derivSurv(survival(iter_time,3), column) += survival(iter_time,2); // derivative regarding S(t+tau)
		  derivSurvD(survival(iter_time,4), column) -= survival(iter_time,1); // derivative regarding dS(t-)
		  derivSurvD(survival(iter_time,5), column) += survival(iter_time,1); // derivative regarding dS(t+)
		}
	  }
	}
	
	// add extra contribution to the bound
	if(stopdSurv){
	  if(lastdSurv > 0){
		integral[1] = integral[0] - survival(nJump-1,1)*lastdSurv; // upper bound (minus because dSurv = (0 - surv(tmax))
		if(returnDeriv){ //	  if(returnDeriv && (survival(nJump-1,1)*lastdSurv) ){
		  derivSurv(survival(nJump-1,3), column) -= lastdSurv; // derivative regarding S(t+tau)
		  derivSurvD(derivSurvD.n_rows-1, column) -= survival(nJump-1,1); // derivative regarding dS(t-)
		}
	  }else{
		integral[1] = integral[0];
	  }
	}
  }
  return(integral);
}



