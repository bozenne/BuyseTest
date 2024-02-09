// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:

inline std::vector< double > calcOnePair_Continuous(double diff, double threshold);
 
double normalCDF(double mu, double sigma, double x);
inline std::vector< double > calcOnePair_Gaussian(double mean_C, double mean_T, double std_C, double std_T, double rho, double threshold);

inline std::vector< double > calcOnePair_TTEgehan(double diff, double status_C, double status_T, double threshold, bool threshold0);

inline std::vector< double > calcOnePair_TTEgehan2(double diff, double status_C, double status_T, double threshold, bool threshold0);
 
inline std::vector< double > calcOnePair_SurvPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold, double restriction,
						   arma::rowvec survTimeC, arma::rowvec survTimeT,
						   const arma::mat& survJumpC, const arma::mat& survJumpT, double lastSurvC, double lastSurvT,
						   arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						   int p_C, int p_T, bool precompute, int returnIID);

inline std::vector< double > calcOnePair_CRPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold,
						 arma::rowvec cifTimeC_vec,  arma::rowvec cifTimeT_vec, const arma::mat& cifJumpC, const arma::mat& cifJumpT,
						 double lastCif1C, double lastCif1T, double lastCif2C, double lastCif2T,
						 arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						 int p_C, int p_T, bool precompute, int returnIID);

std::vector<double> calcIntegralSurv_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
					 bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD);

double calcIntegralCif_cpp(const arma::mat& cifJump, double start_val, double stop_val, arma::rowvec cifTimeT, double lastCIF, int type,
			   bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD);

// * calcOnePair_Continuous
// author Brice Ozenne
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
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}

// * calcOnePair_Gaussian
// author Brice Ozenne
double normalCDF(double mu, double sigma, double x){// Phi(-oo, x) 
  return 0.5 * (1 + std::erf((x-mu)/(sigma * std::sqrt(2))));
}

inline std::vector< double > calcOnePair_Gaussian(double mean_C, double mean_T, double std_C, double std_T, double rho, double threshold){

  // ** initialize
  std::vector< double > score(4,0.0);
  double diffstd = std::sqrt(std::pow(std_C,2.0) + std::pow(std_T,2.0) - 2*rho*std_C*std_T);
  // Rcpp::Rcout << "( )" << std_C << " " << std_T << " " << rho << " " << 2*rho*std_C*std_T << " " << diffstd << std::endl;  
  // ** score
  if(R_IsNA(mean_T) || R_IsNA(mean_C) || R_IsNA(diffstd)){ // missing data: uninformative
    score[3] = 1.0;
  }else{
    // P[Y>=X+t] = 1 - P[Y-X<t] where Y-X is N(muY-muX,\sqrt(sigmaY^2+sigma^X-2 rho sigmaY sigmaX))
    score[0] = 1 - normalCDF(mean_T-mean_C,diffstd,threshold);
    // P[X>=Y+t] = 1 - P[X-Y<t] where X-Y is N(muX-muY,\sqrt(sigmaX^2+sigma^Y-2 rho sigmaX sigmaY))
    score[1] = 1 - normalCDF(mean_C-mean_T,diffstd,threshold);
    score[2] = 1 - (score[0] + score[1]);
  }

  // ** export
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}


// * calcOnePair_TTEgehan
// author Brice Ozenne
inline std::vector< double > calcOnePair_TTEgehan(double diff, double status_C, double status_T, double threshold, bool threshold0){
  
  // ** initialize
  std::vector< double > score(4,0.0);
  // Rcpp::Rcout << diff << " " << status_T << " " << status_C << " " << threshold << std::endl;

  // ** score
  if(status_T==0.5 && status_C==0.5){
    
    score[2] = 1.0; // both restricted: neutral
    
  }else if(status_T==1 || status_T==0.5){
    
    if(status_C==1 || status_C==0.5){ // (treatment event/restricted, control event/restricted)
      
      if(diff >= threshold){         // >= tau    : favorable
        score[0] = 1.0;
      }else if(diff <= -threshold){ // <= -tau    : unfavorable
	score[1] = 1.0;
      }else{                        // ]-tau;tau[ : neutral
	score[2] = 1.0;
      }
      
    }else if(status_C==0){ // (treatment event, control right-censored)
	
      if(diff <= -threshold || (threshold0 && status_T==1 && diff==0)){ // <= -tau   : unfavorable
	// also unfavorable when 0-difference, infinitesimal threshold and observation is censored and the other is observed since event must happen after censoring
	score[1] = 1.0;
      }else{ // otheriwise : uninformative
	score[3] = 1.0;
      }
      
    }else if(status_C==2){ // (treatment event, control competing risk)

      if(status_T==1){
	score[1] = 1.0; //  unfavorable
      }else if(status_T==0.5){
	score[2] = 1.0; //  neutral
      }

    }
    
  }else if(status_T==0){
    
    if(status_C==1 || status_C==0.5){ // (treatment right-censored, control event/restricted)
    
      if(diff >= threshold ||  (threshold0 && status_C==1 && diff==0)){ // >= tau    : favorable
	// also favorable when 0-difference, infinitesimal threshold and observation is censored and the other is observed since event must happen after censoring
	score[0] = 1.0;
      }else{ // otherwise: uninformative
	score[3] = 1.0;
      }
    
    }else{ // otherwise: uninformative
      score[3] = 1.0;
    }
    
  }else if(status_T==2){ 

    if(status_C==1){ // (treatment competing risk, control event): favorable
      score[0] = 1.0;
    }else if(status_C==0){ // (treatment competing risk, control right-censored): uninformative
      score[3] = 1.0;
    }else if(status_C==2 || status_C==0.5){ // (treatment competing risk, control competing risk/restricted): neutral
      score[2] = 1.0;
    }
    
  }
  // ** export
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}

// * calcOnePair_TTEgehan2
// author Brice Ozenne
inline std::vector< double > calcOnePair_TTEgehan2(double diff, double status_C, double status_T, double threshold, bool threshold0){
  
  // ** initialize
  std::vector< double > score(4,0.0);
  // Rcpp::Rcout << diff << " " << status_T << " " << status_C << " " << threshold << std::endl;

  // ** score
  if(status_T==0.5 && status_C==0.5){
    
    score[2] = 1.0; // both restricted: neutral
    
  }else if(status_T==1 || status_T==0.5){
    
    if(status_C==1 || status_C==0.5){ // (treatment event/restricted, control event/restricted)
      
      if(diff >= threshold){         // >= tau    : favorable
        score[0] = 1.0;
      }else if(diff <= -threshold){ // <= -tau    : unfavorable
	score[1] = 1.0;
      }else{                        // ]-tau;tau[ : neutral
	score[2] = 1.0;
      }
      
    }else if(status_C==0){ // (treatment event, control left-censored)
	
      if(diff >= threshold || (threshold0 && status_T==1 && diff==0)){ // >= tau   : favorable
	// also favorable when 0-difference, infinitesimal threshold and observation is censored and the other is observed since event must happen after censoring
	score[0] = 1.0;
      }else{ // otheriwise : uninformative
	score[3] = 1.0;
      }
      
    }else if(status_C==2){ // (treatment event, control competing risk)

      if(status_T==1){
	score[1] = 1.0; //  unfavorable
      }else if(status_T==0.5){
	score[2] = 1.0; //  neutral
      }

    }
    
  }else if(status_T==0){
    
    if(status_C==1 || status_C==0.5){ // (treatment left-censored, control event/restricted)
    
      if(diff <= -threshold ||  (threshold0 && status_C==1 && diff==0)){ // <= tau    : unfavorable
	// also unfavorable when 0-difference, infinitesimal threshold and observation is censored and the other is observed since event must happen after censoring
	score[1] = 1.0;
      }else{ // otherwise: uninformative
	score[3] = 1.0;
      }
    
    }else{ // otherwise: uninformative
      score[3] = 1.0;
    }
    
  }else if(status_T==2){ 

    if(status_C==1){ // (treatment competing risk, control event): favorable
      score[0] = 1.0;
    }else if(status_C==0){ // (treatment competing risk, control censored): uninformative
      score[3] = 1.0;
    }else if(status_C==2 || status_C==0.5){ // (treatment competing risk, control competing risk/restricted): neutral
      score[2] = 1.0;
    }
    
  }
  
  // ** export
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}

// * calcOneScore_SurvPeron
// author Brice Ozenne
inline std::vector< double > calcOnePair_SurvPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold, double restriction,
						   arma::rowvec survTimeC, arma::rowvec survTimeT,
						   const arma::mat& survJumpC, const arma::mat& survJumpT, double lastSurvC, double lastSurvT,
						   arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						   int p_C, int p_T, bool precompute, int returnIID){
  
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
  double upperFavorable = 0.0;
  double upperUnfavorable = 0.0;
  bool testNoRestrictionC = (R_IsNA(restriction) || endpoint_C+threshold < restriction);
  bool testNoRestrictionT = (R_IsNA(restriction) || endpoint_T+threshold < restriction);
  
  if(returnIID > 1){
    Dscore_Dnuisance_C.fill(0.0); // initialized to 0
    Dscore_Dnuisance_T.fill(0.0); // initialized to 0
  }
  
  // ** deal with null survival
  // according to the survival the observation will die immediatly after the observation time.
  // so we treat it as if was an event
  if((status_C==0 || status_C==1/2) && (survTimeC(2) == 0) ){
    status_C = 1;
  }
  if((status_T==0 || status_T==1/2) && (survTimeT(5) == 0) ){
    status_T = 1;
  }
  	
  // Rcpp::Rcout << " (" << status_T << ";" << status_C << ")";
  // ** compute favorable and unfavorable
  if(status_T == 0.5 && status_C==0.5){

    // both restricted so neutral pair
    // score[0] = score[1] = upperFavorable =  upperUnfavorable = 0;      

  }else if(status_T==1 || status_T==0.5){
    if(status_C==1 || status_C==0.5){
     // Rcpp::Rcout << "(1) ";
      
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
	
    }else{ // status_C==0 

      // favorable
      if(diff >= threshold && testNoRestrictionC){
	// Rcpp::Rcout << "(2+a) ";
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
	// score[0] = upperFavorable = 0.0;
      }

      // unfavorable
      if(status_T==1){
	if(diff <= -threshold){ 
	  score[1] = 1.0;
	  upperUnfavorable = score[1];
	}else if(testNoRestrictionT){
	  // Rcpp::Rcout << "(2-b) " << "";
	  if((R_IsNA(survTimeT(3))==false) && (survTimeT(3) > 0)){
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
	}else{ // cannot be unfavorable when the outcome in the treatment group is close to the restriction (i.e. difference below threshold)
	  // score[1] = upperUnfavorable = 0.0;
	}
	
      }else{ // cannot be unfavorable when pair in the treatment group restricted
	// score[1] = upperUnfavorable = 0.0;
      }

    }
    
  }else{  // status_T==0 
    
    if(status_C==1 || status_C==0.5){ 
      
      // favorable
      if(status_C==1){
	if(diff >= threshold){ // 
	  score[0] = 1.0;
	  upperFavorable = score[0];
	}else if(testNoRestrictionC){
	  // Rcpp::Rcout << "(3+b) ";
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
	}else{ // cannot be favorable when the outcome in the control group is close to the restriction (i.e. difference below threshold)
	  // score[0] = upperFavorable = 0.0;
	}
      }else{ // cannot be favorable when pair in the control group restricted
	  // score[0] = upperFavorable = 0.0;
      }

      // unfavorable
      if(diff <= -threshold && testNoRestrictionT){
	// Rcpp::Rcout << "(3-a) ";
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
	// score[1] = upperUnfavorable = 0.0;
      }
      
     }else{ // status_T==0 && status_C==0

      double denom = survTimeT(5)*survTimeC(2);
      std::vector< double > intFavorable(2); 
      std::vector< double > intUnfavorable(2);
      arma::colvec intDscore_Dnuisance_C;
      arma::colvec intDscore_Dnuisance_T;
      if(returnIID>1 && precompute==false){
	intDscore_Dnuisance_C.resize(p_C); // initialized in calcIntegralSurv_cpp
	intDscore_Dnuisance_T.resize(p_T); // initialized in calcIntegralSurv_cpp
      }
	  
      // favorable
      if(diff >= threshold){
	// Rcpp::Rcout << "(4+a) ";

	if(precompute){
	  intFavorable[0] = survTimeT(13);
	  intFavorable[1] = survTimeT(14);	  
	}else{
	  intFavorable = calcIntegralSurv_cpp(survJumpC, endpoint_T-threshold, lastSurvT, lastSurvC,
					      (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	}
		
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
	  if(precompute){		
	    for(int iJump=survTimeT(18); iJump>=survTimeT(17); iJump--){
	      Dscore_Dnuisance_T(survJumpC(iJump,2),0) -= survJumpC(iJump,3)/denom;
	      Dscore_Dnuisance_C(survJumpC(iJump,4),0) -= survJumpC(iJump,5)/denom;
	      Dscore_Dnuisance_C(survJumpC(iJump,6),0) -= survJumpC(iJump,7)/denom;
	    }
	  }else{
	    Dscore_Dnuisance_C.col(0) -= intDscore_Dnuisance_C/denom;
	    Dscore_Dnuisance_T.col(0) -= intDscore_Dnuisance_T/denom;
	  }
	}

      }else if(testNoRestrictionC){
	// Rcpp::Rcout << "(4+b) ";

	if(precompute){
	  intFavorable[0] = survTimeC(15);
	  intFavorable[1] = survTimeC(16);
	}else{
	  intFavorable = calcIntegralSurv_cpp(survJumpC, endpoint_C, lastSurvT, lastSurvC,
					      (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C); 
	}
		
	score[0] = -intFavorable[0] / denom; // (lower bound)
	upperFavorable = -intFavorable[1] / denom; // (upper bound)
	if(returnIID>1){ //		if((returnIID>1) && (intFavorable[0]>n0)){
	  Dscore_Dnuisance_C(survTimeC(8),0) += intFavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
	  Dscore_Dnuisance_T(survTimeT(11),0) += intFavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
	  if(precompute){
	    for(int iJump=survTimeC(20); iJump>=survTimeC(19); iJump--){
	      Dscore_Dnuisance_T(survJumpC(iJump,2),0) -= survJumpC(iJump,3)/denom;
	      Dscore_Dnuisance_C(survJumpC(iJump,4),0) -= survJumpC(iJump,5)/denom;
	      Dscore_Dnuisance_C(survJumpC(iJump,6),0) -= survJumpC(iJump,7)/denom;
	    }
	  }else{
	    Dscore_Dnuisance_C.col(0) -= intDscore_Dnuisance_C/denom;
	    Dscore_Dnuisance_T.col(0) -= intDscore_Dnuisance_T/denom;
	  }
	}
      }else{
	// cannot be favorable when the outcome in the control group is close to the restriction (i.e. difference below threshold)
	// score[0] = upperFavorable = 0.0;
      }
      
      // unfavorable
      if(diff <= -threshold){	
	// Rcpp::Rcout << "(4-a) ";

	if(precompute){
	  intUnfavorable[0] = survTimeC(13);
	  intUnfavorable[1] = survTimeC(14);
	}else{
	  intUnfavorable = calcIntegralSurv_cpp(survJumpT, endpoint_C-threshold, lastSurvC, lastSurvT,
						(returnIID > 1), intDscore_Dnuisance_C, intDscore_Dnuisance_T); // -intUnfavorable is already the lower bound
	}

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
	  if(precompute){
	    for(int iJump=survTimeC(18); iJump>=survTimeC(17); iJump--){
	      Dscore_Dnuisance_C(survJumpT(iJump,2),1) -= survJumpT(iJump,3)/denom;
	      Dscore_Dnuisance_T(survJumpT(iJump,4),1) -= survJumpT(iJump,5)/denom;
	      Dscore_Dnuisance_T(survJumpT(iJump,6),1) -= survJumpT(iJump,7)/denom;
	    }
	  }else{
	    Dscore_Dnuisance_C.col(1) -= intDscore_Dnuisance_C/denom;
	    Dscore_Dnuisance_T.col(1) -= intDscore_Dnuisance_T/denom;
	  }
	}
		
      }else if(testNoRestrictionT){
	// Rcpp::Rcout << "(4-b) ";
	    
	if(precompute){
	  intUnfavorable[0] = survTimeT(15);
	  intUnfavorable[1] = survTimeT(16);
	}else{
	  intUnfavorable = calcIntegralSurv_cpp(survJumpT, endpoint_T, lastSurvC, lastSurvT,
						(returnIID > 1), intDscore_Dnuisance_C, intDscore_Dnuisance_T); // -intUnfavorable is already the lower bound
	}
		
	score[1]= -intUnfavorable[0] / denom; // (lower bound)
	upperUnfavorable = -intUnfavorable[1] / denom;  // (upper bound)
	if(returnIID>1){ //		if((returnIID>1) && (intUnfavorable[0]>0)){
	  Dscore_Dnuisance_C(survTimeC(8),1) += intUnfavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
	  Dscore_Dnuisance_T(survTimeT(11),1) += intUnfavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
	  if(precompute){
	    for(int iJump=survTimeT(20); iJump>=survTimeT(19); iJump--){
	      Dscore_Dnuisance_C(survJumpT(iJump,2),1) -= survJumpT(iJump,3)/denom;
	      Dscore_Dnuisance_T(survJumpT(iJump,4),1) -= survJumpT(iJump,5)/denom;
	      Dscore_Dnuisance_T(survJumpT(iJump,6),1) -= survJumpT(iJump,7)/denom;
	    }
	  }else{
	    Dscore_Dnuisance_C.col(1) -= intDscore_Dnuisance_C/denom;
	    Dscore_Dnuisance_T.col(1) -= intDscore_Dnuisance_T/denom;
	  }
	}
      }else{
	// cannot be unfavorable when the outcome in the treatment group is close to the restriction (i.e. difference below threshold)
	// score[1] = upperUnfavorable = 0.0;
      }
      // Rcpp::Rcout << std::endl;
    }
  }

  // Rcpp::Rcout << score[0] << " " << score[1] <<  " " << upperFavorable <<  " " << upperUnfavorable << std::endl;
  
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
  // Rcpp::Rcout << " (" << score[0] << " " << score[1] << ") "  << std::endl << Dscore_Dnuisance_C << std::endl;
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << " (upper) " << upperFavorable << " " << upperUnfavorable << std::endl;
  return(score);  
}


// * calcOnePair_CRPeron
// author Eva Cantagallo
inline std::vector< double > calcOnePair_CRPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold,
						 arma::rowvec cifTimeC_vec,  arma::rowvec cifTimeT_vec, const arma::mat& cifJumpC, const arma::mat& cifJumpT,
						 double lastCif1C, double lastCif1T, double lastCif2C, double lastCif2T,
						 arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						 int p_C, int p_T, bool precompute, int returnIID) {

  // cifTimeC and cifTimeT: cumulative incidence at control/treatment observation times
  //        [1]    times
  //        [2-4]  cif of event of interest estimated in the control arm: time - threshold, time, time + threshold
  //        [5-7]  cif of event of interest estimated in the treatment arm: time - threshold, time, time + threshold
  //        [8]  cif of competing event estimated in the control arm: time
  //        [9]  cif of competing event estimated in the treatment arm: time

  // cifJumpC: cumulative incidence of event of interest in control group at jump times
  //        [1]  jump times in control group (unique values in ascending order)
  //        [2-3]  cif of the treatment group at times-threshold and times+threshold
  //        [4]  d(cif) of control group estimated at time

  double diff = endpoint_T - endpoint_C;
  std::vector< double > score(4, 0.0); // [0] favorable, [1] unfavorable, [2] neutral, [3] uninformative
  // Rcpp::Rcout << diff << " " << status_T << " " << status_C << " " << threshold << std::endl;
  double denomC = 1 - cifTimeC_vec(2) - cifTimeC_vec(7);
  double denomT = 1 - cifTimeT_vec(5) - cifTimeT_vec(7);
 
  if(returnIID > 1){
    Dscore_Dnuisance_C.fill(0.0); // initialized to 0
    Dscore_Dnuisance_T.fill(0.0); // initialized to 0
  }

  if(status_T == 2) {
    if(status_C == 2) { // (2,2)
      score[2] = 1.0; // systematically neutral competing
    } else if(status_C == 1){ // (2,1)
      score[0] = 1.0; // systematically favorable
    } else if(status_C == 0) { // (2,0)
      score[0] = (lastCif1C - cifTimeC_vec(2))/denomC;
      if(returnIID>1){
	Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -1/denomC + score[0]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	Dscore_Dnuisance_C(cifTimeC_vec(14),0) += score[0]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
      }
      // score[1] = 0
      score[2] = (lastCif2C - cifTimeC_vec(7))/denomC;
      if(returnIID>1){
	Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,2) += 1/denomC; // derivative regarding F_2^C(\infty) (i.e. lastCif2C)
	Dscore_Dnuisance_C(cifTimeC_vec(9),2) += score[2]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	Dscore_Dnuisance_C(cifTimeC_vec(14),2) += -1/denomC + score[2]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
      }
      // score[3] = 0 // since the treated patient had the competing event
    }
  } else if(status_T == 1){
    if(status_C == 2){ // (1,2)
      score[1] = 1.0; // systematically defavorable
    } else if(status_C == 1){ // (1,1)
      if(diff >= threshold) {
        score[0] = 1.0;
      } else if(diff <= -threshold) {
        score[1] = 1.0;
      } else { // |diff| < threshold
        score[2] = 1.0;
      }
    } else if(status_C == 0) { // (1,0)
      if(diff >= threshold) {
        if(R_IsNA(cifTimeT_vec(1)) == false) {
          score[0] = (cifTimeT_vec(1) - cifTimeC_vec(2))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeT_vec(8),0) += 1/denomC; // derivative regarding F_1^C(Tt-\tau) (i.e. cifTimeT_vec(1)
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -1/denomC + score[0]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),0) += score[0]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        } else {
          score[0] = (lastCif1C - cifTimeC_vec(2))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -1/denomC + score[0]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),0) += score[0]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        }
        if(R_IsNA(cifTimeT_vec(3)) == false) {
          score[1] = (lastCif1C - cifTimeT_vec(3) + lastCif2C - cifTimeC_vec(7))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,1) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(9),1) += score[1]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeT_vec(10),1) += -1/denomC; // derivative regarding F_1^C(Tt+\tau) (i.e. cifTimeT_vec(3))
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += 1/denomC; // derivative regarding F_2^C(\infty) (i.e. lastCif2C)
	    Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -1/denomC + score[1]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        } else {
          score[1] = (lastCif2C - cifTimeC_vec(7))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += 1/denomC; // derivative regarding F_2^C(\infty) (i.e. lastCif2C)
	    Dscore_Dnuisance_C(cifTimeC_vec(9),1) += score[1]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -1/denomC + score[1]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        }
        if((R_IsNA(cifTimeT_vec(3)) == false) && (R_IsNA(cifTimeT_vec(1)) == false)) {
          score[2] = (cifTimeT_vec(3) - cifTimeT_vec(1))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeC_vec(9),2) += score[2]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeT_vec(10),2) += 1/denomC; // derivative regarding F_1^C(Tt+\tau) (i.e. cifTimeT_vec(3))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),2) += score[2]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	    Dscore_Dnuisance_C(cifTimeT_vec(8),2) += -1/denomC; // derivative regarding F_1^C(Tt-\tau) (i.e. cifTimeT_vec(1))
	  }

        } else if ((R_IsNA(cifTimeT_vec(3)) == true) && (R_IsNA(cifTimeT_vec(1)) == false)) {
          score[2] = (lastCif1C - cifTimeT_vec(1))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeC_vec(9),2) += score[2]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,2) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(14),2) += score[2]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	    Dscore_Dnuisance_C(cifTimeT_vec(8),2) += -1/denomC; // derivative regarding F_1^C(Tt-\tau) (i.e. cifTimeT_vec(1))
	  }
        } else {
          score[2] = 0.0;
        }
      } else if(diff <= -threshold) {
        score[1] = 1.0;
      } else { // |diff| < threshold
        if(R_IsNA(cifTimeT_vec(3)) == false) {
          score[1] = (lastCif1C - cifTimeT_vec(3) + lastCif2C - cifTimeC_vec(7))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,1) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(9),1) += score[1]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeT_vec(10),1) += -1/denomC; // derivative regarding F_1^C(Tt+\tau) (i.e. cifTimeT_vec(3))
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += 1/denomC; // derivative regarding F_2^C(\infty) (i.e. lastCif2C)
	    Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -1/denomC + score[1]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
          score[2] = (cifTimeT_vec(3) - cifTimeC_vec(2))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeC_vec(9),2) += -1/denomC + score[2]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeT_vec(10),2) += 1/denomC; // derivative regarding F_1^C(Tt+\tau) (i.e. cifTimeT_vec(3))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),2) += score[2]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        } else {
          score[1] = (lastCif2C - cifTimeC_vec(7))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeC_vec(9),1) += score[1]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -1/denomC + score[1]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
          score[2] = (lastCif1C - cifTimeC_vec(2))/denomC;
	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeC_vec(9),2) += -1/denomC + score[2]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,2) += 1/denomC; // derivative regarding F_1^C(\infty) (i.e. lastCif1C)
	    Dscore_Dnuisance_C(cifTimeC_vec(14),2) += score[2]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  }
        }
      }
    }
  } else { // status_T == 0
    if(status_C == 2) { // (0,2)
      score[1] = (lastCif1T - cifTimeT_vec(5))/denomT;
      if(returnIID>1){
	Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,1) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -1/denomT + score[1]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	Dscore_Dnuisance_T(cifTimeT_vec(14),1) += score[1]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
      }

      score[2] = (lastCif2T - cifTimeT_vec(7))/denomT;
      if(returnIID>1){
	Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,2) += 1/denomT; // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	Dscore_Dnuisance_T(cifTimeT_vec(12),2) += score[2]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	Dscore_Dnuisance_T(cifTimeT_vec(14),2) += -1/denomT + score[2]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
      }
    } else if(status_C == 1) { // (0,1)
      if(diff >= threshold) {
        score[0] = 1.0;
      } else if(diff <= -threshold) {
        if(R_IsNA(cifTimeC_vec(6)) == false) {
          score[0] = (lastCif1T - cifTimeC_vec(6) + lastCif2T - cifTimeT_vec(7))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,0) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += 1/denomT; // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeC_vec(13),0) += -1/denomT; // derivative regarding F_1^T(Ct+\tau) (i.e. cifTimeC_vec(6))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -1/denomT + score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        } else {
          score[0] = (lastCif2T - cifTimeT_vec(7))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += 1/denomT; // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -1/denomT + score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        }
        if(R_IsNA(cifTimeC_vec(4)) == false) {
          score[1] = (cifTimeC_vec(4) - cifTimeT_vec(5))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(cifTimeC_vec(11),1) += 1/denomT; // derivative regarding F_1^T(Ct-\tau) (i.e. cifTimeC_vec(4))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -1/denomT + score[1]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),1) += score[1]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        } else {
          score[1] = (lastCif1T - cifTimeT_vec(5))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,0) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -1/denomT + score[1]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),1) += score[1]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        }
        if((R_IsNA(cifTimeC_vec(6)) == false) && (R_IsNA(cifTimeC_vec(4)) == false)) {
          score[2] = (cifTimeC_vec(6) - cifTimeC_vec(4))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(cifTimeC_vec(13),2) += 1/denomT; // derivative regarding F_1^T(Ct+\tau) (i.e. cifTimeT_vec(6))
	    Dscore_Dnuisance_T(cifTimeC_vec(11),2) += -1/denomT; // derivative regarding F_1^T(Ct-\tau) (i.e. cifTimeT_vec(4))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),2) += score[2]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),2) += score[2]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        } else if ((R_IsNA(cifTimeC_vec(6)) == true) && (R_IsNA(cifTimeC_vec(4)) == false)) {
          score[2] = (lastCif1T - cifTimeC_vec(4))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 -1,2) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(cifTimeC_vec(11),2) += -1/denomT; // derivative regarding F_1^T(Ct-\tau) (i.e. cifTimeT_vec(4))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),2) += score[2]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),2) += score[2]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        } else {
          score[2] = 0.0;
        }
      } else { // |diff| < threshold
        if(R_IsNA(cifTimeC_vec(6)) == false) {
          score[0] = (lastCif1T - cifTimeC_vec(6) + lastCif2T - cifTimeT_vec(7))/denomT;
          if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,0) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += 1/denomT; // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeC_vec(13),0) += -1/denomT; // derivative regarding F_1^T(Ct+\tau) (i.e. cifTimeC_vec(6))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -1/denomT + score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
	  score[2] = (cifTimeC_vec(6) - cifTimeT_vec(5))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(cifTimeC_vec(13),2) += 1/denomT; // derivative regarding F_1^T(Ct+\tau) (i.e. cifTimeT_vec(4))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),2) += -1/denomT + score[2]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),2) += score[2]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        } else {
          score[0] = (lastCif2T - cifTimeT_vec(7))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += 1/denomT; // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -1/denomT + score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
          score[2] = (lastCif1T - cifTimeT_vec(5))/denomT;
	  if(returnIID>1){
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 -1,2) += 1/denomT; // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),2) += -1/denomT + score[2]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),2) += score[2]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  }
        }
      }
    } else if (status_C == 0) { // (0,0)
      double prob21 = (lastCif2T - cifTimeT_vec(7))*(lastCif1C - cifTimeC_vec(2))/(denomT*denomC);
      double prob12 = (lastCif2C - cifTimeC_vec(7))*(lastCif1T - cifTimeT_vec(5))/(denomT*denomC);
      double prob22 = (lastCif2C - cifTimeC_vec(7))*(lastCif2T - cifTimeT_vec(7))/(denomT*denomC);
      // lower bound of each integral
      double intFav, intDefav, intNeutralEvent1, intNeutralEvent2;
      arma::colvec intDscore_Dnuisance_C;
      arma::colvec intDscore_Dnuisance_T;
      if(returnIID>1 && precompute==false){
	intDscore_Dnuisance_C.resize(p_C); // initialized in calcIntegralCif_cpp
	intDscore_Dnuisance_T.resize(p_T); // initialized in calcIntegralCif_cpp
      }

      if(diff >= threshold) {
	// favorable
	if (R_IsNA(cifTimeT_vec(1)) == false) {
		  
	  score[0] = (cifTimeT_vec(1) - cifTimeC_vec(2))*(lastCif1T - cifTimeT_vec(5))/(denomC*denomT);

	  if(returnIID>1){
	    Dscore_Dnuisance_C(cifTimeT_vec(8),0) += (lastCif1T - cifTimeT_vec(5))/(denomC*denomT); // derivative regarding F_1^C(Tt-\tau) (i.e. cifTimeT_vec(1))
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -(lastCif1T - cifTimeT_vec(5))/(denomC*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,0) += (cifTimeT_vec(1) - cifTimeC_vec(2))/(denomC*denomT); // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += -(cifTimeT_vec(1) - cifTimeC_vec(2))/(denomC*denomT); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
			
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += score[0]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),0) += score[0]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  }
	}else {
	  score[0] = (lastCif1C - cifTimeC_vec(2))*(lastCif1T - cifTimeT_vec(5))/(denomC*denomT);

	  if(returnIID>1){
	    Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += (lastCif1T - cifTimeT_vec(5))/(denomC*denomT); // derivative regarding F_1^C(Tt-\tau) (i.e. cifTimeT_vec(1))
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -(lastCif1T - cifTimeT_vec(5))/(denomC*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,0) += (lastCif1C - cifTimeC_vec(2))/(denomC*denomT); // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += -(lastCif1C - cifTimeC_vec(2))/(denomC*denomT); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
			
	    Dscore_Dnuisance_C(cifTimeC_vec(9),0) += score[0]/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	    Dscore_Dnuisance_C(cifTimeC_vec(14),0) += score[0]/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	    Dscore_Dnuisance_T(cifTimeT_vec(12),0) += score[0]/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	    Dscore_Dnuisance_T(cifTimeT_vec(14),0) += score[0]/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  }
	}

	intFav = calcIntegralCif_cpp(cifJumpC, endpoint_T - threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 1, // cifTimeT_vec(5);
				     (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	score[0] += intFav/(denomT*denomC) + prob21;
	if(returnIID>1){
	  // intFav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob21
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += (lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -(lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += prob21/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += prob21/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += prob21/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += prob21/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}

	// unfavorable
	intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_T + threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 2,
				       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	score[1] = intDefav/(denomT*denomC) + prob12;
	if(returnIID>1){
	  // intDefav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob12
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += (lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -(lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,1) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += prob12/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += prob12/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += prob12/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += prob12/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}

	// neutral
	intNeutralEvent1 = calcIntegralCif_cpp(cifJumpC, endpoint_T - threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 3,
					       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	score[2] = prob22 += intNeutralEvent1/(denomT*denomC);
	if(returnIID>1){
	  // prob22		 
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,2) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,2) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += prob22/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += prob22/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += prob22/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += prob22/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
		
	  // intNeutralEvent1
	  Dscore_Dnuisance_C.col(2) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(2) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}

        intNeutralEvent2 = calcIntegralCif_cpp(cifJumpC, endpoint_T + threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 4,
					       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[2] += intNeutralEvent2/(denomT*denomC);
	if(returnIID>1){
		
	  // intNeutralEvent2
	  Dscore_Dnuisance_C.col(2) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(2) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += intNeutralEvent2/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += intNeutralEvent2/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += intNeutralEvent2/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += intNeutralEvent2/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
      } else if(diff <= -threshold) {
        intFav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 1,
				     (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[0] = intFav/(denomT*denomC) + prob21;
	if(returnIID>1){
	  // intFav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob21
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += (lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -(lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += prob21/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += prob21/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += prob21/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += prob21/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
		
	intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 2,
				       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[1] = intDefav/(denomT*denomC) + prob12;
	if(returnIID>1){
	  // intDefav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob12
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += (lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -(lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,1) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += prob12/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += prob12/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += prob12/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += prob12/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}

	intNeutralEvent1 = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 4,
					       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[2] = prob22 + intNeutralEvent1/(denomT*denomC);
	if(returnIID>1){
	  // prob22		 
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,2) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,2) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += prob22/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += prob22/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += prob22/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += prob22/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
		
	  // intNeutralEvent1
	  Dscore_Dnuisance_C.col(2) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(2) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
		
      } else { // |diff| < threshold
        intFav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 1,
				     (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[0] = intFav/(denomT*denomC) + prob21;
	if(returnIID>1){
	  // intFav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += intFav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += intFav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob21
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,0) += (lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += -(lastCif1C - cifTimeC_vec(2))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows/2 - 1,0) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),0) += prob21/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),0) += prob21/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),0) += prob21/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),0) += prob21/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
		
        intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_T+threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 2,
				       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
        score[1] = intDefav/(denomT*denomC) + prob12;
	if(returnIID>1){
	  // intDefav
	  Dscore_Dnuisance_C.col(1) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(1) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += intDefav/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += intDefav/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  // prob12
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,1) += (lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += -(lastCif1T - cifTimeT_vec(5))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows/2 - 1,1) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(\infty) (i.e. lastCif1T)
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),1) += prob12/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),1) += prob12/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),1) += prob12/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),1) += prob12/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
		
	intNeutralEvent1 = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 3,
					       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	score[2] = prob22 + intNeutralEvent1/(denomT*denomC);
	if(returnIID>1){
	  // prob22		 
	  Dscore_Dnuisance_C(Dscore_Dnuisance_C.n_rows - 1,2) += (lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(\infty) (i.e. lastCif12C)
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += -(lastCif2T - cifTimeT_vec(7))/(denomT*denomC); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(Dscore_Dnuisance_T.n_rows - 1,2) += (lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(\infty) (i.e. lastCif2T)
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += -(lastCif2C - cifTimeC_vec(7))/(denomT*denomC); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
			
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += prob22/denomC; // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += prob22/denomC; // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += prob22/denomT; // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += prob22/denomT; // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
		
	  // intNeutralEvent1
	  Dscore_Dnuisance_C.col(2) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(2) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += intNeutralEvent1/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += intNeutralEvent1/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
		
	intNeutralEvent2 = calcIntegralCif_cpp(cifJumpC, endpoint_T+threshold, endpoint_T + threshold, cifTimeT_vec, lastCif1T, 4,
					       (returnIID > 1), intDscore_Dnuisance_T, intDscore_Dnuisance_C);
	score[2] += intNeutralEvent2/(denomT*denomC);
	if(returnIID>1){
		
	  // intNeutralEvent2
	  Dscore_Dnuisance_C.col(2) += intDscore_Dnuisance_C/(denomC*denomT);
	  Dscore_Dnuisance_T.col(2) += intDscore_Dnuisance_T/(denomC*denomT);
	  Dscore_Dnuisance_C(cifTimeC_vec(9),2) += intNeutralEvent2/(pow(denomC,2)*denomT); // derivative regarding F_1^C(Ct) (i.e. cifTimeC_vec(2))
	  Dscore_Dnuisance_C(cifTimeC_vec(14),2) += intNeutralEvent2/(pow(denomC,2)*denomT); // derivative regarding F_2^C(Ct) (i.e. cifTimeC_vec(7))
	  Dscore_Dnuisance_T(cifTimeT_vec(12),2) += intNeutralEvent2/(denomC*pow(denomT,2)); // derivative regarding F_1^T(Tt) (i.e. cifTimeT_vec(5))
	  Dscore_Dnuisance_T(cifTimeT_vec(14),2) += intNeutralEvent2/(denomC*pow(denomT,2)); // derivative regarding F_2^T(Tt) (i.e. cifTimeT_vec(7))
	}
      }
    }
  }
  
  score[3] = 1 - (score[0] + score[1] + score[2]);
  if(returnIID>1){
    Dscore_Dnuisance_C.col(3) = - Dscore_Dnuisance_C.col(0) - Dscore_Dnuisance_C.col(1) - Dscore_Dnuisance_C.col(2);
    Dscore_Dnuisance_T.col(3) = - Dscore_Dnuisance_T.col(0) - Dscore_Dnuisance_T.col(1) - Dscore_Dnuisance_T.col(2);
  }
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return score;
}


// * calcIntegralSurv_cpp
//' @title C++ Function Computing the Integral Terms for the Peron Method in the survival case. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored.
//' 
//' @param survival [matrix] Contains the jump times in the first column,
//' the survival in the other arm at times plus threshold in the second column,
//' and the jump in survival in the third column.
//' @param start [integer] time at which to start the integral.
//' @param lastSurv [numeric,>0] last survival value for the survival function in the second column.
//' @param lastdSurv [numeric,>0] last survival value for the survival function in the third column.
//' @param returnDeriv [logical] should the derivative regarding the survival parameters be return. 
//' @param derivSurv [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival.
//' @param derivSurvD [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival used to compute the jumps.
//'
//' @keywords function Cpp internal
//' @author Brice Ozenne
//' @export
// [[Rcpp::export(".calcIntegralSurv_cpp")]]
std::vector< double > calcIntegralSurv_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
					   bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD){

  // WARNING here lastdSurv is actually abs(lastdSurv)
  // computes \int_t>tau S dS

  // survival contains:
  // [0] times
  // [1] S
  // [2] dS
  
  std::vector< double > integral(2,0.0); // lower and upper bound
  int nJump = survival.n_rows;
  if(returnDeriv){
    derivSurv.fill(0.0);
    derivSurvD.fill(0.0);
  }
  
  int iter_time = 0;
  while(iter_time<nJump){
    // Rcpp::Rcout << "Jump: " << iter_time << "/" << nJump << std::endl;
	  
    // compute contribution at the current timepoint
    if(survival(iter_time,0) > start){
      if(R_IsNA(survival(iter_time,1))){ // if we are beyond the last known point of the survival curve
	integral[1] += lastSurv*survival(iter_time,2); // increment upper bound
      }else{  
	integral[0] += survival(iter_time,1)*survival(iter_time,2); // increment lower bound
	integral[1] = integral[0]; // upper bound = lower bound
	if(returnDeriv){
	  derivSurv(survival(iter_time,3)) += survival(iter_time,2); // derivative regarding S(t+tau)
	  derivSurvD(survival(iter_time,4)) -= survival(iter_time,1); // derivative regarding dS(t-)
	  derivSurvD(survival(iter_time,5)) += survival(iter_time,1); // derivative regarding dS(t+)
	}
      }
    }
    iter_time++;
  }
    
  // add extra contribution to the bound (minus because dSurv = (0 - surv(tmax))
  if(lastdSurv > 0){
    if(nJump==0){
      integral[1] -= lastSurv*lastdSurv; 
    }else if(R_IsNA(survival(iter_time-1,1))){
      integral[1] -= lastSurv*lastdSurv; 
    }else{
      integral[1] -= survival(iter_time-1,1)*lastdSurv;
    }
  }
  
  return(integral);
}

// * calcIntegralCif_cpp
//' @title C++ Function Computing the Integral Terms for the Peron Method in the presence of competing risks (CR).
//' @description Compute the integral with respect to the jump in CIF for pairs where both outcomes are censored.
//'
//' @param cifJump [matrix] cif[1] = jump times in control group (event of interest), cif[2-3] = CIF of event of interest in group
//' T at times - tau and times + tau, cif[4] : jump in cif of control group at times (event of interest).
//' @param start_val [numeric] Time at which to start the integral.
//' @param stop_val [numeric] Time at which to stop the integral.
//' @param cifTimeT [numeric] CIF of event of interest in group T evaluated at observed time of treatment patient.
//' @param lastCIF [numeric, >0] last value of CIF of event type 1 in group T.
//' @param type [numeric] Indicates the type of integral to compute (1 for wins, 2 for losses, 3 for neutral pairs with two
//' events of interest - integral with t+tau and xi - and 4 for neutral pairs with two events of interest - integral with
//' t+tau and t-tau).
//' @param returnDeriv [logical] should the derivative regarding the survival parameters be return. 
//' @param derivSurv [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival.
//' @param derivSurvD [matrix] matrix column filled of 0 whose number of rows is the number of parameters of the survival used to compute the jumps.
//'
//' @keywords function Cpp internal
//' @author Eva Cantagallo
//' @export
// [[Rcpp::export(".calcIntegralCif_cpp")]]
double calcIntegralCif_cpp(const arma::mat& cifJump, double start_val, double stop_val, arma::rowvec cifTimeT, double lastCIF, int type,
			   bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD){

  double integral = 0.0;
  int nJump = cifJump.n_rows;
  if(returnDeriv){
    derivSurv.fill(0.0);
    derivSurvD.fill(0.0);
  }
  
  if (nJump > 0) {
    if(type == 1) {
      for(int i = 0; i<nJump; i++){
        if(R_IsNA(cifJump(i,2))) {break;}
        if(cifJump(i,0) > start_val) {
          integral = integral + (lastCIF - cifJump(i, 2))*cifJump(i, 3);
	  if(returnDeriv){
	    derivSurv(derivSurv.n_rows/2) += cifJump(i, 3); // derivative regarding CIF(\infty)
	    derivSurv(cifJump(i, 5)) -= cifJump(i, 3); // derivative regarding CIF(t+\tau)
	    derivSurvD(cifJump(i, 6)) -= (lastCIF - cifJump(i, 2)); // derivative regarding dCIF(t-)
	    derivSurvD(cifJump(i, 7)) += (lastCIF - cifJump(i, 2)); // derivative regarding dCIF(t+)
	  }
        }
      }
    } else if(type == 2) {
      double CIF_t = cifTimeT(5);
      int indexCIF_t = cifTimeT(12);
      for(int i = 0; i<nJump; i++){
	if(cifJump(i,0) > start_val) {
	  if(R_IsNA(cifJump(i,1))){
	    integral = integral + (lastCIF - CIF_t)*cifJump(i, 3);
	    if(returnDeriv){
	      derivSurv(derivSurv.n_rows/2) += cifJump(i, 3); // derivative regarding CIF(\infty)
	      derivSurv(indexCIF_t) -= cifJump(i, 3); // derivative regarding CIF(Tt)
	      derivSurvD(cifJump(i, 6)) -= (lastCIF - CIF_t); // derivative regarding dCIF(t-)
	      derivSurvD(cifJump(i, 7)) += (lastCIF - CIF_t); // derivative regarding dCIF(t+)
	    }
	  }else{
	    integral = integral + (cifJump(i,1) - CIF_t)*cifJump(i, 3);
	    if(returnDeriv){
	      derivSurv(cifJump(i,4)) += cifJump(i, 3); // derivative regarding CIF(t-\tau)
	      derivSurv(indexCIF_t) -= cifJump(i, 3); // derivative regarding CIF(Tt)
	      derivSurvD(cifJump(i, 6)) -= (cifJump(i,1) - CIF_t); // derivative regarding dCIF(t-)
	      derivSurvD(cifJump(i, 7)) += (cifJump(i,1) - CIF_t); // derivative regarding dCIF(t+)
	    }
	  }
	}
      }
    }
  } else if(type == 3) {
    double CIF_t = cifTimeT(5);
    int indexCIF_t = cifTimeT(12);
    for(int i = 0; i<nJump; i++){
      if((cifJump(i,0) > start_val) && (cifJump(i,0) <= stop_val)) {
	if(R_IsNA(cifJump(i,2))){
	  integral = integral + (lastCIF - CIF_t)*cifJump(i, 3);
	  if(returnDeriv){
	    derivSurv(derivSurv.n_rows/2) += cifJump(i, 3); // derivative regarding CIF(\infty)
	    derivSurv(indexCIF_t) -= cifJump(i, 3); // derivative regarding CIF(Tt)
	    derivSurvD(cifJump(i, 6)) -= (lastCIF - CIF_t); // derivative regarding dCIF(t-)
	    derivSurvD(cifJump(i, 7)) += (lastCIF - CIF_t); // derivative regarding dCIF(t+)
	  }
	}else{
	  integral = integral + (cifJump(i,2) - CIF_t)*cifJump(i, 3);
	  if(returnDeriv){
	    derivSurv(cifJump(i,5)) += cifJump(i, 3); // derivative regarding CIF(t+\tau)
	    derivSurv(indexCIF_t) -= cifJump(i, 3); // derivative regarding CIF(Tt)
	    derivSurvD(cifJump(i, 6)) -= (cifJump(i,2) - CIF_t); // derivative regarding dCIF(t-)
	    derivSurvD(cifJump(i, 7)) += (cifJump(i,2) - CIF_t); // derivative regarding dCIF(t+)
	  }
	}
      }
    }
  } else if(type == 4) {
    for(int i = 0; i<nJump; i++){
      if(R_IsNA(cifJump(i,1))) {break;}
      if(cifJump(i,0) > start_val) {
	if(R_IsNA(cifJump(i,2))){
	  integral = integral + (lastCIF - cifJump(i, 1)) * cifJump(i, 3);
	  if(returnDeriv){
	    derivSurv(derivSurv.n_rows/2) += cifJump(i, 3); // derivative regarding CIF(\infty)
	    derivSurv(cifJump(i,4)) -= cifJump(i, 3); // derivative regarding CIF(t-\tau)
	    derivSurvD(cifJump(i, 6)) -= (lastCIF - cifJump(i, 1)); // derivative regarding dCIF(t-)
	    derivSurvD(cifJump(i, 7)) += (lastCIF - cifJump(i, 1)); // derivative regarding dCIF(t+)
	  }
	}else{
	  integral = integral + (cifJump(i, 2) - cifJump(i, 1)) * cifJump(i, 3);
	  if(returnDeriv){
	    derivSurv(cifJump(i,5)) += cifJump(i, 3); // derivative regarding CIF(t+\tau)
	    derivSurv(cifJump(i,4)) -= cifJump(i, 3); // derivative regarding CIF(t-\tau)
	    derivSurvD(cifJump(i, 6)) -= (cifJump(i, 2) - cifJump(i, 1)); // derivative regarding dCIF(t-)
	    derivSurvD(cifJump(i, 7)) += (cifJump(i, 2) - cifJump(i, 1)); // derivative regarding dCIF(t+)
	  }
	}
      }
    }
  }

  return integral;

}
