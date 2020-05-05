// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:

inline std::vector< double > calcOnePair_Continuous(double diff, double threshold);
 
inline std::vector< double > calcOnePair_TTEgehan(double diff, double status_C, double status_T, double threshold);

inline std::vector< double > calcOnePair_TTEgehan2(double diff, double status_C, double status_T, double threshold);
 
inline std::vector< double > calcOneScore_SurvPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold,
						    arma::rowvec survTimeC, arma::rowvec survTimeT,
						    const arma::mat& survJumpC, const arma::mat& survJumpT, double lastSurvC, double lastSurvT,
						    arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
						    int p_C, int p_T, int returnIID);


inline std::vector< double > calcOnePair_CRPeron(double endpoint_C, double endpoint_T,												 
						 double status_C, double status_T, double tau,
						 arma::rowvec cifTimeC_vec,  arma::rowvec cifTimeT_vec, const arma::mat& cifJumpC,
						 double lastCif1C, double lastCif1T, double lastCif2C, double lastCif2T);

std::vector<double> calcIntegralSurv_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
					 bool returnDeriv, int column, arma::mat& derivSurv, arma::mat& derivSurvD);

double calcIntegralCif_cpp(const arma::mat& cif, double start_val, double stop_val, double CIF_t,
			   double lastCIF, int type);

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

// * calcOnePair_TTEgehan
// author Brice Ozenne
inline std::vector< double > calcOnePair_TTEgehan(double diff, double status_C, double status_T, double threshold){
  
  // ** initialize
  std::vector< double > score(4,0.0);
  // Rcpp::Rcout << diff << " " << status_T << " " << status_C << " " << threshold << std::endl;

  // ** score
  if(status_T==1){
    
    if(status_C==1){ // (treatment event, control event)
      
      if(diff >= threshold){         // >= tau    : favorable
        score[0] = 1.0;
      }else if(diff <= -threshold){ // <= -tau    : unfavorable
	score[1] = 1.0;
      }else{                        // ]-tau;tau[ : neutral
	score[2] = 1.0;
      }
      
    }else if(status_C==0){ // (treatment event, control censored)
	
      if(diff <= -threshold){ // <= -tau   : unfavorable
	score[1] = 1.0;
      }else{                  // otherwise : uninformative
	score[3] = 1.0;
      }
      
    }else if(status_C==2){ // (treatment event, control competing risk)
      score[1] = 1.0; //  unfavorable
    }
    
  }else if(status_T==0){
    
    if(status_C==1){ // (treatment censored, control event)
    
      if(diff >= threshold){ // > tau    : favorable
	score[0] = 1.0;
      }else{                 // otherwise: uninformative
	score[3] = 1.0;
      }
    
    }else{ // (treatment censored, control censored/competing risk): uninformative
      score[3] = 1.0;
    }
    
  }else if(status_T==2){ 

    if(status_C==1){ // (treatment competing risk, control event): favorable
      score[0] = 1.0;
    }else if(status_C==2){ // (treatment competing risk, control competing risk): neutral
      score[2] = 1.0;
    }else if(status_C==0){ // (treatment competing risk, control censored): uninformative
      score[3] = 1.0;
    }
    
  }

  // ** export
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}

// * calcOnePair_TTEgehan2
// author Brice Ozenne
inline std::vector< double > calcOnePair_TTEgehan2(double diff, double status_C, double status_T, double threshold){
  
  // ** initialize
  std::vector< double > score(4,0.0);
  // Rcpp::Rcout << diff << " " << status_T << " " << status_C << " " << threshold << std::endl;

  // ** score
  if(status_T==1){
    
    if(status_C==1){ // (treatment event, control event)
      
      if(diff >= threshold){         // >= tau    : favorable
        score[0] = 1.0;
      }else if(diff <= -threshold){ // <= -tau    : unfavorable
	score[1] = 1.0;
      }else{                        // ]-tau;tau[ : neutral
	score[2] = 1.0;
      }
      
    }else if(status_C==0){ // (treatment event, control censored)
	
      if(diff >= threshold){ // <= -tau   : unfavorable
	score[0] = 1.0;
      }else{                  // otherwise : uninformative
	score[3] = 1.0;
      }
      
    }else if(status_C==2){ // (treatment event, control competing risk)
      score[1] = 1.0; //  unfavorable
    }
    
  }else if(status_T==0){
    
    if(status_C==1){ // (treatment censored, control event)
    
      if(diff <= -threshold){ // > tau    : favorable
	score[1] = 1.0;
      }else{                 // otherwise: uninformative
	score[3] = 1.0;
      }
    
    }else{ // (treatment censored, control censored/competing risk): uninformative
      score[3] = 1.0;
    }
    
  }else if(status_T==2){ 

    if(status_C==1){ // (treatment competing risk, control event): favorable
      score[0] = 1.0;
    }else if(status_C==2){ // (treatment competing risk, control competing risk): neutral
      score[2] = 1.0;
    }else if(status_C==0){ // (treatment competing risk, control censored): uninformative
      score[3] = 1.0;
    }
    
  }

  // ** export
  // Rcpp::Rcout << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << std::endl;
  return(score);
  
}

// * calcOneScore_SurvPeron
// author Brice Ozenne
inline std::vector< double > calcOneScore_SurvPeron(double endpoint_C, double endpoint_T, double status_C, double status_T, double threshold,
													arma::rowvec survTimeC, arma::rowvec survTimeT,
													const arma::mat& survJumpC, const arma::mat& survJumpT, double lastSurvC, double lastSurvT,
													arma::mat& Dscore_Dnuisance_C, arma::mat& Dscore_Dnuisance_T,
													int p_C, int p_T, int returnIID){
  
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
  	Dscore_Dnuisance_C.fill(0.0); // initialized to 0 (since it is a sparse matrix)
  	Dscore_Dnuisance_T.fill(0.0); // initialized to 0 (since it is a sparse matrix)
  }
  
  // ** deal with null survival
  // according to the survival the observation will die immediatly after the observation time.
  // so we treat it as if was an event
  if(status_C==0 && (survTimeC(2) == 0) ){
    status_C = 1;
  }
  if(status_T==0 && (survTimeT(5) == 0) ){
    status_T = 1;
  }
  	
  // Rcpp::Rcout << " (" << status_T << ";" << status_C << ")";
  // ** compute favorable and unfavorable
  if(status_T==1){
    if(status_C==1){
      
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
	
    }else{ // statusC[iter_C]==0

      // favorable
      if(diff >= threshold){
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
		// score[0] = 0.0;
		upperFavorable = score[0];
      }

      // unfavorable
      if(diff <= -threshold){ 
		score[1] = 1.0;
		upperUnfavorable = score[1];
      }else {
		// Rcpp::Rcout << "(2-b) " << "";
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
	
  }else{ // statusT[iter_T]==0

    if(status_C==1){ 
      
      // favorable
      if(diff >= threshold){ // 
		score[0] = 1.0;
		upperFavorable = score[0];
      }else {
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
      }

      // unfavorable
      if(diff <= -threshold){
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
		// score[1] = 0.0;
		upperUnfavorable = score[1];
      }
      
    }else{ // status_C==0

      double denom = survTimeT(5)*survTimeC(2);
      std::vector< double > intFavorable; 
      std::vector< double > intUnfavorable;
   	  arma::mat intDscore_Dnuisance_C;
	  arma::mat intDscore_Dnuisance_T;
	  if(returnIID>1){
   	   intDscore_Dnuisance_C.resize(p_C,4); // initialized in calcIntegralSurv_cpp
	   intDscore_Dnuisance_T.resize(p_T,4); // initialized in calcIntegralSurv_cpp
	  }
	  
	  // favorable
	  if(diff >= threshold){
		// Rcpp::Rcout << "(4+a) ";
		intFavorable = calcIntegralSurv_cpp(survJumpC, endpoint_T-threshold, lastSurvT, lastSurvC,
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
		// Rcpp::Rcout << "end " << std::endl;

	  }else{
		// Rcpp::Rcout << "(4+b) ";

		intFavorable = calcIntegralSurv_cpp(survJumpC, endpoint_C, lastSurvT, lastSurvC,
											(returnIID > 1), 0, intDscore_Dnuisance_T, intDscore_Dnuisance_C); // -intFavorable is already the lower bound
		score[0] = -intFavorable[0] / denom; // (lower bound)
		upperFavorable = -intFavorable[1] / denom; // (upper bound)
		if(returnIID>1){ //		if((returnIID>1) && (intFavorable[0]>0)){
		  Dscore_Dnuisance_C(survTimeC(8),0) += intFavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
		  Dscore_Dnuisance_T(survTimeT(11),0) += intFavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
		  Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
		  Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
		}
		// Rcpp::Rcout << "end" << std::endl;
	  }
      
	  // unfavorable
	  if(diff <= -threshold){	
		// Rcpp::Rcout << "(4-a) ";

		intUnfavorable = calcIntegralSurv_cpp(survJumpT, endpoint_C-threshold, lastSurvC, lastSurvT,
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
		
		// Rcpp::Rcout << "end ";

	  }else{
	    // Rcpp::Rcout << "(4-b) ";
	    intUnfavorable = calcIntegralSurv_cpp(survJumpT, endpoint_T, lastSurvC, lastSurvT,
						  (returnIID > 1), 1, intDscore_Dnuisance_C, intDscore_Dnuisance_T); // -intUnfavorable is already the lower bound
		
	    score[1]= -intUnfavorable[0] / denom; // (lower bound)
	    upperUnfavorable = -intUnfavorable[1] / denom;  // (upper bound)
	    if(returnIID>1){ //		if((returnIID>1) && (intUnfavorable[0]>0)){
	      Dscore_Dnuisance_C(survTimeC(8),1) += intUnfavorable[0] / (denom * survTimeC(2)); // derivative regarding Sc(y_j)
	      Dscore_Dnuisance_T(survTimeT(11),1) += intUnfavorable[0] / (denom * survTimeT(5)); // derivative regarding St(x_i)
	      Dscore_Dnuisance_C -= intDscore_Dnuisance_C/denom;
	      Dscore_Dnuisance_T -= intDscore_Dnuisance_T/denom;
	    }
	    // Rcpp::Rcout << "end ";
	  }
	  // Rcpp::Rcout << std::endl;
    }}

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
inline std::vector< double > calcOnePair_CRPeron(double endpoint_C, double endpoint_T,
						 double status_C, double status_T, double tau,
						 arma::rowvec cifTimeC_vec,  arma::rowvec cifTimeT_vec, const arma::mat& cifJumpC,
						 double lastCif1C, double lastCif1T, double lastCif2C, double lastCif2T) {

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
  double Cif1T_t = cifTimeT_vec(5);
  std::vector< double > proba(5, 0.0); // [0] favorable, [1] unfavorable, [2] neutral competing, [3] neutral event,
  // [4] uninformative
  std::vector< double > proba2(4, 0.0); // [0] favorable, [1] unfavorable, [2] neutral, [3] uninformative
  double denomC = 1 - cifTimeC_vec(2) - cifTimeC_vec(7);
  double denomT = 1 - cifTimeT_vec(5) - cifTimeT_vec(7);

  if(status_T == 2) {
    if(status_C == 2) { // (2,2)
      proba[2] = 1.0; // systematically neutral competing
    } else if(status_C == 1){ // (2,1)
      proba[0] = 1.0; // systematically favorable
    } else if(status_C == 0) { // (2,0)
      proba[0] = (lastCif1C - cifTimeC_vec(2))/denomC;
      // proba[1] = 0
      proba[2] = (lastCif2C - cifTimeC_vec(7))/denomC;
      // proba[3] = 0 // since the treated patient had the competing event
      //proba[4] = 1 - (proba[0] + proba[2]);
    }
  } else if(status_T == 1){
    if(status_C == 2){ // (1,2)
      proba[1] = 1.0; // systematically defavorable
    } else if(status_C == 1){ // (1,1)
      if(diff >= tau) {
        proba[0] = 1.0;
      } else if(diff <= -tau) {
        proba[1] = 1.0;
      } else { // |diff| < tau
        proba[3] = 1.0;
      }
    } else if(status_C == 0) { // (1,0)
      if(diff >= tau) {
        if(R_IsNA(cifTimeT_vec(1)) == false) {
          proba[0] = (cifTimeT_vec(1) - cifTimeC_vec(2))/denomC;
        } else {
          proba[0] = (lastCif1C - cifTimeC_vec(2))/denomC;
        }
        if(R_IsNA(cifTimeT_vec(3)) == false) {
          proba[1] = (lastCif1C - cifTimeT_vec(3) + lastCif2C - cifTimeC_vec(7))/denomC;
        } else {
          proba[1] = (lastCif2C - cifTimeC_vec(7))/denomC;
        }
        if((R_IsNA(cifTimeT_vec(3)) == false) & (R_IsNA(cifTimeT_vec(1)) == false)) {
          proba[3] = (cifTimeT_vec(3) - cifTimeT_vec(1))/denomC;
        } else if ((R_IsNA(cifTimeT_vec(3)) == true) & (R_IsNA(cifTimeT_vec(1)) == false)) {
          proba[3] = (lastCif1C - cifTimeT_vec(1))/denomC;
        } else {
          proba[3] = 0.0;
        }
      } else if(diff <= -tau) {
        proba[1] = 1.0;
      } else { // |diff| < tau
        if(R_IsNA(cifTimeT_vec(3)) == false) {
          proba[1] = (lastCif1C - cifTimeT_vec(3) + lastCif2C - cifTimeC_vec(7))/denomC;
          proba[3] = (cifTimeT_vec(3) - cifTimeC_vec(2))/denomC;
        } else {
          proba[1] = (lastCif2C - cifTimeC_vec(7))/denomC;
          proba[3] = (lastCif1C - cifTimeC_vec(2))/denomC;
        }
      }
      //proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    }
  } else { // status_T == 0
    if(status_C == 2) { // (0,2)
      proba[1] = (lastCif1T - cifTimeT_vec(5))/denomT;
      proba[2] = (lastCif2T - cifTimeT_vec(7))/denomT;
      //proba[4] = 1 - (proba[1] + proba[2]);
    } else if(status_C == 1) { // (0,1)
      if(diff >= tau) {
        proba[0] = 1.0;
      } else if(diff <= -tau) {
        if(R_IsNA(cifTimeC_vec(6)) == false) {
          proba[0] = (lastCif1T - cifTimeC_vec(6) + lastCif2T - cifTimeT_vec(7))/denomT;
        } else {
          proba[0] = (lastCif2T - cifTimeT_vec(7))/denomT;
        }
        if(R_IsNA(cifTimeC_vec(4)) == false) {
          proba[1] = (cifTimeC_vec(4) - cifTimeT_vec(5))/denomT;
        } else {
          proba[1] = (lastCif1T - cifTimeT_vec(5))/denomT;
        }
        if((R_IsNA(cifTimeC_vec(6)) == false) & (R_IsNA(cifTimeC_vec(4)) == false)) {
          proba[3] = (cifTimeC_vec(6) - cifTimeC_vec(4))/denomT;
        } else if ((R_IsNA(cifTimeC_vec(6)) == true) & (R_IsNA(cifTimeC_vec(4)) == false)) {
          proba[3] = (lastCif1T - cifTimeC_vec(4))/denomT;
        } else {
          proba[3] = 0.0;
        }
      } else { // |diff| < tau
        if(R_IsNA(cifTimeC_vec(6)) == false) {
          proba[0] = (lastCif1T - cifTimeC_vec(6) + lastCif2T - cifTimeT_vec(7))/denomT;
          proba[3] = (cifTimeC_vec(6) - cifTimeT_vec(5))/denomT;
        } else {
          proba[0] = (lastCif2T - cifTimeT_vec(7))/denomT;
          proba[3] = (lastCif1T - cifTimeT_vec(5))/denomT;
        }
      }
      //proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    } else if (status_C == 0) { // (0,0)
      double prob21 = (lastCif2T - cifTimeT_vec(7))*(lastCif1C - cifTimeC_vec(2))/(denomT*denomC);
      double prob12 = (lastCif2C - cifTimeC_vec(7))*(lastCif1T - cifTimeT_vec(5))/(denomT*denomC);
      double prob22 = (lastCif2C - cifTimeC_vec(7))*(lastCif2T - cifTimeT_vec(7))/(denomT*denomC);
      if(diff >= tau) {
        // lower bound of each integral
        double intFav = calcIntegralCif_cpp(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent1 = calcIntegralCif_cpp(cifJumpC, endpoint_T - tau, endpoint_T + tau, Cif1T_t, lastCif1T, 3)/(denomT*denomC);
        double intNeutralEvent2 = calcIntegralCif_cpp(cifJumpC, endpoint_T + tau, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        if (R_IsNA(cifTimeT_vec(1)) == false) {
          proba[0] = ((cifTimeT_vec(1) - cifTimeC_vec(2))/denomC)*((lastCif1T - cifTimeT_vec(5))/denomT) +
            intFav + prob21;
        } else {
          proba[0] = ((lastCif1C - cifTimeC_vec(2))/denomC)*((lastCif1T - cifTimeT_vec(5))/denomT) +
            intFav + prob21;
        }
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent1 + intNeutralEvent2;
      } else if(diff <= -tau) {
        double intFav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        proba[0] = intFav + prob21;
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent;
      } else { // |diff| < tau
        double intFav = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 1)/(denomT*denomC);
        double intDefav = calcIntegralCif_cpp(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T_t, lastCif1T, 2)/(denomT*denomC);
        double intNeutralEvent1 = calcIntegralCif_cpp(cifJumpC, endpoint_C, endpoint_T + tau, Cif1T_t, lastCif1T, 3)/(denomT*denomC);
        double intNeutralEvent2 = calcIntegralCif_cpp(cifJumpC, endpoint_T+tau, endpoint_T + tau, Cif1T_t, lastCif1T, 4)/(denomT*denomC);
        proba[0] = intFav + prob21;
        proba[1] = intDefav + prob12;
        proba[2] = prob22;
        proba[3] = intNeutralEvent1 + intNeutralEvent2;
      }
      //proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
    }
  }
  
  proba[4] = 1 - (proba[0] + proba[1] + proba[2] + proba[3]);
  
  // merge the two kinds of neutral pairs
  proba2[0] = proba[0];
  proba2[1] = proba[1];
  proba2[2] = proba[2] + proba[3];
  proba2[3] = proba[4];
  
  return proba2;
}


// * calcIntegralSurv_cpp
//' @title C++ Function Computing the Integral Terms for the Peron Method in the survival case. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored.
//' @name calcIntegralSurv_cpp
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
//' @author Brice Ozenne
//' @export
// [[Rcpp::export]]
std::vector< double > calcIntegralSurv_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv,
					   bool returnDeriv, int column, arma::mat& derivSurv, arma::mat& derivSurvD){

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
	  derivSurv(survival(iter_time,3), column) += survival(iter_time,2); // derivative regarding S(t+tau)
	  derivSurvD(survival(iter_time,4), column) -= survival(iter_time,1); // derivative regarding dS(t-)
	  derivSurvD(survival(iter_time,5), column) += survival(iter_time,1); // derivative regarding dS(t+)
	}
      }
    }
    iter_time++;
  }
    
  // Rcpp::Rcout << std::endl << "lastSurv: " << lastSurv << " " << lastdSurv << " | current value " << integral[1] << std::endl;

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
//' @name calcIntegralCif_cpp
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
//' @author Eva Cantagallo
//' @export
// [[Rcpp::export]]
double calcIntegralCif_cpp(const arma::mat& cif, double start_val, double stop_val, double CIF_t,
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
