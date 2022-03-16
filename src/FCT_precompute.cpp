// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// * calcIntegralSurv2_cpp
//' @title C++ Function pre-computing the Integral Terms for the Peron Method in the survival case. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored, i.e. \eqn{\int S1(t+\tau) dS2(t)}.
//' @name calcIntegralSurv2_cpp
//' 
//' @param time [numeric vector] vector of jump time for S2.
//' @param survival [numeric vector] the survival at each jump time: \eqn{S1(t+\tau)}.
//' @param dsurvival [numeric vector] the jump in survival at each jump time: \eqn{S2(t+)-S2(t-)}
//' @param index_survival [numeric vector] the position of survival parameter \eqn{S1(t+\tau)} among all parameters relative to S1.
//' @param index_dSurvival1 [numeric vector] the position of survival parameter \eqn{S2(t-)} among all parameters relative to S2.
//' @param index_dSurvival2 [numeric vector] the position of survival parameter \eqn{S2(t+)} among all parameters relative to S2.
//' @param lastSurv [numeric] the value of S1 at the end of the follow-up.
//' @param lastSurv [numeric] the value of S2 at the end of the follow-up.
//' @param iidNuisance [logical] should the derivative of the integral relative to the S1 and S2 parameter be output.
//' @param nJump [integer] the number of jump times relative to S2.
//'
//' @keywords function Cpp internal
//' @author Brice Ozenne
//' @export
// [[Rcpp::export]]
Rcpp::List calcIntegralSurv2_cpp(const std::vector<double>& time,
				 const std::vector<double>& survival,
				 const std::vector<double>& dSurvival,
				 const std::vector<int>& index_survival,
				 const std::vector<int>& index_dSurvival1,
				 const std::vector<int>& index_dSurvival2,
				 double lastSurv,
				 double lastdSurv,
				 bool iidNuisance,
				 int nJump){


  std::vector<double> intSurv_lower(nJump+1,0.0);
  std::vector<double> intSurv_upper(nJump+1,0.0);

  arma::mat intSurv_deriv;
  if(iidNuisance){
    intSurv_deriv.resize(nJump+1,8);
    intSurv_deriv.fill(0.0);
    intSurv_deriv(nJump,0) = nJump+1;
    if(nJump>0){
      intSurv_deriv(nJump,1) = time[nJump-1]+1e-12;
    }
  }

  if(nJump == 0){
	
	return(Rcpp::List::create(Rcpp::Named("time") = 0,
							  Rcpp::Named("intSurv_lower") = 0,
							  Rcpp::Named("intSurv_upper") = 0,
							  Rcpp::Named("intSurv_deriv") = intSurv_deriv
							  ));
	
  }

  // ** loop over time
  if(lastdSurv > 0){
    // add extra contribution to the bound after the last jump (minus because dSurv = (0 - surv(tmax))
    if(nJump<0 || abs(nJump) > intSurv_upper.size()-1){ // keep abs() for C++ compiler warning: int versus unsigned int
      throw std::runtime_error("Incorrect nJump value (negative or too large)");
    }else if(nJump==0){ // should never be used as by default the time 0 is included so nJump is at least 1
      intSurv_upper[0] = - lastSurv * lastdSurv;
    }else if((survival[nJump-1]>=0) && arma::is_finite(survival[nJump-1])){ // <0 and is_finite test whether it is a missing value (NA in R)
      intSurv_upper[nJump] = - survival[nJump-1] * lastdSurv;
    }else{
      intSurv_upper[nJump] = - lastSurv * lastdSurv;
    }
  }


  if(nJump>0){
    for(int iJump=nJump-1; iJump>=0; iJump--){
      intSurv_lower[iJump] = intSurv_lower[iJump+1];
      intSurv_upper[iJump] = intSurv_upper[iJump+1];
        
      if(survival[iJump]>=0 && arma::is_finite(survival[iJump])){ // <0 and is_finite test whether it is a missing value (NA in R)
		intSurv_lower[iJump] += survival[iJump] * dSurvival[iJump];
		intSurv_upper[iJump] += survival[iJump] * dSurvival[iJump];

		if(iidNuisance){
		  intSurv_deriv(iJump,0) = iJump;
		  intSurv_deriv(iJump,1) = time[iJump];

		  // derivative regarding S(t+\tau)
		  intSurv_deriv(iJump,2) = index_survival[iJump];
		  intSurv_deriv(iJump,3) = dSurvival[iJump];

		  // derivative regarding S(-)
		  intSurv_deriv(iJump,4) = index_dSurvival1[iJump];
		  intSurv_deriv(iJump,5) = -survival[iJump];

		  // derivative regarding S(+)
		  intSurv_deriv(iJump,6) = index_dSurvival2[iJump];
		  intSurv_deriv(iJump,7) = survival[iJump];
		}
      }else{
		intSurv_upper[iJump] += lastSurv * dSurvival[iJump];
      }
    }
  }

  // ** export
  return(Rcpp::List::create(Rcpp::Named("time") = time,
			    Rcpp::Named("intSurv_lower") = intSurv_lower,
			    Rcpp::Named("intSurv_upper") = intSurv_upper,
			    Rcpp::Named("intSurv_deriv") = intSurv_deriv
			    ));
}

