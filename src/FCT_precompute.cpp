// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// * calcIntegralSurv2_cpp
//' @title C++ Function Computing the Integral Terms for the Peron Method in the survival case. 
//' @description Compute the integral with respect to the jump in survival for pairs where both outcomes are censored, i.e. \int S1(t+\tau) dS2(t).
//' @name calcIntegralSurv_cpp
//' 
//' @param time [numeric vector] vector of jump time for S2.
//' @param survival [numeric vector] the survival at each jump time: S1(t+\tau).
//' @param dsurvival [numeric vector] the jump in survival at each jump time: S2(t+)-S2(t-)
//' @param index_survival [numeric vector] the position of survival parameter S1(t+\tau) among all parameters relative to S1.
//' @param index_dSurvival1 [numeric vector] the position of survival parameter S2(t-) among all parameters relative to S2.
//' @param index_dSurvival2 [numeric vector] the position of survival parameter S2(t+) among all parameters relative to S2.
//' @param lastSurv [numeric] the value of S1 at the end of the follow-up.
//' @param lastSurv [numeric] the value of S2 at the end of the follow-up.
//' @param iidNuisance [logical] should the derivative of the integral relative to the S1 and S2 parameter be output.
//' @param p_Surv [integer] the number of parameters relative to S1.
//' @param p_SurvD [integer] the number of parameters relative to S2.
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
				 int p_Surv,
				 int p_SurvD,
				 int nJump){

  std::vector<double> intSurv_lower(nJump+1,0.0);
  std::vector<double> intSurv_upper(nJump+1,0.0);
  
  std::vector<double> jumpDeriv_Surv(0);
  std::vector<double> indexDeriv_Surv(0);
  std::vector<double> valueDeriv_Surv(0);

  std::vector<double> jumpDeriv_SurvD(0);
  std::vector<double> indexDeriv_SurvD(0);
  std::vector<double> valueDeriv_SurvD(0);

  if(iidNuisance){
    jumpDeriv_Surv.reserve(nJump);
    indexDeriv_Surv.reserve(nJump);
    valueDeriv_Surv.reserve(nJump);
  
    jumpDeriv_SurvD.reserve(2*nJump);
    indexDeriv_SurvD.reserve(2*nJump);
    valueDeriv_SurvD.reserve(2*nJump);
  }
  
  // ** loop over time
  if(lastdSurv > 0){
    // add extra contribution to the bound after the last jump (minus because dSurv = (0 - surv(tmax))
    if((nJump > 0) && (survival[nJump-1]>=0) && arma::is_finite(survival[nJump-1])){ // <0 and is_finite test whether it is a missing value (NA in R)
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
	  // derivative regarding S(t+\tau)
	  jumpDeriv_Surv.push_back(iJump);
	  indexDeriv_Surv.push_back(index_survival[iJump]);
	  valueDeriv_Surv.push_back(dSurvival[iJump]);

	  // derivative regarding S(-)
	  jumpDeriv_SurvD.push_back(iJump);
	  indexDeriv_SurvD.push_back(index_dSurvival1[iJump]);
	  valueDeriv_SurvD.push_back(-survival[iJump]);

	  // derivative regarding S(+)
	  jumpDeriv_SurvD.push_back(iJump);
	  indexDeriv_SurvD.push_back(index_dSurvival2[iJump]);
	  valueDeriv_SurvD.push_back(survival[iJump]);
	}
      }else{
	intSurv_upper[iJump] += lastSurv * dSurvival[iJump];
      }
    }
  }

  // ** export
  int n_derivSurv = jumpDeriv_Surv.size();
  arma::mat intSurv_derivSurv(n_derivSurv,3);
  for(int iter=0; iter<n_derivSurv; iter++){
    if(jumpDeriv_Surv[iter] == nJump){
      if(nJump==0){
	intSurv_derivSurv(n_derivSurv-iter-1,0) = NA_REAL;
      }else{
	intSurv_derivSurv(n_derivSurv-iter-1,0) = time[jumpDeriv_Surv[iter]-1]+1e-12;
      }
    }else{
      intSurv_derivSurv(n_derivSurv-iter-1,0) = time[jumpDeriv_Surv[iter]];
    }
    intSurv_derivSurv(n_derivSurv-iter-1,1) = indexDeriv_Surv[iter];
    intSurv_derivSurv(n_derivSurv-iter-1,2) = valueDeriv_Surv[iter];
  }

  int n_derivSurvD = indexDeriv_SurvD.size();
  arma::mat intSurv_derivSurvD(n_derivSurvD,3);
  for(int iter=0; iter<n_derivSurvD; iter++){
    if(jumpDeriv_SurvD[iter] == nJump){
      if(nJump==0){
	intSurv_derivSurvD(n_derivSurvD-iter-1,0) = NA_REAL;
      }else{
	intSurv_derivSurvD(n_derivSurvD-iter-1,0) = time[jumpDeriv_SurvD[iter]-1]+1e-12;
      }
    }else{
      intSurv_derivSurvD(n_derivSurvD-iter-1,0) = time[jumpDeriv_SurvD[iter]];
    }
    intSurv_derivSurvD(n_derivSurvD-iter-1,1) = indexDeriv_SurvD[iter];
    intSurv_derivSurvD(n_derivSurvD-iter-1,2) = valueDeriv_SurvD[iter];
  }

  
  return(Rcpp::List::create(Rcpp::Named("time") = time,
			    Rcpp::Named("intSurv_lower") = intSurv_lower,
			    Rcpp::Named("intSurv_upper") = intSurv_upper,
			    Rcpp::Named("intSurv_derivSurv") = intSurv_derivSurv,
			    Rcpp::Named("intSurv_derivSurvD") = intSurv_derivSurvD
			    ));
}
