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

// * calcStatistic
void calcStatistic(arma::mat& delta_netBenefit, arma::mat& delta_winRatio, vector<double>& Delta_netBenefit, vector<double>& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   const int& D, const int& n_strata, const vector<double>& n_pairs){
  

  // total number of pairs
  double ntot_pair=0;
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	ntot_pair += n_pairs[iter_strata];
  }

  //
  double strata_favorable=0, strata_unfavorable=0;
  for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints
    
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	   // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the number of pair in the strata
      delta_netBenefit(iter_strata,iter_d) = (Mcount_favorable(iter_strata,iter_d)-Mcount_unfavorable(iter_strata,iter_d))/(double)(n_pairs[iter_strata]);
	  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
      delta_winRatio(iter_strata,iter_d) = Mcount_favorable(iter_strata,iter_d)/(double)(Mcount_unfavorable(iter_strata,iter_d));

	  // accumulate number of favorable and unfavorable
      strata_favorable += Mcount_favorable(iter_strata,iter_d);
      strata_unfavorable += Mcount_unfavorable(iter_strata,iter_d);
    }
    
    Delta_winRatio[iter_d] = strata_favorable/(double)(strata_unfavorable);
    Delta_netBenefit[iter_d] = (strata_favorable-strata_unfavorable)/(double)(ntot_pair);
  }
  
  return ;
}
