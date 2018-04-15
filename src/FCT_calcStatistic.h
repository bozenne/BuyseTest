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

// * calcStatistic
void calcStatistic(arma::mat& delta_netChance, arma::mat& delta_winRatio, vector<double>& Delta_netChance, vector<double>& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   const int& D, const int& n_strata, const double& n_pairs){
  
  double strata_favorable=0, strata_unfavorable=0;
  for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints
    
    if(iter_d==0){
      Delta_netChance[iter_d] = 0;
    }else{
      Delta_netChance[iter_d] = Delta_netChance[iter_d-1];
    }
    
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      delta_netChance(iter_strata,iter_d) = (Mcount_favorable(iter_strata,iter_d)-Mcount_unfavorable(iter_strata,iter_d))/(double)(n_pairs); // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the total number of pairs
      Delta_netChance[iter_d] += delta_netChance(iter_strata,iter_d);
      
      delta_winRatio(iter_strata,iter_d) = Mcount_favorable(iter_strata,iter_d)/(double)(Mcount_unfavorable(iter_strata,iter_d)); // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
      strata_favorable += Mcount_favorable(iter_strata,iter_d);
      strata_unfavorable += Mcount_unfavorable(iter_strata,iter_d);
    }
    
    Delta_winRatio[iter_d] = strata_favorable/(double)(strata_unfavorable);
  }
  
  return ;
}
