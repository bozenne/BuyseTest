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
void calcStatistic(arma::mat& delta_netBenefit, arma::mat& delta_winRatio, std::vector< double >& Delta_netBenefit, std::vector< double >& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   arma::mat& iid_favorable, arma::mat& iid_unfavorable, arma::mat& Mvar, bool returnIID,
				   std::vector< arma::uvec >& indexC, std::vector< arma::uvec >& indexT,
                   const unsigned int& D, const int& n_strata, const std::vector< double >& n_pairs,
		           const std::vector< double >& weight){
  
  
  // total number of pairs and patients in each arm
  double ntot_pair=0;
  double ntot_treatment=0;
  double ntot_control=0;
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	ntot_pair += n_pairs[iter_strata];
	ntot_control += indexC[iter_strata].size();
	ntot_treatment += indexT[iter_strata].size();
  }

  // global iid decomposition
  arma::colvec iCumiid_favorable;
  arma::colvec iCumiid_unfavorable;
  arma::mat cumiid_favorable;
  arma::mat cumiid_unfavorable;
  if(returnIID){
	int nObs = iid_favorable.n_rows;
	iCumiid_favorable.resize(nObs);
	iCumiid_favorable.fill(0.0);
	iCumiid_unfavorable.resize(nObs);
	iCumiid_unfavorable.fill(0.0);

	cumiid_favorable.resize(nObs, D);
	cumiid_favorable.fill(0.0);
	cumiid_unfavorable.resize(nObs, D);
	cumiid_unfavorable.fill(0.0);
  }
  
  //
  double iStrata_favorable, iStrata_unfavorable;
  double iFavorable = 0;
  double iUnfavorable = 0;
  uvec iUvec_iter_d;
  
  for(unsigned int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints

    iStrata_favorable=0;
    iStrata_unfavorable=0;
    iUvec_iter_d = {iter_d};
	
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the number of pair in the strata
      delta_netBenefit(iter_strata,iter_d) = (Mcount_favorable(iter_strata,iter_d)-Mcount_unfavorable(iter_strata,iter_d))/(double)(n_pairs[iter_strata]);
      // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
      delta_winRatio(iter_strata,iter_d) = Mcount_favorable(iter_strata,iter_d)/(double)(Mcount_unfavorable(iter_strata,iter_d));
      // accumulate number of favorable and unfavorable
      iStrata_favorable += Mcount_favorable(iter_strata,iter_d);
      iStrata_unfavorable += Mcount_unfavorable(iter_strata,iter_d);
    }


	// normalize iid
	if(returnIID){
	  iid_favorable.col(iter_d) -= iStrata_favorable / (double)(ntot_pair); // center
	  iid_unfavorable.col(iter_d) -= iStrata_unfavorable / (double)(ntot_pair); // center
	  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // divide by the number of terms i.e. we output \psi_i/n in 1/n \sum_i \psi_i
		iid_favorable.submat(indexC[iter_strata], iUvec_iter_d) /= (double) ntot_control;
		iid_unfavorable.submat(indexC[iter_strata], iUvec_iter_d) /= (double) ntot_control;

		iid_favorable.submat(indexT[iter_strata], iUvec_iter_d) /= (double) ntot_treatment;
		iid_unfavorable.submat(indexT[iter_strata], iUvec_iter_d) /= (double) ntot_treatment;
	}
	}
	
	// cumulate statistic over endpoints
	iFavorable += weight[iter_d] * iStrata_favorable;
	iUnfavorable += weight[iter_d] * iStrata_unfavorable;
  
	Delta_winRatio[iter_d] = iFavorable/(double)(iUnfavorable);
	Delta_netBenefit[iter_d] = (iFavorable-iUnfavorable)/(double)(ntot_pair);

	// cumulate iid over endpoints
	if(returnIID){
	  if(iter_d==0){
		cumiid_favorable.col(iter_d) = weight[iter_d] * iid_favorable.col(iter_d);
		cumiid_unfavorable.col(iter_d) = weight[iter_d] * iid_unfavorable.col(iter_d);
	  }else{
		cumiid_favorable.col(iter_d) = weight[iter_d] * iid_favorable.col(iter_d) + cumiid_favorable.col(iter_d-1);
		cumiid_unfavorable.col(iter_d) = weight[iter_d] * iid_unfavorable.col(iter_d) + cumiid_unfavorable.col(iter_d-1);
	  }
	}
    
  }

  // compute variance
  if(returnIID){
	Mvar.col(0) = trans(sum(pow(cumiid_favorable,2), 0));
	Mvar.col(1) = trans(sum(pow(cumiid_unfavorable,2), 0));
	Mvar.col(2) = trans(sum(cumiid_favorable % cumiid_unfavorable, 0));
  }
  
  return ;
}
