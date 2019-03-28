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
void calcStatistic(arma::mat& delta_netBenefit, arma::mat& delta_winRatio, arma::vec & Delta_netBenefit, arma::vec& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   arma::mat& iid_favorable, arma::mat& iid_unfavorable, arma::mat& Mvar, bool returnIID,
				   std::vector< arma::uvec >& posC, std::vector< arma::uvec >& posT,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs,
		           const arma::vec& weight, int hprojection){
  
  // ** total number of pairs and patients in each arm
  double ntot_pair = 0;
  double ntot_treatment = 0;
  double ntot_control = 0;
  for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	ntot_pair += n_pairs[iter_strata];
	ntot_control += posC[iter_strata].size();
	ntot_treatment += posT[iter_strata].size();
  }
    
  // ** net benefit and win ratio
  arma::vec cumWcount_favorable(D);
  arma::vec cumWcount_unfavorable(D);

  // partial sum
  cumWcount_favorable = conv_to<vec>::from(sum(Mcount_favorable,0)); // sum over strata
  cumWcount_favorable %= weight; // weight endpoints
  cumWcount_favorable = cumsum(cumWcount_favorable); // cumulate over endpoints
  
  cumWcount_unfavorable = conv_to<vec>::from(sum(Mcount_unfavorable,0)); // sum over strata
  cumWcount_unfavorable %= weight;  // weight endpoints
  cumWcount_unfavorable = cumsum(cumWcount_unfavorable); // cumulate over endpoints

  // net benefit equals (number of favorable pairs minus number of unfavorable pairs) divided by number of pairs
  Delta_netBenefit = (cumWcount_favorable - cumWcount_unfavorable)/(double)(ntot_pair);
  delta_netBenefit = (Mcount_favorable - Mcount_unfavorable);
  delta_netBenefit.each_col() /= n_pairs;

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
  Delta_winRatio = cumWcount_favorable / cumWcount_unfavorable;
  delta_winRatio = Mcount_favorable / Mcount_unfavorable;

	  
  // ** iid and variance estimation
  if(returnIID){
					   
  arma::vec meani2_favorable = zeros<vec>(D);
  arma::vec meanj2_favorable = zeros<vec>(D);
  arma::vec meani2_unfavorable = zeros<vec>(D);
  arma::vec meanj2_unfavorable = zeros<vec>(D);
  arma::vec meani_mixed = zeros<vec>(D);
  arma::vec meanj_mixed = zeros<vec>(D);

  // cumulative score
  arma::rowvec rowweight = conv_to<rowvec>::from(weight);

  arma::mat cumWiid_favorable = iid_favorable;
  cumWiid_favorable.each_row() %= rowweight;
  cumWiid_favorable = cumsum(cumWiid_favorable,1);
  
  arma::mat cumWiid_unfavorable = iid_unfavorable;
  cumWiid_unfavorable.each_row() %= rowweight;
  cumWiid_unfavorable = cumsum(cumWiid_unfavorable,1);

  // iid and sufficient statistics
  int iN_control, iN_treatment;
  arma::vec delta_favorable = cumWcount_favorable/(double)(ntot_pair);
  arma::vec delta_unfavorable = cumWcount_unfavorable/(double)(ntot_pair);

  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	iN_control = posC[iter_strata].size();
	iN_treatment = posT[iter_strata].size();
	
	iid_favorable.rows(posC[iter_strata]) = cumWiid_favorable.rows(posC[iter_strata]) / iN_treatment;
	iid_favorable.rows(posT[iter_strata]) = cumWiid_favorable.rows(posT[iter_strata]) / iN_control;
	iid_unfavorable.rows(posC[iter_strata]) = cumWiid_unfavorable.rows(posC[iter_strata]) / iN_treatment;
	iid_unfavorable.rows(posT[iter_strata]) = cumWiid_unfavorable.rows(posT[iter_strata]) / iN_control;

	meani2_favorable += conv_to<vec>::from( sum(pow(cumWiid_favorable.rows(posC[iter_strata]),2),0) / iN_treatment );
	meanj2_favorable += conv_to<vec>::from( sum(pow(cumWiid_favorable.rows(posT[iter_strata]),2),0) / iN_control );
	
	meani2_unfavorable += conv_to<vec>::from( sum(pow(cumWiid_unfavorable.rows(posC[iter_strata]),2),0) / iN_treatment );
	meanj2_unfavorable += conv_to<vec>::from( sum(pow(cumWiid_unfavorable.rows(posT[iter_strata]),2),0) / iN_control );
	
	meani_mixed += conv_to<vec>::from( sum(cumWiid_favorable.rows(posC[iter_strata]) % cumWiid_unfavorable.rows(posC[iter_strata]),0) / iN_treatment );
	meanj_mixed += conv_to<vec>::from( sum(cumWiid_favorable.rows(posT[iter_strata]) % cumWiid_unfavorable.rows(posT[iter_strata]),0) / iN_control );
  }

  // center and scale iid
  // Rcout << endl << "iid" << endl;
  iid_favorable.each_row() -= conv_to<rowvec>::from(delta_favorable);
  iid_unfavorable.each_row() -= conv_to<rowvec>::from(delta_unfavorable);
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	iid_favorable.rows(posT[iter_strata]) /= ntot_treatment;
	iid_favorable.rows(posC[iter_strata]) /= ntot_control;
	iid_unfavorable.rows(posT[iter_strata]) /= ntot_treatment;
	iid_unfavorable.rows(posC[iter_strata]) /= ntot_control;
  }
  
  // compute variance from the sufficient statistics
  // Rcout << endl << "variance" << endl;
  arma::vec delta2_favorable = pow(delta_favorable,2);
  arma::vec delta2_unfavorable = pow(delta_unfavorable,2);
  arma::vec delta2_mixed = (cumWcount_favorable % cumWcount_unfavorable)/pow(ntot_pair, 2);

  if(hprojection==1){
	Mvar.col(0) = (meani2_favorable/ntot_pair - delta2_favorable)/ntot_control + (meanj2_favorable/ntot_pair - delta2_favorable)/ntot_treatment;
	Mvar.col(1) = (meani2_unfavorable/ntot_pair - delta2_unfavorable)/ntot_control + (meanj2_unfavorable/ntot_pair - delta2_unfavorable)/ntot_treatment;
	Mvar.col(2) = (meani_mixed/ntot_pair - delta2_mixed)/ntot_control + (meanj_mixed/ntot_pair - delta2_mixed)/ntot_treatment;
  }else if(hprojection==2){
	Mvar.col(0) = (meani2_favorable/ntot_pair - delta2_favorable)*(ntot_treatment-1)/ntot_pair + (meanj2_favorable/ntot_pair - delta2_favorable)*(ntot_control-1)/ntot_pair + delta_favorable*(1-delta_favorable);
	Mvar.col(1) = (meani2_unfavorable/ntot_pair - delta2_unfavorable)*(ntot_treatment-1)/ntot_pair + (meanj2_unfavorable/ntot_pair - delta2_unfavorable)*(ntot_control-1)/ntot_pair + delta_unfavorable*(1-delta_unfavorable);
	Mvar.col(2) = (meani_mixed/ntot_pair - delta2_mixed)*(ntot_treatment-1)/ntot_pair + (meanj_mixed/ntot_pair - delta2_mixed)*(ntot_control-1)/ntot_pair + delta_mixed*(1-delta_mixed);
  }
  // Rcout << endl << "delta method" << endl;  
  // delta method
  // var(A-B) = var(A) + var(B) - 2 * cov(A,B)
  // indeed (A-B)' = A' - B' so (A-B)^'2 = A'A' + B'B'  - 2*A'B'
  Mvar.col(3) = Mvar.col(0) + Mvar.col(1) - 2 * Mvar.col(2);
  // var(A/B) = var(A)/B^2 + var(B)*(A^2/B^4) - 2*cov(A,B)A/B^3
  // indeed (A/B)' = A'/B - B'A/B^2 so (A/B)^'2 = A'A'/B^2 + B'B'A^2/B^2 - 2B'A' A/B^3
  cumWcount_favorable = cumWcount_favorable/(double)(ntot_pair);
  cumWcount_unfavorable = cumWcount_unfavorable/(double)(ntot_pair);
  Mvar.col(4) = Mvar.col(0)/pow(cumWcount_unfavorable, 2) + Mvar.col(1) % pow(cumWcount_favorable,2)/pow(cumWcount_unfavorable,4) - 2 * Mvar.col(2) % cumWcount_favorable/pow(cumWcount_unfavorable, 3);
  }
  
  return ;
}
