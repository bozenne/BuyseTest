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
                   arma::mat& iid_favorable, arma::mat& iid_unfavorable, const arma::mat& iidNuisance_favorable, const arma::mat& iidNuisance_unfavorable,
		   arma::mat& Mvar, int returnIID,
		   std::vector< arma::uvec >& posC, std::vector< arma::uvec >& posT,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weight, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore){
  
  // ** total number of pairs and patients in each arm
  double ntot_pair = 0;
  for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    ntot_pair += n_pairs[iter_strata];
  }
  double ntot_control = sum(n_control);
  double ntot_treatment = sum(n_treatment);
    
  // ** net benefit and win ratio
  arma::vec count_favorable(D);
  arma::vec count_unfavorable(D);
  arma::vec cumWcount_favorable(D);
  arma::vec cumWcount_unfavorable(D);

  // sum over strata
  count_favorable = conv_to<vec>::from(sum(Mcount_favorable,0));
  count_unfavorable = conv_to<vec>::from(sum(Mcount_unfavorable,0)); 

  // weight endpoints and cumulate over endpoints
  cumWcount_favorable = count_favorable;
  cumWcount_favorable %= weight; 
  cumWcount_favorable = cumsum(cumWcount_favorable); 

  cumWcount_unfavorable = count_unfavorable;
  cumWcount_unfavorable %= weight;
  cumWcount_unfavorable = cumsum(cumWcount_unfavorable); 

  // net benefit equals (number of favorable pairs minus number of unfavorable pairs) divided by number of pairs
  Delta_netBenefit = (cumWcount_favorable - cumWcount_unfavorable)/(double)(ntot_pair);
  delta_netBenefit = (Mcount_favorable - Mcount_unfavorable);
  delta_netBenefit.each_col() /= n_pairs;

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
  Delta_winRatio = cumWcount_favorable / cumWcount_unfavorable;
  delta_winRatio = Mcount_favorable / Mcount_unfavorable;

	  
  // ** iid and variance estimation
  if(returnIID > 0){
    arma::vec delta_favorable = count_favorable/(double)(ntot_pair);
    arma::vec delta_unfavorable = count_unfavorable/(double)(ntot_pair);
    arma::vec cumWdelta_favorable = cumWcount_favorable/(double)(ntot_pair);
    arma::vec cumWdelta_unfavorable = cumWcount_unfavorable/(double)(ntot_pair);
					   
    // *** compute expectation from sum
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      iid_favorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iid_favorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iid_unfavorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iid_unfavorable.rows(posT[iter_strata]) /= n_control[iter_strata];
    }
	
    // *** center
    iid_favorable.each_row() -= conv_to<rowvec>::from(delta_favorable);
    iid_unfavorable.each_row() -= conv_to<rowvec>::from(delta_unfavorable);

    // *** rescale
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
      iid_favorable.rows(posC[iter_strata]) /= ntot_control;
      iid_favorable.rows(posT[iter_strata]) /= ntot_treatment;
      iid_unfavorable.rows(posC[iter_strata]) /= ntot_control;
      iid_unfavorable.rows(posT[iter_strata]) /= ntot_treatment;
    }
	
    // *** weight endpoints and cumulate them
    arma::rowvec rowweight = conv_to<rowvec>::from(weight);

    iid_favorable.each_row() %= rowweight;
    iid_favorable = cumsum(iid_favorable,1);
  
    iid_unfavorable.each_row() %= rowweight;
    iid_unfavorable = cumsum(iid_unfavorable,1);

    // *** sufficient statistics
    arma::vec sigmaC_favorable = zeros<vec>(D);
    arma::vec sigmaT_favorable = zeros<vec>(D);
    arma::vec sigmaC_unfavorable = zeros<vec>(D);
    arma::vec sigmaT_unfavorable = zeros<vec>(D);
    arma::vec sigmaC_mixed = zeros<vec>(D);
    arma::vec sigmaT_mixed = zeros<vec>(D);

    arma::mat iidTot_favorable = iid_favorable;
    arma::mat iidTot_unfavorable = iid_unfavorable;
	if(returnIID>1){
      iidTot_favorable += iidNuisance_favorable;
      iidTot_unfavorable += iidNuisance_unfavorable;
    }
  
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      sigmaC_favorable += conv_to<vec>::from( ntot_control * sum(pow(iidTot_favorable.rows(posC[iter_strata]),2),0) );
      sigmaT_favorable += conv_to<vec>::from( ntot_treatment * sum(pow(iidTot_favorable.rows(posT[iter_strata]),2),0) );
	
      sigmaC_unfavorable += conv_to<vec>::from(ntot_control * sum(pow(iidTot_unfavorable.rows(posC[iter_strata]),2),0) );
      sigmaT_unfavorable += conv_to<vec>::from(ntot_treatment * sum(pow(iidTot_unfavorable.rows(posT[iter_strata]),2),0) );
	
      sigmaC_mixed += conv_to<vec>::from(ntot_control * sum(iidTot_favorable.rows(posC[iter_strata]) % iidTot_unfavorable.rows(posC[iter_strata]),0) );
      sigmaT_mixed += conv_to<vec>::from(ntot_treatment * sum(iidTot_favorable.rows(posT[iter_strata]) % iidTot_unfavorable.rows(posT[iter_strata]),0) );
    }

    Mvar.col(0) = sigmaC_favorable/ntot_control + sigmaT_favorable/ntot_treatment;
    Mvar.col(1) = sigmaC_unfavorable/ntot_control + sigmaT_unfavorable/ntot_treatment; 
    Mvar.col(2) = sigmaC_mixed/ntot_control + sigmaT_mixed/ntot_treatment;
    // Mvar.col(0) = trans(sum(pow(iid_favorable,2), 0));
    // Mvar.col(1) = trans(sum(pow(iid_unfavorable,2), 0));
    // Mvar.col(2) = trans(sum(iid_favorable % iid_unfavorable, 0));

    // second order
    if(hprojection==2){
      // compute variance at the pair level
      arma::vec varUijF,varUijUF,covUijFUF;

      if(keepScore){	  
	arma::mat pairScoreF(ntot_pair,D,fill::zeros);
	arma::mat pairScoreUF(ntot_pair,D,fill::zeros);
	arma::uvec indexRemainingPair;
	arma::uvec iUvec_iter_d(1);
	arma::vec n_cumcontrol = cumsum(n_control);
	arma::vec n_cumpairs = cumsum(n_pairs);
	for(unsigned int iter_d=0; iter_d<D; iter_d++){
	  if(iter_d==0){
	    pairScoreF.col(0) = lsScore[0].col(11);
	    pairScoreUF.col(0) = lsScore[0].col(12);
	  }else{
	    iUvec_iter_d = {iter_d};
	    indexRemainingPair = conv_to<uvec>::from(lsScore[iter_d].col(3));
	    pairScoreF.submat(indexRemainingPair, iUvec_iter_d) = lsScore[iter_d].col(11);
	    pairScoreUF.submat(indexRemainingPair, iUvec_iter_d) = lsScore[iter_d].col(12);
	  }
	}
	pairScoreF.each_row() %= rowweight;
	pairScoreUF.each_row() %= rowweight;
	pairScoreF = cumsum(pairScoreF,1);
	pairScoreUF = cumsum(pairScoreUF,1);

	varUijF = conv_to<vec>::from(sum(pow(pairScoreF,2),0)/ntot_pair) - pow(cumWdelta_favorable,2);
	varUijUF = conv_to<vec>::from(sum(pow(pairScoreUF,2),0)/ntot_pair) - pow(cumWdelta_unfavorable,2);
	covUijFUF = conv_to<vec>::from(sum(pairScoreF % pairScoreUF,0)/ntot_pair) - (cumWdelta_favorable % cumWdelta_unfavorable);
      }else{ // only ok for binary scores i.e. win neutral or loss
	varUijF = cumWdelta_favorable % (1-cumWdelta_favorable);
	varUijUF = cumWdelta_unfavorable % (1-cumWdelta_unfavorable);
	covUijFUF = - cumWdelta_favorable % cumWdelta_unfavorable;
      }

      // compute global variance
      Mvar.col(0) += (varUijF - sigmaC_favorable - sigmaT_favorable)/(ntot_pair);
      Mvar.col(1) += (varUijUF - sigmaC_unfavorable - sigmaT_unfavorable)/(ntot_pair);
      Mvar.col(2) += (covUijFUF - sigmaC_mixed - sigmaT_mixed)/(ntot_pair);
    }
    // Rcout << endl << "delta method" << endl;  
    // delta method
    // var(A-B) = var(A) + var(B) - 2 * cov(A,B)
    // indeed (A-B)' = A' - B' so (A-B)^'2 = A'A' + B'B'  - 2*A'B'
    Mvar.col(3) = Mvar.col(0) + Mvar.col(1) - 2 * Mvar.col(2);
    // var(A/B) = var(A)/B^2 + var(B)*(A^2/B^4) - 2*cov(A,B)A/B^3
    // indeed (A/B)' = A'/B - B'A/B^2 so (A/B)^'2 = A'A'/B^2 + B'B'A^2/B^2 - 2B'A' A/B^3
    Mvar.col(4) = Mvar.col(0)/pow(cumWdelta_unfavorable, 2) + Mvar.col(1) % pow(cumWdelta_favorable,2)/pow(cumWdelta_unfavorable,4) - 2 * Mvar.col(2) % cumWdelta_favorable/pow(cumWdelta_unfavorable, 3);
	// check if no variability then set var(win ratio) to 0.
	for(unsigned int iEndpoint = 0; iEndpoint<D; iEndpoint++){
	  if( (Mvar(iEndpoint,0)==0) && (Mvar(iEndpoint,1)==0)){
		Mvar(iEndpoint,4)=0;
	  }
	}
  }
  
  return ;
}
