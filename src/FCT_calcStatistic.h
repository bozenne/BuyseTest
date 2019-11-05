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

void calcStatistic(arma::mat& delta_netBenefit, arma::mat& delta_winRatio, arma::vec & Delta_netBenefit, arma::vec& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   arma::mat& iidAverage_favorable, arma::mat& iid_unfavorable, arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable,
				   arma::mat& Mvar, int returnIID,
				   std::vector< arma::uvec >& posC, std::vector< arma::uvec >& posT,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
				   const arma::vec& weight, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore);

// * calcStatistic
void calcStatistic(arma::mat& delta_netBenefit, arma::mat& delta_winRatio, arma::vec & Delta_netBenefit, arma::vec& Delta_winRatio,
                   const arma::mat& Mcount_favorable, const arma::mat& Mcount_unfavorable, 
                   arma::mat& iidAverage_favorable, arma::mat& iid_unfavorable, arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable,
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
      iidAverage_favorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_favorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iid_unfavorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iid_unfavorable.rows(posT[iter_strata]) /= n_control[iter_strata];
    }
	
    // *** center
    iidAverage_favorable.each_row() -= conv_to<rowvec>::from(delta_favorable);
    iid_unfavorable.each_row() -= conv_to<rowvec>::from(delta_unfavorable);

    // *** rescale
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
      iidAverage_favorable.rows(posC[iter_strata]) /= ntot_control;
      iidAverage_favorable.rows(posT[iter_strata]) /= ntot_treatment;
      iid_unfavorable.rows(posC[iter_strata]) /= ntot_control;
      iid_unfavorable.rows(posT[iter_strata]) /= ntot_treatment;
    }
	
    // *** weight endpoints and cumulate them
    arma::rowvec rowweight = conv_to<rowvec>::from(weight);

    iidAverage_favorable.each_row() %= rowweight;
    iidAverage_favorable = cumsum(iidAverage_favorable,1);
  
    iid_unfavorable.each_row() %= rowweight;
    iid_unfavorable = cumsum(iid_unfavorable,1);

	if(returnIID>1){
	  iidNuisance_favorable.each_row() %= rowweight;
	  iidNuisance_favorable = cumsum(iidNuisance_favorable,1);
  
	  iidNuisance_unfavorable.each_row() %= rowweight;
	  iidNuisance_unfavorable = cumsum(iidNuisance_unfavorable,1);
	}
	
    // *** sufficient statistics
    arma::vec sigmaC_favorable = zeros<vec>(D);
    arma::vec sigmaT_favorable = zeros<vec>(D);
    arma::vec sigmaC_unfavorable = zeros<vec>(D);
    arma::vec sigmaT_unfavorable = zeros<vec>(D);
    arma::vec sigmaC_mixed = zeros<vec>(D);
    arma::vec sigmaT_mixed = zeros<vec>(D);

    arma::mat iidTotal_favorable = iidAverage_favorable;
    arma::mat iidTotal_unfavorable = iid_unfavorable;
	if(returnIID>1){
      iidTotal_favorable += iidNuisance_favorable;
      iidTotal_unfavorable += iidNuisance_unfavorable;
    }
  
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      sigmaC_favorable += conv_to<vec>::from( ntot_control * sum(pow(iidTotal_favorable.rows(posC[iter_strata]),2),0) );
      sigmaT_favorable += conv_to<vec>::from( ntot_treatment * sum(pow(iidTotal_favorable.rows(posT[iter_strata]),2),0) );
	
      sigmaC_unfavorable += conv_to<vec>::from(ntot_control * sum(pow(iidTotal_unfavorable.rows(posC[iter_strata]),2),0) );
      sigmaT_unfavorable += conv_to<vec>::from(ntot_treatment * sum(pow(iidTotal_unfavorable.rows(posT[iter_strata]),2),0) );
	
      sigmaC_mixed += conv_to<vec>::from(ntot_control * sum(iidTotal_favorable.rows(posC[iter_strata]) % iidTotal_unfavorable.rows(posC[iter_strata]),0) );
      sigmaT_mixed += conv_to<vec>::from(ntot_treatment * sum(iidTotal_favorable.rows(posT[iter_strata]) % iidTotal_unfavorable.rows(posT[iter_strata]),0) );
    }

    Mvar.col(0) = sigmaC_favorable/ntot_control + sigmaT_favorable/ntot_treatment;
    Mvar.col(1) = sigmaC_unfavorable/ntot_control + sigmaT_unfavorable/ntot_treatment; 
    Mvar.col(2) = sigmaC_mixed/ntot_control + sigmaT_mixed/ntot_treatment;
    // Mvar.col(0) = trans(sum(pow(iidAverage_favorable,2), 0));
    // Mvar.col(1) = trans(sum(pow(iid_unfavorable,2), 0));
    // Mvar.col(2) = trans(sum(iidAverage_favorable % iid_unfavorable, 0));

    // second order
    if(hprojection==2){

      if(keepScore){

		// reconstruct individual score
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

		// compute second order h-projection and corresponding moments
		rowvec H2_favorable, H2_unfavorable;
		arma::mat H2_moments(5,D,fill::zeros);
		int iter_strata, iter_C, iter_T;
		
		for(int iter_pair=0; iter_pair<ntot_pair ; iter_pair++){ 
		  iter_strata = lsScore[0](iter_pair,0);
		  iter_C = lsScore[0](iter_pair,4); // index within strata
		  iter_T = lsScore[0](iter_pair,5); // index within strata
		  H2_favorable = pairScoreF.row(iter_pair) - ntot_control * iidTotal_favorable.row(posC[iter_strata][iter_C]) - ntot_treatment * iidTotal_favorable.row(posT[iter_strata][iter_T]);
		  H2_unfavorable = pairScoreUF.row(iter_pair) - ntot_control * iidTotal_unfavorable.row(posC[iter_strata][iter_C]) - ntot_treatment * iidTotal_unfavorable.row(posT[iter_strata][iter_T]);

		  H2_moments.row(0) += H2_favorable;
		  H2_moments.row(1) += H2_unfavorable;
		  H2_moments.row(2) += pow(H2_favorable,2.0);
		  H2_moments.row(3) += pow(H2_unfavorable,2.0);
		  H2_moments.row(4) += H2_favorable % H2_unfavorable;		  
		}
		H2_moments /= ntot_pair;

		// update variance:  Var(H2) = Var( 1/mn \sum_i,j H2_ij) = 1/mn Var(H2_ij)
		Mvar.col(0) += trans((H2_moments.row(2) - pow(H2_moments.row(0),2))/(ntot_pair));
		Mvar.col(1) += trans((H2_moments.row(3) - pow(H2_moments.row(1),2))/(ntot_pair));
		Mvar.col(2) += trans((H2_moments.row(4) - H2_moments.row(0) % H2_moments.row(1))/(ntot_pair));
		
	  }else{ // only ok for binary scores i.e. win neutral or loss
		Mvar.col(0) += (cumWdelta_favorable % (1-cumWdelta_favorable) - sigmaC_favorable - sigmaT_favorable)/(ntot_pair);
		Mvar.col(1) += (cumWdelta_unfavorable % (1-cumWdelta_unfavorable) - sigmaC_unfavorable - sigmaT_unfavorable)/(ntot_pair);
		Mvar.col(2) += (- cumWdelta_favorable % cumWdelta_unfavorable - sigmaC_mixed - sigmaT_mixed)/(ntot_pair);
	  }
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
