// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:

void calcStatistic(arma::cube& delta, arma::mat& Delta,
                   arma::mat Mcount_favorable, arma::mat Mcount_unfavorable, arma::mat Mcount_neutral, 
                   arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, arma::mat& iidAverage_neutral,
		   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, arma::mat& iidNuisance_neutral,
		   arma::mat& Mvar, int returnIID,
		   std::vector< arma::uvec >& posC, std::vector< arma::uvec >& posT,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weight, bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore);

// * calcStatistic
void calcStatistic(arma::cube& delta, arma::mat& Delta,
                   arma::mat Mcount_favorable, arma::mat Mcount_unfavorable, arma::mat Mcount_neutral, 
                   arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, arma::mat& iidAverage_neutral,
		   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, arma::mat& iidNuisance_neutral,
		   arma::mat& Mvar, int returnIID,
		   std::vector< arma::uvec >& posC, std::vector< arma::uvec >& posT,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weight, bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore){
  
  // ** total number of pairs and patients in each arm
  double ntot_pair = 0;
  for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    ntot_pair += n_pairs[iter_strata];
  }
    
  // ** add neutral score
  if(addHalfNeutral){
    Mcount_favorable += 0.5*Mcount_neutral;
    Mcount_unfavorable += 0.5*Mcount_neutral;

    if(returnIID>0){
      iidAverage_favorable += 0.5*iidAverage_neutral;
      iidAverage_unfavorable += 0.5*iidAverage_neutral;
    }
    if(returnIID>1){
      iidNuisance_favorable += 0.5*iidNuisance_neutral;
      iidNuisance_unfavorable += 0.5*iidNuisance_neutral;
    }
  }

  // ** Point estimates

  // *** strata specifc
  // Mann Whitney parameter equals number of (un)favorable pairs divided by the number of pairs
  delta.slice(0) = Mcount_favorable;
  delta.slice(0).each_col() /= n_pairs;

  delta.slice(1) = Mcount_unfavorable;
  delta.slice(1).each_col() /= n_pairs;

  // net benefit equals (number of favorable pairs minus number of unfavorable pairs) divided by number of pairs
  delta.slice(2) = delta.slice(0) - delta.slice(1);

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
  delta.slice(3) = delta.slice(0) / delta.slice(1);

  // *** across strata
  // weight and cumulate endpoints
  // cumsum:  For matrix X, return a matrix containing the cumulative sum of elements in each column (dim = 0), or each row (dim = 1)
  arma::rowvec rowweight = arma::conv_to<arma::rowvec>::from(weight);

  arma::mat cumWcount_favorable = Mcount_favorable;
  cumWcount_favorable.each_row() %= rowweight; 
  cumWcount_favorable = arma::cumsum(cumWcount_favorable,1); 

  arma::mat cumWcount_unfavorable = Mcount_unfavorable;
  cumWcount_unfavorable.each_row() %= rowweight;
  cumWcount_unfavorable = arma::cumsum(cumWcount_unfavorable,1); 
  
  // Mann Whitney parameter equals number of favorable pairs divided by the number of pairs (i.e. proportion in favor of the treatment)  
  Delta.col(0) = arma::trans(arma::sum(cumWcount_favorable,0))/(double)(ntot_pair);
  Delta.col(1) = arma::trans(arma::sum(cumWcount_unfavorable,0))/(double)(ntot_pair);

  // net benefit equals (number of favorable pairs minus number of unfavorable pairs) divided by number of pairs
  Delta.col(2) = Delta.col(0) - Delta.col(1);

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
  Delta.col(3) = Delta.col(0) / Delta.col(1);

  // ** iid and variance estimation
  // \hat{U}-U = \frac{1}{m} h_1^{gamma}(i) + ...
  // where h_1^{gamma}(i) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma}
  if(returnIID > 0){

    // *** compute expectation (E[s^{\gamma}_{l,j}|x_l]) from sum 
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      iidAverage_favorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_favorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iidAverage_unfavorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_unfavorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iidAverage_neutral.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_neutral.rows(posT[iter_strata]) /= n_control[iter_strata];
    }

    // *** center to obtain (endpoint specific) first order projection (h_1^{\gamma}(l) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma})
    arma::uvec ivecD(1);
    
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
      for(unsigned int iter_d=0 ; iter_d < D ; iter_d ++){
	ivecD[0] = iter_d;
	
	iidAverage_favorable(posC[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_favorable(posT[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];
      
	iidAverage_unfavorable(posC[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_unfavorable(posT[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];

	iidAverage_neutral(posC[iter_strata],ivecD) -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_neutral(posT[iter_strata],ivecD) -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
      }
    }

    arma::mat h1_favorable = iidAverage_favorable;
    arma::mat h1_unfavorable = iidAverage_unfavorable;
    if(returnIID>1){ // used when computing the variance of H projection order 2 with keepScore
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
	h1_favorable.rows(posC[iter_strata]) += iidNuisance_favorable.rows(posC[iter_strata]) * n_control[iter_strata];
	h1_unfavorable.rows(posT[iter_strata]) += iidNuisance_favorable.rows(posT[iter_strata]) * n_treatment[iter_strata];
	h1_favorable.rows(posC[iter_strata]) += iidNuisance_unfavorable.rows(posC[iter_strata]) * n_control[iter_strata];
	h1_unfavorable.rows(posT[iter_strata]) += iidNuisance_unfavorable.rows(posT[iter_strata]) * n_treatment[iter_strata];
      }
    }

    // *** rescale to move from empirical average to sum, i.e. obtain h_1^{\gamma}(l)/m
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
      iidAverage_favorable.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_favorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];

      iidAverage_unfavorable.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_unfavorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];

      iidAverage_neutral.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_neutral.rows(posT[iter_strata]) /= n_treatment[iter_strata];
    }

    // *** rescale to account for the pooling across strata i.e. T = \sum_s n_pair(s) T(s) / n_tot
    if(n_strata>1){
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 

	iidAverage_favorable.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	iidAverage_favorable.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	
	iidAverage_unfavorable.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	iidAverage_unfavorable.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	// not used here but ensure consistent export with other iid for R object
	iidAverage_neutral.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	iidAverage_neutral.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;

	if(returnIID>1){
	  iidNuisance_favorable.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	  iidNuisance_favorable.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;

	  iidNuisance_unfavorable.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	  iidNuisance_unfavorable.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	  // not used here but ensure consistent export with other iid for R object
	  iidNuisance_neutral.rows(posC[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	  iidNuisance_neutral.rows(posT[iter_strata]) *= n_pairs[iter_strata] / ntot_pair;
	}
      }
    }

    // *** weight endpoints and cumulate them to obtain (cumulative) first order projection
    arma::mat iidTotal_favorable = iidAverage_favorable;
    arma::mat iidTotal_unfavorable = iidAverage_unfavorable;
    if(returnIID>1){
      iidTotal_favorable += iidNuisance_favorable;
      iidTotal_unfavorable += iidNuisance_unfavorable;
    }

    h1_favorable.each_row() %= rowweight;
    h1_favorable = arma::cumsum(h1_favorable,1);
    iidTotal_favorable.each_row() %= rowweight;
    iidTotal_favorable = arma::cumsum(iidTotal_favorable,1);
  
    h1_unfavorable.each_row() %= rowweight;
    h1_unfavorable = arma::cumsum(h1_unfavorable,1);
    iidTotal_unfavorable.each_row() %= rowweight;
    iidTotal_unfavorable = arma::cumsum(iidTotal_unfavorable,1);
	
    // *** estimate variance
    // first order
    Mvar.col(0) = arma::trans(arma::sum(pow(iidTotal_favorable,2), 0));
    Mvar.col(1) = arma::trans(arma::sum(pow(iidTotal_unfavorable,2), 0));
    Mvar.col(2) = arma::trans(arma::sum(iidTotal_favorable % iidTotal_unfavorable, 0));
 	
    // second order
    if(hprojection==2){

      if(keepScore){
	// reconstruct individual score
	arma::mat pairScoreF(ntot_pair,D,arma::fill::zeros);
	arma::mat pairScoreUF(ntot_pair,D,arma::fill::zeros);
	arma::uvec indexRemainingPair;
	arma::uvec iUvec_iter_d(1);
	arma::vec n_cumcontrol = arma::cumsum(n_control);
	arma::vec n_cumpairs = arma::cumsum(n_pairs);
	for(unsigned int iter_d=0; iter_d<D; iter_d++){
	  if(iter_d==0){
	    pairScoreF.col(0) = lsScore[0].col(11);
	    pairScoreUF.col(0) = lsScore[0].col(12);
	  }else{
	    iUvec_iter_d = {iter_d};
	    indexRemainingPair = arma::conv_to<arma::uvec>::from(lsScore[iter_d].col(3));
	    pairScoreF.submat(indexRemainingPair, iUvec_iter_d) = lsScore[iter_d].col(11);
	    pairScoreUF.submat(indexRemainingPair, iUvec_iter_d) = lsScore[iter_d].col(12);
	  }
	}
	pairScoreF.each_row() %= rowweight;
	pairScoreUF.each_row() %= rowweight;
	pairScoreF = arma::cumsum(pairScoreF,1);
	pairScoreUF = arma::cumsum(pairScoreUF,1);
	// compute second order h-projection, its variance and covariance
	// NOTE: when computing h we use iid = E[]-U so h2 = s - E[] - E[] + U = s - iid - iid - U
	arma::rowvec H2_favorable, H2_unfavorable;
	arma::mat H2_moments(5,D,arma::fill::zeros);
	int iter_strata, iter_C, iter_T;
	arma::mat cumdelta_favorable = arma::cumsum(delta.slice(0),1);
	arma::mat cumdelta_unfavorable = arma::cumsum(delta.slice(1),1);

	// endpoint specific    
	for(int iter_pair=0; iter_pair<ntot_pair ; iter_pair++){ 
	  iter_strata = lsScore[0](iter_pair,0);
	  iter_C = lsScore[0](iter_pair,4); // index within strata
	  iter_T = lsScore[0](iter_pair,5); // index within strata
	  H2_favorable = (pairScoreF.row(iter_pair) - h1_favorable.row(posC[iter_strata][iter_C]) - h1_favorable.row(posT[iter_strata][iter_T]) - cumdelta_favorable.row(iter_strata));
	  H2_unfavorable = (pairScoreUF.row(iter_pair) - h1_unfavorable.row(posC[iter_strata][iter_C]) - h1_unfavorable.row(posT[iter_strata][iter_T]) - cumdelta_unfavorable.row(iter_strata));

	  // Rcpp::Rcout << iter_pair << " (" << iter_strata << ":" << iter_C << " " << iter_T << ") " << H2_favorable[0] << " " << H2_unfavorable[0] << std::endl;

	  // var(H2) = \sum_ij H2 / n_pair[iter_strata]
	  // var(1/nm \sum_ij H2) = \sum_ij H2 / n_pair[iter_strata]^2
	  // w_{pooling} var(1/nm \sum_ij H2) = (n_pair[iter_strata]/ntot_pair)^2 \sum_ij H2 / n_pair[iter_strata]^2 = \sum_ij H2/ntot_pair^2
	  Mvar.col(0) += arma::trans(pow(H2_favorable,2)) / pow(ntot_pair,2);
	  Mvar.col(1) += arma::trans(pow(H2_unfavorable,2)) / pow(ntot_pair,2);
	  Mvar.col(2) += arma::trans(H2_favorable % H2_unfavorable) / pow(ntot_pair,2);
	}

	
      }else{ // only ok for binary scores i.e. win neutral or loss
	arma::colvec strataDelta_favorable(D);
	arma::colvec strataDelta_unfavorable(D);

	for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata

	  // strata specific cumulative endpoint
	  strataDelta_favorable = arma::trans(cumWcount_favorable.row(iter_strata))/n_pairs[iter_strata];
	  strataDelta_unfavorable = arma::trans(cumWcount_unfavorable.row(iter_strata))/n_pairs[iter_strata];

	  // (n_pairs[iter_strata] / ntot_pair)^2 * (1/n_pairs[iter_strata]) = n_pairs[iter_strata]/ntot_pair^2
	  Mvar.col(0) += strataDelta_favorable % (1-strataDelta_favorable) * n_pairs[iter_strata]/pow(ntot_pair,2);
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(pow(iidTotal_favorable.rows(posC[iter_strata]),2),0) / n_treatment[iter_strata] );
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(pow(iidTotal_favorable.rows(posT[iter_strata]),2),0) / n_control[iter_strata] );

	  Mvar.col(1) += strataDelta_unfavorable % (1-strataDelta_unfavorable) * n_pairs[iter_strata]/pow(ntot_pair,2);
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(pow(iidTotal_unfavorable.rows(posC[iter_strata]),2),0) / n_treatment[iter_strata] );
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(pow(iidTotal_unfavorable.rows(posT[iter_strata]),2),0) / n_control[iter_strata] );
	
	  Mvar.col(2) -= strataDelta_favorable % strataDelta_unfavorable * n_pairs[iter_strata]/pow(ntot_pair,2);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iidTotal_favorable.rows(posC[iter_strata]) % iidTotal_unfavorable.rows(posC[iter_strata]),0) / n_treatment[iter_strata]);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iidTotal_favorable.rows(posT[iter_strata]) % iidTotal_unfavorable.rows(posT[iter_strata]),0) / n_control[iter_strata]);
	}
      }

    }
    // Rcpp::Rcout << std::endl << "delta method" << std::endl;  
    // delta method
    // var(A-B) = var(A) + var(B) - 2 * cov(A,B)
    // indeed (A-B)' = A' - B' so (A-B)^'2 = A'A' + B'B'  - 2*A'B'
    Mvar.col(3) = Mvar.col(0) + Mvar.col(1) - 2 * Mvar.col(2);
    // var(A/B) = var(A)/B^2 + var(B)*(A^2/B^4) - 2*cov(A,B)A/B^3
    // indeed (A/B)' = A'/B - B'A/B^2 so (A/B)^'2 = A'A'/B^2 + B'B'A^2/B^2 - 2B'A' A/B^3
    Mvar.col(4) = Mvar.col(0)/pow(Delta.col(1), 2) + Mvar.col(1) % pow(Delta.col(0),2)/pow(Delta.col(1),4) - 2 * Mvar.col(2) % Delta.col(0)/pow(Delta.col(1), 3);
    // Mann-Whitney parameter is the same as the proportion in favor of treatment
    // check if no variability then set var(win ratio) to 0.
    for(unsigned int iEndpoint = 0; iEndpoint<D; iEndpoint++){
      if( (Mvar(iEndpoint,0)==0) && (Mvar(iEndpoint,1)==0)){
	Mvar(iEndpoint,4)=0;
      }
    }
  }
  
  return ;
}
