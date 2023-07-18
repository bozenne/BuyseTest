// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>

// :cppFile:{FCT_buyseTest.cpp}:end:

void calcStatistic(arma::cube& delta, arma::mat& Delta,
                   arma::mat Mcount_favorable, arma::mat Mcount_unfavorable, arma::mat Mcount_neutral, arma::mat Mcount_uninf, 
                   arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, arma::mat& iidAverage_neutral,
		   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, arma::mat& iidNuisance_neutral,
		   arma::mat& Mvar, std::vector< int > returnIID,
		   const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT, const arma::vec& weightObs,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weightEndpoint, double pool, arma::vec& weightPool,
		   bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore, bool paired);

// * calcStatistic
void calcStatistic(arma::cube& delta, arma::mat& Delta,
                   arma::mat Mcount_favorable, arma::mat Mcount_unfavorable, arma::mat Mcount_neutral, arma::mat Mcount_uninf,
                   arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, arma::mat& iidAverage_neutral,
		   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, arma::mat& iidNuisance_neutral,
		   arma::mat& Mvar, std::vector< int > returnIID,
		   const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT, const arma::vec& weightObs,
                   const unsigned int& D, const int& n_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weightEndpoint, double pool, arma::vec& weightPool,
		   bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore, bool paired){

  // ** total number of pairs and patients in each arm
  double ntot_pair = 0;
  for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    ntot_pair += n_pairs[iter_strata];
  }

  // ** add neutral score
  if(addHalfNeutral){
    Mcount_favorable += 0.5*Mcount_neutral;
    Mcount_unfavorable += 0.5*Mcount_neutral;

    if(returnIID[0]>0 || pool>=3){
      iidAverage_favorable += 0.5*iidAverage_neutral;
      iidAverage_unfavorable += 0.5*iidAverage_neutral;
    }
    if(returnIID[1]>0){
      iidNuisance_favorable += 0.5*iidNuisance_neutral;
      iidNuisance_unfavorable += 0.5*iidNuisance_neutral;
    }
  }

  // ** Strata and outcome specifc
  // arma::conv_to<arma::rowvec>::from()
  arma::rowvec rowweightEndpoint = weightEndpoint.t();

  // *** Point estimates
  // Mann Whitney parameter equals number of (un)favorable pairs divided by the number of pairs
  delta.slice(0) = Mcount_favorable;
  delta.slice(0).each_col() /= n_pairs;

  delta.slice(1) = Mcount_unfavorable;
  delta.slice(1).each_col() /= n_pairs;

  delta.slice(2) = Mcount_neutral;
  delta.slice(2).each_col() /= n_pairs;

  delta.slice(3) = Mcount_uninf;
  delta.slice(3).each_col() /= n_pairs;

  // net benefit equals (number of favorable pairs minus number of unfavorable pairs) divided by number of pairs
  delta.slice(4) = delta.slice(0) - delta.slice(1);

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs  
  delta.slice(5) = delta.slice(0) / delta.slice(1);


  // *** iid and variance estimation
  // \hat{U}-U = \frac{1}{m} h_1^{gamma}(i) + ...
  // where h_1^{gamma}(i) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma}

  // storeiid average plus iid nuisance
  arma::mat iidTotal_favorable;
  arma::mat iidTotal_unfavorable;
  // for second order projection
  arma::mat h1_favorable; 
  arma::mat h1_unfavorable;
 
  if(returnIID[0]>0 || pool>=3){

    // **** compute expectation (E[s^{\gamma}_{l,j}|x_l]) from sum 
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      iidAverage_favorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_favorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iidAverage_unfavorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_unfavorable.rows(posT[iter_strata]) /= n_control[iter_strata];
      iidAverage_neutral.rows(posC[iter_strata]) /= n_treatment[iter_strata];
      iidAverage_neutral.rows(posT[iter_strata]) /= n_control[iter_strata];
    }

    // **** center to obtain (endpoint specific) first order projection (h_1^{\gamma}(l) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma})
    arma::uvec ivecD(1);
    
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
      for(unsigned int iter_d=0 ; iter_d < D ; iter_d ++){
	ivecD[0] = iter_d;
	
	iidAverage_favorable(posC[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_favorable(posT[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];
      
	iidAverage_unfavorable(posC[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_unfavorable(posT[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];

	iidAverage_neutral(posC[iter_strata],ivecD) -= delta.slice(2)(iter_strata,iter_d); // same as -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
	iidAverage_neutral(posT[iter_strata],ivecD) -= delta.slice(2)(iter_strata,iter_d); // same as -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
      }
    }

    if(returnIID[0]>0 && hprojection==2 && keepScore){ // used when computing the variance of H projection order 2 with keepScore
      h1_favorable = iidAverage_favorable;
      h1_unfavorable = iidAverage_unfavorable;
      if(returnIID[1]>0){ // add iid of the nuisance parameters
	for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
	  h1_favorable.rows(posC[iter_strata]) += iidNuisance_favorable.rows(posC[iter_strata]) * n_control[iter_strata];
	  h1_unfavorable.rows(posT[iter_strata]) += iidNuisance_favorable.rows(posT[iter_strata]) * n_treatment[iter_strata];
	  h1_favorable.rows(posC[iter_strata]) += iidNuisance_unfavorable.rows(posC[iter_strata]) * n_control[iter_strata];
	  h1_unfavorable.rows(posT[iter_strata]) += iidNuisance_unfavorable.rows(posT[iter_strata]) * n_treatment[iter_strata];
	}
      }
    }

    // **** rescale to move from empirical average to sum, i.e. obtain h_1^{\gamma}(l)/m
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
      iidAverage_favorable.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_favorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];

      iidAverage_unfavorable.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_unfavorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];

      iidAverage_neutral.rows(posC[iter_strata]) /= n_control[iter_strata];
      iidAverage_neutral.rows(posT[iter_strata]) /= n_treatment[iter_strata];
    }
      
    // **** add the two source of uncertainty into a single iid
    if(returnIID[1]>0){
      iidTotal_favorable = iidAverage_favorable + iidNuisance_favorable;
      iidTotal_unfavorable = iidAverage_unfavorable + iidNuisance_unfavorable;
    }else{
      iidTotal_favorable = iidAverage_favorable;
      iidTotal_unfavorable = iidAverage_unfavorable;
    }

    // **** weight endpoints and cumulate them to obtain (cumulative) first order projection
    // logically should be in the next subection (Global) but is put here because used for the point estimate with pooling = 2.1-2.4
    iidTotal_favorable.each_row() %= rowweightEndpoint;
    iidTotal_favorable = arma::cumsum(iidTotal_favorable,1);

    iidTotal_unfavorable.each_row() %= rowweightEndpoint;
    iidTotal_unfavorable = arma::cumsum(iidTotal_unfavorable,1);

    if(returnIID[0]>0 && hprojection==2 && keepScore){ // used when computing the variance of H projection order 2 with keepScore
      h1_favorable.each_row() %= rowweightEndpoint;
      h1_favorable = arma::cumsum(h1_favorable,1);

      h1_unfavorable.each_row() %= rowweightEndpoint;
      h1_unfavorable = arma::cumsum(h1_unfavorable,1);
    }

  }

  // ** Global (across strata)

  // *** Strata specific weights
  if(pool==0){ // Buyse
    weightPool = n_pairs;
  }else if(pool==1){ // CMH
    weightPool = n_pairs/(n_control+n_treatment);
  }else if(pool==2){ // equal
    weightPool.fill(1.0);
  }else if(pool >= 3){ // precision
    arma::vec iIID;
    arma::uvec iPos;
    arma::uvec iUvec = {D-1};
    double iDeltaFA, iDeltaUF;
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
      iPos = arma::join_cols(posC[iter_strata], posT[iter_strata]);
      if(pool==3.1){
	iIID = iidTotal_favorable(iPos,iUvec);
      }else if(pool==3.2){
	iIID = iidTotal_unfavorable(iPos,iUvec);
      }else if(pool==3.3){
	iIID = iidTotal_favorable(iPos,iUvec) - iidTotal_unfavorable(iPos,iUvec);
      }else if(pool==3.4){
	iDeltaFA = sum(delta.slice(0).row(iter_strata) % rowweightEndpoint);
	iDeltaUF = sum(delta.slice(1).row(iter_strata) % rowweightEndpoint);
	iIID = iidTotal_favorable(iPos,iUvec)/iDeltaUF - iidTotal_unfavorable(iPos,iUvec)/(pow(iDeltaUF,2)/iDeltaFA);
      }
      weightPool[iter_strata] = 1/sum(iIID % iIID);      
    }
  }
  weightPool /= sum(weightPool); // normalize weights to sum up to 1

  // *** strata-weighted cumulative endpoint estimate
  arma::mat cumdelta_favorable = arma::cumsum(delta.slice(0).each_row() % rowweightEndpoint,1);
  arma::mat wcumdelta_favorable = cumdelta_favorable.each_col() % weightPool; 
  
  arma::mat cumdelta_unfavorable = arma::cumsum(delta.slice(1).each_row() % rowweightEndpoint,1);
  arma::mat wcumdelta_unfavorable = cumdelta_unfavorable.each_col() % weightPool; 

  arma::mat cumdelta_neutral = delta.slice(2).each_row() % rowweightEndpoint;
  arma::mat wcumdelta_neutral = cumdelta_neutral.each_col() % weightPool; 

  arma::mat cumdelta_uninf = delta.slice(3).each_row() % rowweightEndpoint;
  arma::mat wcumdelta_uninf = cumdelta_uninf.each_col() % weightPool; 

  // *** Point estimates
  // \Delta = \sum_endpoint \sum_strata w(strata) w(endpoint) \delta(strata, endpoint)
  Delta.col(0) = arma::trans(arma::sum(wcumdelta_favorable,0));
  Delta.col(1) = arma::trans(arma::sum(wcumdelta_unfavorable,0));
  Delta.col(2) = arma::trans(arma::sum(wcumdelta_neutral,0));
  Delta.col(3) = arma::trans(arma::sum(wcumdelta_uninf,0));
  
  // net benefit 
  Delta.col(4) = arma::trans(arma::sum(wcumdelta_favorable - wcumdelta_unfavorable,0));

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs
  if(pool==3.4){
    arma::mat winRatioStrata = cumdelta_favorable/cumdelta_unfavorable;
    Delta.col(5) = arma::trans(arma::sum(winRatioStrata.each_col() % weightPool,0));
  }else{
    Delta.col(5) = arma::trans(arma::sum(wcumdelta_favorable,0)/arma::sum(wcumdelta_unfavorable,0));
  }

  // *** iid and variance estimation
  if(returnIID[0] > 0){

    // **** rescale to account for the pooling across strata i.e. T = \sum_s n_pair(s) T(s) / n_tot
    if(n_strata>1){
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 

	iidTotal_favorable.rows(posC[iter_strata]) *= weightPool[iter_strata];
	iidTotal_favorable.rows(posT[iter_strata]) *= weightPool[iter_strata];
	
	iidTotal_unfavorable.rows(posC[iter_strata]) *= weightPool[iter_strata];
	iidTotal_unfavorable.rows(posT[iter_strata]) *= weightPool[iter_strata];
	
      }
    }
  
    // **** estimate variance
    // first order
    arma::mat iid2;

    iid2 = pow(iidTotal_favorable,2);
    Mvar.col(0) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));
    iid2 = pow(iidTotal_unfavorable,2);
    Mvar.col(1) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));
    iid2 = iidTotal_favorable % iidTotal_unfavorable;
    Mvar.col(2) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));

    if(paired){
      Mvar.col(0) += arma::var(cumdelta_favorable, 0, 0)/n_strata;
      Mvar.col(1) += arma::var(cumdelta_unfavorable, 0, 0)/n_strata;
      Mvar.col(2) += arma::cov(cumdelta_favorable, cumdelta_unfavorable, 0)/n_strata;
    }
    
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
	pairScoreF.each_row() %= rowweightEndpoint;
	pairScoreUF.each_row() %= rowweightEndpoint;
	pairScoreF = arma::cumsum(pairScoreF,1);
	pairScoreUF = arma::cumsum(pairScoreUF,1);
	// compute second order h-projection, its variance and covariance
	// NOTE: when computing h we use iid = E[]-U so h2 = s - E[] - E[] + U = s - iid - iid - U
	arma::rowvec H2_favorable, H2_unfavorable;
	arma::mat H2_moments(5,D,arma::fill::zeros);
	int iter_strata, iter_C, iter_T;
	double iter_weight;
	
	// endpoint specific    
	for(int iter_pair=0; iter_pair<ntot_pair ; iter_pair++){

	  iter_strata = lsScore[0](iter_pair,0);
	  iter_C = lsScore[0](iter_pair,4); // index within strata
	  iter_T = lsScore[0](iter_pair,5); // index within strata
	  iter_weight = weightObs(posC[iter_strata][iter_C])*weightObs(posT[iter_strata][iter_T]);

	  H2_favorable = (pairScoreF.row(iter_pair) - h1_favorable.row(posC[iter_strata][iter_C]) - h1_favorable.row(posT[iter_strata][iter_T]) - cumdelta_favorable.row(iter_strata));
	  H2_unfavorable = (pairScoreUF.row(iter_pair) - h1_unfavorable.row(posC[iter_strata][iter_C]) - h1_unfavorable.row(posT[iter_strata][iter_T]) - cumdelta_unfavorable.row(iter_strata));
	  
	  // Rcpp::Rcout << iter_pair << " (" << iter_strata << ":" << iter_C << " " << iter_T << ") " << H2_favorable[0] << " " << H2_unfavorable[0] << std::endl;

	  // var(H2) = \sum_ij H2 / n_pair[iter_strata]
	  // var(1/nm \sum_ij H2) = \sum_ij H2 / n_pair[iter_strata]^2
	  // w_{pooling}^2 var(1/nm \sum_ij H2) = (n_pairs[iter_strata]/ntot_pair)^2 \sum_ij H2 / n_pairs[iter_strata]^2 = \sum_ij H2/ntot_pair^2
	  Mvar.col(0) += iter_weight * arma::trans(pow(H2_favorable,2)) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	  Mvar.col(1) += iter_weight * arma::trans(pow(H2_unfavorable,2)) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	  Mvar.col(2) += iter_weight * arma::trans(H2_favorable % H2_unfavorable) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	}

	
      }else{ // only ok for binary scores i.e. win neutral or loss
	arma::colvec strataDelta_favorable(D);
	arma::colvec strataDelta_unfavorable(D);
	arma::vec strataWeightC, strataWeightT;

	for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	  // strata specific cumulative endpoint
	  strataDelta_favorable = arma::trans(cumdelta_favorable.row(iter_strata));
	  strataDelta_unfavorable = arma::trans(cumdelta_unfavorable.row(iter_strata));
	  strataWeightC = weightObs(posC[iter_strata]);
	  strataWeightT = weightObs(posT[iter_strata]);
	  
	  // (n_pairs[iter_strata] / ntot_pair)^2 * (1/n_pairs[iter_strata]) = n_pairs[iter_strata]/ntot_pair^2
	  Mvar.col(0) += strataDelta_favorable % (1-strataDelta_favorable) * pow(weightPool[iter_strata],2)/n_pairs[iter_strata];
	  iid2 = pow(iidTotal_favorable.rows(posC[iter_strata]),2);
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / n_treatment[iter_strata] );
	  iid2 = pow(iidTotal_favorable.rows(posT[iter_strata]),2);
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / n_control[iter_strata] );
	  
	  Mvar.col(1) += strataDelta_unfavorable % (1-strataDelta_unfavorable) * pow(weightPool[iter_strata],2)/n_pairs[iter_strata];
	  iid2 = pow(iidTotal_unfavorable.rows(posC[iter_strata]),2);
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / n_treatment[iter_strata] );
	  iid2 = pow(iidTotal_unfavorable.rows(posT[iter_strata]),2);
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / n_control[iter_strata] );
	  
	  Mvar.col(2) -= strataDelta_favorable % strataDelta_unfavorable * pow(weightPool[iter_strata],2)/n_pairs[iter_strata];
	  iid2 = iidTotal_favorable.rows(posC[iter_strata]) % iidTotal_unfavorable.rows(posC[iter_strata]);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / n_treatment[iter_strata]);
	  iid2 = iidTotal_favorable.rows(posT[iter_strata]) % iidTotal_unfavorable.rows(posT[iter_strata]);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / n_control[iter_strata]);
	  }
      }

    }
    // Rcpp::Rcout << std::endl << "delta method" << std::endl;  
    // delta method
    // var(A-B) = var(A) + var(B) - 2 * cov(A,B)
    // indeed (A-B)' = A' - B' so (A-B)^'2 = A'A' + B'B'  - 2*A'B'
    Mvar.col(3) = Mvar.col(0) + Mvar.col(1) - 2 * Mvar.col(2);
    // var(A/B) = var(A)/B^2 + var(B)*(A^2/B^4) - 2*cov(A,B)A/B^3
    // indeed (A/B)' = A'/B - B'A/B^2 so (A/B)^'2 = A'A'/B^2 + B'B'A^2/B^4 - 2B'A' A/B^3
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
