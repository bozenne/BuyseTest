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
                   const unsigned int& D, const int& n_strata, arma::mat grid_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weightEndpoint, double pool, arma::vec& weightPool,
		   bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore, bool paired);

// * calcStatistic
void calcStatistic(arma::cube& delta, arma::mat& Delta,
                   arma::mat Mcount_favorable, arma::mat Mcount_unfavorable, arma::mat Mcount_neutral, arma::mat Mcount_uninf,
                   arma::mat& iidAverage_favorable, arma::mat& iidAverage_unfavorable, arma::mat& iidAverage_neutral,
		   arma::mat& iidNuisance_favorable, arma::mat& iidNuisance_unfavorable, arma::mat& iidNuisance_neutral,
		   arma::mat& Mvar, std::vector< int > returnIID,
		   const std::vector< arma::uvec >& posC, const std::vector< arma::uvec >& posT, const arma::vec& weightObs,
                   const unsigned int& D, const int& n_strata, arma::mat grid_strata, const arma::vec& n_pairs, const arma::vec& n_control, const arma::vec& n_treatment,
		   const arma::vec& weightEndpoint, double pool, arma::vec& weightPool,
		   bool addHalfNeutral, int hprojection, const std::vector< arma::mat >& lsScore, bool keepScore, bool paired){

  // ** Total number of pairs and patients in each arm
  double ntot_pair = 0;
  int iter_strataC;
  int iter_strataT;
  for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    ntot_pair += n_pairs[iter_strata];
  }

  int ntot_control=0,ntot_treatment=0;
  arma::vec nObs_strata(grid_strata.max()+1); // number of observations per strata
  if(pool == 3){ // standardization
    for(int iter_strata = 0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
      if(grid_strata(iter_strata,0)==grid_strata(iter_strata,1)){
	ntot_control += n_control[iter_strata];
	ntot_treatment += n_treatment[iter_strata];
	nObs_strata[grid_strata(iter_strata,0)] = n_control[iter_strata]+n_treatment[iter_strata];
      }
    }
  }

  // ** fixed weights to pool across strata (optimal weight defined later)
  if(pool==0){ // Buyse
    weightPool = n_pairs/ntot_pair;
  }else if(pool==1){ // CMH
    weightPool = n_pairs/(n_control+n_treatment);
    weightPool /= sum(weightPool); // normalize weights to sum up to 1
  }else if(pool==2){ // equal
    weightPool.fill(1.0/n_strata);
  }else if(pool==3){ // standardization
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
      weightPool[iter_strata] = nObs_strata(grid_strata(iter_strata,0))*nObs_strata(grid_strata(iter_strata,1))/pow(ntot_control+ntot_treatment,2);
    }      
  } 

  // ** proportion of favorable/unfavorable/neutral/uniformative pairs within each strata and priority 
  delta.slice(0) = Mcount_favorable;
  delta.slice(0).each_col() /= n_pairs;

  delta.slice(1) = Mcount_unfavorable;
  delta.slice(1).each_col() /= n_pairs;

  delta.slice(2) = Mcount_neutral;
  delta.slice(2).each_col() /= n_pairs;

  delta.slice(3) = Mcount_uninf;
  delta.slice(3).each_col() /= n_pairs;

  // ** derive iid in case variance estimator is need to obtain strata-specific weights
  // \hat{U}-U = \frac{1}{m} h_1^{gamma}(i) + ...
  // where h_1^{gamma}(i) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma}

  // arma::conv_to<arma::rowvec>::from()
  arma::rowvec rowweightEndpoint = weightEndpoint.t();

  // store iid average plus iid nuisance
  arma::mat iidTotal_favorable;
  arma::mat iidTotal_unfavorable;
  arma::mat iidTotal_neutral;

  if(returnIID[0]>0 || pool>=4){

    // *** scale to move from a sum to an expectation (E[s^{\gamma}_{l,j}|x_l])
    if(pool!=3){
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	iidAverage_favorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
	iidAverage_favorable.rows(posT[iter_strata]) /= n_control[iter_strata];
	iidAverage_unfavorable.rows(posC[iter_strata]) /= n_treatment[iter_strata];
	iidAverage_unfavorable.rows(posT[iter_strata]) /= n_control[iter_strata];
	iidAverage_neutral.rows(posC[iter_strata]) /= n_treatment[iter_strata];
	iidAverage_neutral.rows(posT[iter_strata]) /= n_control[iter_strata];
      }
    }
    
    // *** center to obtain (endpoint specific) first order projection (h_1^{\gamma}(l) = E[s^{\gamma}_{l,j}|x_l]-U^{\gamma})
    arma::uvec ivecD(1);
    
    for(unsigned int iter_d=0 ; iter_d < D ; iter_d ++){
      ivecD[0] = iter_d;

      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
	if(pool==3){
	  // estimate: \Delta = \sum_{s,s*) w(s,s*) \Delta(s,s^*)
	  // Hajek projection for individual i belonging to the control group with m observations
	  // h_1(i) = \sum_{s,s*) w(s,s*) h_1(i,s,s^*) = \sum_{s,s*) w(s,s*) (1/n_s* \sum_{j=1}^n 1[Y_j>X_i] 1[s_i=s,s_j=s*] - \Delta(s,s^*)) / m_s
	  //        = \sum_{s,s*) w(s,s*)/(n_s* m_s) \sum_{j=1}^n 1[Y_j>X_i] 1[s_i=s,s_j=s*] - \sum_{s,s*) w(s,s*) \Delta(s,s^*) / m_s
	  // here we substract the second term, i.e., w(s,s*) \Delta(s,s^*) / m_s      
	  iter_strataC = grid_strata(iter_strata,0); // index of the strata (may differ between treatment groups when using standardization)
	  iter_strataT = grid_strata(iter_strata,1); // 
	 
	  iidAverage_favorable(posC[iter_strataC],ivecD) -= (weightPool[iter_strata]/n_control[iter_strata]) * delta.slice(0)(iter_strata,iter_d);
	  iidAverage_favorable(posT[iter_strataT],ivecD) -= (weightPool[iter_strata]/n_treatment[iter_strata]) * delta.slice(0)(iter_strata,iter_d);
      
	  iidAverage_unfavorable(posC[iter_strataC],ivecD) -= (weightPool[iter_strata]/n_control[iter_strata]) * delta.slice(1)(iter_strata,iter_d);
	  iidAverage_unfavorable(posT[iter_strataT],ivecD) -= (weightPool[iter_strata]/n_treatment[iter_strata]) * delta.slice(1)(iter_strata,iter_d);

	  iidAverage_neutral(posC[iter_strataC],ivecD) -= (weightPool[iter_strata]/n_control[iter_strata]) * delta.slice(2)(iter_strata,iter_d);
	  iidAverage_neutral(posT[iter_strataT],ivecD) -= (weightPool[iter_strata]/n_treatment[iter_strata]) * delta.slice(2)(iter_strata,iter_d);
	}else{
	  iidAverage_favorable(posC[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];
	  iidAverage_favorable(posT[iter_strata],ivecD) -= delta.slice(0)(iter_strata,iter_d); // same as -= Mcount_favorable(iter_strata,iter_d)/n_pairs[iter_strata];      
	  iidAverage_unfavorable(posC[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];
	  iidAverage_unfavorable(posT[iter_strata],ivecD) -= delta.slice(1)(iter_strata,iter_d); // same as -= Mcount_unfavorable(iter_strata,iter_d)/n_pairs[iter_strata];
	  iidAverage_neutral(posC[iter_strata],ivecD) -= delta.slice(2)(iter_strata,iter_d); // same as -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
	  iidAverage_neutral(posT[iter_strata],ivecD) -= delta.slice(2)(iter_strata,iter_d); // same as -= Mcount_neutral(iter_strata,iter_d)/n_pairs[iter_strata];
	}
      }
    }
  
    // *** rescale such that the sum of squares of the projection equals the variance, i.e. obtain h_1^{\gamma}(l)/m
    if(pool!=3){
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
	iidAverage_favorable.rows(posC[iter_strata]) /= n_control[iter_strata];
	iidAverage_favorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];
	iidAverage_unfavorable.rows(posC[iter_strata]) /= n_control[iter_strata];
	iidAverage_unfavorable.rows(posT[iter_strata]) /= n_treatment[iter_strata];
	iidAverage_neutral.rows(posC[iter_strata]) /= n_control[iter_strata];
	iidAverage_neutral.rows(posT[iter_strata]) /= n_treatment[iter_strata];
      }
    }

    // *** modification for paired design
    if(paired){ // sum within pair so the influence function is at the cluster level

      iidAverage_favorable.fill(0.0);
      iidAverage_unfavorable.fill(0.0);
      iidAverage_neutral.fill(0.0);
      arma::mat Wdelta(3,D); // compute endpoint specific pooled estimator
      Wdelta.row(0) = weightPool.t() * delta.slice(0);
      Wdelta.row(1) = weightPool.t() * delta.slice(1);
      Wdelta.row(2) = weightPool.t() * delta.slice(2);
      double iidNuisance_tempo;

      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){

	// each strata contains a single C and a single T
	iidAverage_favorable.row(posC[iter_strata][0]) = delta.slice(0).row(iter_strata) - Wdelta.row(0);
	iidAverage_unfavorable.row(posC[iter_strata][0]) = delta.slice(1).row(iter_strata) - Wdelta.row(1);
	iidAverage_neutral.row(posC[iter_strata][0]) = delta.slice(2).row(iter_strata) - Wdelta.row(2);
	if(hprojection==2){
	  // small sample correction: later division by n_strata: n_strata-1 = sqrt(1-1/n_strata)^2*n_strata	
	  iidAverage_favorable.row(posC[iter_strata][0]) /= sqrt(1.0-1.0/n_strata);
	  iidAverage_unfavorable.row(posC[iter_strata][0]) /= sqrt(1.0-1.0/n_strata);
	  iidAverage_neutral.row(posC[iter_strata][0]) /= sqrt(1.0-1.0/n_strata);
	}

	if(returnIID[1]>0){
	  iidNuisance_tempo = accu(iidNuisance_favorable.rows(posC[iter_strata])) + accu(iidNuisance_favorable.rows(posT[iter_strata]));
	  iidNuisance_favorable.rows(posC[iter_strata]).fill(0.0);
	  iidNuisance_favorable.rows(posT[iter_strata]).fill(0.0);
	  iidNuisance_favorable.row(posC[iter_strata][0]) = iidNuisance_tempo;
	    
	  iidNuisance_tempo = accu(iidNuisance_unfavorable.rows(posC[iter_strata])) + accu(iidNuisance_unfavorable.rows(posT[iter_strata]));
	  iidNuisance_unfavorable.rows(posC[iter_strata]).fill(0.0);
	  iidNuisance_unfavorable.rows(posT[iter_strata]).fill(0.0);
	  iidNuisance_unfavorable.row(posC[iter_strata][0]) = iidNuisance_tempo;

	  iidNuisance_tempo = accu(iidNuisance_neutral.rows(posC[iter_strata])) + accu(iidNuisance_neutral.rows(posT[iter_strata]));
	  iidNuisance_neutral.rows(posC[iter_strata]).fill(0.0);
	  iidNuisance_neutral.rows(posT[iter_strata]).fill(0.0);
	  iidNuisance_neutral.row(posC[iter_strata][0]) = iidNuisance_tempo;
	}
      }
    }

    // *** add the two source of uncertainty into a single iid
    if(returnIID[1]>0){
      iidTotal_favorable = iidAverage_favorable + iidNuisance_favorable;
      iidTotal_unfavorable = iidAverage_unfavorable + iidNuisance_unfavorable;
      iidTotal_neutral = iidAverage_neutral + iidNuisance_neutral;
    }else{
      iidTotal_favorable = iidAverage_favorable;
      iidTotal_unfavorable = iidAverage_unfavorable;
      iidTotal_neutral = iidAverage_neutral;
    }

    // *** weight endpoints and cumulate them to obtain (cumulative) first order projection
    // logically should be in the next subection (Global) but is put here because used for the point estimate with pooling = 4.1-4.4
    iidTotal_favorable.each_row() %= rowweightEndpoint;
    iidTotal_favorable = arma::cumsum(iidTotal_favorable,1);

    iidTotal_unfavorable.each_row() %= rowweightEndpoint;
    iidTotal_unfavorable = arma::cumsum(iidTotal_unfavorable,1);

    if(addHalfNeutral){ // only add neutral relative to the latest endpoint for each priority
      iidTotal_neutral.each_row() %= rowweightEndpoint;
      iidTotal_favorable += 0.5*iidTotal_neutral;
      iidTotal_unfavorable += 0.5*iidTotal_neutral;
    }
  }

  // ** optimal weights to pool across strata
  if(pool >= 4){ // usual weights defined earlier
    arma::vec iIID;
    arma::uvec iPos;
    arma::uvec iUvec = {D-1};
    double iDeltaFA, iDeltaUF;
    for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
      iPos = arma::join_cols(posC[iter_strata], posT[iter_strata]);
      if(pool==4.1){
	iIID = iidTotal_favorable(iPos,iUvec);
      }else if(pool==4.2){
	iIID = iidTotal_unfavorable(iPos,iUvec);
      }else if(pool==4.3){
	iIID = iidTotal_favorable(iPos,iUvec) - iidTotal_unfavorable(iPos,iUvec);
      }else if(pool==4.4){
	iDeltaFA = sum(delta.slice(0).row(iter_strata) % rowweightEndpoint);
	iDeltaUF = sum(delta.slice(1).row(iter_strata) % rowweightEndpoint);
	if(addHalfNeutral){  // only add neutral relative to the latest endpoint for each priority
	  iDeltaFA += 0.5*delta.slice(2)(iter_strata,D-1) * rowweightEndpoint[D-1];
	  iDeltaUF += 0.5*delta.slice(2)(iter_strata,D-1) * rowweightEndpoint[D-1];
	}
	iIID = iidTotal_favorable(iPos,iUvec)/iDeltaUF - iidTotal_unfavorable(iPos,iUvec)/(pow(iDeltaUF,2)/iDeltaFA);
      }
      weightPool[iter_strata] = 1/sum(iIID % iIID);      
    }
    weightPool /= sum(weightPool); // normalize weights to sum up to 1
  }
  

  // ** Summary statistics (probabilistic index, net benefit, and win ratio)

  // *** strata-weighted cumulative endpoint estimate
  arma::mat cumdelta_favorable = arma::cumsum(delta.slice(0).each_row() % rowweightEndpoint,1);
  arma::mat wcumdelta_favorable = cumdelta_favorable.each_col() % weightPool; 
  
  arma::mat cumdelta_unfavorable = arma::cumsum(delta.slice(1).each_row() % rowweightEndpoint,1);
  arma::mat wcumdelta_unfavorable = cumdelta_unfavorable.each_col() % weightPool; 

  arma::mat cumdelta_neutral = delta.slice(2).each_row() % rowweightEndpoint;
  arma::mat wcumdelta_neutral = cumdelta_neutral.each_col() % weightPool; 

  arma::mat cumdelta_uninf = delta.slice(3).each_row() % rowweightEndpoint;
  arma::mat wcumdelta_uninf = cumdelta_uninf.each_col() % weightPool; 

  if(addHalfNeutral){  // only add neutral relative to the latest endpoint for each priority
    delta.slice(0) += 0.5*delta.slice(2);
    delta.slice(1) += 0.5*delta.slice(2);
    if(returnIID[0]>0 || pool>=4){
      iidAverage_favorable += 0.5*iidAverage_neutral;
      iidAverage_unfavorable += 0.5*iidAverage_neutral;
    }
    if(returnIID[1]>0){
      iidNuisance_favorable += 0.5*iidNuisance_neutral;
      iidNuisance_unfavorable += 0.5*iidNuisance_neutral;
    }
    cumdelta_favorable += 0.5*cumdelta_neutral;
    cumdelta_unfavorable += 0.5*cumdelta_neutral;
    wcumdelta_favorable += 0.5*wcumdelta_neutral;
    wcumdelta_unfavorable += 0.5*wcumdelta_neutral;
  }

  // *** Point estimates
  // Probabilistic index (\Delta = \sum_endpoint \sum_strata w(strata) w(endpoint) \delta(strata, endpoint))
  Delta.col(0) = arma::trans(arma::sum(wcumdelta_favorable,0));
  Delta.col(1) = arma::trans(arma::sum(wcumdelta_unfavorable,0));
  Delta.col(2) = arma::trans(arma::sum(wcumdelta_neutral,0));
  Delta.col(3) = arma::trans(arma::sum(wcumdelta_uninf,0));
  
  // net benefit 
  delta.slice(4) = delta.slice(0) - delta.slice(1);
  Delta.col(4) = Delta.col(0) - Delta.col(1);

  // win ratio equals number of favorable pairs divided by the number of favorable plus unfavorable pairs
  delta.slice(5) = delta.slice(0) / delta.slice(1);
  if(pool==4.4){
    arma::mat winRatioStrata = cumdelta_favorable/cumdelta_unfavorable;
    Delta.col(5) = arma::trans(arma::sum(winRatioStrata.each_col() % weightPool,0));
  }else{
    Delta.col(5) = Delta.col(0) / Delta.col(1);
  }

  // ** Variance and iid estimation
  if(returnIID[0] > 0){

    if(pool!=3){
      // *** rescale to account for the pooling across strata i.e. T = \sum_s n_pair(s) T(s) / n_tot
      for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ 
	iidTotal_favorable.rows(posC[iter_strata]) *= weightPool[iter_strata];
	iidTotal_favorable.rows(posT[iter_strata]) *= weightPool[iter_strata];
	
	iidTotal_unfavorable.rows(posC[iter_strata]) *= weightPool[iter_strata];
	iidTotal_unfavorable.rows(posT[iter_strata]) *= weightPool[iter_strata];
      }
    }
    
    // *** first order variance 
    arma::mat iid2;
    iid2 = pow(iidTotal_favorable,2);
    Mvar.col(0) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));
    iid2 = pow(iidTotal_unfavorable,2);
    Mvar.col(1) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));
    iid2 = iidTotal_favorable % iidTotal_unfavorable;
    Mvar.col(2) = arma::trans(arma::sum(iid2.each_col() % weightObs, 0));

    // *** second order variance
    if(hprojection==2 && (paired==false)){
      if(keepScore){

	// **** retrive h1
	arma::mat h1_favorable; 
	arma::mat h1_unfavorable;
	double iter_weight, iter_weight2=1;

	h1_favorable = iidTotal_favorable;
	h1_unfavorable = iidTotal_unfavorable;
	for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
	  if(pool==3){ // standardization
	    iter_strataC = grid_strata(iter_strata,0); // index of the strata (may differ between treatment groups when using standardization)
	    iter_strataT = grid_strata(iter_strata,1); // 
	    if(iter_strataC == iter_strataT){ 
	      h1_favorable.rows(posC[iter_strataC]) *= ntot_control;
	      h1_favorable.rows(posT[iter_strataT]) *= ntot_treatment;

	      h1_unfavorable.rows(posC[iter_strataC]) *= ntot_control;
	      h1_unfavorable.rows(posT[iter_strataT]) *= ntot_treatment;
	    }
	  }else{
	    h1_favorable.rows(posC[iter_strata]) *= (n_control[iter_strata]/weightPool[iter_strata]);
	    h1_favorable.rows(posT[iter_strata]) *= (n_treatment[iter_strata]/weightPool[iter_strata]);

	    h1_unfavorable.rows(posC[iter_strata]) *= (n_control[iter_strata]/weightPool[iter_strata]);
	    h1_unfavorable.rows(posT[iter_strata]) *= (n_treatment[iter_strata]/weightPool[iter_strata]);
	  }
	}

	// **** reconstruct individual score
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

	// **** compute second order h-projection, its variance and covariance
	// NOTE: when computing h we use iid = E[]-U so h2 = s - E[] - E[] + U = s - iid - iid - U
	arma::rowvec H2_favorable, H2_unfavorable;
	arma::mat H2_moments(5,D,arma::fill::zeros);
	int iter_strata, iter_C, iter_T;
	arma::rowvec strataDelta_favorable(D);
	arma::rowvec strataDelta_unfavorable(D);
	
	// endpoint specific    
	for(int iter_pair=0; iter_pair<ntot_pair ; iter_pair++){

	  iter_strata = lsScore[0](iter_pair,0);
	  iter_strataC = grid_strata(iter_strata,0); // index of the strata (may differ between treatment groups when using standardization)
	  iter_strataT = grid_strata(iter_strata,1); //
	  
	  iter_C = lsScore[0](iter_pair,4); // index within strata
	  iter_T = lsScore[0](iter_pair,5); // index within strata

	  iter_weight = weightObs(posC[iter_strataC][iter_C])*weightObs(posT[iter_strataT][iter_T]);
	  if(pool==3){
	    iter_weight2 = weightPool[iter_strata]/((n_control[iter_strata]*n_treatment[iter_strata])/(ntot_control*ntot_treatment));
	    strataDelta_favorable = arma::trans(Delta.col(0));
	    strataDelta_unfavorable = arma::trans(Delta.col(1));
	  }else{
	    strataDelta_favorable = cumdelta_favorable.row(iter_strata);
	    strataDelta_unfavorable = cumdelta_unfavorable.row(iter_strata);
	  }
	  H2_favorable = (iter_weight2 * pairScoreF.row(iter_pair) - h1_favorable.row(posC[iter_strataC][iter_C]) - h1_favorable.row(posT[iter_strataT][iter_T]) - strataDelta_favorable);
	  H2_unfavorable = (iter_weight2 * pairScoreUF.row(iter_pair) - h1_unfavorable.row(posC[iter_strataC][iter_C]) - h1_unfavorable.row(posT[iter_strataT][iter_T]) - strataDelta_unfavorable);
	  // Rcpp::Rcout << iter_pair << " (" << iter_strata << ":" << iter_C << " " << iter_T << ") " << H2_favorable[0] << " " << H2_unfavorable[0] << std::endl;

	  // var(H2) = \sum_ij H2 / n_pair[iter_strata]
	  // var(1/nm \sum_ij H2) = \sum_ij H2 / n_pair[iter_strata]^2
	  // w_{pooling}^2 var(1/nm \sum_ij H2) = (n_pairs[iter_strata]/ntot_pair)^2 \sum_ij H2 / n_pairs[iter_strata]^2 = \sum_ij H2/ntot_pair^2
	  if(pool==3){
	    Mvar.col(0) += iter_weight * arma::trans(pow(H2_favorable,2))/pow(ntot_pair,2);
	    Mvar.col(1) += iter_weight * arma::trans(pow(H2_unfavorable,2))/pow(ntot_pair,2);
	    Mvar.col(2) += iter_weight * arma::trans(H2_favorable % H2_unfavorable)/pow(ntot_pair,2);
	  }else{
	    Mvar.col(0) += iter_weight * arma::trans(pow(H2_favorable,2)) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	    Mvar.col(1) += iter_weight * arma::trans(pow(H2_unfavorable,2)) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	    Mvar.col(2) += iter_weight * arma::trans(H2_favorable % H2_unfavorable) * pow(weightPool[iter_strata],2)/pow(n_pairs[iter_strata],2);
	  }
	}
	
      }else{ // only ok for binary scores i.e. win neutral or loss
	arma::colvec strataDelta_favorable(D);
	arma::colvec strataDelta_unfavorable(D);
	arma::vec strataWeightC, strataWeightT;
	double iN_pairs, iN_treatment, iN_control, iWeightPool;

	for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
	  
	  iter_strataC = grid_strata(iter_strata,0); // index of the strata (may differ between treatment groups when using standardization)
	  iter_strataT = grid_strata(iter_strata,1); //
	    
	  // strata specific cumulative endpoint
	  if(pool==3){ // standardization
	    if(iter_strataC == iter_strataT){
	      if(iter_strataC == 0){ // add global estimate only once
	        strataDelta_favorable = Delta.col(0);	      
	        strataDelta_unfavorable = Delta.col(1);
	      }else{
		strataDelta_favorable.fill(0.0);	      
		strataDelta_unfavorable.fill(0.0);
	      }
	      iN_pairs = ntot_pair;
	      iN_treatment = ntot_treatment;
	      iN_control = ntot_control;
	      iWeightPool = 1;
	    }else{ // ignore 'artifical' strata (M.F. or F.M.)
	      continue;
	    }
	  }else{
	    strataDelta_favorable = arma::trans(cumdelta_favorable.row(iter_strata));
	    strataDelta_unfavorable = arma::trans(cumdelta_unfavorable.row(iter_strata));
	    iN_pairs = n_pairs[iter_strata];
	    iN_treatment = n_treatment[iter_strata];
	    iN_control = n_control[iter_strata];
	    iWeightPool = weightPool[iter_strata];
	  }
	  strataWeightC = weightObs(posC[iter_strataC]);
	  strataWeightT = weightObs(posT[iter_strataT]);

    
	  // (n_pairs[iter_strata] / ntot_pair)^2 * (1/n_pairs[iter_strata]) = n_pairs[iter_strata]/ntot_pair^2
	  Mvar.col(0) += strataDelta_favorable % (1-strataDelta_favorable) * pow(iWeightPool,2)/ iN_pairs;
	  iid2 = pow(iidTotal_favorable.rows(posC[iter_strataC]),2);
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / iN_treatment );
	  iid2 = pow(iidTotal_favorable.rows(posT[iter_strataT]),2);
	  Mvar.col(0) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / iN_control );
	  
	  Mvar.col(1) += strataDelta_unfavorable % (1-strataDelta_unfavorable) * pow(iWeightPool,2)/ iN_pairs;
	  iid2 = pow(iidTotal_unfavorable.rows(posC[iter_strataC]),2);
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / iN_treatment );
	  iid2 = pow(iidTotal_unfavorable.rows(posT[iter_strataT]),2);
	  Mvar.col(1) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / iN_control );
	  
	  Mvar.col(2) -= strataDelta_favorable % strataDelta_unfavorable * pow(iWeightPool,2)/ iN_pairs;
	  iid2 = iidTotal_favorable.rows(posC[iter_strataC]) % iidTotal_unfavorable.rows(posC[iter_strataC]);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightC, 0) / iN_treatment);
	  iid2 = iidTotal_favorable.rows(posT[iter_strataT]) % iidTotal_unfavorable.rows(posT[iter_strataT]);
	  Mvar.col(2) -= arma::conv_to<arma::vec>::from( arma::sum(iid2.each_col() % strataWeightT, 0) / iN_control);
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
