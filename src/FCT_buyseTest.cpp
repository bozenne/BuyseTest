// * Preambule

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "FCT_calcOnePair.h"
#include "FCT_calcAllPairs.h"
#include "FCT_calcStatistic.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// * Documentation GPC_cpp
//' @title C++ function performing the pairwise comparison over several endpoints. 
//' @description \code{GPC_cpp} call for each endpoint and each strata the pairwise comparison function suited to the type of endpoint and store the results. 
//' @name GPC_cpp
//' 
//' @param Control A matrix containing the values of each endpoint (in columns) for the control group (in rows). \emph{const arma::mat&}.
//' @param Treatment A matrix containing the values of each endpoint (in columns) for the treatment group (in rows). \emph{const arma::mat&}.
//' @param threshold Store the thresholds associated to each endpoint. \emph{const NumericVector&}. Must have length D. The threshold is ignored for binary endpoints. Must have D columns.
//' @param survEndpoint Does each endpoint is a time to event. \emph{const LogicalVector&}. Must have length D.
//' @param delta_Control A matrix containing the nature of observations in the control group (in rows) (0 censoring, 1 event) for each TTE endpoint (in columns) . \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.
//' @param delta_Treatment A matrix containing in the type of event (0 censoring, 1 event) for each TTE endpoint (in columns) and treatment observations (in rows). \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.
//' @param D The number of endpoints. Strictly positive \emph{const int}.
//' @param strataC A list containing the indexes of control observations belonging for each strata. \emph{List&}  of vector containing positive integers. 
//' @param strataT A list containing the indexes of treatment observations belonging for each strata. \emph{List&} of vector containing positive integers. 
//' @param n_strata The number of strata . Strictly positive \emph{const int}.
//' @param n_TTE The number of time-to-event endpoints. Positive \emph{const int}.
//' @param Wscheme The matrix describing the weighting strategy. For each endpoint (except the first) in column, weights of each pair are initialized at 1 and multiplied by the weight of the endpoints in rows where there is a 1. \emph{const arma::mat&}. Must have D lines and D columns.
//' @param index_survivalM1 The position, among all the survival endpoints, of the last same endpoint (computed with a different threshold). If it is the first time that the TTE endpoint is used it is set to -1. \emph{const IntegerVector}. Must have length n_TTE.
//' @param threshold_M1 The previous latest threshold of each TTE endpoint. When it is the first time that the TTE endpoint is used it is set to -1. \emph{const NumericVector}. Must have length n_TTE.
//' @param list_survTimeC A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows). \emph{List&}.
//' @param list_survTimeT A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows). \emph{List&}. 
//' @param list_survJumpC A list of matrix containing the survival estimates and survival jumps when the survival for the control arm jumps. \emph{List&}.
//' @param list_survJumpT A list of matrix containing the survival estimates and survival jumps when the survival for the treatment arm jumps. \emph{List&}. 
//' @param list_lastSurv A list of matrix containing the last survival estimate in each strata (rows) and treatment group (columns). \emph{List&}. 
//' @param correctionUninf Should the uninformative weight be re-distributed to favorable and unfavorable?
//' @param methodTTE The type of method used to compare censored pairs (0 Gehan 1 Peron).
//' @param neutralAsUninf Should paired classified as neutral be re-analyzed using endpoints of lower priority?  \emph{logical}.
//' @param keepScore Should the result of each pairwise comparison be kept? \emph{logical}.
//' @keywords function Cpp BuyseTest

// * Function GPC_cpp
//' @name GPC_cpp
//' @export
// [[Rcpp::export]]
List GPC_cpp(const arma::mat& Control,
			 const arma::mat& Treatment,
			 const NumericVector& threshold,
			 const LogicalVector& survEndpoint,
			 const arma::mat& delta_Control,
             const arma::mat& delta_Treatment,
			 const int D,
			 const std::vector< arma::uvec >& strataC,
			 const std::vector< arma::uvec >& strataT,
			 const int n_strata,
			 const int n_TTE, 
             const arma::mat& Wscheme,
			 const IntegerVector index_survivalM1,
			 const NumericVector threshold_M1, 
			 const std::vector< std::vector< arma::mat > >& list_survTimeC,
             const std::vector< std::vector< arma::mat > >& list_survTimeT,
             const std::vector< std::vector< arma::mat > >& list_survJumpC,
			 const std::vector< std::vector< arma::mat > >& list_survJumpT,
			 const std::vector< arma::mat >& list_lastSurv,
			 const int methodTTE,
			 const int correctionUninf,
			 const bool neutralAsUninf,
			 const bool keepScore){

  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.
  
  /// ** initialization
  
  // *** final results
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs by outcome for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs by outcome for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs by outcome for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs by outcome for each strata
  vector<double> n_pairs(n_strata); // number of pairs sumed over the strats
  
  // *** for a given stata [input]
  // position of a patients belonging to a given strata
  arma::uvec index_strataC; // in the control arm
  arma::uvec index_strataT; // in the treatment arm

  // Endpoint(s) restricted to strata k
  arma::mat ControlK; // in the control arm
  arma::mat TreatmentK; // in the treatment arm

  // Right-censoring status for TTE endpoint restricted to strata k
  arma::mat delta_ControlK ; // in the control arm
  arma::mat delta_TreatmentK ; // in the treatment arm

  // Survival values restricted to strata k
  arma::mat iSurvTimeC_M1; // at times corresponding to event in the control group
  arma::mat iSurvTimeT_M1; // at times corresponding to event in the treatment group
  arma::mat iSurvJumpC_M1; // at jump times for the survival in the control group
  arma::mat iSurvJumpT_M1; // at jump times for the survival in the treatment group

  // index of the neutral pairs
  vector<int> iIndex_neutralC; // in the control arm
  vector<int> iIndex_neutralT;  // in the treatment arm
  
  // index of the uninformative pairs
  vector<int> iIndex_uninfC; // in the control arm
  vector<int> iIndex_uninfT; // in the treatment arm

  // weights of the neutral / uninformative pairs
  arma::mat matWeight;  // for all endpoint up to the current endpoint 
  arma::mat matWeight_M1; // for all endpoint up to the previous endpoint
  vector<double> iWeight; // only for the current endpoint

  // index of the pairs corresponding to the new weights (iWeight) in the weight vector of the previous endpoint (iCumWeight_M1)
  vector<int> iIndexWeight_pair; 

  // product of the weights of the neutral / uninformative pairs up to the previous endpoint
  // NOTE: it is a special product because weights related to previous survival endpoint are ignored
  arma::vec iCumWeight_M1;

  // others
  int iter_d; // the index of the endpoints
  int iter_dTTE; // number of time to event endpoints that have been used
  int iter_oldpair; // index of the remaining pair [local use when updating iMweights]
    
  int size_neutral; // number of neutral pairs (temporary)
  int size_uninf; // number of uninformative pairs (temporary)

  vector<int> tempo_index;  // temporary index vector used to convert the index SEXP into vector<int>
  
  // *** keep track of each comparison
  vector<arma::mat> lsScore(D);
  arma::mat iScore;
  arma::mat iMat;
  int iNpairs;
  bool iMoreEndpoint;
      
  // ** loop over strata
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
    
    iter_d = 0;
    iter_dTTE = 0;
    
    index_strataT = strataT[iter_strata];
    index_strataC = strataC[iter_strata];
    
    // *** affect the endpoints corresponding to each strata (conversion from SEXP to uvec object is necessary for .rows)
    TreatmentK = Treatment.rows(index_strataT); // select the rows in the treatment matrix corresponding to the indexes contained in strata
    ControlK = Control.rows(index_strataC); // select the rows in the control matrix corresponding to the indexes contained in strata
    
    if(n_TTE>0){
      delta_TreatmentK = delta_Treatment.rows(index_strataT); // select the rows in the treatment status matrix corresponding to the indexes contained in strata
      delta_ControlK = delta_Control.rows(index_strataC); // select the rows in the control status matrix corresponding to the indexes contained in strata          
    }
    
    // *** first endpoint
    // Rcout << "** endpoint 0 **" << endl;
    iIndex_neutralT.resize(0); iIndex_neutralC.resize(0); iIndex_uninfT.resize(0); iIndex_uninfC.resize(0); iWeight.resize(0);
    iMoreEndpoint = (D>1);
	
    if(survEndpoint[0]){ // time to event endpoint
      
      if(methodTTE == 0){
		iScore = calcAllPairs_TTEgehan(ControlK.col(0), TreatmentK.col(0), threshold[0],
									   delta_ControlK.col(0), delta_TreatmentK.col(0), 
									   correctionUninf,
									   Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
									   iIndex_neutralC, iIndex_neutralT, iIndex_uninfC, iIndex_uninfT, 
									   iWeight,
									   neutralAsUninf, keepScore, iMoreEndpoint); 
      }else{

		iScore = calcAllPairs_TTEperon(ControlK.col(0), TreatmentK.col(0), threshold[0],
									   delta_ControlK.col(0), delta_TreatmentK.col(0),
									   list_survTimeC[0][iter_strata], list_survTimeT[0][iter_strata], list_survJumpC[0][iter_strata], list_survJumpT[0][iter_strata],
									   list_lastSurv[0](iter_strata,0), list_lastSurv[0](iter_strata,1), 
									   correctionUninf,
									   Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
									   iIndex_neutralC, iIndex_neutralT, iIndex_uninfC, iIndex_uninfT, 
									   iWeight,
									   neutralAsUninf, keepScore, iMoreEndpoint);
      }
      iter_dTTE++; // increment the number of time to event endpoints that have been used
      
    }else { // binary or continuous endpoint
      
      iScore = calcAllPairs_Continuous(ControlK.col(0), TreatmentK.col(0), threshold[0],
									   correctionUninf,
									   Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
									   iIndex_neutralC, iIndex_neutralT, iIndex_uninfC, iIndex_uninfT, 
									   iWeight,
									   neutralAsUninf, keepScore, iMoreEndpoint);
      
    }

    // **** add to the total number of pairs the number of pairs founded for this endpoint
    n_pairs[iter_strata] = Mcount_favorable(iter_strata,0) + Mcount_unfavorable(iter_strata,0) + Mcount_neutral(iter_strata,0) + Mcount_uninf(iter_strata,0);
    size_neutral = iIndex_neutralT.size(); // update the number of neutral pairs
    size_uninf = iIndex_uninfT.size(); // update the number of uninformative pairs
    
	// Rcout << "update weights" << endl;
	
	// **** update weights associated to the remaing pairs
	// only relevant if there is one more endpoint and the weights may differ from 1
    if(D>1 && (methodTTE>0 || correctionUninf>0)){ 
      matWeight.resize(size_neutral+size_uninf,1); 

	  // set the value of the weights to the one computed and stored in iWeight
	  for(int iter_neutral=0 ; iter_neutral<size_neutral ; iter_neutral++){
		matWeight(iter_neutral,0) = iWeight[iter_neutral];
	  }
	  for(int iter_uninf=0 ; iter_uninf<size_uninf ; iter_uninf++){
		matWeight(size_neutral+iter_uninf,0) = iWeight[size_neutral + iter_uninf];
	  }
    }

	// Rcout << "update Scores" << endl;
    // **** update all Scores
    if(keepScore){
      iNpairs = iScore.n_rows;
      iMat.resize(iNpairs,3);
      // add original index
      for(int iPair=0 ; iPair < iNpairs ; iPair ++){
		iMat.row(iPair) = rowvec({(double)iter_strata,
			  (double)index_strataC(iScore(iPair,0)),
			  (double)index_strataT(iScore(iPair,1))});
      }
      // merge with current table and store
      if(iter_strata==0){
		lsScore[0] = arma::join_rows(iMat,iScore);
      }else{
        lsScore[0] = arma::join_cols(lsScore[0],arma::join_rows(iMat,iScore));
      }
    }

    // *** following endpoints
    while(D>iter_d+1 && (Mcount_neutral(iter_strata,iter_d)>0 || Mcount_uninf(iter_strata,iter_d)>0)){ // loop over the following endpoints

      // while there are remaining endpoints and remaining neutral or uniformative pairs
      iter_d++; // increment the index of the endpoints
      // Rcout << "** endpoint " << iter_d << " **" << endl;
	  matWeight_M1 = matWeight; // save the current iMweight
      iWeight.resize(0); iIndexWeight_pair.resize(0);
      iMoreEndpoint = (D>(iter_d+1));

	  // **** compute the current weights of the pairs
	  // initialize iCumWeight_M1
	  iCumWeight_M1.resize(size_neutral+size_uninf);
	  iCumWeight_M1.fill(1.0);

	  if(methodTTE>0 || correctionUninf>0){
	  	for(int iter_endpoint=0 ; iter_endpoint<iter_d ; iter_endpoint++){
		  if(Wscheme(iter_endpoint,iter_d)==1){iCumWeight_M1 %= matWeight.col(iter_endpoint);}
		}
	  }
		
	  // **** computes scores
      if(survEndpoint[iter_d]){ // time to event endpoint

		if(methodTTE==0){
		  iScore = calcSubsetPairs_TTEgehan(ControlK.col(iter_d), TreatmentK.col(iter_d), threshold[iter_d],
											delta_ControlK.col(iter_dTTE), delta_TreatmentK.col(iter_dTTE),
											correctionUninf,
											Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d),
											Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
											iIndex_neutralC, iIndex_neutralT, size_neutral,
											iIndex_uninfC, iIndex_uninfT, size_uninf,
											iCumWeight_M1, iWeight, iIndexWeight_pair, // Wpairs wNeutral index_wNeutral
											neutralAsUninf, keepScore, iMoreEndpoint);
		}else{
        
		  if(threshold_M1[iter_dTTE]<0){ // first time to event endpoint (i.e. no previous threshold)
			iSurvTimeC_M1 = arma::mat(0,0);
			iSurvTimeT_M1 = arma::mat(0,0);
			iSurvJumpC_M1 = arma::mat(0,0);
			iSurvJumpT_M1 = arma::mat(0,0);

		  }else{ // following times
			iSurvTimeC_M1 = list_survTimeC[index_survivalM1[iter_dTTE]][iter_strata];
			iSurvTimeT_M1 = list_survTimeT[index_survivalM1[iter_dTTE]][iter_strata];
			iSurvJumpC_M1 = list_survJumpC[index_survivalM1[iter_dTTE]][iter_strata];
			iSurvJumpT_M1 = list_survJumpT[index_survivalM1[iter_dTTE]][iter_strata];
		  }
	
          iScore = calcSubsetPairs_TTEperon(ControlK.col(iter_d), TreatmentK.col(iter_d), threshold[iter_d],
											delta_ControlK.col(iter_dTTE), delta_TreatmentK.col(iter_dTTE),
											list_survTimeC[iter_dTTE][iter_strata], list_survTimeT[iter_dTTE][iter_strata],
											list_survJumpC[iter_dTTE][iter_strata], list_survJumpT[iter_dTTE][iter_strata],
											list_lastSurv[iter_dTTE](iter_strata,0), list_lastSurv[iter_dTTE](iter_strata,1), 
											correctionUninf,
											Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d),
											Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
											iIndex_neutralC, iIndex_neutralT, size_neutral,
											iIndex_uninfC, iIndex_uninfT, size_uninf,
											iCumWeight_M1, threshold_M1[iter_dTTE],
											iSurvTimeC_M1, iSurvTimeT_M1, iSurvJumpC_M1, iSurvJumpT_M1,
											iWeight, iIndexWeight_pair,
											neutralAsUninf, keepScore, iMoreEndpoint);
		}

		iter_dTTE++; // increment the number of time to event endpoints that have been used
	
      }else{ // binary or continuous endpoint
	
        iScore = calcSubsetPairs_Continuous(ControlK.col(iter_d), TreatmentK.col(iter_d), threshold[iter_d],
											correctionUninf,
											Mcount_favorable(iter_strata,iter_d),
											Mcount_unfavorable(iter_strata,iter_d),
											Mcount_neutral(iter_strata,iter_d),
											Mcount_uninf(iter_strata,iter_d),
											iIndex_neutralC, iIndex_neutralT, size_neutral,
											iIndex_uninfC, iIndex_uninfT, size_uninf,
											iCumWeight_M1, iWeight, iIndexWeight_pair,
											neutralAsUninf, keepScore, iMoreEndpoint);
	
      }

      size_neutral = iIndex_neutralT.size(); // update the number of neutral pairs
      size_uninf = iIndex_uninfT.size(); // update the number of uninformative pairs

	  // Rcout << "update weights" << endl;
      // **** update weights associated to the remaing pairs
	  // only relevant if there is one more endpoint and the weights may differ from 1
      if(D>iter_d+1 && (methodTTE>0 || correctionUninf>0)){
        matWeight.resize(size_neutral+size_uninf,iter_d+1); // update the size of iMweight

		tempo_index=iIndexWeight_pair; // store the position of the remaining pairs in the previous iMweight (i.e. iMweight_M1)
		  
		for(size_t iter_pair=0; iter_pair<tempo_index.size(); iter_pair++){
			
		  iter_oldpair = tempo_index[iter_pair]; // position of the pair in iMweight_M1
			
		  for(int iter_endpointTTE=0 ; iter_endpointTTE<iter_dTTE ; iter_endpointTTE++){                        
			if(iter_endpointTTE==(iter_dTTE-1) && survEndpoint[iter_d]){ // for the last endpoint (first test) add the new weights in case of survival endpoint (second test)
			  matWeight(iter_pair,iter_endpointTTE) = iWeight[iter_pair];
			}else{ // transfert the existing weights to the new matrix
			  matWeight(iter_pair,iter_endpointTTE) = matWeight_M1(iter_oldpair,iter_endpointTTE); // store iMweight_M1 in the iMweight restrected to the remaining pairs
			  // only if Wscheme is one in the column of the new endpoint and the line of the previous endpoint.            
			}
		  }
		} 
	  }

      // **** update all Scores
	  // Rcout << "update Scores" << endl;
	  if(keepScore){
		iNpairs = iScore.n_rows;
		iMat.resize(iNpairs,3);

		// add original index
		for(int iPair=0 ; iPair < iNpairs ; iPair ++){
		  iMat.row(iPair) = rowvec({(double)iter_strata,
				(double)index_strataC(iScore(iPair,0)),
				(double)index_strataT(iScore(iPair,1))});
		}
     
		// merge with current table and store
		if(iter_strata==0){
		  lsScore[iter_d] = arma::join_rows(iMat,iScore);
		}else{
		  lsScore[iter_d] = arma::join_cols(lsScore[iter_d], arma::join_rows(iMat,iScore));
	}
      }
    } // end endpoint

  }
  
  // ** proportion in favor of treatment 
  arma::mat delta_netBenefit(n_strata,D), delta_winRatio(n_strata,D); // matrix containing for each strata and each endpoint the statistic
  vector<double> Delta_netBenefit(D), Delta_winRatio(D); // vector containing for each endpoint the overall statistic
  
  calcStatistic(delta_netBenefit, delta_winRatio, Delta_netBenefit, Delta_winRatio,
                Mcount_favorable, Mcount_unfavorable, 
                D, n_strata, n_pairs);

  // ** export
  return(List::create(
		      Named("count_favorable")  = Mcount_favorable,
		      Named("count_unfavorable")  = Mcount_unfavorable,
		      Named("count_neutral")  = Mcount_neutral,           
		      Named("count_uninf")  = Mcount_uninf,           
		      Named("delta_netBenefit")  = delta_netBenefit,
		      Named("delta_winRatio")  = delta_winRatio,
		      Named("Delta_netBenefit")  = Delta_netBenefit,
		      Named("Delta_winRatio")  = Delta_winRatio,
		      Named("n_pairs")  = n_pairs,
		      Named("tableScore")  = lsScore
		      ));
  
}


