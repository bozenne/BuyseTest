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

//' @name BuyseTest_cpp
//' @aliases BuyseTest_Gehan_cpp BuyseTest_PetoEfronPeron_cpp
//' @title C++ function performing the pairwise comparison over several endpoints. 
//' 
//' @description \code{BuyseTest_Gehan_cpp} and \code{BuyseTest_PetoEfron_cpp} functions calls for each endpoint and each strata the pairwise comparison function suited to the type of endpoint and store the resuts. 
//' 
//' @param Treatment A matrix containing the values of each endpoint (in columns) for the treatment observations (in rows). \emph{const arma::mat&}. Must have D columns.
//' @param Control A matrix containing the values of each endpoint (in columns) for the control observations (in rows). \emph{const arma::mat&}.
//' @param threshold Store the thresholds associated to each endpoint. \emph{const NumericVector&}. Must have length D. The threshold is ignored for binary endpoints. Must have D columns.
//' @param survEndpoint Does each endpoint is a time to event. \emph{const LogicalVector&}. Must have length D.
//' @param delta_Treatment A matrix containing in the type of event (0 censoring, 1 event) for each TTE endpoint (in columns) and treatment observations (in rows). \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.
//' @param delta_Control A matrix containing the nature of observations in the control group (in rows) (0 censoring, 1 event) for each TTE endpoint (in columns) . \emph{const arma::mat&} containing binary integers. Must have n_TTE columns. Ignored if n_TTE equals 0.
//' @param D The number of endpoints. Strictly positive \emph{const int}.
//' @param returnIndex Should the indexes of the neutral or uninformative pairs be returned. \emph{const bool}.
//' @param strataT A list containing the indexes of treatment observations belonging for each strata. \emph{List&} of vector containing positive integers. 
//' @param strataC A list containing the indexes of control observations belonging for each strata. \emph{List&}  of vector containing positive integers. 
//' @param n_strata The number of strata . Strictly positive \emph{const int}.
//' @param n_TTE The number of time-to-event endpoints. Positive \emph{const int}.
//' @param Wscheme The matrix describing the weighting strategy. For each endpoint (except the first) in column, weights of each pair are initialized at 1 and multiplyied by the weight of the endpoints in raws where there is a 1. \emph{const arma::mat&}. Must have n_TTE lines and D-1 columns.
//' @param index_survivalM1 The position, among all the survival endpoints, of the last same endpoint (computed with a different threshold). If it is the first time that the TTE endpoint is used it is set to -1. \emph{const IntegerVector}. Must have length n_TTE.
//' @param threshold_TTEM1 The previous latest threshold of each TTE endpoint. When it is the first time that the TTE endpoint is used it is set to -1. \emph{const NumericVector}. Must have length n_TTE.
//' @param list_survivalT A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the treatment group (in rows). \emph{List&}. Must have length n_TTE. Each matrix must have 3 (if method is Peto, only one survival function is computed) or 11 (if method is Efron or Peron, two survival functions are computed) columns. Ignored if method is Gehan.
//' @param list_survivalC A list of matrix containing the survival estimates (-threshold, 0, +threshold ...) for each event of the control group (in rows). \emph{List&}. Must have length n_TTE. Each matrix must have 3 (if method is Peto) or 11 (if method is Efron or Peron) columns. Ignored if method is Gehan.
//' @param methodTTE The type of method used to compare censored pairs (1 Peto, 2 Efron, 3 Peron).\emph{const int}.
//' @param neutralAsUninf Should paired classified as neutral be re-analysed using endpoints of lower priority. \emph{logical}.
//' 
//' @keywords function Cpp BuyseTest
//' @export
// [[Rcpp::export]]
List GPC_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const LogicalVector& survEndpoint, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
             const int D, const bool returnIndex, const std::vector< arma::uvec >& strataT, const std::vector< arma::uvec >& strataC, const int n_strata, const int n_TTE, 
             const arma::mat& Wscheme, const IntegerVector index_survivalM1, const NumericVector threshold_TTEM1, 
             const std::vector< arma::mat >& list_survivalT, const std::vector< arma::mat >& list_survivalC, const int methodTTE, const double neutralAsUninf){
  
  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.
  
  //// initialization ////  
  // final results
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs by outcome for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs by outcome for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs by outcome for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs by outcome for each strata
  double n_pairs=0; // number of pairs sumed over the strats
  
  vector<int> index_neutralT(0) ; // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0) ; // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0) ; // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0) ; // index of the uninformative pairs of the control arm
  
  // for a given stata [input]
  arma::uvec index_strataT; // position of a patients belonging to a given strata in the treatment arm
  arma::uvec index_strataC; // position of a patients belonging to a given strata in the control arm
  arma::mat TreatmentK; // Endpoint(s) for treated restricted to strata k
  arma::mat ControlK; // Endpoint(s) for controls restricted to stata k
  arma::mat delta_TreatmentK ; // statut for TTE endpoints for treated restricted to strata k
  arma::mat delta_ControlK ; // statut for TTE endpoints for controls restricted to stata k
  arma::mat matKMT_K ; // KM estimates for treated restricted to stata k
  arma::mat matKMC_K ; // KM estimates for controls restricted to stata k
  
  int iter_d; // the index of the endpoints
  int iter_dTTE; // number of time to event endpoints that have been used
  int iter_oldpair; // index of the remaining pair [local use when updating Wpairs]
  
  // for a given stata [output]
  vector<int> iIndex_neutralT; // index of the neutral pairs of the treatment arm
  vector<int> iIndex_neutralC; // index of the neutral pairs of the control arm
  vector<int> iIndex_uninfT; // index of the uninformative pairs of the treatment arm
  vector<int> iIndex_uninfC; // index of the uninformative pairs of the control arm
  vector<double> iw; // weights of the uninformative pairs
  vector<int> iIndex_w; // index of the pair number of the uninformative pairs
  
  int size_neutral; // number of neutral pairs (temporary)
  int size_uninf; // number of uninformative pairs (temporary)
  
  arma::mat Wpairs; // history of the weights
  arma::mat Wpairs_sauve; // the previous update of Wpairs
  arma::vec w; // current weigths
  
  vector<int> tempo_index;  // temporary index vector used to convert the index SEXP into vector<int>
  
  // loop over strata
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
    
    iter_d = 0;
    iter_dTTE = 0;
    
    index_strataT = strataT[iter_strata];
    index_strataC = strataC[iter_strata];
    
    // affect the endpoints corresponding to each strata (conversion from SEXP to uvec object is necessary for .rows)
    TreatmentK = Treatment.rows(index_strataT); // select the rows in the treatment matrix corresponding to the indexes contained in strata
    ControlK = Control.rows(index_strataC); // select the rows in the control matrix corresponding to the indexes contained in strata
    
    if(n_TTE>0){
      delta_TreatmentK = delta_Treatment.rows(index_strataT); // select the rows in the treatment status matrix corresponding to the indexes contained in strata
      delta_ControlK = delta_Control.rows(index_strataC); // select the rows in the control status matrix corresponding to the indexes contained in strata          
    }
    
    //// first endpoint
    iIndex_neutralT.resize(0); iIndex_neutralC.resize(0); iIndex_uninfT.resize(0); iIndex_uninfC.resize(0); iw.resize(0);
    
    if(survEndpoint[0]){ // time to event endpoint  
      if(methodTTE>0){
        matKMT_K = list_survivalT[0].rows(index_strataT);
        matKMC_K = list_survivalC[0].rows(index_strataC);
      }
      
      calcAllPairs_TTE(TreatmentK.col(0),ControlK.col(0),threshold[0],
                       delta_TreatmentK.col(0),delta_ControlK.col(0), matKMT_K, matKMC_K, methodTTE,
                       Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
                       iIndex_neutralT, iIndex_neutralC, iIndex_uninfT, iIndex_uninfC, 
                       iw); 
      
      iter_dTTE++; // increment the number of time to event endpoints that have been used   
    }else { // binary or continuous endpoint
      calcAllPairs_Continuous(TreatmentK.col(0),ControlK.col(0),threshold[0],
                              Mcount_favorable(iter_strata,0), Mcount_unfavorable(iter_strata,0), Mcount_neutral(iter_strata,0), Mcount_uninf(iter_strata,0), 
                              iIndex_neutralT, iIndex_neutralC, iIndex_uninfT, iIndex_uninfC);
    }
    
    // add to the total number of pairs the number of pairs founded for this endpoint
    n_pairs += Mcount_favorable(iter_strata,0)+Mcount_unfavorable(iter_strata,0)+Mcount_neutral(iter_strata,0)+Mcount_uninf(iter_strata,0);
    size_neutral = iIndex_neutralT.size(); // update the number of neutral pairs
    size_uninf = iIndex_uninfT.size(); // update the number of uninformative pairs
    
    //// update Wpairs 
    if(D>1){ // if there is more than one endpoint
      Wpairs.resize(size_neutral+size_uninf,1); // temporary matrix containing the weigth of each remaining pair for each outcome
      Wpairs.fill(1.0);
      w.resize(size_neutral+size_uninf); // temporary vector containing the weight of each remaining pair to be used for the next outcome
      w.fill(1);
      
      if(methodTTE>0 && iter_dTTE>0){ // update the weights for the uninformative pairs in Wpairs and w
        for(int iter_uninf=0 ; iter_uninf<size_uninf ; iter_uninf++){ // neutral pairs have a weight of 1 by construction
          Wpairs(size_neutral+iter_uninf,0) = iw[iter_uninf];
          if(Wscheme(0,0)==1){w(size_neutral+iter_uninf) = iw[iter_uninf];}
        }
      }        
    }
    
    
    //// following endpoints
    while(D>iter_d+1 && (Mcount_neutral(iter_strata,iter_d)>0 || Mcount_uninf(iter_strata,iter_d)>0)){ // loop over the following endpoints
      
      // while there are remaining endpoints and remaining neutral or uniformative pairs
      iter_d++; // increment the index of the endpoints
      Wpairs_sauve=Wpairs; // save the current Wpairs
      if(neutralAsUninf==false){size_neutral = 0;}
      iw.resize(0); iIndex_w.resize(0);
      
      if(survEndpoint[iter_d]){ // time to event endpoint 
        if(methodTTE>0){
          matKMT_K = list_survivalT[iter_dTTE].rows(index_strataT);
          matKMC_K = list_survivalC[iter_dTTE].rows(index_strataC);
        }
        
        if(methodTTE==0 || threshold_TTEM1[iter_dTTE]<0){ // first time the endpoint is used [no threshold-1]
          calcSubsetPairs_TTE(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                                     delta_TreatmentK.col(iter_dTTE),delta_ControlK.col(iter_dTTE), matKMT_K, matKMC_K, methodTTE,
                                     Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
                                     iIndex_neutralT, iIndex_neutralC, size_neutral,
                                     iIndex_uninfT, iIndex_uninfC, size_uninf,
                                     w, -1, arma::mat(1,1), arma::mat(1,1),
                                     iw, iIndex_w); 
        }else{ // following times    
          calcSubsetPairs_TTE(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                              delta_TreatmentK.col(iter_dTTE),delta_ControlK.col(iter_dTTE), matKMT_K, matKMC_K, methodTTE,
                              Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
                              iIndex_neutralT, iIndex_neutralC, size_neutral,
                              iIndex_uninfT, iIndex_uninfC, size_uninf,
                              w, threshold_TTEM1[iter_dTTE], list_survivalT[index_survivalM1[iter_dTTE]].rows(index_strataT), list_survivalC[index_survivalM1[iter_dTTE]].rows(index_strataC), 
                              iw, iIndex_w); 
        }
        iter_dTTE++; // increment the number of time to event endpoints that have been used
      }else{ // binary or continuous endpoint
        calcSubsetPairs_Continuous(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                                   Mcount_favorable(iter_strata,iter_d), Mcount_unfavorable(iter_strata,iter_d), Mcount_neutral(iter_strata,iter_d), Mcount_uninf(iter_strata,iter_d), 
                                   iIndex_neutralT, iIndex_neutralC, size_neutral,
                                   iIndex_uninfT, iIndex_uninfC, size_uninf,
                                   w, iIndex_w);   
      }
      
      size_neutral = iIndex_neutralT.size(); // update the number of neutral pairs
      size_uninf = iIndex_uninfT.size(); // update the number of uninformative pairs
      
      // update Wpairs
      if(D>iter_d+1){
        Wpairs.resize(size_neutral+size_uninf,iter_dTTE); // update the size of Wpairs
        w.resize(size_neutral+size_uninf); // update the size of w
        w.fill(1);
        
        if(methodTTE>0 && iter_dTTE>0){
          tempo_index=iIndex_w; // store the position of the remaining pairs in the previous Wpairs (i.e. Wpairs_sauve)
          for(size_t iter_pair=0; iter_pair<tempo_index.size(); iter_pair++){
            
            iter_oldpair = tempo_index[iter_pair]; // position of the pair in Wpairs_sauve
            
            for(int iter_endpointTTE=0 ; iter_endpointTTE<iter_dTTE ; iter_endpointTTE++){          
              
              if(iter_endpointTTE==(iter_dTTE-1) && survEndpoint[iter_d]){ // for the last endpoint (first test) add the new weights in case of survival endpoint (second test)
                Wpairs(iter_pair,iter_endpointTTE) = iw[iter_pair];
                if(Wscheme(iter_endpointTTE,iter_d)==1){w(iter_pair) *= iw[iter_pair];} // iter_d - 1 + 1 because the first column is missing but we are interested in the next endpoint
              }else{ // transfert the existing weights to the new matrix
                Wpairs(iter_pair,iter_endpointTTE) = Wpairs_sauve(iter_oldpair,iter_endpointTTE); // store Wpairs_sauve in the Wpairs restrected to the remaining pairs
                if(Wscheme(iter_endpointTTE,iter_d)==1){w(iter_pair) *= Wpairs_sauve(iter_oldpair,iter_endpointTTE);} // make the product over the endpoints of the weights associated to each remaining pair 
                // only if Wscheme is one in the column of the new endpoint and the line of the previous endpoint.            
              }
            }
          } 
        }
      }
      
    } // end endpoint
    
    // 
    if(returnIndex==true){ // store the neutral and the uninformative pairs 
      index_neutralT.insert(index_neutralT.end(), iIndex_neutralT.begin(), iIndex_neutralT.end()); // insert the neutral pairs of the treatment arm after those already founded
      index_neutralC.insert(index_neutralC.end(), iIndex_neutralC.begin(), iIndex_neutralC.end()); // insert the neutral pairs of the control arm after those already founded
      index_uninfT.insert(index_uninfT.end(), iIndex_uninfT.begin(), iIndex_uninfT.end());
      index_uninfC.insert(index_uninfC.end(), iIndex_uninfC.begin(), iIndex_uninfC.end());    
    }
    
  }
  
  //// proportion in favor of treatment ////
  arma::mat delta_netChance(n_strata,D), delta_winRatio(n_strata,D); // matrix containing for each strata and each endpoint the statistic
  vector<double> Delta_netChance(D), Delta_winRatio(D); // vector containing for each endpoint the overall statistic
  
  calcStatistic(delta_netChance, delta_winRatio, Delta_netChance, Delta_winRatio,
                Mcount_favorable, Mcount_unfavorable, 
                D, n_strata, n_pairs);
  
  //// export ////
  if(returnIndex==true){
    
    return(List::create(
        Named("count_favorable")  = Mcount_favorable,
        Named("count_unfavorable")  = Mcount_unfavorable,
        Named("count_neutral")  = Mcount_neutral,           
        Named("count_uninf")  = Mcount_uninf,           
        Named("delta_netChance")  = delta_netChance,
        Named("delta_winRatio")  = delta_winRatio,
        Named("Delta_netChance")  = Delta_netChance,
        Named("Delta_winRatio")  = Delta_winRatio,
        Named("index_neutralT")  = index_neutralT,
        Named("index_neutralC")  = index_neutralC,
        Named("index_uninfT")  = index_uninfT,
        Named("index_uninfC")  = index_uninfC,
        Named("n_pairs")  = n_pairs   
    ));
  }else{
    return(List::create(
        Named("delta_netChance")  = delta_netChance,
        Named("delta_winRatio")  = delta_winRatio,
        Named("Delta_netChance")  = Delta_netChance,
        Named("Delta_winRatio")  = Delta_winRatio
    ));
  }                   
  
}


