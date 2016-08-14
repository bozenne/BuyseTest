// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "FCT_calcOnePair.h"
#include "FCT_calcAllPairs.h"

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
//' @param type The type of each endpoint (1 binary, 2 continuous, 3 TTE). \emph{const IntegerVector&}. Must have length D.
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
//' @param PEP The type of method used to compare censored pairs (1 Peto, 2 Efron, 3 Peron).\emph{const int}.
//' 
//' @keywords function Cpp BuyseTest

//' @rdname BuyseTest_cpp
//' @export
// [[Rcpp::export]]
List BuyseTest_Gehan_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const IntegerVector& type, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
                         const int D, const bool returnIndex, const std::vector< arma::uvec >& strataT, const std::vector< arma::uvec >& strataC, const int n_strata, const int n_TTE){
  
  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  
  //// initialization ////
  structAllPairs resK; // store the result of the pairwize comparison of each endpoint and each strata in a SEXP list
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the number of favorable pairs by outcome for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the number of unfavorable pairs by outcome for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the number of neutral pairs by outcome for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the number of uninf pairs by outcome for each strata
  double n_pairs=0; // number of pairs sumed over the strats
  
  vector<int> index_neutralT(0) ; // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0) ; // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0) ; // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0) ; // index of the uninformative pairs of the control arm
  
  arma::mat TreatmentK; // Endpoint(s) for treated restrected to strata k
  arma::mat ControlK; // Endpoint(s) for controls restrescted to stata k
  arma::mat delta_TreatmentK ; // statut for TTE endpoints for treated restrected to strata k
  arma::mat delta_ControlK ; // statut for TTE endpoints for controls restrescted to stata k
  
  //// loop over strata ////
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
    // affect the endpoints corresponding to each strata (conversion from SEXP to uvec object is necessary for .rows)
    TreatmentK = Treatment.rows(strataT[iter_strata]); // select the rows in the treatment matrix corresponding to the indexes contained in strata
    ControlK = Control.rows(strataC[iter_strata]);  // select the rows in the control matrix corresponding to the indexes contained in strata
    int iter_dTTE=0; // number of time to event endpoints that have been used
    
    if(n_TTE>0){ // if there is TTE endpoints form the status matrix (conversion from SEXP to uvec object is necessary for .rows)
      delta_TreatmentK = delta_Treatment.rows(strataT[iter_strata]); // select the rows in the treatment status matrix corresponding to the indexes contained in strata
      delta_ControlK = delta_Control.rows(strataC[iter_strata]); // select the rows in the control status matrix corresponding to the indexes contained in strata     
    }
    
    //// first endpoint
    if(type[0]==1){ // binary endpoint
      resK = calcAllPairs_BinaryOutcome_cpp(TreatmentK.col(0),ControlK.col(0));
    }
    if(type[0]==2){ // continuous endpoint
      resK = calcAllPairs_ContinuousOutcome_cpp(TreatmentK.col(0),ControlK.col(0),threshold[0]);
    }
    if(type[0]==3){ // time to event endpoint
      resK = calcAllPairs_TTEOutcome_Gehan_cpp(TreatmentK.col(0),ControlK.col(0),threshold[0],
                                               delta_TreatmentK.col(0),delta_ControlK.col(0));
      iter_dTTE++; // increment the number of time to event endpoints that have been used
    }
    
    // store the number of pairs found in each catergory after conversion from SEXP to int. 
    Mcount_favorable(iter_strata,0) = resK.count_favorable;
    Mcount_unfavorable(iter_strata,0) = resK.count_unfavorable;
    Mcount_neutral(iter_strata,0) = resK.count_neutral;
    Mcount_uninf(iter_strata,0) = resK.count_uninf;
    // add to the total number of pairs the number of pairs founded for this endpoint
    n_pairs += Mcount_favorable(iter_strata,0)+Mcount_unfavorable(iter_strata,0)+Mcount_neutral(iter_strata,0)+Mcount_uninf(iter_strata,0); 
    
    //// following endpoints
    int iter_d = 0; // the index of the endpoints
    
    while(D>iter_d+1 && (Mcount_neutral(iter_strata,iter_d)>0 || Mcount_uninf(iter_strata,iter_d)>0)){ // loop over the following endpoints
      // while there are remaining endpoints and remaining neutral or uniformative pairs
      iter_d++; // increment the index of the endpoints
      
      if(type[iter_d]==1){ // binary endpoint
        resK = calcSubsetPairs_BinaryOutcome_cpp(TreatmentK.col(iter_d), ControlK.col(iter_d), 
                                                 resK.index_neutralT, resK.index_neutralC, resK.count_neutral,
                                                 resK.index_uninfT, resK.index_uninfC, resK.count_uninf);    
      }
      
      if(type[iter_d]==2){ // continuous endpoint
        resK = calcSubsetPairs_ContinuousOutcome_cpp(TreatmentK.col(iter_d), ControlK.col(iter_d), threshold[iter_d],
                                                     resK.index_neutralT, resK.index_neutralC, resK.count_neutral,
                                                     resK.index_uninfT, resK.index_uninfC, resK.count_uninf);   
      }
      if(type[iter_d]==3){ // time to event endpoint
        resK = calcSubsetPairs_TTEOutcome_Gehan_cpp(TreatmentK.col(iter_d), ControlK.col(iter_d), threshold[iter_d],
                                                    delta_TreatmentK.col(iter_dTTE), delta_ControlK.col(iter_dTTE),
                                                    resK.index_neutralT, resK.index_neutralC, resK.count_neutral,
                                                    resK.index_uninfT, resK.index_uninfC, resK.count_uninf);   
        iter_dTTE++; // increment the number of time to event endpoints that have been used
      }
      
      
      // store the number of pairs found in each catergory after conversion from SEXP to int. 
      Mcount_favorable(iter_strata,iter_d) = resK.count_favorable;
      Mcount_unfavorable(iter_strata,iter_d) = resK.count_unfavorable;
      Mcount_neutral(iter_strata,iter_d) = resK.count_neutral;
      Mcount_uninf(iter_strata,iter_d) = resK.count_uninf;    
    } // end endpoints
    
    if(returnIndex==true){ // store the neutral and the uninformative pairs for each strata
      index_neutralT.insert(index_neutralT.end(), resK.index_neutralT.begin(), resK.index_neutralT.end()); // insert the neutral pairs of the treatment arm after those already founded
      index_neutralC.insert(index_neutralC.end(), resK.index_neutralC.begin(), resK.index_neutralC.end()); // insert the neutral pairs of the control arm after those already founded
      index_uninfT.insert(index_uninfT.end(), resK.index_uninfT.begin(), resK.index_uninfT.end()); // insert the uninformative pairs of the treatment arm after those already founded
      index_uninfC.insert(index_uninfC.end(), resK.index_uninfC.begin(), resK.index_uninfC.end()); // insert the uninformative pairs of the control arm after those already founded   
    }
  }
  
  //// proportion in favor of treatment ////
  arma::mat delta(n_strata,D); // matrix containing for each strata and each endpoint the proportion in favor of treatment
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints
      delta(iter_strata,iter_d) = (Mcount_favorable(iter_strata,iter_d)-Mcount_unfavorable(iter_strata,iter_d))/(double)(n_pairs); // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the total number of pairs
    }  
  } 
  
  //// export ////
  if(returnIndex==true){
    return(List::create(
        Named("delta")  = delta,
        Named("count_favorable")  = Mcount_favorable,
        Named("count_unfavorable")  = Mcount_unfavorable,
        Named("count_neutral")  = Mcount_neutral,           
        Named("count_uninf")  = Mcount_uninf,           
        Named("index_neutralT")  = index_neutralT,
        Named("index_neutralC")  = index_neutralC,
        Named("index_uninfT")  = index_uninfT,
        Named("index_uninfC")  = index_uninfC,
        Named("n_pairs")  = n_pairs
    ));
  }else{
    return(List::create(
        Named("delta")  = delta
    ));
  }                   
  
}

//' @rdname BuyseTest_cpp
//' @export
// [[Rcpp::export]]
List BuyseTest_PetoEfronPeron_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const IntegerVector& type, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
                                  const int D, const bool returnIndex, const std::vector< arma::uvec >& strataT, const std::vector< arma::uvec >& strataC, const int n_strata, const int n_TTE, 
                                  const arma::mat& Wscheme, const IntegerVector index_survivalM1, const NumericVector threshold_TTEM1, 
                                  const std::vector< arma::mat >& list_survivalT, const std::vector< arma::mat >& list_survivalC, const int PEP){
  
  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.
  
  //// initialization ////  
  structAllPairs resK; // store the result of the pairwize comparison of each endpoint and each strata in a SEXP list
  arma::mat Mcount_favorable(n_strata,D,fill::zeros); // store the total weight of favorable pairs by outcome for each strata
  arma::mat Mcount_unfavorable(n_strata,D,fill::zeros); // store the total weight of unfavorable pairs by outcome for each strata
  arma::mat Mcount_neutral(n_strata,D,fill::zeros); // store the total weight of neutral pairs by outcome for each strata
  arma::mat Mcount_uninf(n_strata,D,fill::zeros); // store the total weight of uninf pairs by outcome for each strata
  double n_pairs=0; // number of pairs sumed over the strats
  
  int size_neutral; // number of neutral pairs (temporary)
  int size_uninf; // number of uninformative pairs (temporary)
  
  vector<int> index_neutralT(0) ; // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0) ; // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0) ; // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0) ; // index of the uninformative pairs of the control arm
  vector<int> tempo_index;  // temporary index vector used to convert the index SEXP into vector<int>
  
  arma::uvec index_strataT; // position of a patients belonging to a given strata in the treatment arm
  arma::uvec index_strataC; // position of a patients belonging to a given strata in the control arm
  arma::mat TreatmentK; // Endpoint(s) for treated restrected to strata k
  arma::mat ControlK; // Endpoint(s) for controls restrescted to stata k
  arma::mat delta_TreatmentK ; // statut for TTE endpoints for treated restrected to strata k
  arma::mat delta_ControlK ; // statut for TTE endpoints for controls restrescted to stata k
  
  // loop over strata
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){
    
    index_strataT = strataT[iter_strata];
    index_strataC = strataC[iter_strata];
    
    // affect the endpoints corresponding to each strata (conversion from SEXP to uvec object is necessary for .rows)
    TreatmentK = Treatment.rows(index_strataT); // select the rows in the treatment matrix corresponding to the indexes contained in strata
    ControlK = Control.rows(index_strataC); // select the rows in the control matrix corresponding to the indexes contained in strata
    int iter_dTTE=0; // number of time to event endpoints that have been used
    
    delta_TreatmentK = delta_Treatment.rows(index_strataT); // select the rows in the treatment status matrix corresponding to the indexes contained in strata
    delta_ControlK = delta_Control.rows(index_strataC); // select the rows in the control status matrix corresponding to the indexes contained in strata          
    
    //// first endpoint
    if(type[0]==1){ // binary endpoint
      resK = calcAllPairs_BinaryOutcome_cpp(TreatmentK.col(0),ControlK.col(0));      
    }
    if(type[0]==2){ // continuous endpoint
      resK = calcAllPairs_ContinuousOutcome_cpp(TreatmentK.col(0),ControlK.col(0),threshold[0]);
    }
    if(type[0]==3){ // time to event endpoint   
      resK = calcAllPairs_TTEOutcome_PetoEfronPeron_cpp(TreatmentK.col(0),ControlK.col(0),threshold[0],
                                                        delta_TreatmentK.col(0),delta_ControlK.col(0),
                                                        list_survivalT[0].rows(index_strataT),
                                                        list_survivalC[0].rows(index_strataC), 
                                                        PEP); 
      iter_dTTE++; // increment the number of time to event endpoints that have been used   
    }
    
    // store the total weight corresponding to each catergory after conversion from SEXP to double. 
    Mcount_favorable(iter_strata,0) = resK.count_favorable;
    Mcount_unfavorable(iter_strata,0) = resK.count_unfavorable;
    Mcount_neutral(iter_strata,0) = resK.count_neutral;
    Mcount_uninf(iter_strata,0) = resK.count_uninf;
    // add to the total number of pairs the number of pairs founded for this endpoint
    n_pairs += Mcount_favorable(iter_strata,0)+Mcount_unfavorable(iter_strata,0)+Mcount_neutral(iter_strata,0)+Mcount_uninf(iter_strata,0);
    
    //// update Wpairs 
    size_neutral = resK.index_neutralT.size(); // update the number of neutral pairs
    size_uninf = resK.index_uninfT.size(); // update the number of uninformative pairs
    arma::mat Wpairs(size_neutral+size_uninf,1,fill::ones); // temporary matrix containing the weigth of each remaining pair for each outcome
    arma::vec w(size_neutral+size_uninf); // temporary vector containing the weight of each remaining pair to be used for the next outcome
    w.fill(1);
    if(type[0]==3 && D>1){ // update the weights for the uninformative pairs in Wpairs and w
      for(int iter_uninf=0 ; iter_uninf<size_uninf ; iter_uninf++){ // neutral pairs have a weight of 1 by construction
        Wpairs(size_neutral+iter_uninf,0) = resK.w[iter_uninf];
        if(Wscheme(0,0)==1){w(size_neutral+iter_uninf) = resK.w[iter_uninf];}
      }
    }        
    
    //// following endpoints
    int iter_d = 0; // the index of the endpoints
    arma::mat Wpairs_sauve; // the previous update of Wpairs
    
    while(D>iter_d+1 && (Mcount_neutral(iter_strata,iter_d)>0 || Mcount_uninf(iter_strata,iter_d)>0)){ // loop over the following endpoints
      
      // while there are remaining endpoints and remaining neutral or uniformative pairs
      iter_d++; // increment the index of the endpoints
      Wpairs_sauve=Wpairs; // save the current Wpairs
      
      if(type[iter_d]==1){ // binary endpoint        
        resK = calcSubsetPairs_WeightedBinaryOutcome_cpp(TreatmentK.col(iter_d),ControlK.col(iter_d),
                                                         resK.index_neutralT,resK.index_neutralC, size_neutral,
                                                         resK.index_uninfT,resK.index_uninfC, size_uninf,
                                                         w);    
      }
      if(type[iter_d]==2){ // continuous endpoint
        resK = calcSubsetPairs_WeightedContinuousOutcome_cpp(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                                                             resK.index_neutralT,resK.index_neutralC, size_neutral,
                                                             resK.index_uninfT,resK.index_uninfC, size_uninf,
                                                             w);   
      }
      if(type[iter_d]==3){ // time to event endpoint 
        
        if(threshold_TTEM1[iter_dTTE]<0){ // first time the endpoint is used
          resK = calcSubsetPairs_TTEOutcome_PetoEfronPeron_cpp(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                                                               delta_TreatmentK.col(iter_dTTE),delta_ControlK.col(iter_dTTE),
                                                               list_survivalT[iter_dTTE].rows(index_strataT),
                                                               list_survivalC[iter_dTTE].rows(index_strataC),
                                                               resK.index_neutralT,resK.index_neutralC, size_neutral,
                                                               resK.index_uninfT,resK.index_uninfC, size_uninf,
                                                               w, -1, arma::mat(1,1), arma::mat(1,1), PEP); 
        }else{ // following times    
          
          resK = calcSubsetPairs_TTEOutcome_PetoEfronPeron_cpp(TreatmentK.col(iter_d),ControlK.col(iter_d),threshold[iter_d],
                                                               delta_TreatmentK.col(iter_dTTE),delta_ControlK.col(iter_dTTE),
                                                               list_survivalT[iter_dTTE].rows(index_strataT),
                                                               list_survivalC[iter_dTTE].rows(index_strataC), 
                                                               resK.index_neutralT,resK.index_neutralC, size_neutral,
                                                               resK.index_uninfT,resK.index_uninfC, size_uninf,
                                                               w, threshold_TTEM1[iter_dTTE], 
                                                                                 list_survivalT[index_survivalM1[iter_dTTE]].rows(index_strataT),
                                                                                 list_survivalC[index_survivalM1[iter_dTTE]].rows(index_strataC), 
                                                                                 PEP); 
        }
        iter_dTTE++; // increment the number of time to event endpoints that have been used
      }
      
      // store the number of pairs found in each catergory after conversion from SEXP to double. 
      Mcount_favorable(iter_strata,iter_d) = resK.count_favorable;
      Mcount_unfavorable(iter_strata,iter_d) = resK.count_unfavorable;
      Mcount_neutral(iter_strata,iter_d) = resK.count_neutral;
      Mcount_uninf(iter_strata,iter_d) = resK.count_uninf;    
      
      // update Wpairs
      size_neutral = resK.index_neutralT.size(); // update the number of neutral pairs
      size_uninf = resK.index_uninfT.size(); // update the number of uninformative pairs
      Wpairs.resize(size_neutral+size_uninf,max(1,iter_dTTE)); // update the size of Wpairs
      w.resize(size_neutral+size_uninf); // update the size of w
      w.fill(1);
      
      tempo_index=resK.index_w; // store the position of the remaining pairs in the previous Wpairs (i.e. Wpairs_sauve)
      int iter_oldpair; // index of the remaining pair
      
      if(iter_dTTE>0 && D>iter_d+1){
        
        for(size_t iter_pair=0; iter_pair<tempo_index.size(); iter_pair++){
          
          iter_oldpair = tempo_index[iter_pair]; // position of the pair in Wpairs_sauve
          
          for(int iter_endpointTTE=0 ; iter_endpointTTE<iter_dTTE ; iter_endpointTTE++){          
            
            if(iter_endpointTTE==(iter_dTTE-1) && type[iter_d]==3){ // for the last endpoint (first test) add the new weights in case of survival endpoint (second test)
              Wpairs(iter_pair,iter_endpointTTE) = resK.w[iter_pair];
              if(Wscheme(iter_endpointTTE,iter_d)==1){w(iter_pair) *= resK.w[iter_pair];} // iter_d - 1 + 1 because the first column is missing but we are interested in the next endpoint
            }else{ // transfert the existing weights to the new matrix
              Wpairs(iter_pair,iter_endpointTTE) = Wpairs_sauve(iter_oldpair,iter_endpointTTE); // store Wpairs_sauve in the Wpairs restrected to the remaining pairs
              if(Wscheme(iter_endpointTTE,iter_d)==1){w(iter_pair) *= Wpairs_sauve(iter_oldpair,iter_endpointTTE);} // make the product over the endpoints of the weights associated to each remaining pair 
              // only if Wscheme is one in the column of the new endpoint and the line of the previous endpoint.            
            }
          }
          
          
        }
      }
      
    } // end endpoint
    
    
    if(returnIndex==true){ // store the neutral and the uninformative pairs 
      index_neutralT.insert(index_neutralT.end(), resK.index_neutralT.begin(), resK.index_neutralT.end()); // insert the neutral pairs of the treatment arm after those already founded
      index_neutralC.insert(index_neutralC.end(), resK.index_neutralC.begin(), resK.index_neutralC.end()); // insert the neutral pairs of the control arm after those already founded
      index_uninfT.insert(index_uninfT.end(), resK.index_uninfT.begin(), resK.index_uninfT.end());
      index_uninfC.insert(index_uninfC.end(), resK.index_uninfC.begin(), resK.index_uninfC.end());    
    }
    
  }
  
  //// proportion in favor of treatment ////
  arma::mat delta(n_strata,D); // matrix containing for each strata and each endpoint the proportion in favor of treatment
  for(int iter_strata=0 ; iter_strata < n_strata ; iter_strata ++){ // loop over strata
    for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints
      delta(iter_strata,iter_d) = (Mcount_favorable(iter_strata,iter_d)-Mcount_unfavorable(iter_strata,iter_d))/(double)(n_pairs); // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the total number of pairs
    }  
  } 
  
  
  //// export ////
  if(returnIndex==true){
    
    return(List::create(
        Named("delta")  = delta,
        Named("count_favorable")  = Mcount_favorable,
        Named("count_unfavorable")  = Mcount_unfavorable,
        Named("count_neutral")  = Mcount_neutral,           
        Named("count_uninf")  = Mcount_uninf,           
        Named("index_neutralT")  = index_neutralT,
        Named("index_neutralC")  = index_neutralC,
        Named("index_uninfT")  = index_uninfT,
        Named("index_uninfC")  = index_uninfC,
        Named("n_pairs")  = n_pairs   
    ));
  }else{
    return(List::create(
        Named("delta")  = delta
    ));
  }                   
  
}

//////////////////////////////////// PROGRAMME 2 : BuysePower //////////////////////////////////////////////////////////////


//////  fct3 : perform pairwize comparison with imputation for censored data ///////////////////////////////////////////////// 
//// [[Rcpp::export]]
////List BuyseTest_Gehan_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const IntegerVector& type, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
////const int D, const bool returnIndex, List& strataT, List& strataC, const int n_strata, const int n_TTE);
//
//List BuysePower_Gehan_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const IntegerVector& type, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
//const int D, const int n_TTE, 
//const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
//  // numT : index corresponding to each patient in the treatment arm
//  // numC : index corresponding to each patient in the control arm
//  // sizeT : sample sizes of the treatment arm
//  // sizeC : sample sizes of the control arm
//  
//  //// initialization ////
//  List resK; // store the result of the pairwize comparison of each endpoint and each sample size in a SEXP list
//  arma::mat Mcount_favorable(D, n_size,fill::zeros);
//  arma::mat Mcount_unfavorable(D, n_size,fill::zeros);
//  
//  vector<int> n_pairs(n_size,0); // number of pairs sumed over the strats
//
//  vector<int> index_neutralT(0) ; // index of the neutral pairs of the treatment arm
//  vector<int> index_neutralC(0) ; // index of the neutral pairs of the control arm
//  vector<int> index_uninfT(0) ; // index of the uninformative pairs of the treatment arm
//  vector<int> index_uninfC(0) ; // index of the uninformative pairs of the control arm
//  vector<int> tempo_index0; // temporary index vector used to convert the index from SEXP into vector<int>
//  vector<int> tempo_index1; // temporary index vector used to convert the index from SEXP into vector<int>
//  vector<int> tempo_index2; // temporary index vector used to convert the index from SEXP into vector<int>
//  vector<int> tempo_index3; // temporary index vector used to convert the index from SEXP into vector<int>
//  
//  int iter_dTTE=0; // number of time to event endpoints that have been used
//    
//  //// first endpoint
//    if(type[0]==1){ // binary endpoint
//    resK = calcAllPairsPower_BinaryEndpoint_cpp(Treatment.col(0),Control.col(0),
//    numT, numC, sizeT, sizeC, n_size);
//    }
//    
//    if(type[0]==2){ // continuous endpoint
//    resK = calcAllPairsPower_ContinuousEndpoint_cpp(Treatment.col(0),Control.col(0),threshold[0],
//    numT, numC, sizeT, sizeC, n_size);
//    }
//    
//    if(type[0]==3){ // time to event endpoint
//    resK = calcAllPairsPower_TTEEndpoint_Gehan_cpp(Treatment.col(0),Control.col(0),threshold[0],
//    delta_Treatment.col(0),delta_Control.col(0),
//    numT, numC, sizeT, sizeC, n_size);
//    
//    iter_dTTE++; // increment the number of time to event endpoints that have been used
//    }
//    
//    // store the number of pairs found in each catergory after conversion from SEXP to int. 
//    tempo_index0 = resK[0];
//    tempo_index1 = resK[1];
//    tempo_index2 = resK[2];
//    tempo_index3 = resK[3];
//    
//    for(int iter_size=0; iter_size<n_size ; iter_size++){
//      Mcount_favorable(0,iter_size) = tempo_index0[iter_size];
//      Mcount_unfavorable(0,iter_size) = tempo_index1[iter_size];
//      
//      // add to the total number of pairs the number of pairs founded for this endpoint
//      n_pairs[iter_size]=tempo_index0[iter_size]+tempo_index1[iter_size]+tempo_index2[iter_size]+tempo_index3[iter_size]; 
//    }
//   
//    //// following endpoints
//    int iter_d = 0; // the index of the endpoints
//    
//    while(D>iter_d+1 && (tempo_index2[n_size-1]>0 || tempo_index3[n_size-1]>0) ){ // loop over the following endpoints
//    // while there are remaining endpoints and remaining neutral or uniformative pairs
//    iter_d++; // increment the index of the endpoints
//    
//    if(type[iter_d]==1){ // binary endpoint
//    resK = calcSubsetPairsPower_BinaryEndpoint_cpp(Treatment.col(iter_d), Control.col(iter_d),
//    resK[4], resK[5], tempo_index2[n_size-1],
//    resK[6], resK[7], tempo_index3[n_size-1],
//    numT, numC, sizeT, sizeC, n_size);    
//    }
//    if(type[iter_d]==2){ // continuous endpoint
//    resK = calcSubsetPairsPower_ContinuousEndpoint_cpp(Treatment.col(iter_d), Control.col(iter_d), threshold[iter_d],
//    resK[4], resK[5], tempo_index2[n_size-1],
//    resK[6], resK[7], tempo_index3[n_size-1], 
//    numT, numC, sizeT, sizeC, n_size);   
//    }
//    if(type[iter_d]==3){ // time to event endpoint
//    resK = calcSubsetPairsPower_TTEEndpoint_Gehan_cpp(Treatment.col(iter_d), Control.col(iter_d), threshold[iter_d],
//    delta_Treatment.col(iter_dTTE), delta_Control.col(iter_dTTE),
//    resK[4], resK[5], tempo_index2[n_size-1],
//    resK[6], resK[7], tempo_index3[n_size-1], 
//    numT, numC, sizeT, sizeC, n_size);   
//    
//    iter_dTTE++; // increment the number of time to event endpoints that have been used
//    }
//    
//    // store the number of pairs found in each catergory after conversion from SEXP to int. 
//    tempo_index0 = resK[0];
//    tempo_index1 = resK[1];
//    
//    for(int iter_size=0; iter_size<n_size ; iter_size++){
//      Mcount_favorable(iter_d,iter_size) = tempo_index0[iter_size];
//      Mcount_unfavorable(iter_d,iter_size) = tempo_index1[iter_size];
//    }
//    }
//  
//  
//  //// proportion in favor of treatment ////
//  arma::mat delta(D, n_size,fill::zeros); // array containing for each sample size and each endpoint the proportion in favor of treatment
//  for(int iter_size=0 ; iter_size < n_size ; iter_size ++){ // loop over sample size
//         for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints
//           delta(iter_d,iter_size) = (Mcount_favorable(iter_d,iter_size)-Mcount_unfavorable(iter_d,iter_size))/(double)(n_pairs[iter_size]); // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the total number of pairs
//    } 
//  } 
//  
//  //// export ////
//  return(List::create(Named("delta")  = delta));
//    
//}


