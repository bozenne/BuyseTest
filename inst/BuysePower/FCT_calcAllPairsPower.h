#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

//  fct2 : perform pairwise comparisons over all possible pairs for a binary endpoint //////////////////////////////////////////
// [[Rcpp::export]]
List calcAllPairsPower_BinaryEndpoint_cpp(const arma::colvec& Treatment, const arma::colvec& Control, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> count_neutral(n_size,0); // number of neutral pairs for each sample size
  vector<int> count_uninf(n_size,0); // number of uninformative pairs for each sample size
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
    index_uninfT.push_back(iter_T);
    index_uninfC.push_back(iter_C);    
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_uninf[iter_size]++;
       }
      }
  }else if(Treatment[iter_T]==1){
    
    if(Control[iter_C]==1){ 
      index_neutralT.push_back(iter_T); 
      index_neutralC.push_back(iter_C);
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_neutral[iter_size]++;
       }
      }
    }else{ // Control[iter_C]==0
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_favorable[iter_size]++;
      }
    }
    }          
    
  }else{ // Treatment[iter_T]==0
  
  if(Control[iter_C]==1){
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_unfavorable[iter_size]++;
      }
    }
  }else{ // Control[iter_C]==0
  index_neutralT.push_back(iter_T);
  index_neutralC.push_back(iter_C);
  for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_neutral[iter_size]++;
       }
  }
  }
  
  }
  
  }}
  
  //// export ////
  return(List::create(
    Named("count_favorable")  = count_favorable, // 0
    Named("count_unfavorable")  = count_unfavorable, // 1
    Named("count_neutral")  = count_neutral, // 2  
    Named("count_uninf")  = count_uninf, // 3      
    Named("index_neutralT")  = index_neutralT, // 4
    Named("index_neutralC")  = index_neutralC, // 5
    Named("index_uninfT")  = index_uninfT, // 6 
    Named("index_uninfC")  = index_uninfC // 7
    ));
    
}

//  fct2bis : perform pairwise comparisons over a prespecified subset of pairs for a binary endpoint ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List pairSampleBin2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);      
      }else if(Treatment[iter_T]==1){
        
        if(Control[iter_C]==1){
          indexNew_neutralT.push_back(iter_T);
          indexNew_neutralC.push_back(iter_C);
        }else{ // Control[iter_C]==0
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_favorable[iter_size] = count_favorable[iter_size] + 1;
          }
        }
        }        
        
      }else{ // Treatment[iter_T]==0
      
      if(Control[iter_C]==1){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
          }
        }
      }else{ // Control[iter_C]==0
      indexNew_neutralT.push_back(iter_T);
      indexNew_neutralC.push_back(iter_C);
      }
      
      }}}
      
      //// loop over the uninformative pairs ////
      if(nUninf_pairs>0){
        
        for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
          
          iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
          iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
          
          if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
            indexNew_uninfT.push_back(iter_T);
            indexNew_uninfC.push_back(iter_C);      
          }else if(Treatment[iter_T]==1){
            
            if(Control[iter_C]==1){
              indexNew_neutralT.push_back(iter_T);
              indexNew_neutralC.push_back(iter_C);
            }else{ // Control[iter_C]==0
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + 1;
              }
            }
            }        
            
          }else{ // Treatment[iter_T]==0
          
          if(Control[iter_C]==1){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
              }
            }
          }else{ // Control[iter_C]==0
          indexNew_neutralT.push_back(iter_T);
          indexNew_neutralC.push_back(iter_C);
          }
          
          }
          
        }}
        
        //// export ////
        int count_neutral  = indexNew_neutralT.size();           
        int count_uninf  = indexNew_uninfT.size();
        
        return(List::create(
          Named("count_favorable")  = count_favorable, // 0
          Named("count_unfavorable")  = count_unfavorable, // 1
          Named("count_neutral")  = count_neutral, // 2  // indexNew_neutralT.size()        
          Named("count_uninf")  = count_uninf, // 3    // indexNew_uninfT.size()
          Named("index_neutralT")  = indexNew_neutralT, // 4 
          Named("index_neutralC")  = indexNew_neutralC, // 5
          Named("index_uninfT")  = indexNew_uninfT, // 6
          Named("index_uninfC")  = indexNew_uninfC // 7
          ));
          
}

//  fct3 : perform pairwise comparisons over all possible pairs for a continuous endpoint ///////////////////////////////
// [[Rcpp::export]]
List pairSampleCont_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> count_neutral(n_size,0); // number of neutral pairs for each sample size
  vector<int> count_uninf(n_size,0); // number of uninformative pairs for each sample size
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
    index_uninfT.push_back(iter_T);
    index_uninfC.push_back(iter_C);     
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_uninf[iter_size] = count_uninf[iter_size] + 1;
      }
    }
  }else{
    
  diff = Treatment[iter_T]-Control[iter_C];
  
  if(diff>=threshold){
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_favorable[iter_size] = count_favorable[iter_size] + 1;
      }
    }    
  }else if(diff<= -threshold){
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
      }
    }
  }else{
    index_neutralT.push_back(iter_T);
    index_neutralC.push_back(iter_C);
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
      if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
        count_neutral[iter_size] = count_neutral[iter_size] + 1;
      }
    }
  }
  
  }}}
  
  //// export ////
  return(List::create(
    Named("count_favorable")  = count_favorable, // 0
    Named("count_unfavorable")  = count_unfavorable, // 1  
    Named("count_neutral")  = count_neutral, // 2  
    Named("count_uninf")  = count_uninf, // 3      
    Named("index_neutralT")  = index_neutralT, // 4
    Named("index_neutralC")  = index_neutralC, // 5
    Named("index_uninfT")  = index_uninfT, // 6
    Named("index_uninfC")  = index_uninfC // 7
    ));
    
}

//  fct3bis : perform pairwise comparisons over a prespecified subset of pairs for a continuous endpoint//////////////////////
// [[Rcpp::export]]
List pairSampleCont2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs, 
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);      
      }else{
        
      diff = Treatment[iter_T]-Control[iter_C];
      
      if(diff>=threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_favorable[iter_size] = count_favorable[iter_size] + 1;
          }
        }  
      }else if(diff<= -threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
          }
        }        
      }else{
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
      }  
      
    }}}
    
    //// loop over the uninformative pairs ////
    if(nUninf_pairs>0){
      
      for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
        iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
        iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix 
        
        if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C);      
        }else{
        
        diff = Treatment[iter_T]-Control[iter_C];
        
        if(diff>=threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
            if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
              count_favorable[iter_size] = count_favorable[iter_size] + 1;
            }
          }      
        }else if(diff<= -threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
            if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
              count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
            }
          }       
        }else{
          indexNew_neutralT.push_back(iter_T);
          indexNew_neutralC.push_back(iter_C);
        }  
        
      }}}
      
      //// export ////
      int count_neutral  = indexNew_neutralT.size();           
      int count_uninf  = indexNew_uninfT.size();
      
      return(List::create(
        Named("count_favorable")  = count_favorable, // 0
        Named("count_unfavorable")  = count_unfavorable, // 1 
        Named("count_neutral")  = count_neutral, // 2 //          indexNew_neutralT.size()
        Named("count_uninf")  = count_uninf, // 3  // indexNew_uninfT.size()
        Named("index_neutralT")  = indexNew_neutralT, // 4
        Named("index_neutralC")  = indexNew_neutralC, // 5
        Named("index_uninfT")  = indexNew_uninfT, // 6
        Named("index_uninfC")  = indexNew_uninfC // 7
        ));
        
}


//  fct4 : perform pairwise comparisons over all possible pairs for a TTE endpoint ///////////////////////////////////////////
// [[Rcpp::export]]
List pairSampleTTE_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> count_neutral(n_size,0); // number of neutral pairs for each sample size
  vector<int> count_uninf(n_size,0); // number of uninformative pairs for each sample size
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  diff = Treatment[iter_T]-Control[iter_C];
  
  if(deltaT[iter_T]==1){
    if(deltaC[iter_C]==1){
      
      if(diff>=threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_favorable[iter_size] = count_favorable[iter_size] + 1;
          }
        }  
      }else if(diff<= -threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
          }
        }     
      }else{ 
        index_neutralT.push_back(iter_T);
        index_neutralC.push_back(iter_C);
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_neutral[iter_size] = count_neutral[iter_size] + 1;
          }
        }  
      }      
      
    }else{ // deltaC[iter_C]==0
    
    if(diff<= -threshold){
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
        if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
          count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
        }
      }         
    }else{
      index_uninfT.push_back(iter_T);
      index_uninfC.push_back(iter_C);
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_uninf[iter_size] = count_uninf[iter_size] + 1;
          }
        }
    } 
    }
    
  }else{ // deltaT[iter_T]==0
  if(deltaC[iter_C]==1){
    
    if(diff>=threshold){
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
        if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
          count_favorable[iter_size] = count_favorable[iter_size] + 1;
        }
      }    
    }else{
      index_uninfT.push_back(iter_T);
      index_uninfC.push_back(iter_C);
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_uninf[iter_size] = count_uninf[iter_size] + 1;
          }
        }
    }
    
  }else{ // deltaC[iter_C]==0
  index_uninfT.push_back(iter_T);
  index_uninfC.push_back(iter_C); 
  for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_uninf[iter_size] = count_uninf[iter_size] + 1;
          }
        }
  }
  }
  
  }}
  
  //// export ////  
  return(List::create(
    Named("count_favorable")  = count_favorable, // 0
    Named("count_unfavorable")  = count_unfavorable, // 1     
    Named("count_neutral")  = count_neutral, // 2     
    Named("count_uninf")  = count_uninf, // 3    
    Named("index_neutralT")  = index_neutralT, // 4
    Named("index_neutralC")  = index_neutralC, // 5
    Named("index_uninfT")  = index_uninfT, // 6
    Named("index_uninfC")  = index_uninfC // 7 
    ));
    
}

//  fct4bis : perform pairwise comparisons over a prespecified subset of pairs for a TTE endpoint ///////////////////////////
// [[Rcpp::export]]
List pairSampleTTE2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, const arma::colvec& deltaT, const arma::colvec& deltaC, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<int> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<int> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
      
      diff = Treatment[iter_T]-Control[iter_C];
      
      if(deltaT[iter_T]==1){
        if(deltaC[iter_C]==1){
          
          if(diff>=threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + 1;
              }
            }   
          }else if(diff<= -threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
              }
            }   
          }else{
            indexNew_neutralT.push_back(iter_T);
            indexNew_neutralC.push_back(iter_C);
          }      
          
        }else{ // deltaC[iter_C]==0
        
        if(diff<= -threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
            if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
              count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
            }
          }         
        }else{
          indexNew_uninfT.push_back(iter_T);
          indexNew_uninfC.push_back(iter_C);
        } 
        }
        
      }else{ // deltaT[iter_T]==0
      
      if(deltaC[iter_C]==1){
        
        if(diff>=threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
            if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
              count_favorable[iter_size] = count_favorable[iter_size] + 1;
            }
          }   
        }else{
          indexNew_uninfT.push_back(iter_T);
          indexNew_uninfC.push_back(iter_C);
        }     
      }else{ // deltaC[iter_C]==0
      indexNew_uninfT.push_back(iter_T);
      indexNew_uninfC.push_back(iter_C); 
      }
      }
      
    }}
    
    //// loop over the uninformative pairs ////
    if(nUninf_pairs>0){
      
      for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
        iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment and deltaT matrix
        iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control and deltaC matrix
        
        diff = Treatment[iter_T]-Control[iter_C];
        
        if(deltaT[iter_T]==1){
          if(deltaC[iter_C]==1){
            
            if(diff>=threshold){
              for(int iter_size=0 ; iter_size<n_size ; iter_size++){
                if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                  count_favorable[iter_size] = count_favorable[iter_size] + 1;
                }
              }        
            }else if(diff<= -threshold){
              for(int iter_size=0 ; iter_size<n_size ; iter_size++){
                if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                  count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
                }
              }       
            }else{
              indexNew_neutralT.push_back(iter_T);
              indexNew_neutralC.push_back(iter_C);
            }      
            
          }else{ // deltaC[iter_C]==0
          
          if(diff<= -threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
              }
            }         
          }else{
            indexNew_uninfT.push_back(iter_T);
            indexNew_uninfC.push_back(iter_C);
          } 
          }
          
        }else{ // deltaT[iter_T]==0
        
        if(deltaC[iter_C]==1){
          
          if(diff>=threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + 1;
              }
            }    
          }else{
            indexNew_uninfT.push_back(iter_T);
            indexNew_uninfC.push_back(iter_C);
          }     
        }else{ // deltaC[iter_C]==0
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C); 
        }
        }
        
      }}
      
      //// export ////
      int count_neutral  = indexNew_neutralT.size();           
      int count_uninf  = indexNew_uninfT.size();
      
      return(List::create(
        Named("count_favorable")  = count_favorable, // 0
        Named("count_unfavorable")  = count_unfavorable, // 1
        Named("count_neutral")  = count_neutral, // 2     //  indexNew_neutralT.size()     
        Named("count_uninf")  = count_uninf, // 3    //  indexNew_uninfT.size()
        Named("index_neutralT")  = indexNew_neutralT, // 4
        Named("index_neutralC")  = indexNew_neutralC, // 5
        Named("index_uninfT")  = indexNew_uninfT, // 6
        Named("index_uninfC")  = indexNew_uninfC // 7
        ));
        
}



////////////////// END PROGRAM 1 (no KMimputation) //////////////
////////////////// BEGIN PROGRAM 2 (KMimputation) //////////////


////  fct5 : perform pairwize comparison with imputation for censored data ///////////////////////////////////////////////// 
// [[Rcpp::export]]
List pairSampleMainStrataKM_cpp(const arma::mat& Treatment, const arma::mat& Control, const NumericVector& threshold, const IntegerVector& type, const arma::mat& delta_Treatment, const arma::mat& delta_Control,
int D, int n_TTE, const arma::mat& Wscheme, List& listKMT, List& listKMC, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  // WARNING : strataT and strataC should be passed as const argument but it leads to an error in the conversion to arma::uvec.
  // NOTE : each pair has an associated weight initialized at 1. The number of pairs and the total weight are two different things.
  // (ex : 3 pairs with weights 0.5 0.75 0.5 have total weight 1.75). 
  // The first endpoint begin with complete weight for each pair (i.e 1). If the pair is classed favorable or unfavorable the whole weight is affected to this category (i.e count_favorable++ or count_unfavorable++) and the pair is not used for the following outcomes.
  // If the pair is partially or completely classed uninformative or neutral, the remaining weight is used for the following endpoints.
  // EX : pair 1 is 0.15 favorable and 0.45 unfavorable for survival endpoint 1. Then the remaining 0.4 are passed to endpoint 2 that class it into favorable.
  
  //// initialization ////  
  List resK; // store the result of the pairwize comparison of each endpoint and each strata in a SEXP list
  arma::mat Mcount_favorable(D, n_size,fill::zeros); // store the total weight of favorable pairs by outcome for each strata
  arma::mat Mcount_unfavorable(D, n_size,fill::zeros); // store the total weight of unfavorable pairs by outcome for each strata

  vector<double> n_pairs(n_size,0); // number of pairs sumed over the strats
  
  int size_neutral; // number of neutral pairs (temporary)
  int size_uninf; // number of uninformative pairs (temporary)
  
  vector<int> index_neutralT(0) ; // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0) ; // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0) ; // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0) ; // index of the uninformative pairs of the control arm
  vector<double> tempo_index0; // temporary index vector used to convert the index from SEXP into vector<double>
  vector<double> tempo_index1; // temporary index vector used to convert the index from SEXP into vector<double>
  vector<double> tempo_index2; // temporary index vector used to convert the index from SEXP into vector<double>
  vector<double> tempo_index3; // temporary index vector used to convert the index from SEXP into vector<double>
  vector<double> tempo_w;  // temporary index vector used to convert the weights from SEXP into vector<double>
  
  arma::mat Wpairs_return; // matrix containing the weight of each remaining pair for each outcome for the last strata. Will be exported.
  
    int iter_dTTE=0; // number of time to event endpoints that have been used
    
    //// first endpoint
    if(type[0]==1){ // binary endpoint
    resK = pairSampleBin_cpp(Treatment.col(0),Control.col(0),
    numT, numC, sizeT, sizeC, n_size);      
    }
    if(type[0]==2){ // continuous endpoint
    resK = pairSampleCont_cpp(Treatment.col(0),Control.col(0),threshold[0],
    numT, numC, sizeT, sizeC, n_size);
    }
    if(type[0]==3){ // time to event endpoint   
    resK = pairSampleTTEKM_cpp(Treatment.col(0),Control.col(0),threshold[0],
    delta_Treatment.col(0),delta_Control.col(0),listKMT[0],listKMC[0],
    numT, numC, sizeT, sizeC, n_size);      
    iter_dTTE++; // increment the number of time to event endpoints that have been used   
    tempo_w=resK[8]; // store the weight corresponding to each pair
    }
    
    // store the total weight corresponding to each catergory after conversion from SEXP to double. 
    tempo_index0 = resK[0];
    tempo_index1 = resK[1];
    tempo_index2 = resK[2];
    tempo_index3 = resK[3];
    
    for(int iter_size=0; iter_size<n_size ; iter_size++){
      Mcount_favorable(0,iter_size) = tempo_index0[iter_size];
      Mcount_unfavorable(0,iter_size) = tempo_index1[iter_size];
      
      // add to the total number of pairs the number of pairs founded for this endpoint
      n_pairs[iter_size]=tempo_index0[iter_size]+tempo_index1[iter_size]+tempo_index2[iter_size]+tempo_index3[iter_size]; 
    }
    
      resK[2] = tempo_index2[n_size-1];
      resK[3] = tempo_index3[n_size-1];

    //// update Wpairs 
    size_neutral = as<vector<int> >(resK[4]).size(); // update the number of neutral pairs
    size_uninf = as<vector<int> >(resK[6]).size(); // update the number of uninformative pairs
    arma::mat Wpairs(size_neutral+size_uninf,1,fill::ones); // temporary matrix containing the weigth of each remaining pair for each outcome
    arma::vec w(size_neutral+size_uninf); // temporary vector containing the weight of each remaining pair to be used for the next outcome
    w.fill(1);
    if(type[0]==3){ // update the weights for the uninformative pairs in Wpairs and w
    for(int iter_uninf=0 ; iter_uninf<size_uninf ; iter_uninf++){ // neutral pairs have a weight of 1 by construction
    Wpairs(size_neutral+iter_uninf,0) = tempo_w[iter_uninf];
    w(size_neutral+iter_uninf) = tempo_w[iter_uninf];
    }
    }        
      
    //// following endpoints
    int iter_d = 0; // the index of the endpoints
    arma::mat Wpairs_sauve; // the previous update of Wpairs
    while(D>iter_d+1 && (as<double>(resK[2])>0 || as<double>(resK[3])>0)){ // loop over the following endpoints
    
    // while there are remaining endpoints and remaining neutral or uniformative pairs
    iter_d++; // increment the index of the endpoints
    Wpairs_sauve=Wpairs; // save the current Wpairs
    
    if(type[iter_d]==1){ // binary endpoint        
    resK = pairSampleBinKM2_cpp(Treatment.col(iter_d),Control.col(iter_d),
    resK[4],resK[5], size_neutral,
    resK[6],resK[7], size_uninf,
    w, numT, numC, sizeT, sizeC, n_size);    
    }
    if(type[iter_d]==2){ // continuous endpoint
    resK = pairSampleContKM2_cpp(Treatment.col(iter_d),Control.col(iter_d),threshold[iter_d],
    resK[4],resK[5], size_neutral,
    resK[6],resK[7], size_uninf,
    w, numT, numC,  sizeT, sizeC, n_size);   
    }
    if(type[iter_d]==3){ // time to event endpoint
    resK = pairSampleTTEKM2_cpp(Treatment.col(iter_d),Control.col(iter_d),threshold[iter_d],
    delta_Treatment.col(iter_dTTE),delta_Control.col(iter_dTTE),listKMT[iter_dTTE],listKMC[iter_dTTE],
    resK[4],resK[5], size_neutral,
    resK[6],resK[7], size_uninf,
    w, numT, numC, sizeT, sizeC, n_size); 
    iter_dTTE++; // increment the number of time to event endpoints that have been used
    tempo_w=resK[8]; // store the weight corresponding to each pair
    }
    
    // store the number of pairs found in each catergory after conversion from SEXP to int. 
    tempo_index0 = resK[0];
    tempo_index1 = resK[1];
     
    for(int iter_size=0; iter_size<n_size ; iter_size++){
      Mcount_favorable(iter_d,iter_size) = tempo_index0[iter_size];
      Mcount_unfavorable(iter_d,iter_size) = tempo_index1[iter_size];
     }

    // update Wpairs
    size_neutral = as<vector<int> >(resK[4]).size(); // update the number of neutral pairs
    size_uninf = as<vector<int> >(resK[6]).size(); // update the number of uninformative pairs
    Wpairs.resize(size_neutral+size_uninf,max(1,iter_dTTE)); // update the size of Wpairs
    w.resize(size_neutral+size_uninf); // update the size of w
    w.fill(1);
    tempo_index0=resK[9]; // store the position of the remaining pairs in the previous Wpairs (i.e. Wpairs_sauve)
    int iter_oldpair; // index of the remaining pair
    int n_TTE_iter_d; // number of used survival endpoints at the current iteration
    if(type[iter_d]==3){n_TTE_iter_d=iter_dTTE-1;}else{n_TTE_iter_d=iter_dTTE;} // if the last endpoint that have been processed is a TTE endpoint then the number of used survival endpoint is decreased by one to enable specific process for the last endpoint.
    
    if(iter_dTTE>0){
      for(size_t iter_pair=0; iter_pair<tempo_index0.size(); iter_pair++){
        
        iter_oldpair = tempo_index0[iter_pair]; // position of the pair in Wpairs_sauve
        
        if(n_TTE_iter_d>0){ 
          for(int iter_endpointTTE=0 ; iter_endpointTTE<n_TTE_iter_d ; iter_endpointTTE++){          
            Wpairs(iter_pair,iter_endpointTTE) = Wpairs_sauve(iter_oldpair,iter_endpointTTE); // store Wpairs_sauve in the Wpairs restrected to the remaining pairs
            if(Wscheme(iter_endpointTTE,iter_dTTE-1)==1){w(iter_pair) = w(iter_pair)*Wpairs_sauve(iter_oldpair,iter_endpointTTE);} // make the product over the endpoints of the weights associated to each remaining pair 
            // only if Wscheme is one in the column of the new endpoint and the line of the previous endpoint.
          }
        }
        if(type[iter_d]==3){ // if the last endpoint is a TTE endpoints
        Wpairs(iter_pair,iter_dTTE-1) = tempo_w[iter_pair]; // update Wpairs with tempo_w
        w(iter_pair) = w(iter_pair)*tempo_w[iter_pair]; // update w with tempo_w
        }
      }
    }
    }
 
  //// proportion in favor of treatment ////
  arma::mat delta(D, n_size,fill::zeros); // array containing for each sample size and each endpoint the proportion in favor of treatment
  for(int iter_d=0; iter_d<D ; iter_d++){ // loop over endpoints       
       for(int iter_size=0 ; iter_size < n_size ; iter_size ++){ // loop over sample size
           delta(iter_d,iter_size) = (Mcount_favorable(iter_d,iter_size)-Mcount_unfavorable(iter_d,iter_size))/(double)(n_pairs[iter_size]); // proportion in favor of treatment equals number of favorable pairs minus unfavorable pairs divided by the total number of pairs
       }  
  } 
  
  //// export ////
  return(List::create(
    Named("delta")  = delta
    ));
    
}


//  fct6 : perform a weighted pairwise comparisons over a prespecified subset of pairs for a binary endpoint ///////////////////////////
// [[Rcpp::export]]
List pairSampleBinKM2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs,
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C;  // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<double> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<double> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T); 
        indexNew_uninfC.push_back(iter_C);
        index_wUninf.push_back(iter_pairs); // index of the pair relative to Wpairs
        count_uninf=count_uninf+Wpairs(iter_pairs);
      }else if(Treatment[iter_T]==1){
        
        if(Control[iter_C]==1){
          indexNew_neutralT.push_back(iter_T);
          indexNew_neutralC.push_back(iter_C);
          index_wNeutral.push_back(iter_pairs); // index of the pair relative to Wpairs
          count_neutral=count_neutral+Wpairs(iter_pairs);
          
        }else{ // Control[iter_C]==0
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs);
          }
        }
        }        
        
      }else{ // Treatment[iter_T]==0
      
      if(Control[iter_C]==1){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
          if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
            count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs);
          }
        }        
      }else{ // Control[iter_C]==0
      indexNew_neutralT.push_back(iter_T);
      indexNew_neutralC.push_back(iter_C);
      index_wNeutral.push_back(iter_pairs); // index of the pair relative to Wpairs
      count_neutral=count_neutral+Wpairs(iter_pairs);
      }
      
      }}}
      
      //// loop over the uninformative pairs ////
      if(nUninf_pairs>0){
        
        for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
          
          iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
          iter_C = index_uninfC[iter_pairs];  // index of the control patient of the pair in the Control matrix
          
          if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
            indexNew_uninfT.push_back(iter_T); 
            indexNew_uninfC.push_back(iter_C);
            index_wUninf.push_back(nNeutral_pairs+iter_pairs);  // index of the pair relative to Wpairs 
            count_uninf=count_uninf+Wpairs(iter_pairs);
          }else if(Treatment[iter_T]==1){
            
            if(Control[iter_C]==1){
              indexNew_neutralT.push_back(iter_T);
              indexNew_neutralC.push_back(iter_C);
              index_wNeutral.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs
              count_neutral=count_neutral+Wpairs(nNeutral_pairs+iter_pairs);
            }else{ // Control[iter_C]==0
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
              }
            }
           }        
            
          }else{ // Treatment[iter_T]==0
          
          if(Control[iter_C]==1){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
              }
            }
          }else{ // Control[iter_C]==0
            indexNew_neutralT.push_back(iter_T);
            indexNew_neutralC.push_back(iter_C);
            index_wNeutral.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs
            count_neutral=count_neutral+Wpairs(nNeutral_pairs+iter_pairs);
          }
          
          }
          
        }}
        
        //// export ////
         index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // merging

          return(List::create(
          Named("count_favorable")  = count_favorable, // 0
          Named("count_unfavorable")  = count_unfavorable, // 1                   
          Named("count_neutral")  = count_neutral, // 2  
          Named("count_uninf")  = count_uninf, // 3
          Named("index_neutralT")  = indexNew_neutralT, // 4
          Named("index_neutralC")  = indexNew_neutralC, // 5
          Named("index_uninfT")  = indexNew_uninfT, // 6
          Named("index_uninfC")  = indexNew_uninfC, // 7
          Named("w")  = -1, // 8
          Named("index_w")  = index_wNeutral // 9
          ));
          
}

//  fct7 : perform a weighted pairwise comparisons over a prespecified subset of pairs for a continuous endpoint //////////////////
// [[Rcpp::export]]
List pairSampleContKM2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs, 
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<double> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<double> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  vector<int> index_wNeutral(0); // index of the neutral and uninformative pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the neutral and uninformative pairs relative to Wpairs
  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T);
        indexNew_uninfC.push_back(iter_C); 
        index_wUninf.push_back(iter_pairs); // index of the pair relative to Wpairs 
        count_uninf=count_uninf+Wpairs(iter_pairs);
      }else{
        
      diff = Treatment[iter_T]-Control[iter_C];
      
      if(diff>=threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs);
            }
        }        
      }else if(diff<= -threshold){
         for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs);
            }
        }        
      }else{
        indexNew_neutralT.push_back(iter_T);
        indexNew_neutralC.push_back(iter_C);
        index_wNeutral.push_back(iter_pairs); // index of the pair relative to Wpairs
        count_neutral = count_neutral + Wpairs(iter_pairs);
      }  
      
    }}}
    
    //// loop over the uninformative pairs ////
    if(nUninf_pairs>0){
      
      for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
        iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
        iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
        
        if(R_IsNA(Treatment[iter_T]) || R_IsNA(Control[iter_C])){
        indexNew_uninfT.push_back(iter_T);  
        indexNew_uninfC.push_back(iter_C);  
        index_wUninf.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs
        count_uninf=count_uninf+Wpairs(iter_pairs);
        }else{
          
          diff = Treatment[iter_T]-Control[iter_C];
        
        if(diff>=threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
            }
          }      
        }else if(diff<= -threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
            }
          }      
        }else{
          indexNew_neutralT.push_back(iter_T);
          indexNew_neutralC.push_back(iter_C);
          index_wNeutral.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs
          count_neutral = count_neutral + Wpairs(nNeutral_pairs+iter_pairs);
        }  
        
      }}}
      
      //// export ////
       index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end());

      return(List::create(
        Named("count_favorable")  = count_favorable, // 0
        Named("count_unfavorable")  = count_unfavorable, // 1
        Named("count_neutral")  = count_neutral, // 2
        Named("count_uninf")  = count_uninf, // 3
        Named("index_neutralT")  = indexNew_neutralT, // 4
        Named("index_neutralC")  = indexNew_neutralC, // 5
        Named("index_uninfT")  = indexNew_uninfT, // 6
        Named("index_uninfC")  = indexNew_uninfC, // 7
        Named("w")  = -1, // 8
        Named("index_w")  = index_wNeutral // 9
        ));
        
}


//  fct8 : perform a weighted pairwise comparisons over all possible pairs for a TTE endpoint //////////////////////
// [[Rcpp::export]]
List pairSampleTTEKM_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold,
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int n_Treatment=Treatment.size(); // number of patients from the treatment arm
  int n_Control=Control.size(); // number of patients from the control arm
  vector<double> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<double> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  vector<double> count_neutral(n_size,0); // number of neutral pairs for each sample size
  vector<double> count_uninf(n_size,0); // number of uninformative pairs for each sample size
  vector<int> index_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> index_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> index_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> index_uninfC(0); // index of the uninformative pairs of the control arm
  vector<double> wUninf(0); // weights of the uninformative pairs
  
  double proba_favorable; // proba of uninformative pairs to be in favor of the treatment using KM survival
  double proba_unfavorable; // proba of uninformative pairs to be in behind of the treatment using KM survival
  double proba_uninf; // proba of uninformative pairs to be uninformative using KM survival
 
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the pairs ////
  for(int iter_T=0; iter_T<n_Treatment ; iter_T++){ // over treatment patients
  for(int iter_C=0; iter_C<n_Control ; iter_C++){ // over control patients
  
  diff = Treatment[iter_T]-Control[iter_C];
  
  if(deltaT[iter_T]==1){
    if(deltaC[iter_C]==1){
      
      if(diff>=threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + 1;
            }
          }
      }else if(diff<= -threshold){
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
            }
          }
      }else{ 
        index_neutralT.push_back(iter_T);
        index_neutralC.push_back(iter_C);
        for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_neutral[iter_size] = count_neutral[iter_size] + 1;
            }
          }
      }      
      
    }else{ // deltaC[iter_C]==0
    
    if(diff<= -threshold){
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + 1;
            }
          }
    }else{
      index_uninfT.push_back(iter_T);
      index_uninfC.push_back(iter_C);
      
      if(R_IsNA(matKMT(iter_T,2)) ){
        proba_favorable=0;
        proba_unfavorable=0;   
      }else if(diff>threshold){
        proba_favorable=1-matKMT(iter_T,0)/matKMC(iter_C,1);
        proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);              
      }else{
        proba_favorable=0;
        proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);   
      }
      
      proba_uninf = 1-(proba_favorable+proba_unfavorable) ;
      wUninf.push_back(proba_uninf);  // update with the weight of the uninformative pair
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + proba_unfavorable;
                count_uninf[iter_size]=count_uninf[iter_size]+proba_uninf;
           }
      }
      
       } 
    }
    
  }else{ // deltaT[iter_T]==0
  if(deltaC[iter_C]==1){
    
    if(diff>=threshold){
       for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + 1;
            }
          }          
    }else{
      index_uninfT.push_back(iter_T);
      index_uninfC.push_back(iter_C);
      
      if(R_IsNA(matKMC(iter_C,2)) ){
        proba_favorable=0;
        proba_unfavorable=0;   
      }else if(diff< -threshold){
        proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
        proba_unfavorable=1-matKMC(iter_C,0)/matKMT(iter_T,1);
      }else{
        proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
        proba_unfavorable=0;                        
      }
      
       proba_uninf = 1-(proba_favorable+proba_unfavorable) ;
       wUninf.push_back(proba_uninf);  // update with the weight of the uninformative pair
       for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + proba_unfavorable;
                count_uninf[iter_size]=count_uninf[iter_size]+proba_uninf;
           }
      }
           
    } 
    
  }else{ // deltaC[iter_C]==0
  index_uninfT.push_back(iter_T);
  index_uninfC.push_back(iter_C);       
  
  if(diff>threshold && (R_IsNA(matKMT(iter_T,2))==false) ){
    proba_favorable=1-0.5*matKMT(iter_T,0)/matKMC(iter_C,1);
    proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);   
  }else if(diff< -threshold && (R_IsNA(matKMC(iter_C,2))==false) ){
    proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
    proba_unfavorable=1-0.5*matKMC(iter_C,0)/matKMT(iter_T,1);        
  }else if( (R_IsNA(matKMC(iter_C,2))==false) && (R_IsNA(matKMT(iter_T,2))==false) ){ 
    proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
    proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);   
  }else{
    proba_favorable=0;
    proba_unfavorable=0;   
  }
  
   proba_uninf = 1-(proba_favorable+proba_unfavorable) ;
   wUninf.push_back(proba_uninf);   // update with the weight of the uninformative pair
    for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + proba_unfavorable;
                count_uninf[iter_size]=count_uninf[iter_size]+proba_uninf;
           }
      }
   }}
  
  }}
  
  //// export ////
  return(List::create(
    Named("count_favorable")  = count_favorable, // 0
    Named("count_unfavorable")  = count_unfavorable, // 1     
    Named("count_neutral")  = count_neutral, // 2
    Named("count_uninf")  = count_uninf, // 3
    Named("index_neutralT")  = index_neutralT, // 4
    Named("index_neutralC")  = index_neutralC, // 5
    Named("index_uninfT")  = index_uninfT, // 6
    Named("index_uninfC")  = index_uninfC, // 7
    Named("w")  = wUninf // 8 
    ));
    
}

//  fct9 : perform a weighted pairwise comparisons over a prespecified subset of pairs for a TTE endpoint //////////////////
// [[Rcpp::export]]
List pairSampleTTEKM2_cpp( const arma::colvec& Treatment, const arma::colvec& Control, double threshold, 
const arma::colvec& deltaT, const arma::colvec& deltaC, const arma::mat& matKMT, const arma::mat& matKMC, 
const IntegerVector& index_neutralT, const IntegerVector& index_neutralC, int nNeutral_pairs, 
const IntegerVector& index_uninfT, const IntegerVector& index_uninfC, int nUninf_pairs,
const arma::vec& Wpairs, 
const IntegerVector& numT, const IntegerVector& numC, const IntegerVector& sizeT, const IntegerVector& sizeC, int n_size){
  
  int iter_T,iter_C; // index of the treatment / control patient of the pair in the treatment / control arm
  
  vector<double> count_favorable(n_size,0); // number of favorable pairs for each sample size
  vector<double> count_unfavorable(n_size,0); // number of unfavorable pairs for each sample size
  double count_neutral=0; // number of neutral pairs
  double count_uninf=0; // number of uninf pairs
  vector<int> indexNew_neutralT(0); // index of the neutral pairs of the treatment arm
  vector<int> indexNew_neutralC(0); // index of the neutral pairs of the control arm
  vector<int> indexNew_uninfT(0); // index of the uninformative pairs of the treatment arm
  vector<int> indexNew_uninfC(0); // index of the uninformative pairs of the control arm
  
  vector<double> wNeutral(0);  // weigth of the neutral pairs 
  vector<double> wUninf(0);  // weigth of the uninformative pairs 
  vector<int> index_wNeutral(0); // index of the neutral pairs relative to Wpairs
  vector<int> index_wUninf(0); // index of the uninformative pairs relative to Wpairs
  
  double proba_favorable; // proba of uninformative pairs to be in favor of the treatment using KM survival
  double proba_unfavorable; // proba of uninformative pairs to be in behind of the treatment using KM survival  
  double diff; // difference between the endpoints from the treatment and control patients of the pair
  
  //// loop over the neutral pairs ////
  if(nNeutral_pairs>0){
    
    for(int iter_pairs=0; iter_pairs<nNeutral_pairs ; iter_pairs++){
      iter_T = index_neutralT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
      iter_C = index_neutralC[iter_pairs]; // index of the control patient of the pair in the Control matrix
      
      diff = Treatment[iter_T]-Control[iter_C];
      
      if(deltaT[iter_T]==1){
        if(deltaC[iter_C]==1){
          
          if(diff>=threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs);
             }
          }         
          }else if(diff<= -threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs);
             }
          }
          }else{
            indexNew_neutralT.push_back(iter_T);
            indexNew_neutralC.push_back(iter_C);
            wNeutral.push_back(1); // update with the weight of the neutral pair
            index_wNeutral.push_back(iter_pairs); // index of the pair relative to Wpairs
            count_neutral = count_neutral + Wpairs(iter_pairs);
          }      
          
        }else{ // deltaC[iter_C]==0
        
        if(diff<= -threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs);
             }
          }
        }else{
          indexNew_uninfT.push_back(iter_T); // [last pairs are the last]
          indexNew_uninfC.push_back(iter_C); // [last pairs are the last]
          
          if(R_IsNA(matKMT(iter_T,2)) ){
            proba_favorable=0;
            proba_unfavorable=0;   
          }else if(diff>threshold){
            proba_favorable=1-matKMT(iter_T,0)/matKMC(iter_C,1);
            proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);            
          }else{
            proba_favorable=0;
            proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);                    
          }
          
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs)*proba_unfavorable;
             }
          }
          
          wUninf.push_back(1-(proba_favorable+proba_unfavorable)); // update with the weight of the uninformative pair [last pairs are the last]
          index_wUninf.push_back(iter_pairs); // index of the pair relative to Wpairs [uninformative pairs last]
          count_uninf = count_uninf + Wpairs(iter_pairs)*(1-(proba_favorable+proba_unfavorable));
        } 
        }
        
      }else{ // deltaT[iter_T]==0
      
      if(deltaC[iter_C]==1){
        
        if(diff>=threshold){
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs);
             }
          }      
        }else{
          indexNew_uninfT.push_back(iter_T); // [last pairs are the last]
          indexNew_uninfC.push_back(iter_C); // [last pairs are the last]
          
          if(R_IsNA(matKMC(iter_C,2)) ){
            proba_favorable=0;
            proba_unfavorable=0;   
          }else if(diff< -threshold){
            proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
            proba_unfavorable=1-matKMC(iter_C,0)/matKMT(iter_T,1);           
          }else{
            proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
            proba_unfavorable=0;                
          }
          
          for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs)*proba_unfavorable;
             }
          }
          
          wUninf.push_back(1-(proba_favorable+proba_unfavorable)); // update with the weight of the uninformative pair [last pairs are the last]     
          index_wUninf.push_back(iter_pairs); // index of the pair relative to Wpairs [uninformative pairs last]
          count_uninf = count_uninf + Wpairs(iter_pairs)*(1-(proba_favorable+proba_unfavorable));
        }     
      }else{ // deltaC[iter_C]==0
      indexNew_uninfT.push_back(iter_T); // [last pairs are the last] 
      indexNew_uninfC.push_back(iter_C); // [last pairs are the last] 
      
      if(diff>threshold && (R_IsNA(matKMT(iter_T,2))==false) ){
        proba_favorable=1-0.5*matKMT(iter_T,0)/matKMC(iter_C,1);
        proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);        
      }else if(diff< -threshold && (R_IsNA(matKMC(iter_C,2))==false) ){
        proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
        proba_unfavorable=1-0.5*matKMC(iter_C,0)/matKMT(iter_T,1);       
      }else if( (R_IsNA(matKMC(iter_C,2))==false) && (R_IsNA(matKMT(iter_T,2))==false) ){ 
        proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
        proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);        
      }else{
       proba_favorable=0;
       proba_unfavorable=0;   
      }   
      
      for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(iter_pairs)*proba_unfavorable;
             }
        }
          
      wUninf.push_back(1-(proba_favorable+proba_unfavorable)); // update with the weight of the uninformative pair [last pairs are the last]               
      index_wUninf.push_back(iter_pairs); // index of the pair relative to Wpairs [uninformative pairs last]
      count_uninf = count_uninf + Wpairs(iter_pairs)*(1-(proba_favorable+proba_unfavorable));
      }
      }
      
    }}
    
    //// loop over the uninformative pairs ////
    if(nUninf_pairs>0){
      
      for(int iter_pairs=0; iter_pairs<nUninf_pairs ; iter_pairs++){
        iter_T = index_uninfT[iter_pairs]; // index of the treatment patient of the pair in the Treatment matrix
        iter_C = index_uninfC[iter_pairs]; // index of the control patient of the pair in the Control matrix
        
        diff = Treatment[iter_T]-Control[iter_C];
        
        if(deltaT[iter_T]==1){
          if(deltaC[iter_C]==1){
            
            if(diff>=threshold){
              for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
                }
              }
             }else if(diff<= -threshold){
             for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
                }
              } 
            }else{              
              indexNew_neutralT.push_back(iter_T);
              indexNew_neutralC.push_back(iter_C);
              wNeutral.push_back(1); // update with the weight of the neutral pair
              index_wNeutral.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs
              count_neutral = count_neutral + Wpairs(nNeutral_pairs+iter_pairs);
            }      
            
          }else{ // deltaC[iter_C]==0
          
          if(diff<= -threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
                }
              }
          }else{
            indexNew_uninfT.push_back(iter_T); // [last pairs are the last]
            indexNew_uninfC.push_back(iter_C); // [last pairs are the last]
            
            if(R_IsNA(matKMT(iter_T,2)) ){
              proba_favorable=0;
              proba_unfavorable=0;   
            }else if(diff>threshold){
              proba_favorable=1-matKMT(iter_T,0)/matKMC(iter_C,1);
              proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);            
            }else{
              proba_favorable=0;
              proba_unfavorable=matKMT(iter_T,2)/matKMC(iter_C,1);                   
            }
            
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_unfavorable;
             }
        }
            wUninf.push_back(1-(proba_favorable+proba_unfavorable)); // update with the weight of the uninformative pair [uninformative pairs last]    
            index_wUninf.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs [last pairs are the last]    
            count_uninf = count_uninf + Wpairs(nNeutral_pairs+iter_pairs)*(1-(proba_favorable+proba_unfavorable));
          } 
          }
          
        }else{ // deltaT[iter_T]==0
        
        if(deltaC[iter_C]==1){
          
          if(diff>=threshold){
            for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs);
                }
              }         
          }else{
            indexNew_uninfT.push_back(iter_T); // [last pairs are the last]
            indexNew_uninfC.push_back(iter_C); // [last pairs are the last]
            
            if(R_IsNA(matKMC(iter_C,2)) ){
              proba_favorable=0;
              proba_unfavorable=0;   
            }else if(diff< -threshold){
              proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
              proba_unfavorable=1-matKMC(iter_C,0)/matKMT(iter_T,1);         
            }else{
              proba_favorable=matKMC(iter_C,2)/matKMT(iter_T,1);
              proba_unfavorable=0;                        
            }
            
              for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_unfavorable;
             }
        }            
            wUninf.push_back(1-(proba_favorable+proba_unfavorable));    // update with the weight of the uninformative pair   [uninformative pairs last] 
            index_wUninf.push_back(nNeutral_pairs+iter_pairs); // index of the pair relative to Wpairs [last pairs are the last]
            count_uninf = count_uninf + Wpairs(nNeutral_pairs+iter_pairs)*(1-(proba_favorable+proba_unfavorable));
          }     
        }else{ // deltaC[iter_C]==0
        indexNew_uninfT.push_back(iter_T); // [last pairs are the last]
        indexNew_uninfC.push_back(iter_C); // [last pairs are the last]
        
        if(diff>threshold && (R_IsNA(matKMT(iter_T,2))==false) ){
          proba_favorable=1-0.5*matKMT(iter_T,0)/matKMC(iter_C,1);
          proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);         
        }else if(diff< -threshold && (R_IsNA(matKMC(iter_C,2))==false) ){
          proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
          proba_unfavorable=1-0.5*matKMC(iter_C,0)/matKMT(iter_T,1);        
        }else if( (R_IsNA(matKMC(iter_C,2))==false) && (R_IsNA(matKMT(iter_T,2))==false) ){ 
          proba_favorable=0.5*matKMC(iter_C,2)/matKMT(iter_T,1);
          proba_unfavorable=0.5*matKMT(iter_T,2)/matKMC(iter_C,1);         
        }else{
          proba_favorable=0;
          proba_unfavorable=0;   
        }   
        
           for(int iter_size=0 ; iter_size<n_size ; iter_size++){
              if(numT[iter_T] < sizeT[iter_size] && numC[iter_C] < sizeC[iter_size]){
                count_favorable[iter_size] = count_favorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_favorable;
                count_unfavorable[iter_size] = count_unfavorable[iter_size] + Wpairs(nNeutral_pairs+iter_pairs)*proba_unfavorable;
             }
        }
        
        wUninf.push_back(1-(proba_favorable+proba_unfavorable));   // update with the weight of the uninformative pair [uninformative pairs last] 
        index_wUninf.push_back(nNeutral_pairs+iter_pairs);  // index of the pair relative to Wpairs [last pairs are the last]
        count_uninf = count_uninf + Wpairs(nNeutral_pairs+iter_pairs)*(1-(proba_favorable+proba_unfavorable));
        }
        }
        
      }}
      
      //// export ////
      wNeutral.insert(wNeutral.end(),wUninf.begin(),wUninf.end()); // merging
      index_wNeutral.insert(index_wNeutral.end(),index_wUninf.begin(),index_wUninf.end()); // merging

      return(List::create(
        Named("count_favorable")  = count_favorable, // 0
        Named("count_unfavorable")  = count_unfavorable, // 1        
        Named("count_neutral")  = count_neutral, // 2
        Named("count_uninf")  = count_uninf, // 3
        Named("index_neutralT")  = indexNew_neutralT, // 4
        Named("index_neutralC")  = indexNew_neutralC, // 5
        Named("index_uninfT")  = indexNew_uninfT, // 6
        Named("index_uninfC")  = indexNew_uninfC, // 7
        Named("w")  = wNeutral, // 8
        Named("index_w")  = index_wNeutral // 9
        ));
        
}

