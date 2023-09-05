#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>// or <armadillo>
#include <math.h>//for INFINITY, 
#include "backward.h"
#include "eln.h"

void backward
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *log_transition_ptr,
 //inputs above, outputs below.
 double *log_beta_ptr
 ){
  arma::mat log_emission_mat    //copy_aux_mem, strict(no size change)
    (log_emission_ptr, N_data, N_states, false, true);
  arma::mat log_transition_mat
    (log_transition_ptr, N_states, N_states, false, true);
  arma::mat log_beta_mat
    (log_beta_ptr, N_data, N_states, false, true);
  for(int state_i=0; state_i<N_states; state_i++){
    log_beta_mat(N_data-1, state_i) = 0;
  }
  for(int data_t=N_data-2; 0<=data_t; data_t--){
    for(int state_i=0; state_i<N_states; state_i++){
      double log_beta = -INFINITY;
      for(int state_j=0; state_j<N_states; state_j++){
        double log_emission = log_emission_mat(data_t+1, state_j);
        double next_beta = log_beta_mat(data_t+1, state_j);
        double emission_beta = elnproduct(log_emission, next_beta);
        double log_transition = log_transition_mat(state_i, state_j);
        double emission_beta_transition = elnproduct
          (emission_beta, log_transition);
	log_beta = elnsum(log_beta, emission_beta_transition);
      }
      log_beta_mat(data_t,state_i) = log_beta;
    }
  }
}
