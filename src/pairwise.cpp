#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>// or <armadillo>
#include <math.h>//for INFINITY
#include "pairwise.h"
#include "eln.h"
//arma::set_cerr_stream( Rcpp::Rcerr );

void pairwise
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *log_transition_ptr,
 double *log_alpha_ptr,
 double *log_beta_ptr,
 //inputs above, outputs below.
 double *log_xi_ptr
 ){
  arma::mat log_emission_mat    //copy_aux_mem, strict(no size change)
    (log_emission_ptr, N_data, N_states, false, true);
  arma::mat log_transition_mat
    (log_transition_ptr, N_states, N_states, false, true);
  arma::mat log_alpha_mat
    (log_alpha_ptr, N_data, N_states, false, true);
  arma::mat log_beta_mat
    (log_beta_ptr, N_data, N_states, false, true);
  arma::cube log_xi_arr
    (log_xi_ptr, N_states, N_states, N_data, false, true);
  for(int data_t=0; data_t<N_data-1; data_t++){
    double normalizer = -INFINITY;
    for(int state_i=0; state_i<N_states; state_i++){
      for(int state_j=0; state_j<N_states; state_j++){
	double emission_beta = elnproduct
	  (log_emission_mat(data_t+1, state_j),
           log_beta_mat(data_t+1, state_j));
        double emission_beta_transition = elnproduct
          (emission_beta,
           log_transition_mat(state_i, state_j));
        double value = elnproduct
          (emission_beta_transition,
           log_alpha_mat(data_t, state_i));
        normalizer = elnsum(normalizer, value);
        log_xi_arr(state_i, state_j, data_t) = value;
      }
    }
    for(int state_i=0; state_i<N_states; state_i++){
      for(int state_j=0; state_j<N_states; state_j++){
        log_xi_arr(state_i, state_j, data_t) = elnproduct
          (log_xi_arr(state_i, state_j, data_t),
           -normalizer);
      }
    }
  }
}
