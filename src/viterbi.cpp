#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>// or <armadillo>
#include <math.h>//for INFINITY
#include "viterbi.h"
#include "eln.h"

void viterbi
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *log_transition_ptr,
 double *log_initial_prob_ptr,
 //inputs above, outputs below.
 double *log_max_prob_ptr,
 int *best_state_ptr,
 int *state_seq_ptr
 ){
  arma::mat log_emission_mat    //copy_aux_mem, strict(no size change)
    (log_emission_ptr, N_data, N_states, false, true);
  arma::mat log_transition_mat
    (log_transition_ptr, N_states, N_states, false, true);
  arma::vec log_initial_prob_vec
    (log_initial_prob_ptr, N_states, false, true);
  arma::mat log_max_prob_mat
    (log_max_prob_ptr, N_data, N_states, false, true);
  arma::Mat<int> best_state_mat
    (best_state_ptr, N_data, N_states, false, true);
  arma::Col<int> state_seq_vec
    (state_seq_ptr, N_data, false, true);
  double best_log_prob, candidate_log_prob;
  int best_state = -2; //not used but avoids compiler warning.
  for(int data_t=0; data_t<N_data; data_t++){
    for(int state_j=0; state_j<N_states; state_j++){
      if(data_t == 0){
        best_log_prob = log_initial_prob_vec(state_j);
        best_state = -1;
      }else{
        best_log_prob = -INFINITY;
        for(int state_i=0; state_i<N_states; state_i++){
          candidate_log_prob = elnproduct
            (log_max_prob_mat(data_t-1, state_i),
             log_transition_mat(state_i, state_j));
          if(best_log_prob < candidate_log_prob){
            best_log_prob = candidate_log_prob;
            best_state = state_i;
          }
        }
      }
      log_max_prob_mat(data_t, state_j) = elnproduct
        (best_log_prob,
         log_emission_mat(data_t, state_j));
      best_state_mat(data_t, state_j) = best_state;
    }
  }
  best_log_prob = -INFINITY;
  for(int state_i=0; state_i<N_states; state_i++){
    candidate_log_prob = log_max_prob_mat(N_data-1, state_i);
    if(best_log_prob < candidate_log_prob){
      best_log_prob = candidate_log_prob;
      best_state = state_i;
    }
  }
  state_seq_vec(N_data-1) = best_state;
  for(int data_t=N_data-2; 0<=data_t; data_t--){
    state_seq_vec(data_t) = best_state_mat
      (data_t+1, state_seq_vec(data_t+1));
  }
}
