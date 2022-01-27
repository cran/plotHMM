#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>// or <armadillo>
#include <math.h>//for INFINITY
#include "transition.h"
#include "eln.h"

void transition
(int N_transitions,
 int N_states,
 double *log_gamma_ptr,
 double *log_xi_ptr,
 //inputs above, outputs below.
 double *transition_ptr
 ){
  arma::mat log_gamma_mat    //copy_aux_mem, strict(no size change)
    (log_gamma_ptr, N_transitions, N_states, false, true);
  arma::cube log_xi_arr
    (log_xi_ptr, N_states, N_states, N_transitions, false, true);
  arma::mat log_transition_mat
    (transition_ptr, N_states, N_states, false, true);
  for(int state_i=0; state_i<N_states; state_i++){
    for(int state_j=0; state_j<N_states; state_j++){
      double numerator = -INFINITY;
      double denominator = -INFINITY;
      for(int data_t=0; data_t<N_transitions; data_t++){
        numerator = elnsum(numerator, log_xi_arr(state_i, state_j, data_t));
        denominator = elnsum(denominator, log_gamma_mat(data_t, state_i));
      }
      log_transition_mat(state_i, state_j) =
        elnproduct(numerator, -denominator);
    }
  }
}
