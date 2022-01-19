#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include <math.h>
#include "multiply.h"
#include "eln.h"

void multiply 
(int N_data,
 int N_states,
 double *log_alpha_ptr,
 double *log_beta_ptr,
 //inputs above, outputs below.
 double *log_gamma_ptr
 ){
  arma::mat log_alpha_mat    //copy_aux_mem, strict(no size change)
    (log_alpha_ptr, N_data, N_states, false, true);
  arma::mat log_beta_mat
    (log_beta_ptr, N_data, N_states, false, true);
  arma::mat log_gamma_mat
    (log_gamma_ptr, N_data, N_states, false, true);
  for(int data_t=0; data_t<N_data; data_t++){
    double normalizer = -INFINITY;
    for(int state_i=0; state_i<N_states; state_i++){
      log_gamma_mat(data_t, state_i) = elnproduct
        (log_alpha_mat(data_t, state_i),
         log_beta_mat(data_t, state_i));
      normalizer = elnsum(normalizer, log_gamma_mat(data_t, state_i));
    }
    for(int state_i=0; state_i<N_states; state_i++){
      log_gamma_mat(data_t, state_i) = elnproduct
        (log_gamma_mat(data_t, state_i),
         -normalizer);
    }
  }
}
