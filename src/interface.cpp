#include <Rcpp.h>
#include "forward.h"
#include "backward.h"
#include "multiply.h"
#include "pairwise.h"
#include "transition.h"
#include "viterbi.h"

// [[Rcpp::export]]
Rcpp::List forward_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericMatrix log_transition_mat,
 Rcpp::NumericVector log_initial_prob_vec
 ) {
  int N_data = log_emission_mat.nrow();
  int N_states = log_emission_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_emission_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_emission_mat must have at least one col");
  }
  if(log_transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(log_transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(log_initial_prob_vec.length() != N_states){
    Rcpp::stop("length of log_initial_prob_vec must be same as number of columns of log_emission_mat");
  }
  Rcpp::NumericMatrix log_alpha_mat(N_data, N_states);
  Rcpp::NumericVector log_lik_vec(1);
  forward
    (N_data,
     N_states,
     &log_emission_mat[0],
     &log_transition_mat[0],
     &log_initial_prob_vec[0],
     //inputs above, outputs below.
     &log_alpha_mat[0],
     &log_lik_vec[0]
     );
  return Rcpp::List::create
    (Rcpp::Named("log_alpha", log_alpha_mat),
     Rcpp::Named("log_lik", log_lik_vec));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix backward_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericMatrix log_transition_mat
 ) {
  int N_data = log_emission_mat.nrow();
  int N_states = log_emission_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_emission_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_emission_mat must have at least one col");
  }
  if(log_transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(log_transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  Rcpp::NumericMatrix log_beta_mat(N_data, N_states);
  backward
    (N_data,
     N_states,
     &log_emission_mat[0],
     &log_transition_mat[0],
     //inputs above, outputs below.
     &log_beta_mat[0]
     );
  return log_beta_mat;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix multiply_interface
(Rcpp::NumericMatrix log_alpha_mat,
 Rcpp::NumericMatrix log_beta_mat
 ) {
  int N_data = log_alpha_mat.nrow();
  int N_states = log_alpha_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_alpha_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_alpha_mat must have at least one col");
  }
  if(log_beta_mat.nrow() != N_data){
    Rcpp::stop("nrow(log_beta_mat) must be same as nrow(log_alpha_mat)");
  }
  if(log_beta_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_beta_mat) must be same as ncol(log_alpha_mat)");
  }
  Rcpp::NumericMatrix log_gamma_mat(N_data, N_states);
  multiply
    (N_data,
     N_states,
     &log_alpha_mat[0],
     &log_beta_mat[0],
     //inputs above, outputs below.
     &log_gamma_mat[0]
     );
  return log_gamma_mat;
}

// [[Rcpp::export]]
Rcpp::NumericVector pairwise_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericMatrix log_transition_mat,
 Rcpp::NumericMatrix log_alpha_mat,
 Rcpp::NumericMatrix log_beta_mat
 ) {
  int N_data = log_alpha_mat.nrow();
  int N_states = log_alpha_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_alpha_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_alpha_mat must have at least one col");
  }
  if(log_beta_mat.nrow() != N_data){
    Rcpp::stop("nrow(log_beta_mat) must be same as nrow(log_alpha_mat)");
  }
  if(log_beta_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_beta_mat) must be same as ncol(log_alpha_mat)");
  }
  if(log_emission_mat.nrow() != N_data){
    Rcpp::stop("nrow(log_emission_mat) must be same as nrow(log_alpha_mat)");
  }
  if(log_emission_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_emission_mat) must be same as ncol(log_alpha_mat)");
  }
  if(log_transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(log_transition_mat) must be same as ncol(log_alpha_mat)");
  }
  if(log_transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_transition_mat) must be same as ncol(log_alpha_mat)");
  }
  Rcpp::NumericVector log_xi_arr(N_states*N_states*(N_data-1));
  log_xi_arr.attr("dim") = Rcpp::IntegerVector::create
    (N_states, N_states, N_data-1);
  pairwise
    (N_data,
     N_states,
     &log_emission_mat[0],
     &log_transition_mat[0],
     &log_alpha_mat[0],
     &log_beta_mat[0],
     //inputs above, outputs below.
     &log_xi_arr[0]
     );
  return log_xi_arr;
}

// [[Rcpp::export]]
Rcpp::NumericVector transition_interface
(Rcpp::NumericMatrix log_gamma_mat,
 Rcpp::NumericVector log_xi_array
 ) {
  int N_transitions = log_gamma_mat.nrow();
  int N_states = log_gamma_mat.ncol();
  if(N_transitions < 1){
    Rcpp::stop("log_gamma_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_gamma_mat must have at least one col");
  }
  if(log_xi_array.length() != N_states * N_states * N_transitions){
    Rcpp::stop("length(log_xi_array) must be S x S x N where N=nrow(log_gamma_mat) and S=ncol(log_gamma_mat)");
  }
  Rcpp::NumericMatrix log_transition_mat(N_states, N_states);
  transition
    (N_transitions,
     N_states,
     &log_gamma_mat[0],
     &log_xi_array[0],
     //inputs above, outputs below.
     &log_transition_mat[0]
     );
  return log_transition_mat;
}

// [[Rcpp::export]]
Rcpp::List viterbi_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericMatrix log_transition_mat,
 Rcpp::NumericVector log_initial_prob_vec
 ){
  int N_data = log_emission_mat.nrow();
  int N_states = log_emission_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_emission_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_emission_mat must have at least one col");
  }
  if(log_transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(log_transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(log_transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(log_initial_prob_vec.length() != N_states){
    Rcpp::stop("length of log_initial_prob_vec must be same as number of columns of log_emission_mat");
  }
  Rcpp::NumericMatrix log_max_prob_mat(N_data, N_states);
  Rcpp::IntegerMatrix best_state_mat(N_data, N_states);
  Rcpp::IntegerVector state_seq_vec(N_data);
  viterbi
    (N_data,
     N_states,
     &log_emission_mat[0],
     &log_transition_mat[0],
     &log_initial_prob_vec[0],
     //inputs above, outputs below.
     &log_max_prob_mat[0],
     &best_state_mat[0],
     &state_seq_vec[0]
     );
  return Rcpp::List::create
    (Rcpp::Named("log_max_prob", log_max_prob_mat),
     Rcpp::Named("best_state", best_state_mat+1),
     Rcpp::Named("state_seq", state_seq_vec+1));
}
