void viterbi
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *transition_ptr,
 double *initial_prob_vec,
 //inputs above, outputs below.
 double *log_max_prob_ptr,
 int *best_state_ptr,
 int *state_seq_ptr
 );
