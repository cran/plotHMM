\name{viterbi_interface}
\alias{viterbi_interface}
\title{Viterbi algorithm}
\description{Efficient implementation of Viterbi algorithm in C++ code,
  for N data and S states.}
\usage{viterbi_interface(
  log_emission_mat, log_transition_mat, log_initial_prob_vec)}
\arguments{
  \item{log_emission_mat}{N x S numeric matrix of log likelihood of
    observing each data point in each state.}
  \item{log_transition_mat}{S x S numeric matrix; log_transition_mat[i,j] is the
    log probability of going from state i to state j.}
  \item{log_initial_prob_vec}{S numeric vector of log probabilities of observing
    each state at the beginning of the sequence.}
}

\value{list with elements
\item{log_max_prob}{N x S numeric matrix of max log probabilities.}
\item{best_state}{N x S integer matrix. First row is
  fixed at zero, other rows indicate best states (from 1 to S).}
\item{state_seq}{N integer vector, best overall state sequence (entries
  from 1 to S).}
}

\author{Toby Dylan Hocking}

\examples{

##simulated data.
seg.mean.vec <- c(2, 0, -1, 0)
data.mean.vec <- rep(seg.mean.vec, each=2)
N.data <- length(data.mean.vec)
sd.param <- 0.1
set.seed(1)
y.vec <- rnorm(N.data, data.mean.vec, sd.param)
##model.
state.mean.vec <- unique(seg.mean.vec)
n.states <- length(state.mean.vec)
log.A.mat <- log(matrix(1/n.states, n.states, n.states))
log.pi.vec <- log(rep(1/n.states, n.states))
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
plotHMM::viterbi_interface(log.emission.mat, log.A.mat, log.pi.vec)

}
