\name{backward_interface}
\alias{backward_interface}
\title{Backward algorithm}
\description{Efficient implementation of backward algorithm in C++ code,
  for N data and S states.}
\usage{backward_interface(
  log_emission_mat, log_transition_mat)}
\arguments{
  \item{log_emission_mat}{N x S numeric matrix of log likelihood of
    observing each data point in each state.}
  \item{log_transition_mat}{S x S numeric matrix; log_transition_mat[i,j] is the
    log probability of going from state i to state j.}
}

\value{N x S numeric matrix of backward log likelihood.}

\author{Toby Dylan Hocking}

\examples{

##simulated data.
seg.mean.vec <- c(2, 0, -1, 0)
data.mean.vec <- rep(seg.mean.vec, each=10)
set.seed(1)
N.data <- length(data.mean.vec)
y.vec <- rnorm(N.data, data.mean.vec)
##model.
n.states <- 3
log.A.mat <- log(matrix(1/n.states, n.states, n.states))
state.mean.vec <- c(-1, 0, 1)*0.1
sd.param <- 1
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
plotHMM::backward_interface(log.emission.mat, log.A.mat)

}
