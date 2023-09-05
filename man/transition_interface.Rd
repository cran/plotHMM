\name{transition_interface}
\alias{transition_interface}
\title{Transition algorithm}
\description{Efficient implementation of transition algorithm in C++ code,
  for T transitions and S states.}
\usage{transition_interface(
  log_gamma_mat, log_xi_array)}
\arguments{
  \item{log_gamma_mat}{
    T x S numeric matrix, taken by removing the last row from the log
    probabilities from multiply.}
  \item{log_xi_array}{S x S x T numeric array of log probabilities
    from pairwise.}
}

\value{
  S x S numeric array of log probabilities (new transition matrix).
}

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
log.pi.vec <- log(rep(1/n.states, n.states))
f.list <- plotHMM::forward_interface(log.emission.mat, log.A.mat, log.pi.vec)
b.mat <- plotHMM::backward_interface(log.emission.mat, log.A.mat)
log.gamma.mat <- plotHMM::multiply_interface(f.list$log_alpha, b.mat)
prob.mat <- exp(log.gamma.mat)
log.xi.array <- plotHMM::pairwise_interface(
  log.emission.mat, log.A.mat, f.list$log_alpha, b.mat)
plotHMM::transition_interface(log.gamma.mat[-N.data,], log.xi.array)

}
