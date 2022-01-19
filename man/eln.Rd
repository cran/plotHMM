\name{eln}
\alias{elnproduct}
\alias{elnsum}
\alias{logsumexp}
\title{Log probability arithmetic}
\description{Binary operators in log probability space, to avoid
  numerical underflow.}
\usage{
elnproduct(elnx, elny)
elnsum(elnx, elny)
logsumexp(exponents.vec)
}
\arguments{
  \item{elnx, elny, exponents.vec}{numeric vectors of log probabilities.}
}

\value{Numeric vector with one (logsumexp) or more (others) log
  probability value(s).}

\author{Toby Dylan Hocking}

\examples{

px <- c(0.1, 0.5, 0.9)
py <- c(0.001, 0.123, 0.999)
lx <- log(px)
ly <- log(py)
library(plotHMM)
elnproduct(lx, ly)
elnsum(lx, ly)
logsumexp(ly)

}

\references{\url{http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf}}
