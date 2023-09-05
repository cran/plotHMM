## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
data(neuroblastoma, package="neuroblastoma")
library(data.table)
nb.dt <- data.table(neuroblastoma[["profiles"]])
nb.dt[, data.i := rank(position), keyby=.(profile.id, chromosome)]
(pro.dt <- nb.dt[profile.id=="4" & chromosome %in% paste(1:10)])
library(ggplot2)
ggplot()+
  facet_grid(. ~ chromosome, scales="free", space="free")+
  scale_x_continuous(
    "Position/index in data sequence",
    breaks=seq(0, 1000, by=100))+
  scale_y_continuous(
    "logratio (data values)")+
  geom_point(aes(
    data.i, logratio),
    data=pro.dt)

## -----------------------------------------------------------------------------
pro.dt[, row := 1:.N]
all.y.vec <- pro.dt[["logratio"]]
data.list <- split(pro.dt, paste(pro.dt[["chromosome"]]))
first.row.vec <- pro.dt[data.i==1, row]
last.row.vec <- pro.dt[, .SD[data.i==.N], by=chromosome][["row"]]
n.states <- 4
log.A.mat <- log(matrix(1/n.states, n.states, n.states))
set.seed(1)
mean.vec <- rnorm(n.states)
sd.param <- 1
log.pi.vec <- log(rep(1/n.states, n.states))

## -----------------------------------------------------------------------------
all.log.gamma.mat <- matrix(NA, nrow(pro.dt), n.states)
all.log.xi.array <- array(NA, c(n.states, n.states, nrow(pro.dt)))
for(chrom.i in seq_along(data.list)){
  one.chrom <- data.list[[chrom.i]]
  row.vec <- one.chrom[["row"]]
  y.vec <- one.chrom[["logratio"]]
  N.data <- length(y.vec)
  log.emission.mat <- dnorm(
    matrix(y.vec, N.data, n.states, byrow=FALSE),
    matrix(mean.vec, N.data, n.states, byrow=TRUE),
    sd.param,
    log=TRUE)
  fwd.list <- plotHMM::forward_interface(
    log.emission.mat, log.A.mat, log.pi.vec)
  log.alpha.mat <- fwd.list[["log_alpha"]]
  log.beta.mat <- plotHMM::backward_interface(log.emission.mat, log.A.mat)
  all.log.gamma.mat[row.vec, ] <- plotHMM::multiply_interface(
    log.alpha.mat, log.beta.mat)
  all.log.xi.array[,, row.vec[-N.data] ] <- plotHMM::pairwise_interface(
    log.emission.mat, log.A.mat, log.alpha.mat, log.beta.mat)
}

## -----------------------------------------------------------------------------
(new.log.pi.vec <- apply(
  all.log.gamma.mat[first.row.vec,]-log(length(first.row.vec)),
  2, plotHMM::logsumexp))

## -----------------------------------------------------------------------------
prob.mat <- exp(all.log.gamma.mat)
(new.mean.vec <- colSums(all.y.vec*prob.mat)/colSums(prob.mat))
resid.mat <- all.y.vec-matrix(
  new.mean.vec, length(all.y.vec), n.states, byrow=TRUE)
var.est <- sum(prob.mat * resid.mat^2) / sum(prob.mat)
(new.sd.param <- sqrt(var.est))

## -----------------------------------------------------------------------------
new.log.A.mat <- plotHMM::transition_interface(
  all.log.gamma.mat[-last.row.vec,],
  all.log.xi.array[,, -last.row.vec])

