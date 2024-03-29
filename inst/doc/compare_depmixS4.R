## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
if(
  require(data.table) && 
  requireNamespace("neuroblastoma") && 
  require(ggplot2) && 
  requireNamespace("depmixS4")
){
  data(neuroblastoma, package="neuroblastoma")
  nb.dt <- data.table(neuroblastoma$profiles)
  one.pro <- nb.dt[profile.id=="4" & chromosome%in%1:10]
  ntimes <- rle(as.integer(one.pro$chromosome))
  n.states <- 4
  model.spec <- depmixS4::depmix(
    logratio ~ 1, data=one.pro,
    nstates=n.states, ntimes=ntimes$lengths)
  set.seed(1)
  unconstrained.fit <- depmixS4::fit(model.spec)
  param.names <- c(mean="(Intercept)", sd="sd")
  par.vec <- depmixS4::getpars(unconstrained.fit)
  matrix(
    par.vec[names(par.vec) %in% param.names],
    ncol=length(param.names),
    byrow=TRUE,
    dimnames=list(state=1:n.states, parameter=names(param.names)))
  one.pro[, viterbi := factor(unconstrained.fit@posterior[,1]) ]
  ggplot()+
    geom_point(aes(
      position/1e6, logratio, color=viterbi),
      data=one.pro)+
    facet_grid(. ~ chromosome, scales="free", space="free")
}

## -----------------------------------------------------------------------------
if(requireNamespace("depmixS4")){
  par.vec <- depmixS4::getpars(model.spec)
  equal.groups <- rep(1, length(par.vec))
  equal.groups[names(par.vec)=="sd"] <- 2
  if(FALSE){
    constrained.fit <- depmixS4::fit(model.spec, equal=equal.groups)
  }
}

## -----------------------------------------------------------------------------
if(requireNamespace("depmixS4")){
  one.chrom <- nb.dt[profile.id=="4" & chromosome=="2"]
  n.states <- 3
  model.spec <- depmixS4::depmix(
    logratio ~ 1,
    data=one.chrom,
    nstates=n.states)
  log.emission.mat <- log(model.spec@dens[,1,])
  log.transition.mat <- log(model.spec@trDens[1,,])
  log.init.vec <- log(model.spec@init[1,])
  microbenchmark::microbenchmark(depmixS4={
    result <- depmixS4::forwardbackward(model.spec)
  }, plotHMM={
    fwd.list <- plotHMM::forward_interface(
      log.emission.mat, log.transition.mat, log.init.vec)
    back.mat <- plotHMM::backward_interface(
      log.emission.mat, log.transition.mat)
    mult.mat <- plotHMM::multiply_interface(
      fwd.list$log_alpha, back.mat)
    pairwise.array <- plotHMM::pairwise_interface(
      log.emission.mat, log.transition.mat,
      fwd.list$log_alpha, back.mat)
  }, times=5)  
}

## -----------------------------------------------------------------------------
if(requireNamespace("depmixS4")){
  microbenchmark::microbenchmark(depmixS4={
    depmixS4::viterbi(model.spec)
  }, plotHMM={
    plotHMM::viterbi_interface(
      log.emission.mat, log.transition.mat, log.init.vec)
  }, times=5)
}

