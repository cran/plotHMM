logsumexp <- function(exponents.vec){
  m <- max(exponents.vec)
  m + log(sum(exp(exponents.vec-m)))
}

elnsum <- function(elnx, elny){
  ifelse(
    elnx == -Inf, elny, 
  ifelse(
    elny == -Inf, elnx,
  ifelse(
    elny < elnx,
    elnx + log(1+exp(elny-elnx)),
    elny + log(1+exp(elnx-elny)))))
}

elnproduct <- function(elnx, elny){
  ifelse(
    elnx == -Inf | elny == -Inf,
    -Inf,
    elnx+elny)
}
