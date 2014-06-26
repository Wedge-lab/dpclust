log.f.of.y <- function(y1, n1, kappa1, x) {
  #x=1 and kappa=1 causes problems
  x[x>0.999 & kappa1==1] = 0.999
  #allow kappa = 0, for mutations on deleted chromosomes
  if (class(kappa1) == 'numeric') {
    # Case input consists of vectors
    no.kappa.nonzero = length(which(kappa1!=0))
    no.subsamples = length(y1)
  } else {
    # Case input consists of matrices
    no.kappa.nonzero = rowSums(kappa1!=0)
    no.subsamples = ncol(y1)
  }
  kappa1[kappa1==0] = NA
  res = lchoose(n1, y1) + y1 * log(kappa1*x) + (n1-y1) * log(1-kappa1*x)
  if (class(res) == "numeric") { res = matrix(res, nrow=1) }
  resSums = rowSums(res,na.rm=T) * no.subsamples/no.kappa.nonzero
  # Drop the Inf and other NaNs
  return(sum(resSums[!is.nan(resSums)]))
}

calc.new.likelihood2 = function(y, n, kappa, thetas) {
  lfoy = log.f.of.y(y,n,kappa,thetas)
  new.likelihood = sum(lfoy)
  return(new.likelihood)
}

aic <- function(likelihood, num.samples, num.trees) {
  return(- 2 * likelihood + 2 * num.samples * num.trees)
}

bic <- function(likelihood, num.samples, num.trees, num.muts) {
  return(-2 * likelihood + 2 * num.samples * num.trees * log(num.muts))
}

dic <- function(likelihood, likelihood.theta.mean) {
  return(-1*(2*mean(likelihood)-mean(likelihood.theta.mean)))
}
