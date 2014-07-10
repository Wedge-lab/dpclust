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
  return(resSums)
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

dic <- function(y, n, kappa, all.likelihoods, all.thetas) {
  # all.likelihoods: matrix where each row corresponds to all likelihoods for MCMC iterations until now of a single mutation
  # all.thetas: list that contains a matrix for each subsample with each row a mutation and each column the theta's for the subsample of each MCMC iteration
  if (length(dim(all.likelihoods)) == 2) {
    mean.likelihoods = rowMeans(-2 * all.likelihoods, na.rm=T)
    mean.thetas = matrix(unlist(lapply(all.thetas, function(x) rowMeans(x, na.rm=T) )), ncol=length(all.thetas))
  } else {
    # A one dimensional array was given. rowMeans therefore does not work, while the values in the array also represent the mean.
    mean.likelihoods = -2 * all.likelihoods
    mean.thetas = matrix(unlist(lapply(all.thetas, function(x) x )), ncol=length(all.thetas))
  }
  
  likelihoods.mean.thetas = log.f.of.y(y,n,kappa,mean.thetas)
  return(sum(-2 * likelihoods.mean.thetas + 2*log(1) + 2*(mean.likelihoods + 2*likelihoods.mean.thetas)))
}