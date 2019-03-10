negtwo_loglik <- function(vhat, dvec, cvec){
# USAGE: compute -2 * log likelihood for a given K
# INPUT:
#	vhat: vector of length num_snp, vhat := Q'*uhat
#	dvec: vector of length num_snp, eigenvalues of Rhat
#	cvec: vector of length K, positive scalars
# OUTPUT:
#	obj: -2*log likelihood value

  num_c <- length(cvec)

  tmp <- dvec

  for (i in 1:num_c){
    tmp <- tmp + cvec[i] * (dvec^(3-i))
  }

  obj <- sum(log(tmp))
  obj <- obj + sum((vhat^2) / tmp)
  obj <- obj + length(vhat) * log(2*base::pi)
  return(obj)
}
