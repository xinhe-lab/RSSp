#' Function for simulating "true" orthogonalized effects given eigenvalues and PVE parameters
simulate_v <- function(cvec,D){
  
  p <- length(D)
  k <- length(cvec)
  dvar <- numeric(p)
  for (i in 1:k) {
    dvar <- dvar + cvec[i]*D^( - ( i - 1))
  }
  V <- rnorm(n = p,mean = 0,sd = sqrt(dvar))
  return(V)
}

simulate_vhat <- function(cvec,D){
  v <- simulate_v(cvec,D)
  p <- length(v)
  vhat <- rnorm(p,mean = D*v,sd = sqrt(D))
  return(vhat)
}
