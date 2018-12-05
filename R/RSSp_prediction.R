
posterior_mean_u <- function(sigu,dvec,Q,quh){
  varu <- sigu^2
  # p_D <- (dvec*varu)/(a+dvec*dvec*varu+dvec)
  p_D <- 1/(dvec+1/varu)
  return(Q%*%(p_D*quh))
}

posterior_mean_v <- function(sigu,dvec,quh){
  return((1/(dvec+1/(sigu^2)))*quh)
}


posterior_mean_beta <- function(sigu,dvec,se,Q,quh){
  return(se*posterior_mean_u(sigu,dvec,Q,quh))
}


estimate_pve <- function(dvec,quh,cvec,N,n_samples=0){
  
  num_c <- length(cvec)
  
  if(num_c>1){
    stop("multiple cvec terms not yet implemented")
  }
  if(n_samples!=0){
    stop("sampling based pve estimate not yet implemented")
  }
  # lambda_init <- 0
  
  
  tmp <- 1+1/(cvec*dvec)
  rto <- (quh^2)/((dvec)^2)
  rto <- dvec*rto
  rto <- rto/(tmp^2)
  
  pve_mean <- sum(1/tmp)+sum(rto)
  pve_mean <- pve_mean/N
  
  
  return(pve_mean)
  
}


# 
# d_sigma <- function(sigu,dvec){
#   # D_\Sigma from https://crerecombinase.github.io/PolygenicRSS/RSSp_Posterior.html
#   return(((dvec)*sigu^2)/((dvec+confound)+dvec^2*sigu^2))
# }
# 
# posterior_var_u <- function(sigu,vec,Q){
#   pss_Di <- d_sigma(sigu,dvec)
#   return(Q%*%diag(pss_Di)%*%t(Q))
# }
# 
# posterior_var_beta <- function(sigu,confound=0,dvec,se,Q){
#   return(diag(se)*posterior_var_u(sigu,confound,dvec,Q)%*%diag(se))
# }
# 
# 
# posterior_mean_y <- function(sigu,confound=0,dvec,se,Q,quh,x){
#   return(x%*%posterior_mean_beta(sigu,confound,dvec,se,Q,quh))
# }
# 
# posterior_var_y <- function(sigu,confound=0,dvec,se,Q,quh,x){
#   return(x%*%posterior_var_beta(sigu,confound,dvec,se,Q)%*%t(x))
# }

