
posterior_mean_u <- function(sigu,confound=0,dvec,Q,quh){
  varu <- sigu^2
  a <- confound
  p_D <- (dvec*varu)/(a+dvec*dvec*varu+dvec)
  return(Q%*%(p_D*quh))
}



posterior_mean_beta <- function(sigu,confound=0,dvec,se,Q,quh){
  return(se*posterior_mean_u(sigu,confound,dvec,Q,quh))
}

posterior_var_u <- function(sigu,confound=0,dvec,Q){
  varu <- sigu^2
  a <- confound
  # D**^{-1} from https://crerecombinase.github.io/PolygenicRSS/RSSp_Posterior.html
  pss_Di <- varu*(dvec+a)/(a+dvec*dvec*varu+dvec)
  return(Q%*%diag(pss_Di)%*%t(Q))
}

posterior_var_beta <- function(sigu,confound=0,dvec,se,Q){
  return(diag(se)*poster_var_u(sigu,confound,dvec,Q)%*%diag(se))
}


posterior_mean_y <- function(sigu,confound=0,dvec,se,Q,quh,x){
  return(x%*%posterior_mean_beta(sigu,confound,dvec,se,Q,quh))
}

posterior_var_y <- function(sigu,confound=0,dvec,se,Q,quh,x){
  return(x%*%posterior_var_beta(sigu,confound,dvec,se,Q)%*%t(x))
}

