//' evd_dnorm
//' 
//' This function computes the RSSp negative log likelihood
//' @param par a length 2 numeric vector.  `par[1]` corresponds to the sigma_u^2 parameter, and `par[2]` corresponds to the confounding parameter
//' @param dvec. The eigenvalues of the LD matrix
//' @param quh The precomputed matrix vector product `crossprod(Q,u_hat)` (passed as a vector)