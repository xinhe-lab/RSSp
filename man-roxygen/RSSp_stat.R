#' @param par a (up to) length 2 numeric vector.  `par[1]` corresponds to the sigma_u^2 parameter, and `par[2]` corresponds to the confounding parameter. 
#' If the model does not include confounding, `par` can be a length one vector.
#' @param D. The eigenvalues of the LD matrix, passed as a length `p` vector
#' @param quh The precomputed matrix vector product `crossprod(Q,u_hat)` (passed as a  length `p` vector)
#'
