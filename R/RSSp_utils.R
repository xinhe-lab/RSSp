#' Convert summary statistics to PC-space
#' @param uhat a vector of length `p` containing marginal association Z scores.  Alternatively, a  matrix with `p` rows, each column
#' representing marginal association with a different trait
#' @param Q A `q` by `p` matrix consisting of `q` eigenvectors of the LD matrix computed at p loci.
convert_quh <- function(uhat,Q){
    if(NCOL(uhat)>1)
        return(crossprod(Q,uhat))
    return(c(crossprod(Q,uhat)))
}

#' Convert chunks of summary statistics to PC-space
#' @param chunklist a vector or list specifying LD regions. each elemenent of chunklist will be passed to `uhat_uf` and
#' @param uhat_uf a function that takes a single element of `chunklist` and returns either a length p vector, or a matrix with p rows
#' @param Q_uf a function that takes a single element of `chunklist` and returns a `q` by `p` matrix
map_convert_quh <- function(chunklist,uhat_uf,q_uf,map_fun=purrr::map,...){
    pb <- dplyr::progress_estimated(length(chunklist))
    stopifnot(is.function(uhat_uf),is.function(q_uf))
    map_fun(chunklist,~{ pb$tick()$print()
        convert_quh(uhat_uf(.x),q_uf(.x))})
}



#' Function factory for obtaining PC-transformed summary statistics
#' @param uh uhat a vector of length `p` containing marginal association Z scores.  Alternatively, a  matrix with `p` rows, each column
#' representing marginal association with a different trait
#' @param chunklist a vector of length p indicating which chunk each locus is assigned to
gensplit_uh_uf <- function(uh,region_id){

    sf <- ifelse(NCOL(uh)>1,split.data.frame,split)
    uhl <- sf(uh,region_id)

    retfun <- function(idx){
        return(uhl[[as.character(idx)]])
    }
    return(retfun)
}



#' Function for obtaining pve estimate from the RSSp optimizer
#' @param quh uhat a vector of length `p` containing marginal association Z scores transformed into PC space.
#' @param dvec uhat a vector of length `p` containing eigenvalues
#' @param cvec an estimate of the variance
#' @param N the sample size
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
