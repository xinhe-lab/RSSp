




gen_trait_uuid <- function(){
    paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
}


#' Calculate PVE from sd of u and p/n
#'
#' @param sigu sd of u
#' @param p_n p/n
#'
#' @return
#' @export
#'
calc_pve <- function(sigu,p_n){
  return(p_n*sigu^2)
}

#' calculate 
#'
#' @param pve pve
#' @param p_n p/n
#'
#' @export
#'
calc_sigu <- function(pve,p_n){
  return(sqrt(pve/p_n))
}

#' calculate variance of u from pve and p_n
#'
#' @param pve pve
#' @param p_n p/n
#'
#' @export
#'
calc_varu <- function(pve,p_n){
  return(pve/p_n)
}



#' Modify LD matrix 
#'
#' @param R LD matrix 
#' @param z vector of z scores
#' @param ref_sample_size reference panel size
#'
#' @return matrix of same dimension as D
#' @export
#'
#' @examples
modify_R <- function(R,z,ref_sample_size){
  z_ld_weight <- 1/ref_sample_size
  cov2cor((1-z_ld_weight) * R + z_ld_weight * tcrossprod(z))
}

#' Estimate heritability using RSSp
#' 
#' The three most important arguments  
#' 
#' @param quh transformed summary statistics
#' @param D vector of eigenvalues of the reference LD matrix
#' @param sample_size sample size of GWAS
#' @param trait_id identifier for the trait
#' @param pve_bounds boundaries for pve estimation
#' @param nterms number of terms in the model, by default there is one term that corresponds to heritability
#' @param p number of variants in the dataset.  By default (and in almost all scenarios) this is the sum of the eigenvalues.
#' @param useGradient a boolean indicating whether or not to optimize using an analytic form for the gradient.  
#'
#' @export
RSSp_estimate <- function(quh,
                          D,
                          sample_size=NULL,
                          trait_id=gen_trait_uuid(),
                          pve_bounds=c(.Machine$double.eps,1-.Machine$double.eps),
                          nterms=1,
                          p=sum(D),useGradient=T){

  stopifnot(NROW(quh)>0,NROW(quh)==length(D),all(D>0))
  calc_H <- F
  if(length(sample_size)>1){
    warning("only using first element of `sample_size`")
  }
  sample_size <- unique(sample_size)
  stopifnot(length(sample_size)==1)
  if(length(p)>1){
    warning("only using the first value of `p`")
  }
  p <- unique(p)
  p_n <- (p/sample_size)
  stopifnot(length(p_n)==1,!is.null(p_n))
  varu_bound_mat <- replicate(nterms,calc_varu(pve_bounds,p_n))
  stopifnot(dim(varu_bound_mat)==c(2,nterms))
  optim_method <- "L-BFGS-B"
  if(useGradient|nterms==1){
    grf <- evd_dnorm_grad_stan
  }else{
    grf <- NULL
  }
  if(nterms==1){
    optim_method <- "Brent"
  }
  stopifnot(!anyNA(varu_bound_mat))
  par0 <- stats::runif(nterms,min=c(varu_bound_mat[1,]),max = c(varu_bound_mat[2,]))
  ldat  <- stats::optim(par=par0,fn=evd_dnorm,gr=grf,lower=c(varu_bound_mat[1,]),upper=c(varu_bound_mat[2,]),D=D,quh=quh,method=optim_method,hessian=calc_H)
  par_ret <- ldat$par
  varuv <- par_ret[1]
  siguv <- sqrt(varuv)
  lnzv <- ldat$value
  if(ldat$convergence!=0){
    warning(paste0("RSSp did not converge on trait:",trait_id," message:",ldat$message))
  }
 
  pve=estimate_pve(cvec=par_ret,D = D,quh=quh,sample_size = sample_size)
  
  
  retdf <- tibble::tibble(sigu=siguv,bias=list(tibble::tibble(term_no=seq_along(par_ret[-1]),value=par_ret[-1])),lnZ=lnzv,
                              convergence=ldat$convergence==0,
                              trait_id=as.character(trait_id),
                              nterms=nterms,
                              pve=pve)
  return(retdf)
}




