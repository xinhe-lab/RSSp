



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param D_df
##' @param normD
##' @param add_region_id
##' @return
##' @author
##' @export
filter_normD <- function(D_df,normD=0,add_region_id=F){
  if(pvv==0){
    return(D_df)
  }
  stopifnot(!is.null(D_df$D))
  if(is.null(D_df[["norm_D"]])){
    if(is.null(D_df$region_id)){
      if(add_region_id){
        D_df <- dplyr::mutate(D_df,region_id=1)
      }else{
        stop("column `region_id` missing from D_df, (use `add_region_id=TRUE` to use only one block)")
      }
    }
    stopifnot(!is.null(D_df$region_id))
    D_df <- dplyr::group_by(D_df,region_id) %>% mutate(norm_D=D/mean(D)) %>% ungroup()
  }
  stopifnot(!is.null(D_df[["norm_D"]]))
  ret_df <- D_df %>% dplyr::filter(norm_D>=norm_D) %>%
    dplyr::ungroup()
  stopifnot(nrow(ret_df)>0)
  return(ret_df)
}




gen_trait_uuid <- function(){
    paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
}


#' Estimate heritability using RSSp
#' @param quh
#' @param D
#' @param sample_size
#' @param trait_id
#' @param pve_bounds
#' @param nterms
#' @param eigenvalue_cutoff
#' @param calc_H
#' @param alt_pve
#' @param useGradient
#' @param save_data
#' @export
RSSp_estimate <- function(quh,
                          D,
                          sample_size=NULL,
                          trait_id=gen_trait_uuid(),
                          pve_bounds=c(.Machine$double.eps,1-.Machine$double.eps),
                          nterms=1,
                          p=sum(D),
                          eigenvalue_cutoff=1e-3,
                          calc_H=F,
                          alt_pve=F,useGradient=T){
  #useGradient <- T
  if(NCOL(quh)>1){
    quh <- quh[D>eigenvalue_cutoff,]
  }else{
    quh <- quh[D>eigenvalue_cutoff]
  }
  D <- D[D>eigenvalue_cutoff]
  #data_df <- data_frame(quh=quh,D=D) %>% dplyr::filter(D>eigenvalue_cutoff)
  stopifnot(NROW(quh)>0,NROW(quh)==length(D),all(D>0))
  #stopifnot(all(data_df$D>0))
  sample_size <- unique(sample_size)
  stopifnot(length(sample_size)==1)
  p <- unique(p)
  p_n <- (p/sample_size)
   stopifnot(length(p_n)==1,!is.null(p_n))
  varu_bound_mat <- replicate(nterms,calc_varu(pve_bounds,p_n))
  #sigu_bound_mat <- calc_sigu(pve_bounds,p_n)
  stopifnot(dim(varu_bound_mat)==c(2,nterms))
  # varu_bounds <- sigu_bounds^2
  optim_method <- "L-BFGS-B"
  if(useGradient){
    grf <- evd_dnorm_grad_stan
    if(nterms==1){
      optim_method <- "Brent"
    }
  }else{
    grf <- NULL
  }
  stopifnot(!anyNA(varu_bound_mat))
  par0 <- stats::runif(nterms,min=c(varu_bound_mat[1,]),max = c(varu_bound_mat[2,]))
  ldat  <- stats::optim(par=par0,fn=evd_dnorm,gr=grf,lower=c(varu_bound_mat[1,]),upper=c(varu_bound_mat[2,]),D=D,quh=quh,method=optim_method,hessian=calc_H)
  par_ret <- ldat$par
  varuv <- par_ret[1]
  siguv <- sqrt(varuv)
  lnzv <- ldat$value
  conv <- ldat$convergence
  convm <- ldat$message
  if(is.null(convm)){
      convm <- "Converged"
  }
  pve=estimate_pve(cvec=par_ret,D = D,quh=quh,sample_size = sample_size)
  
  
  retdf <- tibble::tibble(sigu=siguv,bias=list(tibble::data_frame(term_no=seq_along(par_ret[-1]),value=par_ret[-1])),lnZ=lnzv,
                              convergence=conv,
                              trait_id=as.character(trait_id),
                              nterms=nterms,
                              pve=pve)
  return(retdf)
}

calc_pve <- function(sigu,p_n){
  return(p_n*sigu^2)
}

calc_sigu <- function(pve,p_n){
  return(sqrt(pve/p_n))
}

calc_varu <- function(pve,p_n){
  return(pve/p_n)
}


quh_mat <- function(Ql,uhmat){
  ret_quh_mat <- block_mat_mul(Ql,uhmat,transpose_mat_l = T)
  colnames(ret_quh_mat) <- colnames(uhmat)
  return(ret_quh_mat)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param Ql
##' @param uhmat
##' @param indl
##' @param Dl
##' @return
##' @author
gen_quh <- function(Ql,uhmat,indl,Dl){
  uh_l <- lapply(indl,function(x,mat)return(mat[x,,drop=F]),mat=uhmat)
  buhl<- tibble::data_frame(
    region_id=names(Ql),
    quh=purrr::pmap(
      list(xm=Ql,ym=uh_l,im=indl,rd=Dl),
      function(xm,ym,im,rd){
        tibble::as_data_frame(crossprod(xm,ym)) %>%
          dplyr::mutate(snp_index=im,D=rd) %>%
          tidyr::gather(fgeneid,quh,-snp_index,-D)
      }
    )
  )
  return(buhl)
}


