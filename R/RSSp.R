RSSp_evd_mvd <- function(par,dvec,quh){
  return(-evd_dnorm(par=par,dvec=dvec,quh=quh))
}

RSSp_evd_lmvd <- function(par,dvec,quh){
  return(-evd_dnorm(par=exp(par),dvec=dvec,quh=quh))
}

RSSp_evd_mvd_grad <- function(par,dvec,quh){
  return(-evd_dnorm_grad(par=par,dvec=dvec,quh=quh))
}
RSSp_evd_lmvd_grad <- function(par,dvec,quh){
  return(-evd_dnorm_grad(par=exp(par),dvec=dvec,quh=quh))
}

RSSp_evd_lmvd_prec <- function(par,dvec,quh){
  return(-evd_dnorm_prec(par=exp(par),dvec=dvec,quh=quh))
}

RSSp_evd_mvd_prec <- function(par,dvec,quh){
  return(-evd_dnorm_prec(par=par,dvec=dvec,quh=quh))
}


filter_pvv <- function(D_df,pvv=0,add_region_id=F){
  if(pvv==0){
    return(D_df)
  }
  stopifnot(!is.null(D_df$D))
  if(is.null(D_df[["rel_D"]])){
    if(is.null(D_df$region_id)){
        if(add_region_id){
            D_df <- dplyr::mutate(D_df,region_id=1)
        }else{
            stop("column `region_id` missing from D_df, (use `add_region_id=TRUE` to use only one block)")
        }
    }
    stopifnot(!is.null(D_df$region_id))
    D_df <- dplyr::group_by(D_df,region_id) %>% arrange(desc(D)) %>% 
    dplyr::mutate(rel_D=D/sum(D)) %>% ungroup()
  }
  stopifnot(!is.null(D_df[["rel_D"]]))
  ret_df <- D_df %>% dplyr::filter(rel_D>=pvv) %>%
    dplyr::ungroup()
  stopifnot(nrow(ret_df)>0)
  return(ret_df)
}



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



RSSp_hess <- function(par,dvec,quh){
  s <- par[1]
  a <- ifelse(length(par)==2,par[2],0)
  Hmat <- matrix(0,2,2)
  pmult <- (a+dvec^2*s+dvec-2*quh^2)/(2*(a+dvec^2*s+dvec)^3)
  
  
  Hmat[1,1] <- sum(pmult*dvec^4)
  Hmat[1,2]  <-sum(pmult*dvec^2)
  Hmat[2,1] <- Hmat[1,2]
  Hmat[2,2] <- sum(pmult)
  return(Hmat)
}

RSSp_hessi <-function(par,dvec,quh){
  s <- par[1]
  a <- ifelse(length(par)==2,par[2],0)
 
  Hmati <- matrix(0,2,2)
  
  
  
  pmult <- (2*(a+dvec^2*s+dvec)^3)/(a+dvec^2*s+dvec-2*quh^2)
  
  
  Hmati[1,1] <- sum(pmult*dvec^4)
  Hmati[1,2]  <-sum(pmult*dvec^2)
  Hmati[2,1] <- Hmati[1,2]
  Hmati[2,2] <- sum(pmult)

  return(Hmati)
}



#' Estimate the quantile of the maximum likelihood estimate

# RSSp_cov <- function(par0,dvec,quh,perc=0.95){
#
#   stopifnot(length(dvec)==length(quh))
#   Imat <- 1/(RSSp_hess(par = par0,dvec = dvec,quh = quh))
#   return(mvtnorm::qmvnorm(p = perc,mean = par0,sigma = Imat,tail = "both")$quantile)
# }

#' Title
#'
#' @param par
#' @param dvec
#' @param quh
#'
#' @return
#' @export
#'
#' @examples
RSSp_information_nc <- function(par,dvec){
  var <- par[1]
  return(matrix(0.5*sum(dvec^4/(var*dvec^2+dvec)^2),1,1))
}

RSSp_information <- function(par,dvec){
  var <- par[1]
  a <- ifelse(length(par)==2,par[2],0)
  Im <- matrix(0,2,2)
  Im[1,1] <- 0.5*sum(dvec^4/(var*dvec^2+dvec+a)^2)
  Im[2,2] <- 0.5*sum(1/(var*dvec^2+dvec+a)^2)
  Im[1,2] <- 0.5*sum(dvec^2/(var*dvec^2+dvec+a)^2)
  Im[2,1] <- Im[1,2]
  return(Im)
  
}



#' Prepare data for RSSp by transforming GWAS summary statistics based on the precomputed eigenvalue decomposition of the LD matrix
#' @param Ql A list of square matrices where each column of each matrix is an eigenvector of the LD matrix
#' @param Dl a list of vectors where each element of each vector is an eigenvalue
#' @param U A vector or matrix of GWAS standardized effect sizes. If `U` is a matrix, there should be one column for each trait.
#' It is generally a good idea to keep track of the traits by naming the columns
#' @return
#'  Returns a dataframe with columns `quh`, which is the transformed effect size,`D` which are the eigenvalues and `fgeneid`
#' which is the trait name (taken from the column names of `U`)
#'
prep_RSSp_evd <- function(Ql,Dl,U){
  tot_p <- sum(sapply(Ql,nrow))
  if(is.null(dim(U))){
    stopifnot(length(U)==tot_p)
    U <- matrix(U,nrow = length(U))
  }
  indf <- unlist(purrr::imap(Ql,function(x,i){rep(i,nrow(x))}))
  indl <- split(1:tot_p,indf)
  stopifnot(nrow(U)==tot_p)
  stopifnot(all(lengths(indl)==lengths(Dl)))

  return(gen_quh(Ql = Ql,uhmat = U,indl=indl,Dl=Dl) %>% tidyr::unnest() %>% dplyr::select(-region_id))
}


#' Prepare data for RSSp by performing eigenvalue decomposition on a list of LD matrices,
#' and transforming GWAS summary statistics based on the decomposition
#' @param Rl A list of symmetric positive-semidefinite matrices, with one matrix for each LD region
#' @param U A vector or matrix of GWAS standardized effect sizes. If `U` is a matrix, there should be one column for each trait.
#' It is generally a good idea to keep track of the traits by naming the columns
#' @return
#'  Returns a dataframe with columns `quh`, which is the transformed effect size,`D` which are the eigenvalues and `fgeneid`
#' which is the trait name (taken from the column names of `U`)
#'
prep_RSSp <- function(Rl,U,N){
  tot_p <- sum(sapply(Rl,nrow))
  if(is.null(dim(U))){
    stopifnot(length(U)==tot_p)
    U <- matrix(U,nrow = length(U))
  }
  evdRl <- purrr::transpose(lapply(Rl,eigen))
  return(prep_RSSp_evd(Ql = evdRl$vectors,Dl=evdRl$values,U=U))
}


#' RSSp
#'
#' Run RSSp on one trait

RSSp <- function(v_hat,
                 eigenvalues,
                 p_n,
                 fgeneid="1",
                 eigenvalue_cutoff=1e-3,
                 pve_bounds=c(.Machine$double.eps,1- .Machine$double.eps),
                 confound_bounds=c(0,0)){
  RSSp_estimate(
    tibble::data_frame(
      fgeneid=fgeneid,D=eigenvalues,quh=v_hat,p_n=p_n
    ),
    pve_bounds=pve_bounds,a_bounds=confound_bounds,
    p_n=p_n,eigenvalue_cutoff = eigenvalue_cutoff,
    doConfound = any(confound_bounds!=0)
    )
}



RSSp_estimate_grid <- function(v_hat,eigenvalues,
                               p_n,
                               grid_points=10,
                               sigu_grid=calc_sigu(seq(from = .Machine$double.eps,
                                                                                  to = 1- .Machine$double.eps,length.out = grid_points),p_n),
                               bias_grid=c(0),...){
  p_n <- unique(p_n)
  stopifnot(length(p_n)==1)
  varu_grid <- sigu_grid^2
  fnf <- RSSp_evd_mvd

  stopifnot(!anyNA(varu_grid))
  stopifnot(!anyNA(bias_grid))
  stopifnot(all(eigenvalues>0))
  #v_seq <- unique(seq(varu_bounds[1],varu_bounds[2],length.out = grid_points))
  #cf_seq <- unique(seq(bias_bounds[1],bias_bounds[2],length.out = grid_points))
  ret_df <- purrr::cross2(varu_grid,bias_grid) %>% 
    purrr::map(purrr::as_vector) %>% 
    purrr::map_df(~tibble::data_frame(sigu=sqrt(.x[1]),
                       pve=calc_pve(sigu,p_n),
                       varu=.x[1],
                       bias=.x[2],
                       lnZ=-evd_dnorm(par=.x,dvec=eigenvalues,quh=v_hat))) %>% return()
}



RSSp_mom <- function(data_df,p_n=NULL,doConfound=F){
  if(!doConfound){
    tlm <- lm(I(quh^2)~I(D^2)+offset(D)+0,data=data_df)
  }else{
    tlm <- lm(I(quh^2)~I(D^2)+offset(D)+0,data=data_df)
  }
  if(!is.null(data_df$p_n)){
    p_n <- unique(data_df$p_n)
  }
  stopifnot(length(p_n)==1)
  sigu <- sqrt(coef(tlm)[1])
  varu <- sigu^2
  av <- ifelse(doConfound,tlm[2],0)
  lnZ <- RSSp_evd_mvd(par = c(varu,av),dvec = data_df$D,quh = data_df$quh)
  retdf <- tibble::data_frame(sigu=sigu,bias=av,lnZ=lnZ,
                              fgeneid=as.character(data_df[["fgeneid"]][1]),
                              log_params=F,
                              method=ifelse(doConfound,"Confound","NoConfound"),
                              optim=F,
                              pve=p_n*varu)
              
}



RSSp_pvv <- function(data_df,pvv=1.0){
  
}

RSSp_estimate <- function(data_df,pve_bounds=c(.Machine$double.eps,1-.Machine$double.eps),p_n=NULL,doConfound=F,a_bounds=c(0,0),eigenvalue_cutoff=1e-3,calc_H=F,calc_var=F){
  truncate_data <- T
  useGradient=F
  log_params <- F
  
  data_df <- dplyr::filter(data_df,D>eigenvalue_cutoff)
  stopifnot(nrow(data_df)>0)
  stopifnot(all(data_df$D>0))
  stopifnot(length(unique(data_df$fgeneid))==1)
  stopifnot(length(a_bounds)==2)
  if(!is.null(data_df[["p_n"]])){
    p_n <- unique(data_df$p_n)
  }
  stopifnot(length(p_n)==1)
  sigu_bounds <- calc_sigu(pve_bounds,p_n)

  stopifnot(length(sigu_bounds)==2)
  varu_bounds <- sigu_bounds^2
  optim_method <- "L-BFGS-B"
  if(useGradient){
    if(log_params){
      grf <- RSSp_evd_lmvd_grad
    }else{
      grf <- RSSp_evd_mvd_grad
    }
  }else{
    grf <- NULL
    if(!doConfound){
      optim_method <- "Brent"
    }
  }

  if(log_params){
    fnf <- RSSp_evd_lmvd
    sigu_bounds <- log(sigu_bounds)
    varu_bounds <- log(varu_bounds)
    a_bounds <- log(a_bounds+.Machine$double.eps)
  }else{
    fnf <- RSSp_evd_mvd
  }
  stopifnot(!anyNA(sigu_bounds))
  stopifnot(!anyNA(varu_bounds))
  stopifnot(!anyNA(a_bounds))
  if(doConfound){
    minb <- c(varu_bounds[1],
              a_bounds[1])
    maxb <- c(varu_bounds[2],
              a_bounds[2])
  }else{
    minb <- c(varu_bounds[1])
    maxb <- c(varu_bounds[2])
  }
  par0 <- stats::runif(length(minb),min=minb,max = maxb)
  dvec <- data_df$D
  quh <- data_df$quh
  stopifnot(!is.null(dvec),
            !is.null(quh))
  ldat  <- stats::optim(par=par0,fn=fnf,gr=grf,lower=minb,upper=maxb,dvec=dvec,quh=quh,method=optim_method,hessian=calc_H)
  par_ret <- ldat$par
  if(log_params){
    par_ret <- exp(par_ret)
    sigu_bounds <- exp(sigu_bounds)
  }
  varuv <- par_ret[1]
  siguv <- sqrt(varuv)
  av <- ifelse(doConfound,par_ret[2],0)
  lnzv <- ldat$value
  conv <- ldat$convergence
  convm <- ldat$message
  if(is.null(convm)){
    convm <- "Converged"
  }
  retdf <- tibble::data_frame(sigu=siguv,bias=av,lnZ=lnzv,
                              convergence=conv,
                              message=convm,
                              fgeneid=as.character(data_df[["fgeneid"]][1]),
                              log_params=log_params,
                              useGradient=useGradient,
                              optim=T,
                              method=ifelse(doConfound,"Confound","NoConfound"),
                              pve=p_n*sigu^2)
  
  if(calc_H){
    Hmat <- ldat$hessian
    if(doConfound){
      Hvar <- purrr::possibly(solve,otherwise = matrix(NA_real_,2,2))(Hmat)/sqrt(length(dvec))
      sigu_var_h <- Hvar[1,1]
      bias_var_h <- Hvar[2,2]
      sigu_bias_cov_h <- Hvar[1,2]
    }else{
      Hvar <- (1/(Hmat))/sqrt(length(dvec))
      sigu_var_h <- Hvar[1,1]
      bias_var_h <- 0
      sigu_bias_cov_h <- 0
    }
    H_det <- det(Hmat)
    retdf <- dplyr::mutate(retdf,
                           sigu_var_h=sigu_var_h,
                           bias_var_h=bias_var_h,
                           sigu_bias_cov_h=sigu_bias_cov_h,
                           H_det=H_det)
  }
    if(calc_var){
      if(doConfound){
        Imat <- RSSp_information(par = ldat$par,dvec = dvec)
        sigu_var <- 1/Imat[1,1]
        bias_var <- 1/Imat[2,2]
        sigu_bias_cov <- 1/Imat[1,2]
      }else{
        Imat <- RSSp_information_nc(par=ldat$par,dvec=dvec)
        sigu_var <- 1/Imat[1,1]
        bias_var <- 0
        sigu_bias_cov <- 0
      }
      Fisher_det <- det(Imat)
      retdf <- dplyr::mutate(retdf,
                             sigu_var=sigu_var,
                             bias_var=bias_var,
                             sigu_bias_cov=sigu_bias_cov,
                             fisher_det=Fisher_det)
  }
  
  

  retdf <- dplyr::mutate(retdf,
                         sigu_min_bound=sigu_bounds[1],
                         sigu_max_bound=sigu_bounds[2],
                         confound_min_bound=ifelse(doConfound,minb[2],0),
                         confound_max_bound=ifelse(doConfound,maxb[2],0))


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


