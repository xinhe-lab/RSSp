RSSp_evd_mvd <- function(par,dvec,quh){
  return(-evd_dnorm(par=par,dvec=dvec,quh=quh))
  #  return(-sum(mapply(evd_dnorm,dvec=dvecl,quh=quhl,MoreArgs=list(par=par),SIMPLIFY=T)))
}
RSSp_evd_lmvd <- function(par,dvec,quh){
  return(-evd_dnorm(par=exp(par),dvec=dvec,quh=quh))
  #  return(-sum(mapply(evd_dnorm,dvec=dvecl,quh=quhl,MoreArgs=list(par=par),SIMPLIFY=T)))
}

RSSp_evd_mvd_grad <- function(par,dvec,quh){
  return(-evd_dnorm_grad(par=par,dvec=dvec,quh=quh))
}
RSSp_evd_lmvd_grad <- function(par,dvec,quh){
  return(-evd_dnorm_grad(par=exp(par),dvec=dvec,quh=quh))
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
  return(prep_RSSp_evd(Ql = evdRl$vectors,Dl=evdRl$values,U=U,N=N))
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





#' RSSp
#'
#' Run RSSp on one or several traits
#'
#' @param fgeneid a character vector specifying the name of the traits being analyzed
#' @param D A vector of length `p` indicating the eigenvalues obtained from an EVD of the LD matrix.
#' @param quh a `p` by `ntrait` matrix with transformed, standardized effect sizes obtained from running `prep_RSSp`
#'  (If only one trait is being analyzed, `quh` can also be a length `p` vector)
#' @param doConfound a logical indicating whether to use the two parameter model (`doConfound=T`) or the one parameter model without confounding
#' @param log_params a logical indicating whether to optimize parameters in log space or not (not recommended)
#' @param a_bounds a length two vector indicating the bounds of the parameter `a` to use when optimizing, `a` is a measure of confounding
#' @param sigu_bounds a length two vector indicating bounds in the parameter `sigu` to use when optimizing `sigu` is directly related to PVE
RSSp <- function(fgeneid=NULL,D,quh,sigu_bounds=c(1e-7,2.6),a_bounds=c(0,1),p_n,doConfound=T,log_params=F,useGradient=T){
  RSSp_estimate(
    tibble::data_frame(
      fgeneid=fgeneid,D=D,quh=quh
    ),
    sigu_bounds=sigu_bounds,a_bounds=a_bounds,
    p_n=p_n,doConfound = doConfound,
    log_params=log_params,
    useGradient=useGradient
    )
}



RSSp_run_mat_quh <- function(quh_mat_d,D,n,sigu_bounds=c(1e-07,2.6),a_bounds=c(0,1),doConfound=T,log_params=F,useGradient=T){
  fgeneids <- colnames(quh_mat_d)
  stopifnot(nrow(quh_mat_d)==length(D))
  p <- length(D)
  p_n <- p/n
  quhl <- purrr::array_branch(quh_mat_d,2)
  names(quhl) <- fgeneids
  return(purrr::imap_dfr(
    quhl,
    function(quh,fgeneid,D,p_n){
      RSSp(
        fgeneid = fgeneid,
        D = D,
        quh = quh,
        p_n = p_n,
        doConfound=doConfound,
        log_params=log_params,
        useGradient=useGradient,
        a_bounds=a_bounds,sigu_bounds=sigu_bounds
      )
    },
    D=D,
    p_n=p/n
  )
  )
}

RSSp_estimate <- function(data_df,sigu_bounds= c(1e-7,2.6),a_bounds=c(0,1),p_n,doConfound=T,log_params=F,useGradient=T){

  stopifnot(length(unique(data_df$fgeneid))==1)
  if(log_params){
    sigu_bounds <- log(sigu_bounds)
    a_bounds <- log(a_bounds+.Machine$double.eps)
    if(useGradient){
      grf <- RSSp_evd_lmvd_grad
    }else{
      grf <- NULL
    }
  }else{
    if(useGradient){
      grf <- RSSp_evd_mvd_grad
    }else{
      grf <- NULL
    }
  }
  stopifnot(!anyNA(sigu_bounds))
  stopifnot(!anyNA(a_bounds))
  minb <- c(sigu_bounds[1],
            a_bounds[1])
  maxb <- c(sigu_bounds[2],
            a_bounds[2])
  
  ctrl <- list(lmm=7)
  par0 <- runif(2,min=minb,max = maxb)
  # dvecl <- lapply(data_df,"[[","rd")
  dvec <- data_df$D
  quh <- data_df$quh
  stopifnot(!is.null(dvec),
            !is.null(quh))
  # quhl <- lapply(data_df,"[[","quh")
  if(log_params){
    if(!doConfound){
      ldat  <- optimise(f=RSSp_evd_lmvd,interval=sigu_bounds,dvec=dvec,quh=quh)
      siguv <- exp(ldat$minimum)
      av <- 0
      lnzv <- ldat$objective
      conv <-0
      convm <- "Converged"
    }else{
      ldat <- optim(par0,fn=RSSp_evd_lmvd,gr=grf,dvec=dvec,quh=quh)
      siguv <- exp(ldat$par[1])
      av <- exp(ldat$par[2])
      lnzv <- ldat$value

      conv <- ldat$convergence
      convm <- ldat$message
      if(is.null(convm)){
        convm <- "Converged"
      }

    }
  }else{
    if(!doConfound){
      ldat  <- optimise(f=RSSp_evd_mvd,interval=sigu_bounds,dvec=dvec,quh=quh)
      siguv <- ldat$minimum
      av <- 0
      lnzv <- ldat$objective

      conv <-0
      convm <- "Converged"
    }else{
      ldat  <- optim(par0,fn=RSSp_evd_mvd,gr=grf,dvec=dvec,quh=quh,lower = minb,upper = maxb,method = "L-BFGS-B",control=ctrl)
      siguv <- ldat$par[1]
      av <- ldat$par[2]
      lnzv <- ldat$value
      conv <- ldat$convergence
      convm <- ldat$message
      if(is.null(convm)){
        convm <- "Converged"
      }
    }
  }

  return(tibble::data_frame(sigu=siguv,a_hat=av,lnZ=lnzv,
                    convergence=conv,
                    message=convm,
                    fgeneid=as.character(unique(data_df[["fgeneid"]])),
                    log_params=log_params,
                    useGradient=useGradient,
                    method=ifelse(doConfound,"Confound","NoConfound"),
                    pve=p_n*sigu^2))
}

calc_pve <- function(sigu,p_n){
  return(p_n*sigu^2)
}

calc_sigu <- function(pve,p_n){
  return(sqrt(pve/p_n))
}


gen_quh_chunk_mat <- function(Zmat,evdf,gw_snpi){
  ld_snpi <- RcppEigenH5::read_ivec(evdf,"LDinfo","snp_id")
  p <-length(gw_snpi)
  subset_cols <- RcppEigenH5::match_sorted(gw_snpi,target = ld_snpi)+1
  if(length(subset_cols)==length(ld_snpi)){
    Dvec <- RcppEigenH5::read_dvec(evdf,"EVD","D")
    Q <- RcppEigenH5::read_2d_mat_h5(h5file = evdf,groupname = "EVD",dataname = "Q")
  }else{
    tR <- RcppEigenH5::read_2d_index_h5(
      h5file = evdf,
      groupname = "LD",
      dataname = "R",
      indvec = subset_cols)[subset_cols,]
    evdR <- eigen(tR)
    Dvec <- evdR$values
    Q <- evdR$vectors
  }
  stopifnot(all.equal(sum(Dvec),p))
  quh <- crossprod(Q,Zmat)
  colnames(quh) <- colnames(Zmat)
  return(list(quh=quh,D=Dvec,snp_id=gw_snpi))
}

gen_quh_chunk <- function(gwas_df,evdf){

  ld_snpi <- RcppEigenH5::read_ivec(evdf,"LDinfo","snp_id")
  gw_snpi <- gwas_df$snp_id
  p <-length(gw_snpi)
  subset_cols <- RcppEigenH5::match_sorted(gw_snpi,target = ld_snpi)+1
  if(length(subset_cols)==length(ld_snpi)){
    Dvec <- RcppEigenH5::read_dvec(evdf,"EVD","D")
    Q <- RcppEigenH5::read_2d_mat_h5(h5file = evdf,groupname = "EVD",dataname = "Q")
  }else{
    tR <- RcppEigenH5::read_2d_index_h5(
      h5file = evdf,
      groupname = "LD",
      dataname = "R",
      indvec = subset_cols)[subset_cols,]
    evdR <- eigen(tR)
    Dvec <- evdR$values
    Q <- evdR$vectors
  }
  stopifnot(all.equal(sum(Dvec),p))
  LD_info <- mutate(gwas_df,quh=crossprod(Q,Z),D=Dvec)
  return(LD_info)
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


