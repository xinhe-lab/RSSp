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

#' prep_RSSp
#' 
#' Prepare data for RSSp by performing eigenvalue decomposition on a list of LD matrices,
#' and transforming GWAS summary statistics based on the decomposition
#' @param Rl A list of symmetric positive-semidefinite matrices, with one matrix for each LD region 
#' @param U A vector or matrix of GWAS standardized effect sizes. If `U` is a matrix, there should be one column for each trait. 
#' It is generally a good idea to keep track of the traits by naming the columns
#' @return 
#'  Returns a dataframe with columns `quh`, which is the transformed effect size,`D` which are the eigenvalues and `fgeneid` 
#' which is the trait name (taken from the column names of `U`)
#'  
#prep_RSSp <- function(Rl,U,N)

#' RSSp
#' 
#' Run RSSp on one or several traits
#' 
#' @param fgeneid a character vector specifying the name of the traits being analyzed
#' @param D A vector of length `p` indicating the eigenvalues obtained from an EVD of the LD matrix.
#' @param quh a `p` by `ntrait` matrix with transformed, standardized effect sizes obtained from running `prep_RSSp`
#'  (If only one trait is being analyzed, `quh` can also be a length `p` vector)
#'  @param sigu_bounds 
RSSp <- function(fgeneid=NULL,D,quh,sigu_bounds=c(1e-7,2.6),a_bounds=c(0,1),p_n,doConfound=T,log_params=F,itermax=5,useGrad=T){
  RSSp_estimate(data_frame(fgeneid=fgeneid,rd=D,quh=quh),
                sigu_bounds=sigu_bounds,a_bounds=a_bounds,
                p_n=p_n,doConfound = doConfound,
                log_params=log_params,
                itermax=itermax,useGrad=useGrad)
  
}

RSSp_run_mat <- function(Ql_df,uh_mat,n){
  library(dplyr)
  library(tidyr)
  library(SeqArray)
  library(LDshrink)
  library(purrr)
  fgeneids <- colnames(uh_mat)
  Ql <- Ql_df$Ql
  p <- sum(sapply(Ql,nrow))

  Dl <- Ql_df$Dl
  stopifnot(sum(lengths(Dl))==p)
  p_n <- p/n
  Ql_df <- unnest(Ql_df$LDinfo) %>% mutate(snp_id=1:n())
  indl <- split(Ql_df$snp_id,Ql_df$region_id)
  quh <- gen_quh(Ql = Ql, uhmat = uh_mat ,indl = indl, rdl = Dl) %>% unnest()
  quhl <- split(quh$quh,quh$fgeneid)
  rdl <- split(quh$rd,quh$fgeneid)
  fgeneid=names(quhl)
  return(pmap(.l = list(fgeneid=names(quhl),
                        D=rdl,
                        quh=quhl),RSSp,p_n=p_n))
  
  
  
}


RSSp_run <- function(Ql_df,sumstat_df,fgeneid=NULL,n=NULL){
  library(dplyr)
  library(tidyr)
  library(SeqArray)
  library(LDshrink)
  library(purrr)
  stopifnot(!is.null(n))
  if(is.null(fgeneid)){
    fgeneid<-"1"
  }
  Ql <- Ql_df[["Ql"]]
  rdl <- Ql_df[["Dl"]]
  Ql_df <- Ql_df[["LDinfo"]] %>% unnest() %>% mutate(snp_id=1:n())
  chrom <- unique(Ql_df$chr)
  indl <- split(Ql_df$snp_id,Ql_df$region_id)  
  stopifnot(all(names(Ql)==names(rdl)))
  stopifnot(sapply(Ql,nrow)==lengths(rdl))
  uh_mat <- matrix(sumstat_df$Z,nrow = length(sumstat_df$Z))
  stopifnot(sum(lengths(rdl))==nrow(uh_mat))
  p <- nrow(uh_mat)
  p_n <- p/n
  colnames(uh_mat) <- as.character(fgeneid)
  quh <- gen_quh(Ql = Ql,uhmat = uh_mat,indl = indl,rdl = rdl) %>% unnest()
  
  return(RSSp(fgeneid = fgeneid,D = quh$rd,quh = quh$quh,p_n = p_n))

}
  

RSSp_estimate <- function(data_df,sigu_bounds= c(1e-7,2.6),a_bounds=c(0,1),p_n,doConfound=T,log_params=F,itermax=5,useGrad=T){
  library(purrr)
  stopifnot(length(unique(data_df$fgeneid))==1)
  if(log_params){
    sigu_bounds <- log(sigu_bounds)
    a_bounds <- log(a_bounds+.Machine$double.eps)
    if(useGrad){
      grf <- RSSp_evd_lmvd_grad
    }else{
      grf <- NULL
    }
  }else{
    if(useGrad){
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
  
  ctrl <- list(lmm=itermax)
  par0 <- runif(2,min=minb,max = maxb)
  # dvecl <- lapply(data_df,"[[","rd")
  dvec <- data_df$rd
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

  return(data_frame(sigu=siguv,a_hat=av,lnZ=lnzv,
                    convergence=conv,
                    message=convm,
                    fgeneid=unique(data_df[["fgeneid"]]),
                    method=ifelse(doConfound,"Confound","NoConfound"),pve=p_n*sigu^2))
}




chunk_mat_rows <- function(mat,factor,index){
  stopifnot(length(factor)==length(index))
  stopifnot(max(factor)<=nrow(mat))
  indexl <- split(index,factor)
  retmatl <- lapply(indexl,function(x,mat)return(mat[x,,drop=F]),mat=mat)
  names(retmatl) <- names(indexl)
  return(retmatl)
}

gen_quh <- function(Ql,uhmat,indl,rdl){
  ind_df <- data_frame(snp_index=indl) %>% mutate(chunki=1:n()) %>% unnest()
  uh_l <- chunk_mat_rows(uhmat,ind_df$chunki,ind_df$snp_index)
  buhl<- data_frame(chunki=names(Ql),quh=pmap(list(xm=Ql,ym=uh_l,im=indl,rd=rdl),function(xm,ym,im,rd){
    as_data_frame(crossprod(xm,ym)) %>% mutate(snp_index=im,rd=rd) %>% gather(fgeneid,quh,-snp_index,-rd)
  }))
  return(buhl)
}


gen_rss_est_df <- function(Ql,rdl,indl,uhmat){
  library(dplyr)
  library(tidyr)
  library(purrr)
  if(is.null(colnames(uhmat))){
    colnames(uhmat) <- as.character(1:ncol(uhmat))
  }
  if(is.null(names(Ql))){
    names(Ql) <- as.character(1:length(Ql))
    names(rdl) <- as.character(1:length(rdl))
    names(indl) <- as.character(1:length(indl))
  }
  ind_df <- data_frame(snp_index=indl) %>% mutate(chunki=1:n()) %>% unnest()
  uh_l <- chunk_mat_rows(uhmat,ind_df$chunki,ind_df$snp_index)
  quhl <- gen_quh(Ql,uhmat,indl,rdl)
  #  buhl<- lapply(map2(Ql,uh_l,function(x,y){crossprod(data.matrix(x),data.matrix(y))}),as_data_frame)
  #buh_df <- data_frame(buhdf=quhl,chunki=as.character(1:length(quhl))) %>% inner_join(rd_df) %>% unnest() %>% mutate(snp_index=1:n()) %>% gather(key = fgeneid,value = buh,-chunki,-snp_index,-rdl)
#  tbuh_l <-buh_df %>% split(.$fgeneid) %>% map(function(x)split(x,x$chunki))
  # res_f <- bind_rows(lapply(tbuh_l,evd_df,sigu_bounds=sigu_bounds,p_n=p/n))
  # tresults <-mutate(tparamdf,fgeneid=as.character(fgeneid)) %>% inner_join(res_f)
  return(unnest(quhl))
}
