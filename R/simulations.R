#Code for generating simulations


gen_tparamdf_norm <- function(pve, bias = 0, nreps, n, p,fgeneid=NULL){
  library(dplyr)
  tfgeneid <- fgeneid
  stopifnot(nreps>0)
  if(length(nreps)>1){
    warning("length(nreps)>1, only using first element")
    nreps <- nreps[1]
  }
  tfp <- list(tpve=pve, tbias=bias) %>% purrr::cross_df() %>%
    distinct()
  rfp <- bind_rows(replicate(nreps, tfp, simplify = F)) %>%
    group_by(tpve, tbias) %>%
    mutate(replicate=1:n()) %>%
    ungroup() %>%
    mutate(fgeneid=1:n(), tsigu=sqrt(n/p*tpve),tbias=tpve*tbias) %>% 
    select(-replicate)
  if(!is.null(tfgeneid)){
    stopifnot(nrow(rfp)==length(tfgeneid))
    rfp <- mutate(rfp,fgeneid=tfgeneid)
  }
  return(rfp)
}



gen_ti <- function(vy, PVE){
  return(vy*(1/PVE-1))
}

parse_ldsc_h2log <- function(h2lf){
  library(purrr)
  library(tidyr)
  h2_dat <- scan(h2lf,what=character(),sep = "\n")
  h2_rowi <- grep("Total Observed scale",h2_dat)
  h2_row <- h2_dat[h2_rowi]
  h2_data <- h2_dat[h2_rowi:length(h2_dat)]
  h2_data <- h2_data[grep("^[^:]+:[^:]+$",h2_data)]
  h2_datd <- transpose(strsplit(h2_data,split=":"))
  names(h2_datd) <- c("Variable","Value")
  h2_datdf <- unnest(as_data_frame(h2_datd)) %>% mutate(Variable=chartr(" ^","__",Variable),Value=trimws(Value)) %>% separate(Value,c("Est","SD"),sep = "[ s\\(\\)]+",remove=T,fill="right",convert=T)
  return(h2_datdf)
}



sim_uh_A <- function(sigu,bias,Q,D,fgeneids,usim){
  nd <- sigu^2*D^2+D+bias
  A <- Q %*% (t(Q) * sqrt(pmax(nd, 0)))
  uhmat <- t(usim%*%A)
  colnames(uhmat) <- fgeneids
  return(uhmat)
}


simuh_dir <-  function(sigu,bias,Q,D,fgeneids,seed=NULL,gen_quh=F){
  p <- nrow(Q)
  if(!is.null(seed)){
    stopifnot(is.integer(seed))
    set.seed(seed)
  }
  nreps <- length(fgeneids[[1]])
  fgeneids <- as.character(fgeneids[[1]])


  usim <- matrix(rnorm(n = p*nreps),nrow=nreps,byrow = T)
  uhmat <- sim_uh_A(sigu,bias,Q,D,fgeneids,usim)
  if(gen_quh){
    quhmat <- crossprod(Q,uhmat)
  }
  
  colnames(uhmat) <- fgeneids
  if(gen_quh){
    colnames(quhmat) <- fgeneids
  }
  retdf <- as_data_frame(uhmat) %>% mutate(snp_index=as.character(1:n())) %>% gather(fgeneid,uh,-snp_index) %>% mutate(tsigu=sigu,tbias=bias)
  if(gen_quh){  
    retdf <- as_data_frame(quhmat) %>% mutate(snp_index=as.character(1:n())) %>% gather(fgeneid,quh,-snp_index) %>% mutate(tsigu=sigu,tbias=bias) %>% inner_join(retdf,by=c("snp_index","fgeneid","tsigu","tbias"))
  }
  return(retdf)
}


estimate_RSSp_files <- function(evdf,simf,genof,chunksize,R_dataname="R",doConfound=T,doNoConfound=F,doLog=F,useGradient=T,result_dir=tempdir(),use_ldetect=F){
  library(purrr)
  library(dplyr)
  library(RcppEigenH5)
  library(tidyr)

  stopifnot((doConfound+doNoConfound)>0)
  library(LDshrink)
  if(!file.exists(evdf)){
    datal <- write_eigen_chunks(genof,evdf,chunksize,dataname=R_dataname)
  }else{
    if(!group_exists(evdf,as.character(chunksize))){
      datal <- write_eigen_chunks(genof,evdf,chunksize,dataname=R_dataname)
    }else{
      datal <- read_eigen_chunks(evdf,chunksize)
    }
  }
  base_simf <- basename(simf)

  asimd <- readRDS(simf)
  sim_md5 <- asimd[["md5"]]
  uhmat <- asimd[["bias_uh_mat"]]
  if(is.null(sim_md5)){
    asimd[["md5"]] <- PKI::PKI.digest(as.raw(asimd[["bias_uh_mat"]]))
    sim_md5 <- asimd[["md5"]]
    saveRDS(asimd,simf)
  }
  if(!dir.exists(result_dir)){
    dir.create(result_dir,recursive=T)
  }
  result_f <- file.path(result_dir,base_simf)
  # if(file.exists(result_f)){
  #   # resultl <- readRDS(result_f)
  #   # res_md5 <- resultl[["md5"]]
  #   # if(res_md5==sim_md5){
  #   #   t_res <- resultl[["result"]]
  #   # }
  # }else{
    p <- asimd$p
    n <- asimd$n
    
    prep_dfl <- prep_RSSp_evd(Ql=datal[["Ql"]],
                              Dl=datal[["Dl"]],
                              U=uhmat,N=n)
    df_l <- split(prep_dfl,prep_dfl$fgeneid) 
    cf_res <- NULL
    ncf_res <- NULL
    
    if(doConfound){
      cf_res <- bind_rows(lapply(df_l,RSSp_estimate,p_n=p/n,doConfound=T,log_params=doLog,useGradient=useGradient))
    }
    if(doNoConfound){
      ncf_res <- bind_rows(lapply(df_l,RSSp_estimate,p_n=p/n,doConfound=F,log_params=doLog,useGradient=useGradient))
    }
    f_res <- bind_rows(cf_res,ncf_res)
    stopifnot(!is.null(f_res))
    t_res <- mutate(asimd$tparam_df,fgeneid=as.character(fgeneid),chunksize=chunksize) %>% inner_join(f_res)
    resultl <- list(result=t_res,md5=sim_md5)
    saveRDS(resultl,result_f)
  # }
  return(t_res)
}
# 
# est_sim_evdf <- function(evdf,resl){
#   library(LDshrink)
#   library(RcppEigen)
#   resl$quh_mat <- sparse_matmul()
#   
# }



est_sim <- function(resl,Ql=NULL,D=NULL,doConfound=T,log_params=F,useGradient=T,a_bounds=c(0,.3),sigu_bounds=c(1e-4,1.5)){
  library(purrr)
  library(dplyr)
  stopifnot(!is.null(D))  
  if(is.null(Ql)){
    stopifnot(!is.null(resl$quh_mat))
  }
  if(is.null(resl$quh_mat)){
    resl$quh_mat <- quh_mat(Ql,resl$bias_uh_mat)
  }
  rss_res <- cross(list(doConfound=doConfound,log_params=log_params,useGradient=useGradient)) %>% 
    invoke_map_dfr("RSSp_run_mat_quh",.,quh_mat_d=resl$quh_mat,D=D,n=resl$n) %>% inner_join(resl$tparam_df)
  return(rss_res)
}

    
gen_sim_direct_evd <- function(Q,D,pve,bias=0,nreps,n,p,fgeneid=NULL,gen_quh=F){
  library(dplyr)
  library(purrr)
  library(tidyr)
  tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n,p,fgeneid) %>%  
    group_by(tpve,tbias,tsigu) %>% 
    mutate(fgeneid=as.character(fgeneid)) %>% 
    nest(fgeneid)
  asims <- tparam_df %>% rowwise() %>% do(simuh_dir(sigu = .$tsigu,bias = .$tbias,Q = Q,D=D,fgeneids = .$data,gen_quh=gen_quh))
  if(gen_quh){
    quh_mat <- select(asims,snp_index,fgeneid,quh) %>% spread(key = fgeneid,value = quh) %>% select(-snp_index) %>% data.matrix
  }
  bias_uh_mat <- select(asims,snp_index,fgeneid,uh) %>% spread(key = fgeneid,value = uh) %>% select(-snp_index) %>% data.matrix
  tparam_df <- unnest(tparam_df)
  retl <- list(bias_uh_mat=bias_uh_mat,
               tparam_df=tparam_df, n=n, p=p)
  if(gen_quh){
    retl[["quh_mat"]] <- quh_mat
  }
  return(retl)
}
    

gen_sim_direct <- function(R,pve,bias=0,nreps,n=NULL,gen_quh=F){
  stopifnot(isSymmetric(R),!is.null(n))
  library(dplyr)
  library(purrr)
  library(tidyr)
  evdR <- eigen(R)
  Q <- evdR$vectors
  D <- evdR$values
  p <- nrow(R)
  return(gen_sim_direct_evd(Q,D,pve,bias,nreps,n,p,gen_quh))
}


gen_sim_gds <- function(gds,pve,bias=0,nreps=1,seed=NULL,chunksize=10000,fgeneid=NULL,evdf=NULL,cores=1){
  library(SeqArray)
  library(LDshrink)
  library(dplyr)
  library(purrr)
  if(!is.null(seed)){
    set.seed(seed)
  }
  # gds <- seqOpen(gds_file,readonly = T)
  isHaplo <- LDshrink::is_haplo(gds)
  stopifnot(!isHaplo)
  n <- length(seqGetData(gds,"sample.id")) 
  p <- length(seqGetData(gds,"variant.id"))
  
  
  tparam_df <- gen_tparamdf_norm(pve, bias, nreps, n, p,fgeneid=fgeneid)

  sim_S <- function(index,x,sigu){
    n <- nrow(x)
    p <- ncol(x)
    u_mat <- sapply(sigu,function(tsigu,p){rnorm(n=p,mean=0,sd=tsigu)},p=p)
    S <- 1/sqrt(n)*1/apply(scale(x,center = T,scale = F),2,sd)
    beta <- u_mat*S
    return(list(u=u_mat,beta=beta,ty=scale(x,center=T,scale=F)%*%beta))
  }
  
  S_U <- transpose(seqBlockApply(gds,c(x="$dosage"),
                       sim_S,
                       margin="by.variant",
                       as.is = "list",
                       .progress = T,var.index="relative",
                       sigu=tparam_df$tsigu,bsize = chunksize,parallel=cores))

  ty <- Reduce("+",S_U$ty)
  vy <- apply(ty, 2, var)
  residvec <- gen_ti(vy, tparam_df$tpve)
  residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
  ymat <- scale(ty+residmat, center=T, scale=F)
  beta_se <- transpose(seqBlockApply(gds,c(x="$dosage"),
                       function(x,ymat){
                         sx <- scale(x,center=T,scale=F)
                         betahat <- RSSReQTL::map_beta_exp(sx, ymat)
                         se_mat <- RSSReQTL::map_se_exp(sx, ymat, betahat)
                         return(list(betahat=betahat,se_mat=se_mat))
                       },ymat=ymat,
                       margin="by.variant",
                       as.is = "list",
                       .progress = T,bsize = chunksize,parallel=cores))
  betahat_mat <- do.call("rbind",beta_se$betahat)
  se_mat <- do.call("rbind",beta_se$se_mat)
  colnames(betahat_mat) <- as.character(tparam_df$fgeneid)

  colnames(se_mat) <- as.character(tparam_df$fgeneid)
  
  bias_mat <- sapply(tparam_df$tbias, function(ta, p){rnorm(n=p, mean=0, sd=sqrt(ta))}, p=p)
  colnames(bias_mat) <- as.character(tparam_df$fgeneid)
  uh_mat <- betahat_mat/se_mat
  colnames(uh_mat) <- as.character(tparam_df$fgeneid)
  #  u_mat <- betamat/se_mat
  bias_uh_mat <- uh_mat+bias_mat
  colnames(bias_uh_mat) <- as.character(tparam_df$fgeneid)
  retl <- list(bias_uh_mat=bias_uh_mat,
               tparam_df=tparam_df, n=n, p=p)
  if(!is.null(evdf)){
    ql_df <- LDshrink::read_si_ql(evdf)
    # quhm <- gen_quh(ql_df$Ql,uhmat <- bias_
    retl[["quh_mat"]] <- quh_mat(ql_df$Ql,uhmat = bias_uh_mat)
  }
  return(retl)
}



read_SNPinfo_ldsc_gwas <- function(gds,zmat,N=NULL){
  library(tidyr)
  if(is.null(N)){
    N <- length(seqGetData(gds,"sample.id"))
    if(LDshrink::is_haplo(gds)){
      N <- N/2
    }
  }
  tdf <- tibble::data_frame(SNP=seqGetData(gds,var.name="annotation/id"),
                            allele=seqGetData(gds,var.name="allele")) %>% 
    separate(allele,into = c("A1","A2"),sep = ",") %>% 
    mutate(N=N,snp_id=1:n())
  tdf<- as_data_frame(zmat) %>% mutate(snp_id=1:n()) %>% 
    gather(fgeneid,Z,-snp_id) %>% 
    inner_join(tdf) %>% 
    select(SNP,A1,A2,Z,N,fgeneid)
  return(split(select(tdf,-fgeneid),tdf$fgeneid))  
}





gen_sim_gds_direct_ldsc <- function(Ql,Dl,gds,pve,bias=0,nreps=1,seed=NULL,fgeneid=NULL,gen_quh=F){
  library(LDshrink)
  library(readr)
  library(purrr)
  library(SeqArray)
  if(!is.null(seed)){
    set.seed(seed)
  }
  p <- sum(lengths(Dl))
  stopifnot(sum(sapply(Ql,nrow))==p)
  isHaplo <- is_haplo(gds)
  
  
  n <- length(seqGetData(gds,"sample.id"))
  if(isHaplo){
    n <- n/2
  }
  sim_l<- transpose(map2(Ql, Dl,gen_sim_direct_evd,
                         pve = pve, bias = bias, nreps = nreps,
                         n = n, p = p, fgeneid = fgeneid))
  resl <- list()
  resl[["bias_uhmat"]] <- do.call("rbind",sim_l$bias_uh_mat)
  resl[["p"]] <- unique(unlist(sim_l$p))
  resl[["n"]] <- unique(unlist(sim_l$n))
  resl[["ldsc_df_l"]] <- read_SNPinfo_ldsc_gwas(gds,resl$bias_uhmat,N=resl$n)
  return(resl)
}
  
  
gen_sim_gds_ldsc <- function(gds,pve,bias=0,nreps=1,seed=NULL,chunksize=10000,fgeneid=NULL,evdf=NULL){
  library(LDshrink)
  library(readr)
  
  sim_l <- gen_sim_gds(gds,pve=pve,bias=bias,nreps=nreps,seed=seed,chunksize=chunksize,fgeneid=fgeneid)
  
  sim_l[["ldsc_df_l"]] <- read_SNPinfo_ldsc_gwas(gds,sim_l$bias_uh_mat,N=sim_l$n)
  return(sim_l)
}




gen_sim_list <- function(SNP, pve, bias=0, nreps,seed=NULL){
  library(dplyr)
  if(!is.null(seed)){
    set.seed(seed)
  }
  n <- nrow(SNP)
  p <- ncol(SNP)

  stopifnot(all.equal(colMeans(SNP),rep(0,p),tolerance=1e-5))
  tparam_df <- gen_tparamdf_norm(pve, bias, nreps, n, p) %>% mutate(fgeneid=as.character(fgeneid))
  
  u_mat <- sapply(tparam_df$tsigu, function(sigu, p){rnorm(n = p, mean = 0, sd = sigu)}, p=p)
  
  S <- 1/sqrt(n)*(1/apply(SNP, 2, sd))
  betamat <- u_mat*S
  vy <- apply(SNP%*%betamat, 2, var)
  residvec <- gen_ti(vy, tparam_df$tpve)
  residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
  ymat <- scale(SNP%*%betamat+residmat, center=T, scale=F)
  betahat_mat <- RSSReQTL::map_beta_exp(SNP, ymat)
  colnames(betahat_mat) <- as.character(tparam_df$fgeneid)
  se_mat <- RSSReQTL::map_se_exp(SNP, ymat, betahat_mat)
  colnames(se_mat) <- as.character(tparam_df$fgeneid)
  
  bias_mat <- sapply(tparam_df$tbias, function(ta, p){rnorm(n=p, mean=0, sd=sqrt(ta))}, p=p)
  colnames(bias_mat) <- as.character(tparam_df$fgeneid)
  uh_mat <- betahat_mat/se_mat
  colnames(uh_mat) <- as.character(tparam_df$fgeneid)
  #  u_mat <- betamat/se_mat
  bias_uh_mat <- uh_mat+bias_mat
  colnames(bias_uh_mat) <- as.character(tparam_df$fgeneid)
  


  retl <- list(bias_uh_mat=bias_uh_mat,
               tparam_df=tparam_df, n=n, p=p)
  
  return(retl)
}

gen_sim <- function(SNP, pve, bias=0, nreps, outdir, collapse_all=T, overwrite=F,debug=F){
  library(dplyr)
  n <- nrow(SNP)
  p <- ncol(SNP)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  stopifnot(all.equal(colMeans(SNP),rep(0,p),tolerance=1e-5))
  tparam_df <- gen_tparamdf_norm(pve, bias, nreps, n, p)
  
  u_mat <- sapply(tparam_df$tsigu, function(sigu, p){rnorm(n = p, mean = 0, sd = sigu)}, p=p)
  
  S <- 1/sqrt(n)*(1/apply(SNP, 2, sd))
  betamat <- u_mat*S
  vy <- apply(SNP%*%betamat, 2, var)
  residvec <- gen_ti(vy, tparam_df$tpve)
  residmat <- sapply(residvec, function(ti, n){rnorm(n=n, mean=0, sd=sqrt(ti))}, n=n)
  ymat <- scale(SNP%*%betamat+residmat, center=T, scale=F)
  betahat_mat <- RSSReQTL::map_beta_exp(SNP, ymat)
  colnames(betahat_mat) <- as.character(tparam_df$fgeneid)
  se_mat <- RSSReQTL::map_se_exp(SNP, ymat, betahat_mat)
  colnames(se_mat) <- as.character(tparam_df$fgeneid)
  
  bias_mat <- sapply(tparam_df$tbias, function(ta, p){rnorm(n=p, mean=0, sd=sqrt(ta))}, p=p)
  colnames(bias_mat) <- as.character(tparam_df$fgeneid)
  uh_mat <- betahat_mat/se_mat
  colnames(uh_mat) <- as.character(tparam_df$fgeneid)
  #  u_mat <- betamat/se_mat
  bias_uh_mat <- uh_mat+bias_mat
  colnames(bias_uh_mat) <- as.character(tparam_df$fgeneid)
  
  if(!collapse_all){
    if(debug){
      saveRDS(betamat, file.path(outdir, "betamat.RDS"))
      saveRDS(residvec, file.path(outdir, "residvec.RDS"))
      saveRDS(residmat, file.path(outdir, "residmat.RDS"))
      saveRDS(ymat, file.path(outdir, "ymat.RDS"))
      saveRDS(betahat, file.path(outdir, "betahat_mat.RDS"))
      saveRDS(se_mat, file.path(outdir, "se_mat.RDS"))
      saveRDS(bias_mat, file.path(outdir, "bias_mat.RDS"))
      saveRDS(u_mat, file.path(outdir, "u_mat.RDS"))
      saveRDS(uh_mat, file.path(outdir, "uh_mat.RDS"))
    }
    saveRDS(bias_uh_mat, file.path(outdir, "bias_uh_mat.RDS"))
    saveRDS(tparam_df, file.path(outdir, "tparam_df.RDS"))
  }else{
    outfilep <- file.path(outdir, "simulation.RDS")
    i <- 0
    if(!overwrite){
      while(file.exists(outfilep)){
        i <- i+1
        outfilep <- file.path(outdir, paste0("simulation", i, ".RDS"))
      }
    }
    if(debug){
      retl <- list(betamat=betamat,
                   residsd=residvec,
                   residmat=residmat,
                   ymat=ymat,
                   betahat=betahat_mat,
                   se_mat=se_mat,
                   bias_mat=bias_mat,
                   u_mat=u_mat,
                   uh_mat=uh_mat,
                   bias_uh_mat=bias_uh_mat,
                   tparam_df=tparam_df, n=n, p=p)
      retl[["md5"]] <- PKI::PKI.digest(as.raw(retl[["bias_uh_mat"]]))
      saveRDS(retl, file = outfilep)
    }else{
      retl <- list(bias_uh_mat=bias_uh_mat,
                   tparam_df=tparam_df, n=n, p=p)
      retl[["md5"]] <- PKI::PKI.digest(as.raw(retl[["bias_uh_mat"]]))
      saveRDS(retl, file = outfilep)
    }
    return(outfilep)
  }
}

# gen_sim_ldsc <- function(df,tpve=0.56,){
#   
#   
#   
# }
