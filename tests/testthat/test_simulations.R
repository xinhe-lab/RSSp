context("Checking Simulations")

test_that("We can simulate data",{
  data("cohort_SNP")
  
  
  pve <- as.numeric(seq(0.01,0.9,length.out = 2))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 2
  SNP <- cohort_SNP
  p <- ncol(cohort_SNP)
  temp_dir <- tempdir()
  file.remove(dir(temp_dir,full.names = T,pattern="RDS"))
  sim_file <- gen_sim(SNP,pve,bias,nreps,temp_dir,overwrite = T)
  
  expect_true(file.exists(sim_file))
  expect_equal(file.path(temp_dir,"simulation.RDS"),sim_file)
  expect_equal(length(dir(temp_dir,pattern="RDS")),1)
  file.remove(sim_file)
})






test_that("We can simulate data directly from LD",{
  data("cohort_SNP")
  data("shrink_R")
  
  pve <- as.numeric(seq(0.01,0.9,length.out = 2))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 2
  n <- nrow(cohort_SNP)
  p <- ncol(cohort_SNP)
  temp_dir <- tempdir()
  file.remove(dir(temp_dir,full.names = T,pattern="RDS"))
  sim_file <- gen_sim_direct(shrink_R,pve,bias,nreps,temp_dir,n = n)
  
  expect_true(file.exists(sim_file))
  expect_equal(file.path(temp_dir,"simulation.RDS"),sim_file)
  file.remove(sim_file)
  # expect_equal(length(dir(temp_dir,pattern="RDS")),1)
})




test_that("We can estimate parameters",{
  data("shrink_R")
  data("cohort_SNP")
  library(LDshrink)
  pve <- as.numeric(seq(0.7,0.9,length.out = 10))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 3
  p <- ncol(cohort_SNP)
  n <- nrow(cohort_SNP)
  chunksize <- p/1
  SNP <- cohort_SNP
  temp_dir <- tempdir()
  sim_file <- gen_sim(cohort_SNP,pve,bias,nreps,temp_dir)
  asimd <- readRDS(sim_file)
  Rl <- write_eigen_chunks(chunksize=chunksize,R = shrink_R,write=F)
  names(Rl[["Ql"]]) <- 1:length(Rl$Ql)
  names(Rl$rdl) <- 1:length(Rl$rdl)
  names(Rl$indl) <- 1:length(Rl$indl)
  quh <- gen_quh(Rl$Ql,uhmat = asimd$bias_uh_mat,indl = Rl$indl,Dl = Rl$rdl) %>% unnest()

  qhl <-  split(quh$quh,quh$fgeneid)
  rdl <- split(quh$rd,quh$fgeneid)
  fgeneid <- names(quh$rd)
  t_res_df <- pmap_dfr(.l = list(fgeneid=names(qhl),
                        D=rdl,
                        quh=qhl),RSSp,p_n=p/n)
  tval <- inner_join(t_res_df,mutate(asimd$tparam_df,fgeneid=as.character(fgeneid)))

})


# 
test_that("Both methods of simulating MVN work",{
  data("shrink_R")
  data("cohort_SNP")
  library(tidyverse)
  pve <- as.numeric(seq(0.01,0.9,length.out = 8))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 4
  n <- nrow(cohort_SNP)

  test_simuh_dir <- function(sigu,bias,R,Rsq,fgeneids,seed=NULL){
    p <- nrow(R)
    if(!is.null(seed)){
      stopifnot(is.integer(seed))
      set.seed(seed)
    }
    nreps <- length(fgeneids[[1]])
    fgeneids <- as.character(fgeneids[[1]])
    mcov <- sigu^2*Rsq+R+bias*diag(p)
    uhmat <- t(mvtnorm::rmvnorm(n = nreps,mean = rep(0,p),sigma = mcov))
    colnames(uhmat) <- fgeneids
    retdf <- as_data_frame(uhmat) %>% mutate(snp_index=1:n()) %>% gather(fgeneid,uh,-snp_index) %>% mutate(tsigu=sigu,tbias=bias)
    return(retdf)
  }



  p <- ncol(cohort_SNP)
  tparam_df <- gen_tparamdf_norm(pve, bias, nreps=nreps, n, p) %>%
    group_by(tpve,tbias,tsigu) %>% mutate(fgeneid=as.character(fgeneid)) %>% nest(fgeneid) %>% slice(1)

  eigenR <- eigen(shrink_R)
  Rsq <- shrink_R%*%shrink_R
  asims <- tparam_df %>% rowwise() %>% do(test_simuh_dir(sigu = .$tsigu,bias = .$tbias,R = shrink_R,Rsq=Rsq,fgeneids = .$data,seed=123L)) %>%
    ungroup() %>% mutate(sim="Old")
  aasims <- tparam_df %>% rowwise() %>% do(test_simuh_dir(sigu = .$tsigu,bias = .$tbias,R = shrink_R,Rsq=Rsq,fgeneids = .$data,seed=123L)) %>%
    ungroup() %>% mutate(sim="Old")
  expect_identical(asims,aasims)

  Q <- eigenR$vectors
  D <- eigenR$values
  bsims <- tparam_df %>% rowwise() %>% do(simuh_dir(sigu = .$tsigu,bias = .$tbias,Q = Q,D=D,fgeneids = .$data,seed=123L)) %>%
    ungroup() %>% mutate(sim="New")
  bbsims <- tparam_df %>% rowwise() %>% do(simuh_dir(sigu = .$tsigu,bias = .$tbias,Q = Q,D=D,fgeneids = .$data,seed=123L)) %>%
    ungroup() %>% mutate(sim="New")
  expect_identical(bsims,bbsims)
  expect_equal(asims$uh,bsims$uh)
  })



test_that("Multivariate density is straightforward to compute when the covariance is diagonal and there's no bias",{
  
  data("shrink_R")
  p <- nrow(shrink_R)
  R <- diag(p)
  Revd <- eigen(R)
  
  sigu <- 0.5
  bias <- 0.0
  Rsq <- R%*%R
  tcov <- sigu^2*Rsq+R+bias*diag(p)
  rdat <- c(mvtnorm::rmvnorm(n=1,mean=rep(0,p),sigma=tcov))
  qrdat <- c(crossprod(Revd$vectors,rdat))
  mrnorm <- sum(dnorm(rdat,mean=0,sd=sqrt(diag(tcov)),log=T))
  R_dens <- mvtnorm::dmvnorm(x = rdat,mean = rep(0,p),sigma = tcov,log = T)
  expect_equal(mrnorm,R_dens)
  cpp_dens <- evd_dnorm(c(sigu,bias),dvec=Revd$values,quh = qrdat)-0.5*p*log(2*pi)
  expect_equal(cpp_dens,R_dens)
  expect_equal(optimise(RSSp_evd_mvd,interval = c(0,1),dvec=Revd$values,quh=qrdat)$minimum,sigu,tolerance=0.2)
})


test_that("Multivariate density is computed correctly without bias",{
  
  data("shrink_R")
  R <- shrink_R
  Revd <- eigen(R)
  p <- nrow(R)
  sigu <- 0.5
  bias <- 0.0
  Rsq <- R%*%R
  tcov <- sigu^2*Rsq+R+bias*diag(p)
  rdat <- c(mvtnorm::rmvnorm(n=1,mean=rep(0,p),sigma=tcov))
  qrdat <- c(crossprod(Revd$vectors,rdat))
  R_dens <- mvtnorm::dmvnorm(x = rdat,mean = rep(0,p),sigma = tcov,log = T)
  cpp_dens <- evd_dnorm(c(sigu,bias),dvec=Revd$values,quh = qrdat)-0.5*p*log(2*pi)
  expect_equal(cpp_dens,R_dens)
})


test_that("Multivariate density is computed correctly with bias",{
  
  data("shrink_R")
  R <- shrink_R
  Revd <- eigen(R)
  p <- nrow(R)
  sigu <- 0.5
  bias <- 0.1
  Rsq <- R%*%R
  tcov <- sigu^2*Rsq+R+bias*diag(p)
  rdat <- c(mvtnorm::rmvnorm(n=1,mean=rep(0,p),sigma=tcov))
  qrdat <- c(crossprod(Revd$vectors,rdat))
  R_dens <- mvtnorm::dmvnorm(x = rdat,mean = rep(0,p),sigma = tcov,log = T)
  cpp_dens <- evd_dnorm(c(sigu,bias),dvec=Revd$values,quh = qrdat)-0.5*p*log(2*pi)
  expect_equal(cpp_dens,R_dens)
})



test_that("We can estimate parameters with direct simulations really well when using the identity as an LD matrix( and not estimating confounding",{
  library(dplyr)
  data("shrink_R")
  data("cohort_SNP")
  library(LDshrink)
  data("map_dat")
  n <- 1000
  p <- nrow(shrink_R)
  R <- diag(p)
  #generate PSD LD matrix
  pve <- 0.5
  bias <- 0
  nreps <- 5
  diag_sim <- gen_sim_direct_evd(diag(p),D = rep(1,p),pve = pve,bias = bias,nreps = nreps,n=n,p = p,gen_quh = T)
  rss_res <-  imap_dfr(as_data_frame(diag_sim$quh_mat),function(quh,fgeneid,D,p_n){
    RSSp(fgeneid = fgeneid,D = D,quh = quh,p_n =p_n ,doConfound = F)
  },p_n=p/n,D=Dl) %>% inner_join(diag_sim$tparam_df)
  expect_equal(rss_res$sigu,rss_res$tsigu,tolerance=0.2)
})

test_that("We can estimate parameters with direct simulations really well without using the identity",{
  library(dplyr)
  data("shrink_R")
  data("cohort_SNP")
  library(LDshrink)
  data("map_dat")
  n <- 1000
  p <- nrow(shrink_R)
  R <- shrink_R
  #generate PSD LD matrix
  pve <- 0.5
  bias <- 0
  nreps <- 5
  evdR <- eigen(R)
  diag_sim <- gen_sim_direct_evd(Q = evdR$vectors  ,D = evdR$values,pve = pve,bias = bias,nreps = nreps,n=n,p = p,gen_quh = T)
  
  # t_sigu <- diag_sim$tparam_df$tsigu[1]
  # sigu_range <- seq(t_sigu*0.5,t_sigu*1.5,length.out = 1000)
  # bias_range <- 0
  # lnzvec <- numeric(length(sigu_range))
  # for(i in 1:length(lnzvec)){
  #   lnzvec[i] <- RSSp_evd_mvd(par = c(sigu_range[i],bias),dvec = D,quh = quh)
  # }
  # plot(sigu_range,lnzvec)
  # outer(X = pve_range,Y = bias_range,FUN = function(pve,bias,D,quh){
  #   },D=evdR$values,quh=diag_sim$bias_uh_mat[,1])
  
  rss_res <-  imap_dfr(as_data_frame(diag_sim$quh_mat),function(quh,fgeneid,D,p_n){
    RSSp(fgeneid = fgeneid,D = D,quh = quh,p_n =p_n ,doConfound = F)
  },p_n=p/n,D=evdR$values) %>% inner_join(diag_sim$tparam_df)
  expect_equal(rss_res$sigu,rss_res$tsigu,tolerance=0.2)
})





test_that("My gradient is the same as the stan gradient",{
  
  tparam <- runif(2)
  dvec <- runif(5)
  quh <- rnorm(5)
  expect_equal(evd_dnorm_grad_stan(par=tparam,dvec = dvec,quh = quh)[-1],evd_dnorm_grad(par = tparam,dvec = dvec,quh = quh))
  
  
  
})

test_that("My cpp mvd transformation is the sam as the R version",{
  
  
  data("shrink_R")
  data("cohort_SNP")
  library(LDshrink)
  t_evd <- eigen(shrink_R)
  Q <- t_evd$vectors
  p <- nrow(Q)
  n <- 4
  D <- t_evd$values
  sigu <- 0.03
  bias <- 0.001
  nreps <- 2
  usim <- t(matrix(rnorm(n = p*nreps),nrow=nreps,byrow = T))
  
    nd <- sigu^2*D^2+D+bias
    A <- (t(Q) * sqrt(pmax(nd, 0)))
    uhmat <- A%*%usim
    colnames(uhmat) <- fgeneids
  #   return(uhmat)
  # }

  tsA =sim_uh_A(sigu,bias,Q,D,NULL,usim)
  cppA = simuh_dir_cpp(sigu,bias,nreps,Q,D,"",usim)
  expect_equal(cppA,tsA)
    
})

# test_that("Gradient works as expected (similar to numeric gradient)",{
#   
#   data("shrink_R")
#   data("cohort_SNP")
#   
#   p <- 3
#   dvec <- runif(p)
#   quh <- runif(p)
#   tsigu <- runif(1)
#   ta <- 0.01
#   par0 <- c(tsigu,ta)
#   evd_dnorm_grad_stan(par0,dvec,quh)
#   
#   ## First reconstruct the density function
#   evd_dnorm_R <- function(par,dvec,quh){
#     sigu <- par[1];
#     varu=sigu*sigu;
#     a <- par[2]
#     tsum <-  sum(log(dvec*dvec*varu+dvec+a))
#     tprod <- sum(quh^2/(dvec^2*varu+dvec+a))
# #    double tprod = ((quh*(1/(dvec*dvec*varu+dvec+a)))*(quh)).sum();
#     return( -0.5*(tsum+tprod))
#   }
#   
#   #My first shot at the gradient function
#   evd_dnorm_grad_R <- function(par,dvec,quh){
#     sigu <- par[1]
#     varu=sigu*sigu
#     a <- par[2]
#     grad_tsum <-evd_dnorm_grad_tsum(par,dvec,quh)
#     grad_tprod <- evd_dnorm_grad_tprod(par,dvec,quh)
#     sgrad <- -sum((0.5*(a+dvec^2*varu+dvec-quh^2))/(a+dvec^2*varu+dvec)^2)
#     return(c(grad_tsum+grad_tprod,sgrad))
#   
#   }
#   evd_dnorm_grad_a <- function(par,dvec,quh){
#     return(evd_dnorm_grad_stan(par,dvec,quh)[-1])
#   }
#   
#   expect_equal(evd_dnorm(par0,dvec,quh),evd_dnorm_grad_stan(par0,dvec,quh)[1])
#   expect_equal(evd_dnorm_R(par0,dvec,quh),evd_dnorm(par0,dvec,quh))
#   expect_equal(evd_dnorm_grad(par,dvec,quh),evd_dnorm_grad_stan(par,dvec,quh)[-1])
#   #First term in the sigu density
#   evd_dnorm_tsum <- function(par,dvec,quh){
#     sigu <- par[1];
#     varu=sigu*sigu;
#     a <- par[2]
#     return(-0.5*sum(log(dvec*dvec*varu+dvec+a)))
#   }
#   #Second term in the sigu density
#   evd_dnorm_tprod <- function(par,dvec,quh){
#     sigu <- par[1];
#     varu=sigu*sigu;
#     a <- par[2]
#     return(-0.5*sum(quh^2/(dvec^2*varu+dvec+a)))
#   }
#   #Gradient of second term in sigu density
#   evd_dnorm_grad_tprod <- function(par,dvec,quh){
#     sigu <- par[1];
#     varu=sigu*sigu;
#     a <- par[2]
#     return(sum((dvec^2*sigu*quh^2)/(a+dvec^2*varu+dvec)^2))
#   }
#   #Gradient of first term in sigu density
#   evd_dnorm_grad_tsum<- function(par,dvec,quh){
#     sigu <- par[1];
#     varu=sigu*sigu;
#     a <- par[2]
#     return(sum(-(dvec^2*sigu)/(a+dvec^2*varu+dvec)))
#   }
# 
#   expect_equal(evd_dnorm_R(par0,dvec,quh),(evd_dnorm_tsum(par0,dvec,quh)+evd_dnorm_tprod(par0,dvec,quh)))
#   expect_equal(evd_dnorm_grad_R(par0,dvec,quh)[1],(evd_dnorm_grad_tsum(par0,dvec,quh)+evd_dnorm_grad_tprod(par0,dvec,quh)))
#   expect_equal(evd_dnorm_grad_R(par0,dvec,quh),evd_dnorm_grad_a(par,dvec,quh))
# 
# })



test_that("We can use the gradient to maximize",{
  data("shrink_R")
  data("cohort_SNP")
  
  pve <- as.numeric(seq(0.01,0.9,length.out = 8))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 4
  n <- nrow(cohort_SNP)
  p <- ncol(cohort_SNP)
  temp_dir <- tempdir()
  sim_file <- gen_sim_direct(shrink_R,pve,bias,nreps,temp_dir,n = n)
  # sim_file <- gen_sim(SNP,pve,bias,nreps,temp_dir)
  chunksize <- p
  asimd <- readRDS(sim_file)
  Rl <- LDshrink::write_eigen_chunks(chunksize=chunksize,R = shrink_R,write=F)
  m_est_df <- gen_rss_est_df(Rl$Ql,Rl$rdl,Rl$indl,asimd$bias_uh_mat)
  splitl <- split(m_est_df,m_est_df$fgeneid) %>% map(function(x)split(x,x$chunki))
  m_res_grad <- bind_rows(lapply(splitl,RSSp_estimate,p_n=p/n))
  t_res_grad <- inner_join(asimd$tparam_df,m_res_grad) %>% mutate(useGrad=T)
  
  m_res_ngrad <- bind_rows(lapply(splitl,RSSp_estimate,p_n=p/n,useGrad=F))
  t_res_ngrad <- inner_join(asimd$tparam_df,m_res_ngrad)%>% mutate(useGrad=F)
  expect_lte(sum(abs(t_res_grad$pve-t_res_grad$tpve)),sum(abs(t_res_ngrad$pve-t_res_ngrad$tpve)))
  
})




# test_that("Fit improves with increasing the number of starts",{
#   data("shrink_R")
#   data("cohort_SNP")
#   pve <- as.numeric(seq(0.01,0.9,length.out = 10))
#   bias <- as.numeric(seq(0.0,0.1,length.out = 2))
#   nreps <- 3
#   SNP <- cohort_SNP
#   p <- ncol(SNP)
#   n <- nrow(SNP)
#   chunksize <- p/1
#   
#   temp_dir <- tempdir()
#   sim_file <- gen_sim(SNP,pve,bias,nreps,temp_dir)
#   asimd <- readRDS(sim_file)
#   Rl <- LDshrink::write_eigen_chunks(chunksize=chunksize,R = shrink_R,write=F)
#   m_est_df <- gen_rss_est_df(Rl$Ql,Rl$rdl,Rl$indl,asimd$bias_uh_mat)
#   
#   # p <- asimd$p
#   # n <- nrow(asimd$ymat)
# 
#   df_l <- split(m_est_df,m_est_df$fgeneid) %>% map(function(x)split(x,x$chunki))
#   niter <- 1:10
#   df_liter <- cross2(df_l,niter) %>% map(function(x){RSSp_estimate(x[[1]],niter=x[[2]],p_n =p/n)})
#   f_res <- bind_rows(df_liter)
# 
#   tparamdf <- asimd$tparam_df %>% mutate(fgeneid=as.character(fgeneid))
# 
#   ffres <- inner_join(f_res,tparamdf) %>% mutate(rmse=abs(pve-tpve))
#   ggplot(ffres,aes(x=niter,y=lnZ))+geom_point()
#   f_res <- bind_rows(apply(df_l,RSSp_estimate,p_n=p/n))
#   t_res <- mutate(asimd$tparam_df,fgeneid=as.character(fgeneid),chunksize=chunksize) %>% inner_join(f_res)
#   
#   
#   expect_equal(t_res$tpve,t_res$pve,tolerance=0.2)
# })



