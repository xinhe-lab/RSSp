context("Checking Simulations")


rssp_dir <- system.file("gds",package="RSSp")
haplogdsf <- file.path(rssp_dir,"sub_19_haplo.gds")

genogdsf <- file.path(rssp_dir,"sub_19_geno.gds")




test_that("We can estimate parameters",{
  

  
  genogds <- SeqArray::seqOpen(genogdsf)
  haplogds <- SeqArray::seqOpen(haplogdsf)
  
  
  # data("shrink_R")
  # data("cohort_SNP")
  library(LDshrink)
  pve <- as.numeric(seq(0.7,0.9,length.out = 10))
  bias <- as.numeric(seq(0.0,0.1,length.out = 2))
  nreps <- 3
  
  p <- calc_p(genogds)
  n <- calc_N(genogds)
  chunksize <- p/1
  # SNP <- cohort_SNP
  # temp_dir <- tempdir()
  msim <- gen_sim_gds(gds,pve=pve,bias=bias,nreps=nreps)
  R <- LDshrink::calc_LD_gds(haplogds)
  evdR <- eigen(R)
  Ql <- list("1"=evdR$vectors)
  D <- evdR$values
  
 res <- est_sim(resl = msim,Ql = Ql,D = D,doConfound = T) %>% select(-log_params,-useGradient)
 expect_true(all(res$pve<=1))
 expect_true(all(res$a_hat>=0))
 expect_equal(res$sigu,res$tsigu,tolerance=0.1)
 
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
  Dl <- rep(1,p)
  tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n,p)
  diag_sim <- sim_quh_dir(tsigu = tparam_df$tsigu[1],tbias = tparam_df$tbias[1],Q = evdR$vectors,D = evdR$values,fgeneid = tparam_df$fgeneid,snp_id = 1:p)
  rss_res <- RSSp_run_mat_quh(quh_mat_d = diag_sim$quh,D=diag_sim$D,n = n) %>% inner_join(tparam_df)
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
  tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n,p)
  diag_sim <- sim_quh_dir(tsigu = tparam_df$tsigu[1],tbias = tparam_df$tbias[1],Q = evdR$vectors,D = evdR$values,fgeneid = tparam_df$fgeneid,snp_id = 1:p)
  
  rss_res <- RSSp_run_mat_quh(quh_mat_d = diag_sim$quh,D=diag_sim$D,n = n) %>% inner_join(tparam_df)
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
  tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n,p)
  diag_sim <- sim_quh_dir(tsigu = tparam_df$tsigu[1],tbias = tparam_df$tbias[1],Q = evdR$vectors,D = evdR$values,fgeneid = tparam_df$fgeneid,snp_id = 1:p)
  
  rss_res <- RSSp_run_mat_quh(quh_mat_d = diag_sim$quh,D=diag_sim$D,n = n) %>% inner_join(tparam_df)
  expect_equal(rss_res$sigu,rss_res$tsigu,tolerance=0.2)
})





test_that("My gradient is the same as the stan gradient",{
  
  tparam <- runif(2)
  dvec <- runif(5)
  quh <- rnorm(5)
  expect_equal(evd_dnorm_grad_stan(par=tparam,dvec = dvec,quh = quh)[-1],evd_dnorm_grad(par = tparam,dvec = dvec,quh = quh))
  
  
  
})



