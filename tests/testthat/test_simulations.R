context("Checking Simulated Data")








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
  cpp_dens <- evd_dnorm(c(sigu^2),dvec=Revd$values,quh = qrdat)
  expect_equal(cpp_dens,R_dens)
})


test_that("Can calculate likelihood at a grid of points",{
  
  p <- 10
  n <- 1000
  v_hat <- rnorm(p)
  eigenvals <- runif(p)
  sigu_grid <- seq(0.01,0.1,length.out = 3)
  cpp_dvec <- -(sigu_grid^2 %>% purrr::map_dbl(evd_dnorm,dvec=eigenvals,quh=v_hat))
  dens_res <- RSSp_estimate_grid(v_hat = v_hat,eigenvalues = eigenvals,p_n = p/n,grid_points = 10,sigu_grid = sigu_grid)
  expect_equal(dens_res$lnZ,cpp_dvec)
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
  cpp_dens <- evd_dnorm(c(sigu^2,bias),dvec=Revd$values,quh = qrdat)
  expect_equal(cpp_dens,R_dens)
})





test_that("My gradient is the same as the stan gradient",{
  
  
  tparam <- runif(2)
  # ttparam <- tparam
  # ttparam[2] <- tparam[2]*tparam[1]
  dvec <- runif(5)
  quh <- rnorm(5)
  expect_equal(evd_dnorm_grad_stan(par=tparam,dvec = dvec,quh = quh),
               evd_dnorm_grad(par = tparam,dvec = dvec,quh = quh))
})


test_that("My hessian is the same as the stan gradient",{
  
  tparam <- runif(2)
  # ttparam <- tparam
  # ttparam[2] <- prod(tparam)
  dvec <- runif(5)
  quh <- rnorm(5)
  stan_hess <- evd_dnorm_hess_stan(par=tparam,dvec=dvec,quh=quh)
  nRh <- RSSp_hess(par = tparam,dvec = dvec,quh = quh)
  expect_equal(stan_hess,nRh)
})










