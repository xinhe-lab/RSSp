context("Checking Simulated Data")



test_that("We can estimate parameters",{
  
    data("simulation_tparam_df")
    data("simulation_R")
    data("simulation_uh_mat")
    
                                        #Use LDshrink to calculate R

                                        #use base R to perform EVD 
    evdR <- eigen(simulation_R)
                                        #est_sim is a helper function for estimating parameters from simulations. 
                                        # instead of taking a single matrix of eigenvectors, it takes a list, with each element representing each LD block

    D <- evdR$values
    n <- unique(simulation_tparam_df$n)
    p <- unique(simulation_tparam_df$p)
    
    quh_mat <- crossprod(evdR$vectors,simulation_uh_mat)
    colnames(quh_mat) <- colnames(simulation_uh_mat)
    
    res <- RSSp_run_mat_quh(quh_mat = quh_mat, D = D,n = n,doConfound = T) %>% dplyr::inner_join(simulation_tparam_df)

    
    expect_true(all(res$pve<=1))
    expect_true(all(res$bias>=0))
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




test_that("Can use Xiang's method of computing the RSSp likelihood",{
  
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
  D <- Revd$values
  x_dens <- function(par,dvec,quh){
    sigu <- par[1]
    bias <- par[2]
    return(sum(dnorm(quh,mean = 0,sd = sqrt((sigu^2*dvec^2+dvec+bias)),log = T)))
  }
  x_res <- x_dens(par = c(sigu,bias),
                  dvec=D,quh = qrdat)
  cpp_dens <- evd_dnorm(c(sigu,bias),dvec=D,quh = qrdat)-0.5*p*log(2*pi)
  
  
  expect_equal(cpp_dens,x_res)
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




test_that("density performed with precision is equal to density using sd",{
  
  tparam <- runif(2)
  tinv <- tparam
  tinv[1] <- 1/tinv[1]
  dvec <- runif(5)
  quh <- rnorm(5)
  expect_equal(evd_dnorm_prec(tinv,dvec,quh),evd_dnorm(tparam,dvec,quh))
  
})

test_that("variance is nonzero when estimating parameters ",{
  
  
  tparam <- c(.Machine$double.eps,0)
  dvec <- runif(5)
  quh <- rnorm(5)
  nRh <- -1/RSSp_hess(par = tparam,dvec = dvec,quh = quh)
  expect_equal(evd_dnorm_grad_stan(par=tparam,dvec = dvec,quh = quh),
               evd_dnorm_grad(par = tparam,dvec = dvec,quh = quh))
  
  
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








