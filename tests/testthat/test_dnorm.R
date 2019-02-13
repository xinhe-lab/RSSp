context("evd_dnorm")

test_that("fancy templated likelihood works like toy R version with confounding",{
  
  quh <- rnorm(1e2)
  D <- runif(1e2)/runif(1e2)  
  par <- c(runif(4))
  gen_D <- function(cvec,D){
    rD <- D
    for(i in 1:length(cvec)){
      rD <- rD+cvec[i]*D^(2-(i-1))
    }
    return(rD)
  }
  mpar <- par
  
  test_D <- gen_D(mpar,D)
  oD <- D+mpar[1]*D*D+mpar[2]*D^1+mpar[3]+mpar[4]/D
  expect_equal(test_D,oD)
  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar,D)),log=T)),
               evd_dnorm(mpar,D,quh)
  )
  expect_equal(  
      -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar,D)),log=T)),
      evd_dnorm(mpar,D,quh),
  )
  expect_equal(
      -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1:3],D)),log=T)),
      evd_dnorm(mpar[1:3],D,quh)
  )
  expect_equal(
      -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1:2],D)),log=T)),
      evd_dnorm(mpar[1:2],D,quh),
  )
  expect_equal(
      -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1],D)),log=T)),
      evd_dnorm(mpar[1],D,quh)
  )
  
 
  
  
  
})


test_that("We can estimate the right answer with no confounding and direct simulation",{
  
  p <- 1e3
  g <- 1000
  D <- runif(p)/runif(p)
  D <- D/(sum(D)/p)
  N <- 1e6
  pvevec <- 0.5
  parvec <- RSSp:::calc_varu(pve = pvevec,p_n=sum(D)/N)
  rvec <- numeric(g)
  pb <- progress::progress_bar$new(total=g)
  # oquh <- rnorm(n = p)
  for(i in 1:g){
    quh <- rnorm(p,mean=0,sd=sqrt(pvevec*(D*D)+D))
    result <- RSSp_estimate(quh,D,sample_size = N,nterms = 2,calc_H = T)
    rvec[i] <- result$pve
    pb$tick()
  }
  
  expect_equal(mean(rvec-pvevec),0,tolerance=1e-1)
  expect_equal(mean(abs(rvec-pvevec)/pvevec),0,tolerance=2e-1)
  
  data_frame(me=rvec,truth=pvevec) %>% ggplot(aes(x=truth,y=(me-truth)/truth))+geom_point()+geom_smooth(method="lm")
  # data_frame(mr=rpar*(D*D)+D,tr=par*(D*D)+D) %>%
  #   ggplot(aes(x=tr,y=mr))+
  #   geom_hex()+
  #   geom_smooth(method="lm")+
  #   geom_abline(intercept=0,slope=1)
  expect_equal((rpar-par)/par,0,tolerance=1e-1)
  
  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(ts),log=T)),
               simple_log_likelihood(quh,ts))
  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(D^2*par+D)),log=T),
               evd_dnorm_t(par,D,quh))
})


testthat::test_that("stan likelihood works like R likelihood",{
  
  quh <- rnorm(1e2)
  D <- runif(1e2)/runif(1e2)  
  par <- c(runif(4))
  gen_D <- function(cvec,D){
    rD <- D
    for(i in 1:length(cvec)){
      rD <- rD+cvec[i]*D^(2-(i-1))
    }
    return(rD)
  }
  mpar <- par
  
  test_D <- gen_D(mpar,D)
  oD <- D+mpar[1]*D*D+mpar[2]*D^1+mpar[3]+mpar[4]/D
  expect_equal(test_D,oD)
  expect_equal(evd_dnorm(mpar,D,quh),
               evd_dnorm_stan(mpar,D,quh)
  )
  expect_equal(  
    -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar,D)),log=T)),
    evd_dnorm(mpar,D,quh),
  )
  expect_equal(
    -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1:3],D)),log=T)),
    evd_dnorm(mpar[1:3],D,quh)
  )
  expect_equal(
    -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1:2],D)),log=T)),
    evd_dnorm(mpar[1:2],D,quh),
  )
  expect_equal(
    -sum(dnorm(quh,mean=0,sd=sqrt(gen_D(mpar[1],D)),log=T)),
    evd_dnorm(mpar[1],D,quh)
  )
  
})



test_that("c++ based l-bfgs works like R version",{
  
  quh <- rnorm(1e2)
  D <- runif(1e2)/runif(1e2)  
  par <- c(runif(4))
  gen_D <- function(cvec,D){
    rD <- D
    for(i in 1:length(cvec)){
      rD <- rD+cvec[i]*D^(2-(i-1))
    }
    return(rD)
  }
  mpar <- par
  
  test_D <- gen_D(mpar,D)
  test_quh <- rnorm(length(D),mean=0,sd = sqrt(test_D))
  test_N <- 10000
  tpve <- estimate_pve(mpar,D,quh,test_N)
  np <- length(mpar)
  mres <- RSSp_estimate(quh = test_quh,D = test_D,nterms=1,sample_size = test_N)
  nmpar <- runif(1)
  cppres <- rssp_cpp(par = nmpar,D = D,quh = quh)
})
