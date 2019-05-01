context("evd_dnorm")

test_that("likelihood is approximately equal to transliterated MATLAB",{
  
  lik_fun <- system.file("m_code/negtwo_loglik.R",package = "RSSp")
  source(lik_fun)
  p <- 100
  quh <- rnorm(p)
  D <- runif(p)

  for(i in 1:4){
    par <- runif(n = i)
    retR <- evd_dnorm(par = par,D = D,quh = quh)
    retM <- negtwo_loglik(vhat = quh,dvec = D,cvec = par)
    expect_equal(retR,retM/2)
  }

})

test_that("fancy templated likelihood works like toy R version with confounding",{
  
  data("reference_genotype",package="ldshrink")
  data("reference_map",package="ldshrink")
  Rl <- ldshrink::ldshrink_evd(reference_genotype,reference_map)
  p <- length(Rl$D)
  # quh <- t(Rl$Q)%*%rnorm(p)
  D <-Rl$D  
  par <- 1
  N <- 100
  quh <- rnorm(p)
  # D <- runif(1e2)/runif(1e2)  
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
  
  data("reference_genotype",package="ldshrink")
  data("reference_map",package="ldshrink")
  Rl <- ldshrink::ldshrink_evd(reference_genotype,reference_map)
  p <- length(Rl$D)
  # quh <- t(Rl$Q)%*%rnorm(p)
  D <-Rl$D  
  N <- 1e6
  pvevec <- 0.5
  parvec <- RSSp:::calc_varu(pve = pvevec,p_n=sum(D)/N)
  g <- 1000
  rvec <- numeric(g)
  pb <- progress::progress_bar$new(total=g)
  # oquh <- rnorm(n = p)
  for(i in 1:g){
    quh <- rnorm(p,mean=0,sd=sqrt(parvec*(D*D)+D))
    result <- RSSp_estimate(quh,D,sample_size = N,nterms = 1)
    rvec[i] <- result$pve
    pb$tick()
  }
  
  expect_equal(mean(rvec-pvevec),0,tolerance=1e-1)
  #expect_equal(mean(abs(rvec-pvevec)/pvevec),0,tolerance=2e-1)
    # data_frame(mr=rpar*(D*D)+D,tr=par*(D*D)+D) %>%
  #   ggplot(aes(x=tr,y=mr))+
  #   geom_hex()+
  #   geom_smooth(method="lm")+
  #   geom_abline(intercept=0,slope=1)
})
# testthat::test_that("rssp is consistent with Xiang's implementation",{
#  res_data <- readRDS("inst/rand_seed_312_k_1_rep_100.RDS")
# tot_res <- map2_dfr(array_branch(res_data$vhat,margin = 2),
#          array_branch(res_data$dvec,margin = 2),~RSSp::RSSp_estimate(.x,.y,sample_size = 1000,pve_bounds = c(.Machine$double.eps,8000)))      
# expect_equal(c(res_data$obj*2),tot_res$lnZ)
#   res_data$vhat
#   res_data$dv2
# par_l <- transpose(list(quh=array_branch(res_data$vhat,margin = 2),
#          dvec=array_branch(res_data$dvec,margin = 2),
#      par=res_data$par))
# pm_db <- map_dbl(par_l,~estimate_pve(cvec=.x$par,D=.x$dvec,quh = .x$quh,sample_size = 1000))
# expect_equal(pm_db,c(res_data$pve))
# })


testthat::test_that("stan likelihood works like R likelihood",{
  
  data("reference_genotype",package="ldshrink")
  data("reference_map",package="ldshrink")
  Rl <- ldshrink::ldshrink_evd(reference_genotype,reference_map)
  p <- length(Rl$D)
  # quh <- t(Rl$Q)%*%rnorm(p)
  D <-Rl$D  
  quh <- rnorm(length(D))
  #D <- runif(1e2)/runif(1e2)  
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

