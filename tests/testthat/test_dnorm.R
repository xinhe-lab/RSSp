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


test_that("fancy templated likelihood works like toy R version without confounding",{
  
  p <- 4
  quh <- rnorm(p)
  D <- runif(p)/runif(p)  
  par <- c(runif(1))

  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(par*(D*D)+D),log=T)),
               evd_dnorm(par,D = D,quh = quh))

  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(ts),log=T)),
               simple_log_likelihood(quh,ts))
  expect_equal(-sum(dnorm(quh,mean=0,sd=sqrt(D^2*par+D)),log=T),
               evd_dnorm_t(par,D,quh))
  

})