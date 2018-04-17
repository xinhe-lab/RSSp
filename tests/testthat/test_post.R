context("test posterior estimation")



test_that("Matrix form for computing D**",{
  
  g <- 10
  p <- 20
  sigu <- runif(g)
  confound <- runif( g)
  dvec <- runif(p)
  Rd <- matrix(0,p,g)

  for(i in 1:g){
    Rd[,i] <- posterior_mean_d(sigu[i],confound[i],dvec)
  }
  expect_equal(
    dim(Rd),
    c(p,g)
  )
  cd <- posterior_mean_D(
    sigu,
    confound,
    dvec)
  expect_equal(dim(Rd),dim(cd))
  expect_equal(cd,Rd)
})

test_that("matrix form computing U",{
  g <- 10
  p <- 20
  sigu <- runif(g)
  confound <- runif( g)
  dvec <- runif(p)
  quh <- matrix(runif(p*g),p,g)
  
  Q <- eigen(cor(matrix(rnorm(p*(p+g)),p+g,p)))$vectors
  
  gen_pm <- function(sigu,confound,dvec,Q,quh){
    p <- length(dvec)
    g <- ncol(quh)
    Rd <- matrix(0,p,g)
    for(i in 1:g){
      Rd[,i] <- posterior_mean_u(sigu[i],confound[i],dvec,Q=Q,quh = quh[,i])
    }
    return(Rd)
  }
Rd <- gen_pm(sigu,confound,dvec,Q,quh)
  cd <- posterior_mean_U(sigu,confound,dvec,quh,Q)
  expect_equal(Rd,cd)  
})

test_that("matrix form computing E[beta]",{
  
  g <- 10
  p <- 20
  sigu <- runif(g)
  confound <- runif( g)
  dvec <- runif(p)
  quh <- matrix(runif(p*g),p,g)
  Q <- eigen(cor(matrix(rnorm(p*(p+g)),p+g,p)))$vectors
  se <- matrix(rnorm(p*g),p,g)
  Rd <- matrix(0,p,g)
  for(i in 1:g){
    Rd[,i] <- posterior_mean_beta(sigu[i],confound[i],dvec,se[,i],Q,quh[,i])
  }
  cd <- posterior_mean_Beta(sigu,confound,dvec,quh,Q,se)
  expect_equal(Rd,cd)  
})

# test_that("matrix form computing Var[beta]",{
#   
#   g <- 10
#   p <- 20
#   sigu <- runif(g)
#   confound <- runif( g)
#   dvec <- runif(p)
#   quh <- matrix(runif(p*g),p,g)
#   Q <- eigen(cor(matrix(rnorm(p*(p+g)),p+g,p)))$vectors
#   se <- matrix(rnorm(p*g),p,g)
#   Rd <- matrix(0,p,g)
#   tRd <- posterior_var_beta(sigu = sigu[1],
#                             confound = confound[1],dvec = dvec,
#                             se,Q)
#   for(i in 1:g){
#     Rd[,i] <- posterior_mean_beta(sigu[i],confound[i],dvec,se[,i],Q,quh[,i])
#   }
#   cd <- posterior_mean_Beta(sigu,confound,dvec,quh,Q,se)
#   expect_equal(Rd,cd)  
# })


test_that("matrix form computing Y",{
  g <- 10
  p <- 20
  n <- 5
  sigu <- runif(g)
  confound <- runif( g)
  dvec <- runif(p)
  quh <- matrix(runif(p*g),p,g)
  Q <- eigen(cor(matrix(rnorm(p*(p+g)),p+g,p)))$vectors
  se <- matrix(rnorm(p*g),p,g)
  x <- matrix(rnorm(n*p),n,p)
  # posterior_mean_y(sigu[i],confound[i],dvec,se[,i],Q,quh[,i],x[i,])
  Rd <- matrix(0,n,g)
  for(i in 1:g){
    for(j in 1:n){
      Rd[j,i] <- posterior_mean_y(sigu[i],confound[i],dvec,se[,i],Q,quh[,i],x[j,])
    }
  }
  cd <- posterior_mean_Y(sigu,confound,dvec,quh,Q,se,x)
  expect_equal(Rd,cd)  
})


test_that("matrix form computing beta with too many estimates gives an error",{
  
  g <- 10
  p <- 20
  sigu <- runif(g*2)
  confound <- runif( g*2)
  dvec <- runif(p)
  quh <- matrix(runif(p*g),p,g)
  Q <- eigen(cor(matrix(rnorm(p*(p+g)),p+g,p)))$vectors
  se <- matrix(rnorm(p*g),p,g)
  Rd <- matrix(0,p,g)
  for(i in 1:g){
    Rd[,i] <- posterior_mean_beta(sigu[i],confound[i],dvec,se[,i],Q,quh[,i])
  }
  expect_error(posterior_mean_Beta(sigu,confound,dvec,quh,Q,se))
})




