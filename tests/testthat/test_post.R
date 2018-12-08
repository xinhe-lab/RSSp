context("test posterior estimation")


test_that("C++ implementation of SPVE",{
  
  
  p <- 4
  quh <- rnorm(p)
  D <- runif(p)/runif(p)  
  par <- c(runif(1))
  N <- 100
  
  R_estimate_pve <- function(cvec,D,quh,N,n_samples=0){
    
    num_c <- length(cvec)
    
    if(num_c>1){
      stop("multiple cvec terms not yet implemented")
    }
    if(n_samples!=0){
      stop("sampling based pve estimate not yet implemented")
    }
    # lambda_init <- 0
    tmp <- 1+1/(cvec*D)
    rto <- (quh^2)/((D)^2)
    rto <- D*rto
    rto <- rto/(tmp^2)
    
    pve_mean <- sum(1/tmp)+sum(rto)
    pve_mean <- pve_mean/N
    
    
    return(pve_mean)
    
  }
  
  
  
  expect_equal(R_estimate_pve(par,D = D,quh = quh,N = N,n_samples = 0),
               estimate_pve(cvec = par,D = D,quh = quh,sample_size =N))
  
  # 
  # spve_R <- function(R,b,bh,se,n){
  #   tbh <- b/sqrt(n*se+bh^2)
  #   return(sum((tb%o%tb)*R))
  # }
  
  
})



