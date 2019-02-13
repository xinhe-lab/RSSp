context("simulating data")

test_that("we can simulate from  2-parameter model",{
  
  data("panel_eigenvalues")
  D <-  panel_eigenvalues
  p <- length(D)
  nterms <- 2
  cvec <- runif(nterms)
  quh_1 <- purrr::flatten_dbl(purrr::rerun(100,RSSp:::simulate_vhat(cvec,D)))
  quh_2 <- rnorm(p*100,mean = 0,sd = sqrt(cvec[1]*D^2 + cvec[2]*D + D))
  expect_equal(mean(quh_1), mean(quh_2), tolerance = 1e-1)
  expect_equal(var(quh_1),var(quh_2), tolerance = 1e-1)
  
})


test_that("we can simulate from and estimate multiparameter model",{
  
  data("panel_eigenvalues")
  D <-  panel_eigenvalues
  p <- length(D)
  nterms <- 3
  N <- 5000
  cvec <- runif(nterms)
  quh_l <- purrr::rerun(100,RSSp:::simulate_vhat(cvec,D))
  pvevec <- purrr::map_dbl(quh_l,~estimate_pve(cvec = cvec,D = D,quh = .x,sample_size = N))  
  resvec <- purrr::map_dbl(quh_l,~RSSp_estimate(.x,D,sample_size=N,nterms = 3)$pve)
  expect_equal(mean(resvec),mean(pvevec),tolerance=1e-2)
})

test_that("we can simulate from and estimate multiparameter model without a gradient",{
  
  data("panel_eigenvalues")
  D <-  panel_eigenvalues
  p <- length(D)
  nterms <- 3
  N <- 5000
  cvec <- runif(nterms)
  quh_l <- purrr::rerun(100,RSSp:::simulate_vhat(cvec,D))
  pvevec <- purrr::map_dbl(quh_l,~estimate_pve(cvec = cvec,D = D,quh = .x,sample_size = N))  
  resvec <- purrr::map_dbl(quh_l,~RSSp_estimate(.x,D,sample_size=N,nterms = 3,useGradient = F)$pve)
  expect_equal(mean(resvec),mean(pvevec),tolerance=1e-2)
})



