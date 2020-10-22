context("simulating data")



test_that("we can simulate from and estimate multiparameter model",{
  
  data("panel_eigenvalues")
  D <-  panel_eigenvalues
  p <- length(D)
  nterms <- 1
  N <- 5000
  cvec <- runif(nterms)
  quh_l <- purrr::rerun(100,RSSp:::simulate_vhat(cvec,D))
  pvevec <- purrr::map_dbl(quh_l,~estimate_pve(cvec = cvec,D = D,quh = .x,sample_size = N))  
  resvec <- purrr::map_dbl(quh_l,~RSSp_estimate(.x,D,sample_size=N,nterms = 1)$pve)
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



