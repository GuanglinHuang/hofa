library(testthat)


context("hofa.sim test")

test_that("hofa.DGP1 test", {
  set.seed(1234)
  n = 100
  t = 200
  k = 2
  par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
  par_e = list(1,0,2,Inf)
  rho_f = c(0.5,0.2)
  alpha = 0
  data = hofa::hofa.DGP1(n,t,k,par_f,par_e,alpha,rho_f,rho_e = 0.2)

  expect_equal(svd(data$X)$d[1],211.711820268,tolerance = 10^-6)

})

test_that("hofa.DGP2 test", {
  set.seed(1234)
  n = 100
  t = 200
  k = 2
  par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
  par_e = list(1,0,2,Inf)
  rho_f = c(0.5,0.2)
  par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))

  data = hofa::hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)

  expect_equal(svd(data$X)$d[1],165.56921837,tolerance = 10^-6)

})


