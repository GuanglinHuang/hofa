library(testthat)
set.seed(1234)

context("hofa.sim test")
test_that("hofa.sim test", {
  n = 100
  t = 200
  k = 2
  par_f = c(0.5,0,2,2,Inf,Inf)
  par_e = c(1,0,2,Inf,0,0,0)
  rho_ar = c(0.5,0.2)

  data = hofa.sim(n,t,k,par_f,par_e,rho_ar)

  expect_equal(svd(data$X)$d[1],164.2617,tolerance = 10^-6)

})
