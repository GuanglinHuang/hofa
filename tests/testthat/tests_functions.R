library(testthat)
set.seed(1)

test_that("M3M works correctly", {
  t = 100;n = 10;
  x = matrix(rnorm(1000), nrow = t, ncol = n)
  m3 = PerformanceAnalytics::M3.MM(x)
  m3m3 = m3%*%t(m3)
  m3m = hofa::M3M(scale(x,scale = F))/t^2
  expect_equal(m3m, m3m3)

})

test_that("M4M works correctly", {
  t = 100;n = 10;
  x = matrix(rnorm(1000), nrow = t, ncol = n)
  m4 = PerformanceAnalytics::M4.MM(x)
  m4m4 = m4%*%t(m4)
  m4m = hofa::M4M(scale(x,scale = F))/t^2
  expect_equal(m4m, m4m4)

})

test_that("C4M works correctly", {
  t = 100;n = 10;
  x = matrix(rnorm(t*n), nrow = t, ncol = n)
  x = scale(x,scale = F)
  c4 = hofa::CUM(x)
  c4c4 = c4%*%t(c4)
  c4m = hofa::C4M(x,cov(x)*(t-1)/t)
  expect_equal(c4m, c4c4)
})

test_that("Portfolio.Moments works correctly", {
  t = 100;n = 10;k = 2;
  w = rep(1/n,n)
  B = matrix(rnorm(k*n),n,k)
  m2f = runif(k,0.5,2)
  m3f = runif(k,-0.5,0.5)*m2f^1.5
  m4f = runif(k,3,6)*m2f^2
  c4f = m4f - 3*m2f^2
  m2e = runif(n,0.5,2)
  m3e = runif(n,-0.5,0.5)*m2e^1.5
  m4e = runif(n,3,6)*m2e^2
  c4e = m4e - 3*m2e^2
  mm_factor = list(m2f,m3f,m4f,c4f)
  mm_eps = list(m2e,m3e,m4e,c4e)

  m2x = B%*%diag(m2f,k,k)%*%t(B) + diag(m2e)
  m3x = B%*%SuperDiag(m3f,3)%*%(t(B)%x%t(B))+ SuperDiag(m3e,3)
  c4x = B%*%SuperDiag(c4f,4)%*%(t(B)%x%t(B)%x%t(B))+ SuperDiag(c4e,4)
  m4x = hofa::C4toM4(c4x,m2x)

  m2p = t(w)%*%m2x%*%w
  m3p = t(w)%*%m3x%*%(w%x%w)
  m4p = t(w)%*%m4x%*%(w%x%w%x%w)

  expect_equal(Portfolio.Moments(w,mm_factor,mm_eps,B), c(m2p,m3p,m4p))


})








