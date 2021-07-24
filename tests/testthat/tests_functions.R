library(testthat)
set.seed(1)

test_that("colCumsum works correctly", {
  t = 100;n = 10;
  x = matrix(rnorm(1000), nrow = t, ncol = n)
  m3 = PerformanceAnalytics::M3.MM(x)
  m3m3 = m3%*%t(m3)
  m3m = hofa::M3M(scale(x,scale = F))/t^2
  expect_equal(m3m, m3m3)

})
