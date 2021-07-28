library(testthat)

context("M2.select test")
test_that("M2.select sim test", {
  data("sp500")
  X = sp500[,1:50]
  ER     = M2.select(X,method = "ER")$R
  GR     = M2.select(X,method = "GR")$R
  BN_IC3 = M2.select(X,method = "BN-IC3")$R
  BN_PC3 = M2.select(X,method = "BN-PC3")$R
  BIC3   = M2.select(X,method = "BIC3")$R
  ON     = M2.select(X,method = "ON")$R
  ACT    = M2.select(X,method = "ACT")$R

  expect_equal(c(ER,GR,BN_IC3,BN_PC3,BIC3,ON,ACT),c(1,1,2,4,1,2,5))

})

context("M3.select test")
test_that("M3.select sim test", {
  data("sp500")
  X = sp500[,1:50]
  GER3  = M3.select(X,method = "GER3")
  GGR3  = M3.select(X,method = "GGR3")

  expect_equal(c(GER3$R,GGR3$R),c(1,1))
  expect_equal(c(GER3$Rh,GGR3$Rh),c(1,1))
  expect_equal(c(GER3$Rg,GGR3$Rg),c(0,0))
})


context("M4.select test")
test_that("M4.select sim test", {
  data("sp500")
  X = sp500[,1:50]
  GER4  = M4.select(X,method = "GER4")
  GGR4  = M4.select(X,method = "GGR4")

  expect_equal(c(GER4$R,GGR4$R),c(1,1))
  expect_equal(c(GER4$Rh,GGR4$Rh),c(1,1))
  expect_equal(c(GER4$Rg,GGR4$Rg),c(0,0))
})
