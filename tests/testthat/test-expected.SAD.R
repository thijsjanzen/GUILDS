context("expected.SAD")

test_that("expected.SAD: use", {
  J <- 10000
  testthat::expect_silent(
    SAD <- expected.SAD(theta = 200, m = 0.09, J)
  )
})

test_that("expected.SAD: abuse", {
  J <- 10000
  expect_error(
      SAD <- expected.SAD(theta = -20, m = 0.09, J),
      "theta can not be below one")

  expect_error(
    SAD <- expected.SAD(theta = 20, m = -0.09, J),
    "m can not be below zero")

  expect_error(
    SAD <- expected.SAD(theta = 20, m = 0.09, J = -10),
    "J can not be below zero")

  expect_error(
    SAD <- expected.SAD(theta = 20, m = 1.09, J),
    "m can not be above 1")
})
