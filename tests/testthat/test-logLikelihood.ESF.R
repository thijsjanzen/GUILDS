context("logLikelihood.ESF")

test_that("logLikelihood.ESF: use", {
  set.seed(42)
  J <- 1000
  theta <- 100
  m <- 0.1
  I <- m * (J - 1) / (1 - m)

  abund <- generate.ESF(theta, I, J)
  LL1 <- logLikelihood.ESF(theta, m, abund)
  LL2 <- logLikelihood.ESF(theta * 2, m * 2, abund)

  testthat::expect_gt(LL1, LL2)
})

test_that("logLikelihood.ESF: abuse", {
  set.seed(42)
  J <- 1000
  theta <- 100
  m <- 0.1
  I <- m * (J - 1) / (1 - m)

  abund <- generate.ESF(theta, I, J)
  LL1 <- logLikelihood.ESF(-10, m, abund)
  testthat::expect_equal(is.infinite(LL1[[1]]), TRUE)
  LL1 <- logLikelihood.ESF(theta, 2, abund)
  testthat::expect_equal(is.infinite(LL1[[1]]), TRUE)

  # enter NA values and expect error:

  testthat::expect_true(is.infinite(logLikelihood.ESF(5, NA, abund)))
  testthat::expect_true(is.infinite(logLikelihood.ESF(NA, 0.1, abund)))



})