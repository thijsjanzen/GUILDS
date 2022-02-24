context("preston_plot")

test_that("preston_plot use", {
  theta <- 100
  m <- 0.1
  J <- 10000
  I <- m * (J - 1) / (1 - m)

  abund <- generate.ESF(theta, I, J)
  abund.expect <- expected.SAD(theta, m, J)

  old_par <- par(no.readonly = TRUE)
  testthat::expect_silent(

    par(mfrow = c(1, 2) )
  )
  testthat::expect_silent(
    preston_plot(abund)
  )
  testthat::expect_silent(
    preston_plot(abund, abund.expect)
  )
  par(old_par)
})
