context("simplex")

test_that("simplex use", {
  a <- 3
  xvals <- 1:10
  yvals <- a * xvals + rnorm(10, mean = 0, sd = 0.1)

  evalfunc <- function(x) {
    sim_y <- xvals * x
    diff <- yvals - sim_y
    ll <- sum(dnorm(diff, mean = 0, sd = 0.1, log = TRUE))
    return(-ll)
  }

testthat::expect_silent(
  v1 <- simplex(initpars = c(2), evalfunc,
               verbose  = FALSE,
               maxiter  = 1000)
  )
testthat::expect_output(
  v2 <- simplex(initpars = c(2), evalfunc,
               verbose  = TRUE,
               maxiter  = 10)
)
  testthat::expect_equal(v1$par, a, tolerance = 1e-2)
  testthat::expect_equal(v1$par, v2$par, tolerance = 1e-3)
})