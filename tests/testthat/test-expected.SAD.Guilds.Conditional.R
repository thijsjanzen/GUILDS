context("expected.SAD.Guilds.Conditional")

test_that("expected.SAD.Guilds.Conditional: use", {
  SAD <- expected.SAD.Guilds.Conditional(theta = 200,
                                         alpha_x = 0.1,
                                         alpha_y = 0.01,
                                         Jx = 200,
                                         Jy = 200,
                                         n_replicates = 2)
  S1 <- sum(SAD$guildX)
  S2 <- sum(SAD$guildY)
  testthat::expect_gt(S1, S2) # because alpha_x > alpha_y

  SAD <- expected.SAD.Guilds.Conditional(theta = 200,
                                         alpha_x = 0.1,
                                         alpha_y = 0.1,
                                         Jx = 200,
                                         Jy = 100,
                                         n_replicates = 2)
  testthat::expect_gt(S1, S2) # because Jx > Jy
})

test_that("expected.SAD.Guilds.Conditional: abuse", {
  n_replicates <- 2
  Jx <- 10
  Jy <- 20

  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 200,
                               alpha_x = 1.0,
                               alpha_y = 1.0,
                               Jx, Jy,
                               n_replicates),
    "alpha_x and alpha_y are both one, leading to I_x = I_y = Inf"
  )

  J <- 100
  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = -20,
                                           alpha_x = 0.09,
                                           alpha_y = 0.5,
                                           Jx, Jy, n_replicates),
    "theta can not be below one")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = -0.09,
                                           alpha_y = 0.5,
                                           Jx, Jy, n_replicates),
    "alpha_x can not be below zero")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = 0.09,
                                           alpha_y = -0.5,
                                           Jx, Jy, n_replicates),
    "alpha_y can not be below zero")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = 1.09,
                                           alpha_y = 0.5,
                                           Jx, Jy, n_replicates),
    "alpha_x can not be above 1")
  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = 0.09,
                                           alpha_y = 1.5,
                                           Jx, Jy, n_replicates),
    "alpha_y can not be above 1")
  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = 0.09,
                                           alpha_y = 0.5,
                                           Jx = -1, Jy, n_replicates),
    "Jx can not be below one")
  testthat::expect_error(
    SAD <- expected.SAD.Guilds.Conditional(theta = 20,
                                           alpha_x = 0.09,
                                           alpha_y = 0.5,
                                           Jx, Jy = -1, n_replicates),
    "Jy can not be below one")
})
