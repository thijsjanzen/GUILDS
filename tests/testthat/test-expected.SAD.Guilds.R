context("expected.SAD.Guilds")

test_that("expected.SAD.Guilds: use", {
  SAD <- expected.SAD.Guilds(theta = 200,
                              alpha_x = 0.1,
                              alpha_y = 0.01,
                              J = 1000,
                              n_replicates = 3)

  S1 <- sum(SAD$guildX)
  S2 <- sum(SAD$guildY)
  testthat::expect_gt(S1, S2) #because alpha_x > alpha_y

  a <- pm_sadaux(x = 1, I = 10, th = 100, j = 1000, k = 100)

  testthat::expect_equal(is.infinite(a),TRUE)
})

test_that("expected.SAD.Guilds: abuse", {
  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 200,
                               alpha_x = 1.0,
                               alpha_y = 1.0,
                               J = 1000,
                               n_replicates = 10),
    "alpha_x and alpha_y are both one"
  )

  J <- 100
  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = -20, alpha_x = 0.09, alpha_y = 0.5, J),
    "theta can not be below one")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 20, alpha_x = -0.09, alpha_y = 0.5, J),
    "alpha_x can not be below zero")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 20, alpha_x = 0.09, alpha_y = -0.5, J),
    "alpha_y can not be below zero")

  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 20, alpha_x = 1.09, alpha_y = 0.5, J),
    "alpha_x can not be above 1")
  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 20, alpha_x = 0.09, alpha_y = 1.5, J),
    "alpha_y can not be above 1")
  testthat::expect_error(
    SAD <- expected.SAD.Guilds(theta = 20, alpha_x = 0.09, alpha_y = 0.5, J = -1),
    "J can not be below one")
})
