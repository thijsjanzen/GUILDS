context("logLikelihood.Guilds.Conditional")

test_that("logLikelihood.Guilds.Conditional: use", {
#  skip_on_cran()
  set.seed(666)
  J <- 1000
  theta <- 100
  alpha_x <- 0.1
  alpha_y <- alpha_x

  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_y,
                                     JX = J / 2, JY = J / 2)
  LL1 <- logLikelihood.Guilds.Conditional(parameters = c(theta, alpha_x),
                              model = "D0",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)
  LL2 <- logLikelihood.Guilds.Conditional(parameters = c(3, 0.00001),
                              model = "D0",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)

  testthat::expect_gt(LL1, LL2)

  set.seed(42 + 666)
  alpha_y <- 0.01

  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  LL1 <- logLikelihood.Guilds.Conditional(parameters =
                                             c(theta, alpha_x, alpha_y),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)
  LL2 <- logLikelihood.Guilds.Conditional(parameters = c(3, 0.99, 0.5),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)

  testthat::expect_gt(LL1, LL2)

  # to test verbose == TRUE
  testthat::expect_output(
   LL2 <- logLikelihood.Guilds.Conditional(parameters = c(3, 0.99, 0.5),
                                          model = "D1",
                                          simul_data$guildX,
                                          simul_data$guildY,
                                          verbose = TRUE)
  )
})

test_that("logLikelihood.Guilds.Conditional: abuse", {
  #skip_on_cran()
  set.seed(666)
  J <- 1000
  theta <- 100
  alpha_x <- 0.1
  alpha_y <- alpha_x

  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_y,
                                     JX = J / 2,
                                     JY = J / 2)
  testthat::expect_error(
    logLikelihood.Guilds.Conditional(parameters = c(theta, alpha_x, alpha_x),
                                     model = "D0",
                                     simul_data$guildX, simul_data$guildY,
                                     verbose = FALSE),
    "Input vector is of incorrect length"
  )
})
