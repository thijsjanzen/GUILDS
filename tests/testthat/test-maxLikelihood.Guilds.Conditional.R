context("maxLikelihood.Guilds.Conditional")

test_that("maxLikelihood.GuildsConditional: use", {
#  skip_on_cran()
  set.seed(42)

  theta <- 100
  alpha_x <- 0.1

  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_x,
                                     JX = 100, JY = 100)

  #initial parameters for the D0 model c(theta,alpha)
  LL <- maxLikelihood.Guilds.Conditional( init_vals = c(theta, alpha_x),
                              model="D0",
                              sadx = simul_data$guildX,
                              sady = simul_data$guildY, verbose = FALSE)
  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.1, scale = 1
  )

  set.seed(666)
  theta <- 1000
  alpha_x <- 0.1
  alpha_y <- 0.001

  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_y,
                                     JX = 100, JY = 100)

  #initial parameters for the D1 model c(theta, alpha_x, alpha_y)
  LL <- maxLikelihood.Guilds.Conditional( init_vals =
                                            c(theta, alpha_x, alpha_y),
                                          model = "D1",
                                          sadx  = simul_data$guildX,
                                          sady  = simul_data$guildY,
                                          verbose = FALSE)

  testthat::expect_equal(
    theta,
    LL$par[1],
    tolerance = 15, scale = 1
  )
  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.01, scale = 1
  )
  testthat::expect_equal(
    alpha_y,
    LL$par[3],
    tolerance = 0.005, scale = 1
  )

  # to test for alpha_x = alpha_y = 1, which leads to I_X = Inf
  testthat::expect_silent(
    simul_data <- generate.Guilds.Cond(theta = 10,
                                     alpha_x = 1.0,
                                     alpha_y = 1.0,
                                     JX = 100, JY = 100)

  )
})


test_that("maxLikelihood.Guilds: abuse", {
#  skip_on_cran()
  set.seed(42)
  J <- 200

  theta <- 20
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(-50, 0.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial theta can not be below one"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, -0.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be below zero"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, 1.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be above 1"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, 0.1, -0.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be below 0"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, 0.1, 1.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be above 1"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, 0.1, 1.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )

  expect_error(
    maxLikelihood.Guilds.Conditional( init_vals = c(50, 0.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )
})