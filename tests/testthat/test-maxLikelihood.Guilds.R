context("maxLikelihood.Guilds")

test_that("maxLikelihood.Guilds: use", {
  skip_on_cran() # takes too long
  testthat::skip("takes too long")
  set.seed(42)
  J <- 1000

  theta <- 1000
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
#  testthat::expect_output(
  LL <- GUILDS::maxLikelihood.Guilds(init_vals = c(theta, alpha_x),
                             model = "D0",
                             sadx  = simul_data$guildX,
                             sady  = simul_data$guildY, verbose = FALSE)
 # )
  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )

  J <- 1000

  theta <- 100
  alpha_x <- 0.1
  alpha_y <- 0.001

  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)

  #testthat::expect_output(
  #initial parameters for the D1 model c(theta, alpha_x, alpha_y)
  LL <- maxLikelihood.Guilds( init_vals = c(theta, alpha_x, alpha_y),
                              model = "D1",
                              sadx  = simul_data$guildX,
                              sady  = simul_data$guildY, verbose = FALSE)
  #)

  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )
  testthat::expect_equal(
    alpha_y,
    LL$par[3],
    tolerance = 0.01, scale = 1
  )

  testthat::expect_gt(LL$par[2], LL$par[3])

  testthat::expect_output(
  LL1 <- GUILDS::maxLikelihood.Guilds(init_vals = c(50, 0.1),
                              model = "D0",
                              sadx  = 1:3,
                              sady  = 1:3, verbose = TRUE)
  )
})


test_that("maxLikelihood.Guilds: abuse", {
  set.seed(42)
  J <- 100

  theta <- 50
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(-50, 0.1),
                              model = "D0",
                              sadx  = simul_data$guildX,
                              sady  = simul_data$guildY, verbose = FALSE),
    "initial theta can not be below one"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, -0.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be below zero"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 1.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be above 1"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, -0.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be below 0"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, 1.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be above 1"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, 1.1),
                          model = "D0",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1),
                          model = "D1",
                          sadx  = simul_data$guildX,
                          sady  = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )
})
