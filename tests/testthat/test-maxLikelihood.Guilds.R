context("maxLikelihood.Guilds")

test_that("maxLikelihood.Guilds: use", {
  skip_on_cran()
  set.seed(42)
  J <- 1000

  theta <- 1000
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
  LL <- maxLikelihood.Guilds(init_vals = c(theta, alpha_x),
                             model="D0",
                             sadx = simul_data$guildX,
                             sady = simul_data$guildY,
                             verbose = FALSE)

  testthat::expect_gte(
    LL$par[1],
    theta
  )
  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )

  set.seed(1)
  J <- 1000

  theta <- 100
  alpha_x <- 0.1
  alpha_y <- 0.001

  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)

  #initial parameters for the D1 model c(theta, alpha_x, alpha_y)
    LL <- maxLikelihood.Guilds( init_vals = c(theta, alpha_x, alpha_y),
                              model="D1",
                              sadx = simul_data$guildX,
                              sady = simul_data$guildY, verbose = FALSE)
  testthat::expect_equal(
    theta,
    LL$par[1],
    tolerance = 55, scale = 1
  )
  testthat::expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )
  testthat::expect_equal(
    alpha_y,
    LL$par[3],
    tolerance = 0.001, scale = 1
  )

   LL1 <- maxLikelihood.Guilds( init_vals = c(50, 0.1),
                              model = "D0",
                              sadx = 1:20,
                              sady = 1:20, verbose = FALSE)

  LL2 <- maxLikelihood.Guilds( init_vals = c(50, 0.1),
                              model = "D0",  #subplex before
                              sadx = 1:20,
                              sady = 1:20, verbose = FALSE)
  LL3 <- maxLikelihood.Guilds( init_vals = c(50, 0.1),
                              model = "D0",
                              sadx = 1:20,
                              sady = 1:20, verbose = FALSE)


  a <- LL1$par
  b <- LL2$par
  c <- LL3$par

  testthat::expect_equal(a, b, tolerance = 0.1, scale = 1)
  testthat::expect_equal(a, c, tolerance = 0.1, scale = 1)
})


test_that("maxLikelihood.Guilds: abuse", {
  #skip_on_cran()
  cat("maxLikelihood.Guilds: abuse\n")
  set.seed(42)
  J <- 20000

  theta <- 50
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(-50, 0.1),
                              model="D0",
                              sadx = simul_data$guildX,
                              sady = simul_data$guildY, verbose = FALSE),
    "initial theta can not be below one"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, -0.1),
                          model="D0",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be below zero"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 1.1),
                          model="D0",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "initial alpha can not be above 1"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, -0.1),
                          model="D1",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be below 0"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, 1.1),
                          model="D1",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "initial alpha_y can not be above 1"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1, 1.1),
                          model="D0",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )

  testthat::expect_error(
    maxLikelihood.Guilds( init_vals = c(50, 0.1),
                          model="D1",
                          sadx = simul_data$guildX,
                          sady = simul_data$guildY, verbose = FALSE),
    "Input vector is of incorrect length"
  )
})
