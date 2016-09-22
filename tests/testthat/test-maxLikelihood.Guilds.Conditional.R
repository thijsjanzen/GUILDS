context("maxLikelihood.Guilds.Conditional")

test_that("maxLikelihood.GuildsConditional: use", {
  set.seed(42)
  initParams <- c(20, 0.1); #Initial parameters for the D0 model, c(theta,alpha)
  maxLikelihood.Guilds.Conditional(initParams, model ="D0", 
                                   method = "subplex",
                                   sadx = 1:20, sady = 1:20, verbose = FALSE)
  
  
  
  
  set.seed(42)
  
  theta = 50
  alpha_x = 0.1
  
  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_x, JX = 10000, JY = 10000)
  
  #initial parameters for the D0 model c(theta,alpha)
  LL <- maxLikelihood.Guilds.Conditional( init_vals = c(theta, alpha_x), 
                              model="D0", method="simplex",
                              sadx = simul_data$guildX, 
                              sady = simul_data$guildY, verbose = FALSE)
  expect_equal(
    theta, 
    LL$par[1],
    tolerance = 10, scale = 1
  )
  expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )
  
  set.seed(666)
  theta = 30
  alpha_x = 0.01
  alpha_y = 0.001
  
  simul_data <- generate.Guilds.Cond(theta, alpha_x, alpha_y, JX = 10000, JY = 10000)
  
  #initial parameters for the D1 model c(theta, alpha_x, alpha_y)
  LL <- maxLikelihood.Guilds.Conditional( init_vals = c(30, 0.1, 0.001), 
                                          model="D1", method="simplex",
                                          sadx = simul_data$guildX, 
                                          sady = simul_data$guildY, verbose = FALSE)
  
  expect_equal(
    theta, 
    LL$par[1],
    tolerance = 7, scale = 1
  )
  expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.01, scale = 1
  )
  expect_equal(
    alpha_y,
    LL$par[3],
    tolerance = 0.0005, scale = 1
  )
})