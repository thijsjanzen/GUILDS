context("maxLikelihood.Guilds")

test_that("maxLikelihood.Guilds: use", {
  set.seed(42)
  J = 20000
  
  theta = 50
  alpha_x = 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)
  
  #initial parameters for the D0 model c(theta,alpha)
  LL <- maxLikelihood.Guilds( initVals = c(50, 0.1), 
                             model="D0", method="simplex",
                             SADX = simul_data$guildX, 
                             SADY = simul_data$guildY, verbose = FALSE)
   
  expect_equal(
    theta, 
    LL$par[1],
    tolerance = 5, scale = 1
  )
  expect_equal(
    alpha_x,
    LL$par[2],
    tolerance = 0.05, scale = 1
  )
  
  
  J = 20000
  
  theta = 50
  alpha_x = 0.1
  alpha_y = 0.001
  
  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  
  #initial parameters for the D1 model c(theta, alpha_x, alpha_y)
  LL <- maxLikelihood.Guilds( initVals = c(50, 0.1, 0.001), 
                              model="D1", method="simplex",
                              SADX = simul_data$guildX, 
                              SADY = simul_data$guildY, verbose = FALSE)
  
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
  expect_equal(
    alpha_y,
    LL$par[3],
    tolerance = 0.0005, scale = 1
  )
  
  
  
})