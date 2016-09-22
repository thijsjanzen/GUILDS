context("logLikelihood.Guilds")

test_that("logLikelihood.Guilds: use", {
  skip("WIP")
  set.seed(42)
  J = 10000
  theta = 100
  alpha_x = 0.1
  alpha_y = alpha_x
  
  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  LL1 <- logLikelihood.Guilds(parameters=c(theta, alpha_x, alpha_y), 
                              model = "D0", 
                              simul_data$guildX, simul_data$guildY)
  LL2 <- logLikelihood.Guilds(parameters=c(theta * 2, alpha_x * 2, alpha_y * 2), 
                              model = "D0", simul_data$guildX, simul_data$guildY)
  
  LL3 <- logLikelihood.Guilds(parameters=c(theta, alpha_x, alpha_y), 
                              model = "D1", simul_data$guildX, simul_data$guildY)
  
  
  
  a <- LL1[[1]] > LL2[[1]]
  expect_equal(a,TRUE)
  
  alpha_y = 0.01
  
  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  LL1 <- logLikelihood.Guilds(parameters=c(theta, alpha_x, alpha_y), 
                              model = "D1", 
                              simul_data$guildX, simul_data$guildY)
  LL2 <- logLikelihood.Guilds(parameters=c(theta * 2, alpha_x * 2, alpha_y * 2), 
                              model = "D1", simul_data$guildX, simul_data$guildY)
  
  a <- LL1[[1]] > LL2[[1]]
  expect_equal(a,TRUE)
  
  
  
  
})

