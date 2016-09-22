context("expected.SAD.Guilds")

test_that("expected.SAD.Guilds: use", {
  SAD <- expected.SAD.Guilds(theta = 200, 
                              alpha_x = 0.1,
                              alpha_y = 0.01,
                              J = 1000,
                              n_replicates = 10)
  
  S1 <- sum(SAD$guildX)
  S2 <- sum(SAD$guildY)
  a <- S1 > S2
  expect_equal(a,TRUE) #because alpha_x > alpha_y

  SAD <- expected.SAD.Guilds(theta = 200, 
                             alpha_x = 1,
                             alpha_y = 1,
                             J = 1000,
                             n_replicates = 10)
  
})




