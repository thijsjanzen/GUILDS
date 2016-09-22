context("expected.SAD.Guilds.Conditional")

test_that("expected.SAD.Guilds.Conditional: use", {
  SAD <- expected.SAD.Guilds.Conditional(theta = 200, 
                                         alpha_x = 0.1,
                                         alpha_y = 0.01,
                                         Jx = 1000,
                                         Jy = 1000,
                                         n_replicates = 10)
  
  S1 <- sum(SAD$guildX)
  S2 <- sum(SAD$guildY)
  a <- S1 > S2
  expect_equal(a,TRUE) #because alpha_x > alpha_y
  
  SAD <- expected.SAD.Guilds.Conditional(theta = 200, 
                                         alpha_x = 0.1,
                                         alpha_y = 0.1,
                                         Jx = 3000,
                                         Jy = 1000,
                                         n_replicates = 10)
  
  S1 <- sum(SAD$guildX)
  S2 <- sum(SAD$guildY)
  a <- S1 > S2
  expect_equal(a,TRUE) #because Jx > Jy
})




