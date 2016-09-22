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
  
  
  a <- pm_sadaux(x=1, I = 10, th = 100, j = 1000, k = 100)
  expect_equal(is.infinite(a),TRUE)
})

test_that("expected.SAD.Guilds: abuse", {
  expect_error(
    SAD <- expected.SAD.Guilds(theta = 200, 
                               alpha_x = 1.0,
                               alpha_y = 1.0,
                               J = 1000,
                               n_replicates = 10),
    "alpha_x and alpha_y are both one"
  )
})