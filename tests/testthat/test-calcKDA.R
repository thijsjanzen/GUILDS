context("calcKDA")

test_that("calcKDA use", {
  sad <- 1:20
  kda_1 <- calcKDA(sad)

  for(i in 1:20) {
    testthat::expect_equal(kda_1[i], 10361.632918)
  }

  testthat::expect_silent( calcKDA(c(sad,1000))   )
  testthat::expect_silent( calcKDA(c(sad,10000))  )
  testthat::expect_silent( calcKDA(c(sad,25000))  )
})