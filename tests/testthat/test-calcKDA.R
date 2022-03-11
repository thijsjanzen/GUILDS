context("calcKDA")

test_that("calcKDA use", {
  sad <- 1:20
  kda_1 <- GUILDS::calcKDA(sad)

  for (i in 1:20) {
    testthat::expect_equal(kda_1[i], 10361.632918)
  }
  testthat::expect_equal(kda_1[21], 0)

  kda_2 <- calcKDA(c(sad, 1000))
  kda_3 <- calcKDA(c(sad, 10000))
  kda_4 <- calcKDA(c(sad, 25000))


  kda_old <- calcKDA_old(sad)
  kda_new <- calcKDA_new(sad)
  diff <- kda_old - kda_new
  testthat::expect_equal(mean(diff), 0.0, tolerance = 0.1)


})