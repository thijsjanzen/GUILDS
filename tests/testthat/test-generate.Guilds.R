context("generate.Guilds")

test_that("generate.Guilds: use", {
  J = 10000
  v <- generate.ESF(theta = 100, I = 10, J)
  expect_equal(
    sum(v),
    J
  )
})

test_that("generate.Guilds: abuse", {
  expect_error(
    generate.ESF(theta = -1, I = 10, J = 100),
    "theta can not be below zero"
  )
  expect_error(
    generate.ESF(theta = 100, I = -10, J = 100),
    "I can not be below zero"
  )
  expect_error(
    generate.ESF(theta = 100, I = 10, J = -100),
    "J can not be below zero"
  )
})