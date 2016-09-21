context("maxLikelihood.ESF")

test_that("maxLikelihood.ESF: use", {
 #first we generate data 
  set.seed(42)
  J = 10000
  theta = 100
  m <- 0.1
  I = m * (J-1) / (1 - m)
  
  v <- generate.ESF(theta, I, J)  
  LL <- maxLikelihood.ESF(initVals = c(100, 0.1),
                          Abund = c(v), verbose = FALSE)
  m_est <- LL$par[2]
  
  expect_equal(m, 
               m_est, 
               tolerance = 0.05)
  expect_equal(theta, 
               LL$par[1], 
               tolerance = 5)
})

test_that("maxLikelihood.ESF: abuse", {
  set.seed(42)
  J = 100
  theta = 100
  m <- 0.1
  I = m * (J-1) / (1 - m)
  
  v <- generate.ESF(theta, I, J)  
  expect_error(
    LL <- maxLikelihood.ESF(initVals = c(-100, 0.1),
                            Abund = v , verbose = FALSE),
    "initial theta can not be below one"
  )
  expect_error(
    LL <- maxLikelihood.ESF(initVals = c(100, -0.1),
                            Abund = v, verbose = FALSE),
    "initial m can not be below zero"
  )
  expect_error(
    LL <- maxLikelihood.ESF(initVals = c(100, 100),
                            Abund = v, verbose = FALSE),
    "initial m can not be above 1 (did you mean to enter I?)"
  )
  
  expect_error(
    LL <- maxLikelihood.ESF(initVals = c(100, 0.1),
                            Abund = c(1), verbose = FALSE),
    "Need more than 1 species in the dataset"
  )
})
  
  
  
  