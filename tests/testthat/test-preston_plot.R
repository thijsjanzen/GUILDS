context("preston_plot")

test_that("preston_plot use", {
  theta = 100
  m = 0.1
  I = m * (J - 1) / (1 - m)
  J = 10000
  
  abund <- generate.ESF(theta, I, J)
  par(mfrow = c(1,2))
  preston_plot(abund)
  abund.expect <- expected.SAD(theta, m, J)
  preston_plot(abund, abund.expect)
})