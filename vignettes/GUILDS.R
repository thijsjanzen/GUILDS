## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)

## ----start--------------------------------------------------------------------
library(GUILDS)
set.seed(42)
theta = 100
m <- 0.1
J <- 10000
I <- m * (J - 1) / (1 - m)
abund <- generate.ESF(theta, I, J)
abund

## -----------------------------------------------------------------------------
preston_plot(abund)

## -----------------------------------------------------------------------------
abund.expect <- expected.SAD(theta, m, J)
preston_plot(abund, abund.expect)

## ----maxlik_ESF---------------------------------------------------------------
LL1 <- maxLikelihood.ESF(init_vals = c(theta, m), abund)
LL2 <- maxLikelihood.ESF(init_vals = c(100, 0.01), abund)
LL3 <- maxLikelihood.ESF(init_vals = c(10, 0.5), abund)
LL1$par
LL2$par
LL3$par

## -----------------------------------------------------------------------------
set.seed(666)
theta <- 200
alpha_x <- 0.1
J <- 10000
J_x <- 8000
J_y <- J - J_x

abund <- generate.Guilds.Cond(theta, alpha_x, alpha_x, J_x, J_y)

## -----------------------------------------------------------------------------
abund.expected <- expected.SAD.Guilds.Conditional(theta, alpha_x, alpha_x,
                                                  J_x, J_y, n_replicates = 5)

par(mfrow = c(1, 2))
preston_plot(abund$guildX, abund.expected$guildX, main = "Guild X")
preston_plot(abund$guildY, abund.expected$guildY, main = "Guild Y")

## ----maxlik_cond_first--------------------------------------------------------
LL <- maxLikelihood.Guilds.Conditional(init_vals = c(theta, alpha_x), 
                                       model = "D0",
                                       sadx = abund$guildX, 
                                       sady = abund$guildY,
                                       verbose = FALSE)
LL$par

## ----maxlik_ESF_single--------------------------------------------------------
LL1 <- maxLikelihood.ESF(init_vals = c(theta, alpha_x),
                         abund = c(abund$guildX, abund$guildY), 
                         verbose = FALSE)
LL1$par

## -----------------------------------------------------------------------------
set.seed(666 + 42)
theta <- 200
alpha_x <- 0.01
alpha_y <- 0.1
J <-  1000
J_x <- 800
J_y <- J - J_x

abund <- generate.Guilds.Cond(theta, alpha_x, alpha_y, J_x, J_y)

## ----expected_SAD_cond--------------------------------------------------------
abund.expected <- expected.SAD.Guilds.Conditional(theta, alpha_x, alpha_y,
                                                  J_x, J_y, n_replicates = 5)

par(mfrow = c(1, 2))
preston_plot(abund$guildX, abund.expected$guildX, main = "Guild X")
preston_plot(abund$guildY, abund.expected$guildY, main = "Guild Y")

## ----maxlik_cond_D1-----------------------------------------------------------
ML1 <- maxLikelihood.Guilds.Conditional(init_vals = c(theta, alpha_x, alpha_y), 
                                       model = "D1",
                                       sadx = abund$guildX, sady = abund$guildY,
                                       verbose = FALSE)

## ----maxlik_cond_D0-----------------------------------------------------------
ML2 <- maxLikelihood.Guilds.Conditional(init_vals = c(theta, alpha_x), 
                                       model = "D0",
                                       sadx = abund$guildX, sady = abund$guildY,
                                       verbose = FALSE)

## ----AIC----------------------------------------------------------------------
AIC_D1 = 2 * 3 - 2 * -ML1$value
AIC_D0 = 2 * 2 - 2 * -ML2$value
AIC_D1
AIC_D0

