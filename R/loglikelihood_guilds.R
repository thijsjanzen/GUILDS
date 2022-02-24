logLikelihood.Guilds <- function(parameters, model,
                                 sadx, sady,
                                 verbose = TRUE) {
  kda_x <- calcKDA(sadx)
  kda_y <- calcKDA(sady)

  x <- c(table(sadx))
  freq_x <- c()
  for (i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor1 <- -1 * (sum(log(sadx)) + sum(lgamma(1 + freq_x)) )

  x2 <- c(table(sady))
  freq_y <- c()
  for (i in seq_along(x2)) freq_y[i] <- x2[[i]]
  prefactor2 <- -1 * (sum(log(sady)) + sum(lgamma(1 + freq_y)) )

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J <- Nx + Ny


  if (model == "D0") {
    #because theta_x = theta_y = theta/2
    theta_x <- parameters[1] * 2
    theta_y <- parameters[1] * 2
    alpha_x <- parameters[2]
    alpha_y <- parameters[2]
  }
  if (model == "D1") {
    theta_x <- parameters[1] * 2
    theta_y <- parameters[1] * 2
    alpha_x <- parameters[2]
    alpha_y <- parameters[3]
  }

  ll <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y, J,
                     Sx, Sy, Nx, Ny, kda_x, kda_y, prefactor1,
                     prefactor2, verbose)

  if ( verbose == TRUE ) {
    cat("Likelihood is ", ll, "\n")
  }

  return(ll)
}