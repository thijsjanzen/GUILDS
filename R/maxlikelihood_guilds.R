maxLikelihood.Guilds <- function(init_vals,
                                 model = "D0",
                                 sadx,
                                 sady,
                                 verbose = FALSE) {
  incorrectlength <- 0
  if (model == "D0" && length(init_vals) != 2) incorrectlength <- 1
  if (model == "D1" && length(init_vals) != 3) incorrectlength <- 1

  if (incorrectlength == 1) {
    stop("maxLikelihood.Guilds: ",
         "Input vector is of incorrect length\n")
  }

  if (init_vals[1] < 1) {
    stop("maxLikelihood.Guilds: ",
         "initial theta can not be below one")
  }
  if (init_vals[2] < 0) {
    stop("maxLikelihood.Guilds: ",
         "initial alpha can not be below zero")
  }
  if (init_vals[2] > 1) {
    stop("maxLikelihood.Guilds: ",
         "initial alpha can not be above 1")
  }
  if (model == "D1") {
    if (init_vals[3] < 0) {
      stop("maxLikelihood.Guilds: ",
           "initial alpha_y can not be below 0")
    }
    if (init_vals[3] > 1) {
      stop("maxLikelihood.Guilds: ",
           "initial alpha_y can not be above 1")
    }
  }

  kda_x <- calcKDA(sadx)
  kda_y <- calcKDA(sady)

  x <- c(table(sadx))
  freq_x <- c()
  for (i in seq_along(x) ) freq_x[i] <- x[[i]]

  prefactor1 <- -(sum(log(sadx)) + sum(lgamma(1 + freq_x)) )

  x2 <- c(table(sady))
  freq_y <- c()
  for (i in seq_along(x2)) freq_y[i] <- x2[[i]]
  prefactor2 <- -(sum(log(sady)) + sum(lgamma(1 + freq_y)) )

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J <- Nx + Ny

  g <- function(v) {
    theta_x <- v[1] * 2
    theta_y <- v[1] * 2
    alpha_x <- v[2]
    alpha_y <- v[2]

    if (model == "D1") {
      alpha_y <- v[3]
    }

    y <- 0
    if (alpha_x < 0 ||
        alpha_y < 0 ||
        theta_x < 1 ||
        theta_y < 1  ||
        alpha_x > (1 - (1e-8)) ||
        alpha_y > (1 - (1e-8))
    ) {
      y <- -Inf
    }
    if (!is.infinite(y)) {
      y <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y,
                        J, Sx, Sy, Nx, Ny, kda_x, kda_y,
                        prefactor1, prefactor2, verbose)
    }
    if (verbose) {
      cat(theta_x, theta_y, alpha_x, alpha_y, y, "\n")
    }
    return(-y)
  }

  x <- subplex::subplex(init_vals, g)
  return(x)
}