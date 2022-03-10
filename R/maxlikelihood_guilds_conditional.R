maxLikelihood.Guilds.Conditional <- function(init_vals, model,
                                             sadx, sady,
                                             verbose = TRUE) {
  incorrectlength <- 0
  if (model == "D0" && length(init_vals) != 2) incorrectlength <- 1
  if (model == "D1" && length(init_vals) != 3) incorrectlength <- 1

  if (incorrectlength == 1) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "Input vector is of incorrect length\n")
  }

  if (init_vals[1] < 1) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial theta can not be below one")
  }
  if (init_vals[2] < 0) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial alpha can not be below zero")
  }
  if (init_vals[2] > 1) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial alpha can not be above 1")
  }
  if (model == "D1") {
    if (init_vals[3] < 0) {
      stop("maxLikelihood.Guilds.Conditional: ",
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
  for (i in seq_along(x)) freq_x[i] <- x[[i]]

  prefactor1 <- -1 * (sum(log(sadx)) + sum(lgamma(1 + freq_x)))

  x2 <- c(table(sady))
  freq_y <- c()
  for (i in seq_along(x2)) freq_y[i] <- x2[[i]]

  prefactor2 <- -1 * (sum(log(sady)) + sum(lgamma(1 + freq_y)))

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J  <- Nx + Ny

  g <- function(x) {
    out <- -1 * conditional.LogLik(x, model, J, Sx, Sy, Nx, Ny,
                                   kda_x, kda_y, prefactor1,
                                   prefactor2, verbose)
    return(out)
  }

  x <- nloptr::nloptr(x0 = init_vals,
                     eval_f = g,
                     opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                 xtol_rel = 1e-4))
  out <- list()
  out$par <- x$solution
  out$value <- x$objective
  return(out)
}