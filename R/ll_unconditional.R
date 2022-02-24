logLikguilds <- function(theta_x, theta_y,
                         alpha_x, alpha_y,
                         J, Sx, Sy, Nx, Ny,
                         KDA_X, KDA_Y,
                         prefactor1, prefactor2,
                         verbose = TRUE) {
  thrs <- 10

  f <- function(x) {
    return(-1 * evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                               J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y))
  }

  maxes <- pracma::fminbnd(f, 0, 1,
                           maxiter = 500, tol = 1e-4)

  ymax <- -1 * maxes$fmin
  xmax <- maxes$xmin
  xlft <- 0
  xrgt <- 1
  eps <- .Machine$double.eps

  check_left_x  <- evaluateLogLik(eps, theta_x, theta_y, alpha_x, alpha_y,
                                  J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs
  check_right_x <- evaluateLogLik(1 - eps, theta_x, theta_y, alpha_x, alpha_y,
                                  J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs

  if (check_left_x < 0) {
    g <- function(x) {
      return(evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                            J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs)
    }
    xlft <- (uniroot(g, c(eps, xmax)))$root
  }

  if (check_right_x < 0) {
    h <- function(x) {
      return(evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                            J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs)
    }
    xrgt <- (uniroot(h, c(xmax, 1 - eps)))$root
  }

  calc_ll_exp <- function(x) {
    out <- exp( evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                               J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax)
    return(out)
  }

  aux <- integrate(f = calc_ll_exp, lower = xlft, upper = xrgt, abs.tol = 1e-9)

  y <- ymax + log(aux$value) + prefactor1 + Sx * log(theta_x / 2) +
    prefactor2 + Sy * log(theta_y / 2) + lgamma((theta_x / 2) + (theta_y / 2)) -
    (lgamma(theta_x / 2) + lgamma(theta_y / 2))

  return(y)
}