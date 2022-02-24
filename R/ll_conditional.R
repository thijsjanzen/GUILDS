evaluate_cond_lik <- function(v, theta_x, theta_y, alpha_x, alpha_y, Nx, Ny) {
  nx <- v
  ny <- 1 - nx
  J <- Nx + Ny
  I_X <- alpha_x * nx * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)
  I_Y <- alpha_y * ny * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)

  a <- lgamma(J + 1)
  b <- rep(0, length(I_X))
  poch_X <- rep(0, length(I_X))
  poch_Y <- rep(0, length(I_X))

  for (cnt in seq_along(I_X)) {
    b[cnt] <- lgamma(I_X[cnt] + I_Y[cnt] + J) - lgamma(I_X[cnt] + I_Y[cnt])
    poch_X[cnt] <- lgamma(I_X[cnt] + Nx) - lgamma(I_X[cnt])
    poch_Y[cnt] <- lgamma(I_Y[cnt] + Ny) - lgamma(I_Y[cnt])
  }

  c <- a - b

  h <- poch_X + poch_Y - (lgamma(Nx + 1) + lgamma(Ny + 1))

  k  <-  lgamma((theta_x / 2) + (theta_y / 2)) -
        (lgamma(theta_x / 2) + lgamma(theta_y / 2))
  l  <- ((theta_x / 2) - 1) * log(nx) + ((theta_y / 2) - 1) * log(ny)

  result <- c + h + k + l
  return(result)
}

calc_conditional <- function(v, model, Nx, Ny) {
  incorrectlength <- FALSE
  if (model == "D0" && length(v) != 2) incorrectlength <- TRUE
  if (model == "D1" && length(v) != 3) incorrectlength <- TRUE

  if (incorrectlength == TRUE) {
    stop("calcConditional:",
         "Input vector is of incorrect length\n")
  }

  theta_x <- v[1] * 2
  theta_y <- v[1] * 2
  alpha_x <- v[2]
  alpha_y <- v[2]

  if (model == "D1") {
    alpha_y <- v[3]
  }

  if (alpha_x < 0 ||
      alpha_y < 0 ||
      theta_x < 1 ||
      theta_y < 1 ||
      alpha_x > (1 - (1e-8)) ||
      alpha_y > (1 - (1e-8))
     ) {
    return(-Inf)
  }

  f <- function(x) {
    return(-1 * evaluate_cond_lik(x, theta_x, theta_y,
                                 alpha_x, alpha_y,
                                 Nx, Ny))
  }

 maxes <- pracma::fminbnd(f, 0, 1, maxiter = 500, tol = 1e-4)

 ymax <- -1 * maxes$fmin
 xmax <- maxes$xmin
 xlft <- 0
 xrgt <- 1
 eps <- .Machine$double.eps
 thrs <- 10

 check_left_x  <- evaluate_cond_lik(eps, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny)  - ymax + thrs
 check_right_x <- evaluate_cond_lik(1 - eps, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny)  - ymax + thrs

 if (check_left_x < 0) {
     g <- function(x) {
        return(evaluate_cond_lik(x, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny) - ymax + thrs)
     }
     xlft <- (uniroot(g, c(eps, xmax)))$root
 }

 if (check_right_x < 0) {
     h <- function(x) {
        return(evaluate_cond_lik(x, theta_x, theta_y,
                               alpha_x, alpha_y,
                               Nx, Ny) - ymax + thrs)
     }
     xrgt <- (uniroot(h, c(xmax, 1 - eps)))$root
 }

 calc_ll_exp <- function(x) {
   out <- exp(evaluate_cond_lik(x, theta_x, theta_y,
                              alpha_x, alpha_y,
                              Nx, Ny) - ymax)
   return(out)
 }

  aux <- integrate(f = calc_ll_exp,
                   lower = xlft,
                   upper = xrgt,
                   abs.tol = 1e-9)

  LL <- log(aux$value) + ymax

  return(LL)
}

conditional.LogLik <- function(v, model, J, Sx, Sy, Nx, Ny, kda_x, kda_y,
                               prefactor1, prefactor2, verbose = TRUE) {

  theta_x <- v[1] * 2
  theta_y <- v[1] * 2
  alpha_x <- v[2]
  alpha_y <- v[2]

  if (model == "D1") {
    alpha_y <- v[3]
  }

  if (alpha_x < 0 ||
     alpha_y < 0 ||
     theta_x < 1 ||
     theta_y < 1 ||
     alpha_x > (1 - (1e-8)) ||
     alpha_y > (1 - (1e-8))
     ) return(-Inf)

  if (is.na(alpha_x) ||
     is.na(alpha_y) ||
     is.na(theta_x) ||
     is.na(theta_y)) {
    cat("warnings! one of the parameters is somehow NA\n")
    cat("displaying alpha_x, alpha_y, theta_x, theta_y\n")
    cat(alpha_x, alpha_y, theta_x, theta_y, "\n")
    return(-Inf)
  }


  ll <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y, J,
                     Sx, Sy, Nx, Ny, kda_x, kda_y,
                     prefactor1, prefactor2, verbose)

  cond_ll <- calc_conditional(v, model, Nx, Ny)
  out <- ll - cond_ll
  return(out)
}