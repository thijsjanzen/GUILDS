#' loglikelihood following the Etienne Sampling formula
#' @title Likelihood of the Etienne sampling formula
#' @description This function calculates the likelihood of the Etienne Sampling
#' Formula, provided abundance data and parameter values.
#' @param theta Parameter value for the fundamental biodiversity number theta
#' @param m Parameter value for migration
#' @param abund Vector containing abundance data
#' @return loglikelihood
#' @references Etienne, R.S. (2005). A new sampling formula for neutral
#' biodiversity. Ecology Letters, 8(3), 253-260.
#' @author Thijs Janzen
#' @examples
#' A <- c(1, 1, 1, 3, 5, 8)   # Artificial abundance dataset
#' LL <- logLikelihood.ESF(theta = 7, m = 0.1, abund = A)
logLikelihood.ESF <- function(theta, m, abund) {
  if (theta < 1 ||
     m > (1 - .Machine$double.eps) ||
     m <= 0) {
     return(-Inf)
  }

  J <- sum(abund)
  S <- length(abund)
  I <- m * (J - 1) / (1 - m)

  kda <- calcKDA(abund)  #confirmed in PARI

  sumkda <- calc_sum_kda(S, J, I, theta, kda)

  x <- c(table(abund))
  freq_x <- c()
  for (i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor1 <- -1 * (sum(log(abund)) + sum(lgamma(1 + freq_x)))

  #J!/[prod(n1)prod(Sx!)]  #confirmed in PARI
  factor1 <- lgamma(J + 1) + prefactor1

  factor2 <- S * log(theta)   -  (lgamma(I + J) - lgamma(I))

  ll <- factor1 + factor2 + sumkda
  return(ll)
}

#' @keywords internal
esf_local <- function(v, abund, prefactor, kda) {
  theta <- v[1]
  m <- v[2]

  if (is.na(theta) ||
     is.na(m)) {
    cat("m is", m, " theta is ", theta, "one of them is NA\n")
    cat(v, "\n")
    return(-Inf)
  }

  if (theta < 1 ||
     m <= 0 ||
     m > (1 - .Machine$double.eps)) {
    return(-Inf)
  }

  J <- sum(abund)
  S <- length(abund)
  I <- m * (J - 1) / (1 - m)

  sumkda <- calc_sum_kda(S, J, I, theta, kda)

  factor2 <- S * log(theta)   -  (lgamma(I + J) - lgamma(I))

  ll <- prefactor + factor2 + sumkda
  return(ll)
}

#' maxmimum likelihood following the Etienne Sampling formula
#' @title Maximization of the loglikelihood given the standard Neutral Model,
#' using the Etienne Sampling Formula
#' @description This function computes the maximum likelihood estimates of the
#' parameters of the Neutral model, using the Etienne Sampling Formula
#' @param init_vals A vector of initial starting values, of the
#' format c(theta, m)
#' @param abund Vector containing a record of the number of individuals per
#' species
#' @param verbose TRUE/FALSE flag, indicates whether intermediate output is
#' shown on screen
#' @return the output is a list containing the following:
#' \item{par}{ a vector containing the parameter values at the maximum
#' likelihood c(theta, m)}
#' \item{fvalues}{ the likelihood at the corresponding parameter values}
#' \item{conv}{ gives a message on convergence of optimization; conv = 0
#' means convergence}
#' @references Etienne, R.S. (2005). A new sampling formula for neutral
#' biodiversity. Ecology Letters, 8(3), 253-260.
#' @author Thijs Janzen
#' examples
#' A <- c(1, 1, 1, 3, 5, 8)
#' maxLikelihood.ESF( c(7, 0.1), abund = A)
#' @export
maxLikelihood.ESF <- function(init_vals, abund, verbose = FALSE) {
  if (init_vals[1] < 1) {
     stop("maxLikelihood.ESF: ",
          "initial theta can not be below one")
   }
  if (init_vals[2] < 0) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be below zero")
  }
  if (init_vals[2] > 1) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be above 1 (did you mean to enter I?)")
  }
  if (length(abund) < 2) {
    stop("maxLikelihood.ESF: ",
         "Need more than 1 species in the dataset")
  }

  kda <- calcKDA(abund)  #confirmed in PARI

  J <- sum(abund)
  x <- c(table(abund))
  freq_x <- c()
  for (i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor <- lgamma(J + 1) - (sum(log(abund)) + sum(lgamma(1 + freq_x)))

  g <- function(x) {
	  out <- -1 * esf_local(x, abund, prefactor, kda)
	  if (verbose) cat(x, out, "\n")
		return(out)
  }

  x <- nloptr::nloptr(x0 = init_vals,
                      eval_f = g,
                      opts = list("algorithm" = "NLOPT_LN_SBPLX",
                                  xtol_rel = 1e-4))
  out <- list()
  out$par <- x$solution
  out$value <- x$objective
  return(out)
}