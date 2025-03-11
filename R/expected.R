#' @keywords internal
sort_aux <- function(A) {
  output <- rep(0, 13)

  for (k in seq_along(output)) {
    start <- 2 ^ (k - 1)
    end <- -1 + 2 ^ (k)
    if (end > length(A)) {
      end <- length(A)
      X <- sum(A[start:end])
      output[k] <- X
      break
    }
    X <- sum(A[start:end])
    output[k] <- X
  }
  return(output)
}

#' expected distribution
#' @title Calculate the expected species abundance distribution of the standard
#' neutral model, given theta, m and J
#' @description This function calculates the expected species abundance
#' distribution of the standard neutral model given theta, m and J, sensu
#' equation 6 from Etienne and Alonso (2005).
#' @param theta Fundamental biodiversity number theta
#' @param m migration parameter
#' @param J Total number of individuals in the local community
#' @return A vector containing the abundances binned into log2 bins
#' (sensu Preston)
#' @references Etienne, R.S., & Alonso, D. (2005). A dispersal-limited sampling
#' theory for species and alleles. Ecology Letters, 8(100), 1147-1156.
#' @author Thijs Janzen & Bart Haegeman
#' @export
#' @examples
#' SAD <- expected.SAD(theta = 42, m = 0.1, J = 200)
#' barplot(SAD,
#'         names.arg = 0:(length(SAD) - 1),
#'         xlab = "Number of individuals (log2)",
#'         ylab = "Number of Species" )
expected.SAD <- function(theta, m, J) {
  if (theta < 1) {
    stop("expected.SAD: ",
         "theta can not be below one")
  }
  if (m < 0) {
    stop("expected.SAD: ",
         "m can not be below zero")
  }
  if (m > 1) {
    stop("expected.SAD: ",
         "m can not be above 1")
  }
  if (J < 0) {
    stop("expected.SAD: ",
         "J can not be below zero")
  }

  I <- (J - 1) * m / (1 - m)
  aux <- pm_sad(theta, I, J)
  sad <- sort_aux(aux)
  return(sad)
}


#' expected distribution under the guilds model
#' @title Estimate the expected species abundance distribution of both guilds
#' using the guilds model, provided theta, alpha_x, alpha_y and J.
#' @description This function estimates the expected species abundance
#' distribution of both guilds using the guilds model, provided theta, alpha_x,
#' alpha_y and J. The expected species abundance distribution is approximated
#' by first drawing px from a beta distribution (equation 4 in
#' Janzen et al. 2015). Then, guild sizes are drawn using equation 3 in
#' Janzen et al. 2015. Because the abundance distributions of the two guilds
#' are independent, the distributions can now be obtained using equation 6 in
#' Etienne and Alonso 2005. Because drawing from the beta distribution and
#' equation 3 is inherently stochastic, this function returns the average
#' over a specified number of replicates.
#' @param theta Fundamental biodiversity number theta
#' @param alpha_x Dispersal ability of guild X
#' @param alpha_y Dispersal ability of guild Y
#' @param J Total number of individuals in the local community, e.g. J = Jx + Jy
#' @param n_replicates Number of replicates to use to estimate the abundance
#' distributions.
#' @returns \item{guildX}{Vector containing the mean abundances of species in
#' Guild X, binned into log2 bins}
#' \item{guildY}{Vector containing the mean abundances of species in Guild Y,
#' inned into log2 bins}
#' @export
#' @references Etienne, R.S., & Alonso, D. (2005). A dispersal-limited sampling
#' theory for species and alleles. Ecology Letters, 8(100), 1147-1156.
#' Janzen, T., Haegeman B., Etienne, R.S. (2015) A sampling formula for
#' communities with multiple dispersal syndromes. Journal of Theoretical
#' Biology 374: 94-106
#' @author Thijs Janzen and Bart Haegeman
#' @examples
#'  SADs <- expected.SAD.Guilds(theta = 42,
#'                              alpha_x = 0.01,
#'                              alpha_y = 0.1,
#'                              J = 1000,
#'                              n_replicates = 3)
#'  par(mfrow=c(1,2));
#'  barplot(SADs$guildX, names.arg = 0:(length(SADs$guildX) - 1),
#'          xlab = "Number of individuals (log2)",
#'          ylab = "Number of Species", main = "Guild X" )
#'  barplot(SADs$guildY, names.arg = 0:(length(SADs$guildY) - 1),
#'          xlab = "Number of individuals (log2)",
#'          ylab = "Number of Species", main = "Guild Y" )
expected.SAD.Guilds <- function(theta, alpha_x, alpha_y,
                                J, n_replicates = 100) {
  if (theta < 1) {
    stop("expected.SAD.Guilds: ",
         "theta can not be below one")
  }
  if (alpha_x < 0) {
    stop("expected.SAD.Guilds: ",
         "alpha_x can not be below zero")
  }
  if (alpha_x > 1) {
    stop("expected.SAD.Guilds: ",
         "alpha_x can not be above 1")
  }
  if (alpha_y < 0) {
    stop("expected.SAD.Guilds: ",
         "alpha_y can not be below zero")
  }
  if (alpha_y > 1) {
    stop("expected.SAD.Guilds: ",
         "alpha_y can not be above 1")
  }
  if (J < 1) {
    stop("expected.SAD: ",
         "J can not be below one")
  }

  meanx <- rep(0, J)
  meany <- rep(0, J)

  for (r in seq_len(n_replicates)) {
    M <- draw_local(theta, alpha_x, alpha_y, J)
    for (m in seq_along(M$guildX)) {
      meanx[m] <- meanx[m] + M$guildX[m]
    }
    for (m in seq_along(M$guildY)) {
      meany[m] <- meany[m] + M$guildY[m]
    }
  }
  meanx <- meanx / n_replicates
  meany <- meany / n_replicates

  gx <- sort_aux(meanx)
  gy <- sort_aux(meany)

  output <- list(guildX = gx, guildY = gy)
  return(output)
}

#' expected distribution under the guilds model, conditional on guild sizes
#' @title EEstimate the expected species abundance distribution of both guilds
#' using the guilds model, provided theta, alpha_x, alpha_y, conditional on the
#' size of guild X, Jx and the size of guild Y, Jy.
#' @description This function estimates the expected species abundance
#' distribution of both guilds using the guilds model, provided theta, alpha_x,
#' alpha_y and J. The expected species abundance distribution is approximated by
#' first drawing px from equation 9. Because the abundance distributions of the
#' two guilds are independent, the distributions can now be obtained using
#' equation 6 in Etienne and Alonso 2005. Because drawing from the beta
#' distribution and equation 3 is inherently stochastic, this function returns
#' the average over a specified number of replicates.
#' @param theta Fundamental biodiversity number theta
#' @param alpha_x Dispersal ability of guild X
#' @param alpha_y Dispersal ability of guild Y
#' @param Jx Total number of individuals in guild X
#' @param Jy Total number of individuals in guild Y
#' @param n_replicates Number of replicates to use to estimate the abundance
#' distributions.
#' @returns \item{guildX}{Vector containing the mean abundances of species in
#' Guild X, binned into log2 bins}
#' \item{guildY}{Vector containing the mean abundances of species in Guild Y,
#' binned into log2 bins}
#' @references Etienne, R.S., & Alonso, D. (2005). A dispersal-limited sampling
#' theory for species and alleles. Ecology Letters, 8(100), 1147-1156.
#' @export
#' @author Thijs Janzen and Bart Haegeman
#' @examples
#'  SADs <- expected.SAD.Guilds.Conditional(theta = 42,
#'                                          alpha_x = 0.01,
#'                                          alpha_y = 0.1,
#'                                          Jx = 100,
#'                                          Jy = 200,
#'                                          n_replicates = 3)
#' par(mfrow=c(1,2))
#' barplot(SADs$guildX, names.arg = 0:(length(SADs$guildX) - 1),
#'         xlab = "Number of individuals (log2)",
#'         ylab = "Number of Species", main = "Guild X" )
#' barplot(SADs$guildY, names.arg = 0:(length(SADs$guildY) - 1),
#'         xlab = "Number of individuals (log2)",
#'         ylab = "Number of Species", main = "Guild Y" )
#' expected.SAD.Guilds.Conditional <- function(theta,
#'                                             alpha_x,
#'                                             alpha_y,
#'                                             Jx,
#'                                             Jy,
#'                                             n_replicates = 100)
expected.SAD.Guilds.Conditional <- function(theta,
                                            alpha_x,
                                            alpha_y,
                                            Jx,
                                            Jy,
                                            n_replicates = 100) {
  if (theta < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "theta can not be below one")
  }
  if (alpha_x < 0) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_x can not be below zero")
  }
  if (alpha_x > 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_x can not be above 1")
  }
  if (alpha_y < 0) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_y can not be below zero")
  }
  if (alpha_y > 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_y can not be above 1")
  }
  if (Jx < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "Jx can not be below one")
  }
  if (Jy < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "Jy can not be below one")
  }

  meanx <- rep(0, Jx)
  meany <- rep(0, Jy)

  for (r in seq_len(n_replicates)) {
    M <- draw_local_cond(theta, alpha_x, alpha_y, Jx, Jy)
    for (m in seq_along(M$guildX)) {
      meanx[m] <- meanx[m] + M$guildX[m]
    }
    for (m in seq_along(M$guildY)) {
      meany[m] <- meany[m] + M$guildY[m]
    }
  }
  meanx <- meanx / n_replicates
  meany <- meany / n_replicates

  gx <- sort_aux(meanx)
  gy <- sort_aux(meany)

  output <- list(guildX = gx, guildY = gy)
  return(output)
}