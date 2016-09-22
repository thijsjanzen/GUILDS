sort_aux <- function(A) {
  output <- rep(0,13)

  for (k in seq_along(output)) {
    start <- 2^(k - 1)
    end <- -1 + 2^(k)
    if(end > length(A)) {
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

expected.SAD.Guilds <- function(theta, alpha_x, alpha_y,
                                J, n_replicates = 100) {
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
  
  output <- list( guildX = gx, guildY = gy)
  return(output)
}

expected.SAD.Guilds.Conditional <- function(theta,
                                            alpha_x,
                                            alpha_y,
                                            Jx,
                                            Jy,
                                            n_replicates = 100) {
  meanx <- rep(0,Jx)
  meany <- rep(0,Jy)

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
  
  output <- list( guildX = gx, guildY = gy)
  return(output)
}