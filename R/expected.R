prestonSort <- function(A) {
  output <- rep(0,13);
  
  for(k in 1:length(output)) {
    start <- 2^(k-1);
    end <- -1 + 2^(k);
    if(end > length(A)) {
      end <- length(A);
      X <- sum(A[start:end]);
      output[k] <- X; 
      break;
    }
    X <- sum(A[start:end]);
    output[k] <- X;
  }
  return(output);
}


expected.SAD <- function(theta, m, J) {
  if(theta < 1) {
    stop("expected.SAD: ",
         "theta can not be below one")
  }
  if(m < 0) {
    stop("expected.SAD: ",
         "m can not be below zero")
  }
  if(m > 1) {
    stop("expected.SAD: ",
         "m can not be above 1")
  }
  if(J < 0) {
    stop("expected.SAD: ",
         "J can not be below zero")
  }
  
  I = (J-1)* m / (1-m)
  aux <- pm_sad(theta, I, J)
  SAD <- prestonSort(aux)
  return(SAD); 
}

expected.SAD.Guilds <- function(theta, alpha_x, alpha_y,
                                J, n_replicates = 100) {
  meanX <- rep(0, J)
  meanY <- rep(0, J)

  for (r in 1:n_replicates) {
		M <- drawLocal(theta, alpha_x, alpha_y, J)
		for (m in 1:length(M$guildX)) {
			meanX[m] <- meanX[m] + M$guildX[m]
		}
		for (m in 1:length(M$guildY)) {
			meanY[m] <- meanY[m] + M$guildY[m]
		}
  }
  meanX <- meanX / n_replicates
  meanY <- meanY / n_replicates

  gX <- prestonSort(meanX)
  gY <- prestonSort(meanY)
  
  output <- list( guildX = gX, guildY = gY)
  return(output);
}

expected.SAD.Guilds.Conditional <- function(theta,
                                            alpha_x,
                                            alpha_y,
                                            Jx,
                                            Jy,
                                            n_replicates = 100) {
  meanX <- rep(0,Jx)
  meanY <- rep(0,Jy)

  for (r in 1:n_replicates) {
		M <- drawLocalCond(theta, alpha_x, alpha_y, Jx, Jy)
		for (m in 1:length(M$guildX)) {
			meanX[m] <- meanX[m] + M$guildX[m]
		}
		for (m in 1:length(M$guildY)) {
			meanY[m] <- meanY[m] + M$guildY[m]
		}
  }
  meanX <- meanX / n_replicates
  meanY <- meanY / n_replicates

  gX <- prestonSort(meanX)
  gY <- prestonSort(meanY)
  
  output <- list( guildX = gX, guildY = gY)
  return(output)
}