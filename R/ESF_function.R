logLikelihood.ESF <- function(theta, m, Abund) {
  if(theta < 1 || 
     m > (1-.Machine$double.eps) ||
     m <= 0
     ) {
     return(-Inf)
  }

  J = sum(Abund)
  S = length(Abund)
  I = m * (J-1) / (1 - m)

  KDA = calcKDA(Abund)  #confirmed in PARI

  sumKDA = calcSumKDA2(S, J, I, theta, KDA)  
 
  x <- c(table(Abund))
  freq_x <- c()
  for(i in 1:length(x)) freq_x[i] <- x[[i]]
  prefactor1 = -( sum(log(Abund)) + sum(lgamma(1+freq_x)) )
  
  factor1 <- lgamma(J+1) + prefactor1  #J!/[prod(n1)prod(Sx!)]  #confirmed in PARI
  
  factor2 = S*log(theta)   -  (lgamma(I + J) - lgamma(I))
  
  LogLikelihood <- factor1 + factor2 + sumKDA
  return(LogLikelihood)
}

ESF_local <- function(v, Abund, prefactor, KDA) {
  theta = v[1]
  m = v[2]
  if(theta < 1 || 
     m <= 0 || 
     m > (1-.Machine$double.eps)) {
    return(-Inf)
  }
 
  J = sum(Abund)
  S = length(Abund)
  I = m * (J-1) / (1 - m)

  sumKDA = calcSumKDA2(S, J, I, theta, KDA);  
 
  factor2 = S*log(theta)   -  (lgamma(I + J) - lgamma(I));
  
  LogLikelihood <- prefactor + factor2 + sumKDA;
  return(LogLikelihood);
}

maxLikelihood.ESF <- function(initVals, Abund, verbose=TRUE) {
  if (initVals[1] < 1) {
     stop("maxLikelihood.ESF: ",
          "initial theta can not be below one")
   }
  if (initVals[2] < 0) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be below zero")
  }
  if (initVals[2] > 1) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be above 1 (did you mean to enter I?)")
  }
  if (length(Abund) < 2) {
    stop("maxLikelihood.ESF: ",
         "Need more than 1 species in the dataset")
  }

  KDA = calcKDA(Abund)  #confirmed in PARI

  J = sum( Abund)
  S = length( Abund)
  x <- c( table( Abund))
  freq_x <- c()
  for (i in 1:length(x)) freq_x[i] <- x[[i]]
  prefactor = lgamma(J+1) -( sum(log(Abund)) + sum(lgamma(1+freq_x)) )
   
  g <- function(x) {
	  out <- -1 * ESF_local(x, Abund, prefactor, KDA)
		return(out)
  }  
  
  optimum <- simplex(initVals, g, verbose)
  return(optimum)
}
