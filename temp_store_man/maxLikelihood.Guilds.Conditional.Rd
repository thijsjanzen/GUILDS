\name{maxLikelihood.Guilds.Conditional}
\alias{maxLikelihood.Guilds.Conditional}

\title{
Maximization of the loglikelihood under the Guilds Model, conditioned on guild size.
}
\description{
 This function computes the maximum likelihood estimates of the parameters of the guilds model, conditioned on guild size.
}
\usage{
maxLikelihood.Guilds.Conditional(init_vals, model, sadx, sady, verbose = TRUE)
}

\arguments{
  \item{init_vals}{\code{init_vals} corresponds to a vector of parameter values in which to start the Maxmimum Likelihood algorithm, depending on the provided model: \cr
		- model: "D0" \code{parameters} = c(theta, alpha) \cr
		- model: "D1" \code{parameters} = c(theta, alpha X, alpha Y) \cr
	}

 \item{model}{ The chosen model to calculate the maximum likelihood for, please note that the vector of parameters should contain the corresponding parameters in the right order. The user can pick one of these models:\cr
		-  "D0" \cr
		-  "D1" \cr
	}

  \item{sadx}{ The Species Abundance Distribution of guild X}
  \item{sady}{ The Species Abundance Distribution of guild Y}
  \item{verbose}{TRUE/FALSE flag, indicates whether intermediate output is shown on screen}
}

\value{
  The output is a list containing the following:
   \item{par}{a vector containing the parameter values at the maximum likelihood}
   \item{value}{the likelihood at the corresponding parameter values}
   \item{counts}{Number of function evaluations required}
   \item{convergence}{
			-2:  invalid input\cr
			-1:  number of maximum function evaluations exceeded \cr
			 0:  success: convergence \cr
			 1:  limit of machine precision reached \cr
	}
   \item{message}{A character string giving a diagnostic message from the optimizer,}
   \item{hessian}{Hessian matrix (not implemented for this package)}
}

\author{
Thijs Janzen
}

\examples{
theta = 20
alpha = 0.1
initParams <- c(theta, alpha)
maxLikelihood.Guilds.Conditional(initParams,
                                 model = "D0",
                                 sadx  = 1:20,
                                 sady  = 1:20,
                                 verbose = TRUE)
}