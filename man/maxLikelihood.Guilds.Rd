\name{maxLikelihood.Guilds}
\alias{maxLikelihood.Guilds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Maximization of the loglikelihood under the Guilds Model.
}
\description{
 This function computes the maximum likelihood estimates of the parameters of the guilds model.
}
\usage{
maxLikelihood.Guilds(init_vals, model = "D0",
                     sadx, sady, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
\dontrun{
  J <- 10000

  theta <- 100
  alpha_x <- 0.1

  simul_data <- generate.Guilds(theta, alpha_x, alpha_x, J)

  #initial parameters for the D0 model c(theta,alpha)
  LL <- maxLikelihood.Guilds(init_vals = c(theta, alpha_x),
                             model = "D0",
                             sadx  = simul_data$guildX,
                             sady  = simul_data$guildY)
}
}