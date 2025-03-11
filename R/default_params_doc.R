#' Default parameter documentation
#' @param parameters \code{parameters} corresponds to a vector of parameter
#' values depending on the provided model: \cr
#' - model: 'D0' \code{parameters} = c(theta, alpha) \cr
#' - model: "D1" \code{parameters} = c(theta, alpha X, alpha Y) \cr
#' @param init_vals \code{init_vals} corresponds to a vector of parameter values
#' in which to start the Maxmimum Likelihood algorithm, depending on the provided
#' model: \cr
#' - model: "D0" \code{parameters} = c(theta, alpha) \cr
#' - model: "D1" \code{parameters} = c(theta, alpha X, alpha Y) \cr
#' @param model The chosen model to calculate the likelihood for,
#' please note that the vector of parameters should contain the corresponding
#' parameters in the right order. The user can pick one of these models:\cr
#' -  "D0" \cr
#' -  "D1" \cr
#' @param sadx The Species Abundance Distribution of guild X
#' @param sady The Species Abundance Distribution of guild Y
#' @param verbose TRUE/FALSE flag, indicates whether intermediate output is
#' shown on screen
#' @param theta Fundamental biodiversity number theta
#' @param alpha_x Dispersal Ability of Guild X
#' @param alpha_y Dispersal Ability of Guild Y
#' @rawNamespace useDynLib(GUILDS)
#' @rawNamespace import(Rcpp)
#' @keywords internal
default_params_doc <- function(parameters,
                               init_vals,
                               model,
                               sadx,
                               sady,
                               verbose,
                               theta,
                               alpha_x,
                               alpha_y) {
  # nothing
}