#' @title Likelihood of the Guilds sampling formula, conditional on guild size
#' @description This function calculates the likelihood of the guilds model,
#' conditional on guild size; provided abundance data and parameter values.
#' @inheritParams default_params_doc
#' @return loglikelihood
#' @author Thijs Janzen
#' @examples
#' exampleData <- generate.Guilds.Cond(theta = 200,
#'                                     alpha_x = 0.005,
#'                                     alpha_y = 0.001,
#'                                     JX = 1000,
#'                                     JY = 2000)
#'            #theta = 200, alpha X = 0.005, alpha Y = 0.001
#'  parametervals <- c(200, 0.005, 0.001)
#'  LL = logLikelihood.Guilds.Conditional(parametervals,
#'                                        model = "D1",
#'                                        exampleData$guildX,
#'                                        exampleData$guildY,
#'                                        verbose = TRUE)
#'  @export
logLikelihood.Guilds.Conditional <- function(parameters, model,
                                             sadx, sady,
                                             verbose = TRUE) {
  Nx <- sum(sadx)
  Ny <- sum(sady)

  if (verbose == TRUE) {
    cat("Chosen model: ", model, "\n")
    cat("Now starting to calculate likelihood of: \n")
    x2 <- parameters
    if (model == "D0") cat("Theta X =", x2[1],
                           " Theta Y =", "Theta X",
                           "\t Alpha X =", x2[2],
                           " Alpha Y =", "Alpha X", "\n")
    if (model == "D1") cat("Theta X =", x2[1],
                           " Theta Y =", "Theta X",
                           "\t Alpha X =", x2[2],
                           " Alpha Y =", x2[3], "\n")

    flush.console()
  }

  ll <- logLikelihood.Guilds(parameters, model, sadx, sady, verbose)

  #conditional part:
  conditional_part <- calc_conditional(parameters, model, Nx, Ny)

  output <- ll - conditional_part
  return(output)
}