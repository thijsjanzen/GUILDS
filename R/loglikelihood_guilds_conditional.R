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