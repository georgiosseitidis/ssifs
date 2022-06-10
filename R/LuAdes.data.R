LuAdes.data <- function(data, ref, Multiarm_studies) {

  ##
  # Specification of comparisons where inconsistency should be added
  ##

  LA <- LuAdes(data, ref, Multiarm_studies)
  IF <- LA$comparison

  ##
  # Generate NMA MCMC data
  ##

  X <- design.matrix(LA) # Structure of the design matrix

  # Generate the data for the NMA model
  data_nma <- NMAdata(data, X, Multiarm_studies, LA)

  Z_MCMC <- data_nma$Z_MCMC
  H <- data_nma$H
  PREC <- data_nma$PREC
  Z <- as.matrix(Z_MCMC[, (dim(Z_MCMC)[2] - length(IF) + 1):dim(Z_MCMC)[2]])
  # ZTZ <- solve(t(Z) %*% Z)
  ZTZ <- (t(Z) %*% Z)
  p <- length(IF)
  y <- Z_MCMC$TE
  N <- dim(Z_MCMC)[1]
  NHtH <- sum(table(Z_MCMC$studlab) == 1)

  res <- list(
    "Z_MCMC" = Z_MCMC,
    "Z" = Z,
    "ZTZ" = ZTZ,
    "p" = p,
    "H" = H,
    "PREC" = PREC,
    "y" = y,
    "IF" = IF,
    "N" = N,
    "NHtH" = NHtH
  )

  res
}
