DBT.data <- function(data, ref, Multiarm_studies, method) {

  ##
  # Find networks designs
  ##

  nmades <- NMAdesigns(data, Multiarm_studies)
  ##
  design <- nmades$designs # Network's designs
  design_multi <- nmades$multi_design # Designs from multi-arm studies
  studies_multi <- nmades$studies_multi # Multi-arm studies for each multi-arm design

  ##
  # Calculate the number of inconsistency factors
  ##

  n_nodes <- length(unique(c(data$treat1, data$treat2))) # Number of network's nodes

  # List with the number of inconsistency factors for the Design-by-treatment methods
  p_DBT <- numberIF.DBT(n_nodes, Multiarm_studies, design, design_multi, method)
  ##
  p_Jackson <- p_DBT$Jackson # Number of inconsistency factor for the Jackson approach
  p_Higgins <- p_DBT$Higgins # Number of inconsistency factor for the Higgins approach
  IFmulti <- p_DBT$IFmulti # Number of inconsistency factors added for each multi-arm design

  ##
  # Generate NMA MCMC data
  ##

  p <- ifelse(method == "DBT", p_Higgins, p_Jackson)

  # List with the NMA data and the treatment's effect
  zy <- dbtdesign(data, ref, Multiarm_studies, studies_multi, design, IFmulti, p_Jackson, method, p_Higgins)
  ##
  Z_MCMC <- zy$Z
  y <- zy$y
  ##
  Z <- as.matrix(Z_MCMC[, 4:dim(Z_MCMC)[2]])
  colnames(Z) <- names(Z_MCMC)[4:dim(Z_MCMC)[2]]
  ##
  # ZTZ <- solve(t(Z) %*% Z)
  ZTZ <- t(Z) %*% Z
  IF <- colnames(Z)
  N <- dim(Z_MCMC)[1]
  NHtH <- sum(table(Z_MCMC$studlab) == 1)

  # Specify covariance matrices
  cov_matrices <- hetmat(data, Z_MCMC, Multiarm_studies)
  ##
  H <- cov_matrices$H
  PREC <- cov_matrices$PREC

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
