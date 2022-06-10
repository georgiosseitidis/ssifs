NMAdata <- function(data, X, Multiarm_studies, LA) {

  ##
  # Specifying matrices X and Z
  ##

  Z_MCMC <- data.frame(data, matrix(nrow = dim(data)[1], ncol = dim(X)[2] - 3), stringsAsFactors = F)
  colnames(Z_MCMC)[6:dim(Z_MCMC)[2]] <- colnames(X)[4:dim(X)[2]]

  for (i in 1:dim(Z_MCMC)[1]) {
    for (j in 6:dim(Z_MCMC)[2]) {
      Z_MCMC[i, j] <- X[which(X$V1 %in% c(Z_MCMC$treat1[i], Z_MCMC$treat2[i]) & X$V2 %in% c(Z_MCMC$treat1[i], Z_MCMC$treat2[i])), j - 2]
    }
  }

  ##
  # Direction of the inconsistency factors
  ##

  for (i in (dim(Z_MCMC)[2] - length(LA$comparison) + 1):dim(Z_MCMC)[2]) {
    Z_MCMC[which(Z_MCMC[, i] == 1 & paste("W", paste(Z_MCMC$treat2, Z_MCMC$treat1, sep = ",")) == colnames(Z_MCMC)[i]), i] <- -1
  }

  ##
  # Keep only the necessary comparisons
  ##

  if (length(Multiarm_studies) > 0) {

    # Comparisons where inconsistency factor is added
    if (length(LA$comparison) == 1) { # if one inconsistency factor is added in the NMA model
      ##
      IF_arm <- Z_MCMC[, (5 + length(LA$basic$V2) + 1):dim(Z_MCMC)[2]]
      ##
    } else { # more than one inconsistency factors in the NMA model
      ##
      IF_arm <- apply(Z_MCMC[, (5 + length(LA$basic$V2) + 1):dim(Z_MCMC)[2]], 1, sum)
      ##
    }
    ##
    Z_MCMC$IF_arm <- IF_arm

    # Exclude one comparison from the multi-arm studies
    for (i in Multiarm_studies) {
      ##
      multi_comp <- which(Z_MCMC$studlab == i) # Rows of the comparisons of the i multi-arm study
      multi_comp_IF <- which(Z_MCMC$studlab == i & Z_MCMC$IF_arm == 1) # Which of the multi_comp has inconsistency factor
      basic_multi <- apply(Z_MCMC[multi_comp, (5 + 1):(5 + length(LA$basic$V2))], 1, sum) # Which of the multi_comp are basic parameters


      if (length(multi_comp_IF) == 0) {
        # None of the multi-arm comparisons have inconsistency factor

        no_basic <- multi_comp[which(basic_multi != 1)] # Which comparisons are not basics

        # Exclude the last comparison
        if (length(no_basic) == 0) {
          Z_MCMC <- Z_MCMC[-multi_comp[length(multi_comp)], ]
        } else {
          Z_MCMC <- Z_MCMC[-no_basic[length(no_basic)], ]
        }

        ##
      } else if (length(multi_comp_IF) > 0) {
        # At least one of the multi-arm comparisons have inconsistency factor

        # Exclude the last comparison
        if (length(multi_comp_IF) == length(multi_comp)) {
          # All multi-arm comparisons have inconsistency factor

          Z_MCMC <- Z_MCMC[-multi_comp[length(multi_comp)], ]

          ##
        } else {
          # At least one multi-arm comparisons have inconsistency factor

          keep <- which(!multi_comp %in% multi_comp_IF)
          multi_comp <- multi_comp[keep]
          basic_multi <- basic_multi[keep]
          no_basic <- multi_comp[which(basic_multi != 1)] # Non-basic comparisons

          if (length(no_basic) == 0) {
            Z_MCMC <- Z_MCMC[-multi_comp[length(multi_comp)], ]
          } else {
            Z_MCMC <- Z_MCMC[-no_basic[length(no_basic)], ]
          }
        }
      }
    }

    # Exclude the IF_arm column
    Z_MCMC <- Z_MCMC[, -dim(Z_MCMC)[2]]
  }

  ##
  # Calculate the variance covariance matrix between random effects (H) and for the sampling errors (PREC)
  ##

  hetmatrix <- hetmat(data, Z_MCMC, Multiarm_studies)
  H <- hetmatrix$H
  PREC <- hetmatrix$PREC

  res <- list(
    "Z_MCMC" = Z_MCMC,
    "H" = H,
    "PREC" = PREC
  )

  res
}
