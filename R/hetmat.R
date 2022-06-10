hetmat <- function(data, Z_MCMC, Multiarm_studies) {
  if (!length(Multiarm_studies) > 0) {
    H <- PREC <- NULL
  } else {
    n_multi_comp <- length(Z_MCMC$studlab[which(Z_MCMC$studlab %in% Multiarm_studies)]) # Number of multi-arm comparisons

    # Matrix H
    H <- diag(x = 1, nrow = n_multi_comp, ncol = n_multi_comp)
    H[which(H != 1)] <- 0.5
    colnames(H) <- row.names(H) <- Z_MCMC$studlab[which(Z_MCMC$studlab %in% Multiarm_studies)]

    # Matrix PREC
    PREC <- matrix(ncol = n_multi_comp, nrow = n_multi_comp)
    ##
    colnames(PREC) <- row.names(PREC) <- paste(
      "Comp", Z_MCMC$treat1[which(Z_MCMC$studlab %in% Multiarm_studies)],
      Z_MCMC$treat2[which(Z_MCMC$studlab %in% Multiarm_studies)], "Study",
      Z_MCMC$studlab[which(Z_MCMC$studlab %in% Multiarm_studies)]
    )
    ##
    t1_E <- Z_MCMC$treat1[which(Z_MCMC$studlab %in% Multiarm_studies)]
    t2_E <- Z_MCMC$treat2[which(Z_MCMC$studlab %in% Multiarm_studies)]
    ##
    stud_E <- Z_MCMC$studlab[which(Z_MCMC$studlab %in% Multiarm_studies)]
    comp_E <- paste(
      Z_MCMC$treat1[which(Z_MCMC$studlab %in% Multiarm_studies)], "vs",
      Z_MCMC$treat2[which(Z_MCMC$studlab %in% Multiarm_studies)]
    )


    for (i in 1:n_multi_comp) {
      for (j in 1:n_multi_comp) {
        ##
        if (i == j) { # Variance of the comparison
          ##
          PREC[i, j] <- data$seTE[which(data$studlab == stud_E[i] & data$treat1 %in% c(t1_E[i], t2_E[i]) & data$treat2 %in% c(t1_E[i], t2_E[i]))]^2
          ##
        } else {
          ##
          if (colnames(H)[i] != row.names(H)[j]) { # Elements of matrix H that do not refer to the same study
            H[i, j] <- 0
          }

          if (stud_E[i] == stud_E[j]) { # Same study
            ##
            if (sum(strsplit(comp_E[i], " vs ")[[1]] %in% strsplit(comp_E[j], " vs ")[[1]]) > 0) { # Common treatment

              # cov(XY, XZ) = (-V(XZ) + V(XY) + V(YZ) )/2

              cmp.1 <- strsplit(comp_E[i], " vs ")[[1]]
              cmp.2 <- strsplit(comp_E[j], " vs ")[[1]]
              f1 <- cmp.1[which(!cmp.1 %in% cmp.2)]
              f2 <- cmp.2[which(!cmp.2 %in% cmp.1)]
              ##
              PREC[i, j] <- (data$seTE[which(data$studlab == stud_E[i] & data$treat1 %in% c(t1_E[i], t2_E[i]) & data$treat2 %in% c(t1_E[i], t2_E[i]))]^2 +
                data$seTE[which(data$studlab == stud_E[i] & data$treat1 %in% c(t1_E[j], t2_E[j]) & data$treat2 %in% c(t1_E[j], t2_E[j]))]^2 -
                data$seTE[which(data$studlab == stud_E[i] & data$treat1 %in% c(f1, f2) & data$treat2 %in% c(f1, f2))]^2) / 2
              ##
            } else { # No common treat
              ##
              PREC[i, j] <- 0
              ##
            }
          } else { # Not in the same study
            ##
            PREC[i, j] <- 0
            ##
          }
        }
      }
    }

    colnames(H) <- row.names(H) <- paste(
      "Comp", Z_MCMC$treat1[which(Z_MCMC$studlab %in% Multiarm_studies)],
      Z_MCMC$treat2[which(Z_MCMC$studlab %in% Multiarm_studies)], "Study",
      Z_MCMC$studlab[which(Z_MCMC$studlab %in% Multiarm_studies)]
    )

    # Invert matrices H and E
    H <- solve(H, tol = 1.1594e-21)
    PREC <- tryCatch(solve(PREC, tol = 1.1594e-21))
  }

  res <- list(
    "H" = H,
    "PREC" = PREC
  )

  res
}
