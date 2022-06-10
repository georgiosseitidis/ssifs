design.matrix <- function(IF) {
  obs_comparisons <- IF$observed # observed comparisons of the network
  obs_comparisons$AvsB <- paste(obs_comparisons$V1, obs_comparisons$V2, sep = "vs")

  # Basic parameters of the network
  basic_parameters <- IF$basic
  basic_parameters$AvsB <- paste(basic_parameters$V1, basic_parameters$V2, sep = "vs")

  # Write the comparisons as a function of basic parameters
  X <- data.frame(obs_comparisons, matrix(nrow = dim(obs_comparisons)[1], ncol = dim(basic_parameters)[1]))
  colnames(X)[which(colnames(X) == "AvsB")] <- c("comparison")
  colnames(X)[(which(colnames(X) == "comparison") + 1):dim(X)[2]] <- as.character(basic_parameters$V2)

  # Set values in the design matrix
  for (i in 1:dim(X)[1]) {
    if (X$comparison[i] %in% basic_parameters$AvsB) { # basic parameter
      X[i, which(colnames(X) == X$V2[i])] <- 1
    } else { # functional parameter

      X[i, which(colnames(X) == X$V1[i])] <- -1
      X[i, which(colnames(X) == X$V2[i])] <- 1
    }
  }

  # Add columns to the design matrix for the inconsistency factors
  IF_comparisons <- IF$comparison
  for (i in 1:length(IF_comparisons)) {
    X[, paste("W", IF_comparisons[i])] <- 0
    X[which(paste(X$V1, X$V2, sep = ",") %in% IF_comparisons[i] |
      paste(X$V2, X$V1, sep = ",") %in% IF_comparisons[i]), paste("W", IF_comparisons[i])] <- 1
  }

  # Replace NAs with zero
  X[is.na(X)] <- 0

  X
}
