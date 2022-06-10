NMAdesigns <- function(data, Multiarm_studies) {

  ##
  # Find all two-arms designs
  ##

  design_2arm <- unique(paste(data$treat1[which(!data$studlab %in% Multiarm_studies)], data$treat2[which(!data$studlab %in% Multiarm_studies)], sep = " ; "))
  # exclude duplicates e.g A | B is the same with B | A
  design_2arm_V2 <- unique(paste(data$treat2[which(!data$studlab %in% Multiarm_studies)], data$treat1[which(!data$studlab %in% Multiarm_studies)], sep = " ; "))
  design <- design_2arm[1]
  i <- 2
  while (i <= length(design_2arm)) {
    if (!design_2arm[i] %in% design_2arm_V2[1:i]) {
      design <- c(design, design_2arm[i])
    }
    i <- i + 1
  }

  ##
  # Find all multi-arms designs
  ##

  if (length(Multiarm_studies) > 0) {
    design_multi <- list()
    studies_multi <- list()

    # Comparisons in every multi-arm study
    for (i in 1:length(Multiarm_studies)) {
      design_multi[[i]] <- unique(c(data$treat1[which(data$studlab %in% Multiarm_studies[i])], data$treat2[which(data$studlab %in% Multiarm_studies[i])]))
      studies_multi[[paste(design_multi[[i]], collapse = " ; ")]] <- data$studlab[which(data$studlab == Multiarm_studies[i])][1]
    }

    # Keep unique designs from multi-arms designs
    design_multi_uniq <- paste(design_multi[[1]], collapse = " ; ")
    i <- 2

    while (i <= length(design_multi)) {
      included <- FALSE

      for (j in 1:(i - 1)) {
        n_str <- length(design_multi[[j]])

        # Check if this design is included in at least one combination of the arms of the included designs
        if (paste(design_multi[[i]], collapse = " ") %in% apply(gtools::permutations(n_str, n_str, design_multi[[j]]), 1, paste, collapse = " ")) {
          included <- TRUE
        }
      }

      # If the design is not included, then add it to the vector with the designs
      if (included == FALSE) {
        design_multi_uniq <- c(design_multi_uniq, paste(design_multi[[i]], collapse = " ; "))
      }

      i <- i + 1
    }

    # Designs of the network
    design <- c(design, design_multi_uniq)
  } else {
    design_multi_uniq <- studies_multi <- NULL
  }

  res <- list(
    "designs" = design,
    "multi_design" = design_multi_uniq,
    "studies_multi" = studies_multi
  )

  res
}
