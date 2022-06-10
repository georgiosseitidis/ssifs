numberIF.DBT <- function(n_nodes, Multiarm_studies, design, design_multi, method) {
  IFmulti <- NULL # Inconsistency factors from multi-arm designs
  p_DBT <- NULL # Number of inconsistency factors for the Higgins approach

  if (length(Multiarm_studies) > 0) {

    ##
    # Number of inconsistency factors for the Jackson approach
    ##

    p_Jackson <- length(design) - length(design_multi)

    # Number of inconsistency factors from multi-arm designs.
    z <- 0
    for (i in (length(design) - length(design_multi) + 1):length(design)) {
      z <- z + length(strsplit(design[i], split = " ; ")[[1]]) - 1
      IFmulti <- c(IFmulti, length(strsplit(design[i], split = " ; ")[[1]]) - 1)
    }

    # Total number of inconsistency factors for the Jackson model
    p_Jackson <- p_Jackson + z

    ##
    # Number of inconsistency factors for the Higgin's approach
    ##

    if (method == "DBT") {
      p_DBT <- p_Jackson - (n_nodes - 1)
      if (p_DBT == 0) {
        stop("DBT approach cannot applied in this network", call. = FALSE)
      }
    }
  } else {

    ##
    # Number of inconsistency factors for the Jackson approach
    ##

    p_Jackson <- length(design)

    ##
    # Number of inconsistency factors for the Higgin's approach
    ##

    if (method == "DBT") {
      p_DBT <- p_Jackson - (n_nodes - 1)
      if (p_DBT == 0) {
        stop("DBT approach cannot applied in this network", call. = FALSE)
      }
    }
  }

  res <- list(
    "Jackson" = p_Jackson,
    "Higgins" = p_DBT,
    "IFmulti" = IFmulti
  )

  res
}
