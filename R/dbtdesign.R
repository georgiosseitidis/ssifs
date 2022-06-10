dbtdesign <- function(data, ref, Multiarm_studies, studies_multi, design, IFmulti, p_Jackson, method, p_Higgins) {
  y <- data$TE # Treatment effect

  # Structure of the design matrix
  Z <- data.frame(data$treat1, data$treat2, data$studlab, matrix(nrow = length(data$treat1), ncol = p_Jackson), stringsAsFactors = F)


  if (length(Multiarm_studies) > 0) {
    colnames(Z) <- c("treat1", "treat2", "studlab", design[1:(length(design) - length(IFmulti))])

    exclude <- list() # List with the comparisons to be excluded

    clname <- NULL # Names of the multi-arm designs
    m2 <- 1

    ##
    # Find comparisons to be excluded and the inconsistency factor's names for the design matrix
    ##

    for (i in (length(design) - length(IFmulti) + 1):length(design)) {

      # All possible comparisons
      zz <- gtools::permutations(IFmulti[m2] + 1, 2, strsplit(design[i], split = " ; ")[[1]])
      m2 <- m2 + 1

      # Keep only the functional parameters
      zz_fun <- which(!zz[, 1] %in% ref & !zz[, 2] %in% ref)
      # Exclude the first functional parameter
      exclude[[design[i]]] <- zz[zz_fun[1], ]

      zz <- zz[-which(apply(zz, 1, paste, collapse = " ") %in% apply(gtools::permutations(2, 2, zz[zz_fun[1], ]), 1, paste, collapse = " ")), ]
      # Keep the unique combinations (e.g AB and BA are the same. Keep one)
      zz_1 <- paste(zz[, 1], zz[, 2], sep = " ; ")
      zz_2 <- paste(zz[, 2], zz[, 1], sep = " ; ")
      zz <- zz_1[1]
      m <- 2

      while (m <= length(zz_1)) {
        if (!zz_1[m] %in% zz_2[1:m]) {
          zz <- c(zz, zz_1[m])
        }
        m <- m + 1
      }
      zz <- paste(zz, "_", paste(strsplit(design[i], split = " ; ")[[1]], collapse = "|"), sep = "")
      clname <- c(clname, zz)
    }

    # Colnames for the multi-arm design IFs
    colnames(Z)[(dim(Z)[2] - length(clname) + 1):dim(Z)[2]] <- clname

    ##
    # Exclude unnecessary comparisons
    ##

    studies_multi_des <- list()

    for (i in Multiarm_studies) {
      ex <- which(Z[, 3] == i) # Comparisons in i study

      # Find the design of the study
      des_z <- unique(c(Z[ex, 1], Z[ex, 2]))
      ##
      ex_comp <- exclude[[names(exclude)[which(names(exclude) %in% apply(gtools::permutations(length(des_z), length(des_z), des_z), 1, paste, collapse = " ; "))]]]
      ##
      studies_multi_des[[i]] <- names(exclude)[which(names(exclude) %in% apply(gtools::permutations(length(des_z), length(des_z), des_z), 1, paste, collapse = " ; "))]
      ##
      ex <- ex[which(Z[ex, 1] %in% ex_comp & Z[ex, 2] %in% ex_comp)]
      ##
      Z <- Z[-ex, ]
      y <- y[-ex]
    }
  } else {
    colnames(Z) <- c("treat1", "treat2", "studlab", design[1:(length(design))])
  }

  ##
  # Generate matrix Z
  ##

  if (method == "Jackson") {
    Z <- jackdes(Z, Multiarm_studies, studies_multi_des)
  } else {
    Z <- higgdes(Z, Multiarm_studies, studies_multi, p_Jackson, p_Higgins)
  }

  res <- list(
    "Z" = Z,
    "y" = y
  )

  res
}
