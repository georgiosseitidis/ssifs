#' @importFrom netmeta netmeta

higgdes <- function(Z, Multiarm_studies, studies_multi, p_Jackson, p_Higgins) {
  IF_parameters <- data.frame(design = NA, IF = colnames(Z)[4:dim(Z)[2]], add = NA, stringsAsFactors = FALSE)
  row_mult_desg <- NULL

  # Find the design of each inconsistency factor
  for (i in 1:dim(IF_parameters)[1]) {
    zzz <- strsplit(IF_parameters$IF[i], split = "_")[[1]] # design i

    if (length(zzz) == 1) {
      IF_parameters$design[i] <- zzz[1]
    } else {
      IF_parameters$design[i] <- zzz[2]
      row_mult_desg <- c(row_mult_desg, i)
    }
  }
  ##
  IF_parameters$IForig <- IF_parameters$IF
  for (i in row_mult_desg) {
    IF_parameters$IF[i] <- strsplit(IF_parameters$IF[i], split = "_")[[1]][1]
  }

  ##
  # Find which inconsistency factor should be included in the NMA model
  ##

  i <- 1
  while (sum(IF_parameters$add, na.rm = TRUE) < p_Higgins) {
    if (i %in% c(1)) {
      ##
      IF_parameters$add[i] <- 0
      ##
    } else if (length(which(IF_parameters$add == 0)) == p_Jackson - p_Higgins) {
      # If number of parameters that cannot have inconsistency factor is reached
      IF_parameters$add[which(is.na(IF_parameters$add))] <- 1
      ##
    } else {
      ##
      z_par <- IF_parameters$IF[1:(i - 1)] # Inconsistency factors names

      # Directly compared
      z_par_list <- strsplit(z_par, split = " ; ")
      z_par_inv <- NULL
      for (j in 1:length(z_par)) {
        z_par_inv <- c(z_par_inv, paste(z_par_list[[j]][2], z_par_list[[j]][1], sep = " ; "))
      }
      ##
      dir_comp <- IF_parameters$IF[i] %in% c(z_par, z_par_inv)

      # Indirectly compared
      ind_comp <- FALSE

      if (dir_comp == FALSE) {
        if (length(z_par) > 1) {
          z_1 <- strsplit(IF_parameters$IF[i], split = " ; ")[[1]][1]
          z_2 <- strsplit(IF_parameters$IF[i], split = " ; ")[[1]][2]

          if ((!z_1 %in% unlist(strsplit(z_par, split = " ; "))) | (!z_2 %in% unlist(strsplit(z_par, split = " ; ")))) {
            # Treatment was not observed in the previous designs
            ind_comp <- FALSE
          } else {

            # Find if can be obtained indirectly with a common treatment

            z_list <- strsplit(z_par, split = " ; ")
            z_1_comp <- NULL
            z_2_comp <- NULL
            for (j in 1:length(z_list)) {
              if (z_1 %in% z_list[[j]]) {
                z_1_comp <- c(z_1_comp, z_par[j])
              }
              if (z_2 %in% z_list[[j]]) {
                z_2_comp <- c(z_2_comp, z_par[j])
              }
            }
            ##
            ##
            for (j in 1:length(z_1_comp)) {
              for (l in 1:length(z_2_comp)) {
                if (sum(table(unlist(strsplit(c(z_1_comp[j], z_2_comp[l]), split = " ; "))) == 2) > 0) { # If a treatment has freq 2 it indicates that they have a common comparator
                  ind_comp <- TRUE
                  break
                }
              }
            }


            # Find if can be obtained indirectly with many ways
            if (ind_comp == FALSE) {
              z_comp <- unlist(strsplit(z_par, split = " ; "))
              data_nma <- data.frame(
                TE = rep(1, length(z_comp) / 2), seTE = rep(0.5, length(z_comp) / 2),
                treat1 = z_comp[seq(1, length(z_comp), 2)],
                treat2 = z_comp[seq(0, length(z_comp), 2)[-1]],
                studlab = 1:(length(z_comp) / 2), stringsAsFactors = FALSE
              )
              ##
              Discon_NMA <- subnet(data_nma)
              if (Discon_NMA$nsub == 1) {
                # Connected network

                m <- netmeta(TE = data_nma$TE, seTE = data_nma$seTE, treat1 = data_nma$treat1, treat2 = data_nma$treat2, studlab = data_nma$studlab, sm = "MD")
                ind_tab <- m$TE.indirect.random

                if (!is.na(ind_tab[which(row.names(ind_tab) == z_1), which(colnames(ind_tab) == z_2)])) {
                  ind_comp <- TRUE
                }
              } else {
                # Disconnected network

                # Find which sub-network have the comparison that is required

                m.comp <- NULL
                for (r in 1:length(Discon_NMA$subnetworks)) {
                  if (sum(c(z_1, z_2) %in% Discon_NMA$subnetworks[[r]]) == 2) {
                    data_nma <- data_nma[which(data_nma$treat1 %in% Discon_NMA$subnetworks[[r]] | data_nma$treat2 %in% Discon_NMA$subnetworks[[r]]), ]
                    m <- netmeta(
                      TE = data_nma$TE, seTE = data_nma$seTE, treat1 = data_nma$treat1, treat2 = data_nma$treat2,
                      studlab = data_nma$studlab, sm = "MD", tol.multiarm.se = 10000
                    )
                    ind_tab <- m$TE.indirect.random

                    if (!is.na(ind_tab[which(row.names(ind_tab) == z_1), which(colnames(ind_tab) == z_2)])) {
                      ind_comp <- TRUE
                      break
                    }
                  }
                }
              }
            }
          }
        }
      }

      # Decision about if inconsistency factor should be added

      if (dir_comp | (ind_comp)) {
        IF_parameters$add[i] <- 1
      } else {
        IF_parameters$add[i] <- 0
      }
    }

    i <- i + 1
  }

  ##
  # Set values to the design matrix Z
  ##

  if (sum(is.na(IF_parameters$add)) > 0) {
    IF_parameters$add[is.na(IF_parameters$add)] <- 0
  }
  ##
  ##
  for (i in which(IF_parameters$add == 1)) {
    if (i %in% row_mult_desg) {
      ##
      z_des_nma <- as.character(sapply(names(studies_multi), FUN = function(x) {
        paste(unlist(strsplit(x, split = " ; ")), sep = "", collapse = "|")
      }, simplify = TRUE))
      Mult_stud_des_i <- which(z_des_nma == IF_parameters$design[i])
      ##
      pos <- which(paste(Z[, 1], Z[, 2], sep = " ; ") == IF_parameters$IF[i] & Z[, 3] %in% studies_multi[[Mult_stud_des_i]])
      pos_inv <- which(paste(Z[, 2], Z[, 1], sep = " ; ") == IF_parameters$IF[i] & Z[, 3] %in% studies_multi[[Mult_stud_des_i]])
      ##
      Z[pos, IF_parameters$IForig[i]] <- 1
      Z[pos_inv, IF_parameters$IForig[i]] <- -1

      ##
    } else {
      ##
      pos <- which(paste(Z[, 1], Z[, 2], sep = " ; ") == IF_parameters$IF[i] & !Z[, 3] %in% Multiarm_studies)
      pos_inv <- which(paste(Z[, 2], Z[, 1], sep = " ; ") == IF_parameters$IF[i] & !Z[, 3] %in% Multiarm_studies)

      Z[pos, IF_parameters$IForig[i]] <- 1
      Z[pos_inv, IF_parameters$IForig[i]] <- -1
      ##
    }
  }
  ##
  Z[is.na(Z)] <- 0
  Z <- Z[, -which(colnames(Z) %in% IF_parameters$IForig[which(IF_parameters$add == 0)])]

  Z
}
