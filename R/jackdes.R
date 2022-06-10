jackdes <- function(X, Multiarm_studies, studies_multi_des) {

  ##
  # Set values to the design matrix
  ##

  for (i in 1:dim(X)[1]) {
    if (!X[i, 3] %in% Multiarm_studies) { # if it is a two-arm studies

      if (length(which(colnames(X) == paste(X[i, 1], X[i, 2], sep = " ; "))) > 0) {
        X[i, which(colnames(X) == paste(X[i, 1], X[i, 2], sep = " ; "))] <- 1
      } else {
        X[i, which(colnames(X) == paste(X[i, 2], X[i, 1], sep = " ; "))] <- -1
      }
    } else { # if it is a multi-arm study

      des_name <- paste(strsplit(studies_multi_des[[which(names(studies_multi_des) == X[i, 3])]], " ; ")[[1]], collapse = "|")

      if (length(which(colnames(X) == paste(paste(X[i, 1], X[i, 2], sep = " ; "), "_", des_name, sep = ""))) > 0) {
        X[i, which(colnames(X) == paste(paste(X[i, 1], X[i, 2], sep = " ; "), "_", des_name, sep = ""))] <- 1
      } else {
        X[i, which(colnames(X) == paste(paste(X[i, 2], X[i, 1], sep = " ; "), "_", des_name, sep = ""))] <- -1
      }
    }
  }

  # Replace NAs with zero
  X[is.na(X)] <- 0

  X
}
