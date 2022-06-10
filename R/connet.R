connet <- function(data) {
  sub_networks <- subnet(data)
  exstud <- NULL # Excluded studies
  omit <- unlist(sub_networks$disnet) # Nodes to be excluded

  # Exclude studies that do not belong in the largest sub-network
  if (!is.null(omit)) {
    exstud <- (data[which(data$treat1 %in% omit | data$treat2 %in% omit), grep("studlab", names(data))])
    data <- data[-(which(data$treat1 %in% omit | data$treat2 %in% omit)), ]
  }

  res <- list(
    "data" = data,
    "excluded_studies" = exstud,
    "n_sub" = sub_networks$nsub,
    "subnet" = sub_networks$subnetworks
  )

  res
}
