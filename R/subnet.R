subnet <- function(data) {
  t <- igraph::graph.data.frame(cbind(as.character(data$treat1), as.character(data$treat2)), directed = FALSE)

  # Find the Strong Connected Components using Kosaraju algorithm
  oldopt <- options()
  on.exit(options(oldopt))

  options(listexpressions = 500000, warn = -1)
  SCC <- RevEcoR::KosarajuSCC(t)

  ##
  subnetworks <- lapply(SCC, names)
  names(subnetworks) <- 1:length(subnetworks)

  nsub <- length(SCC) # Number of sub-networks

  dis_nodes <- list() # Disconnected nodes
  ndis <- 0

  for (i in 1:nsub) {
    if (length(SCC[[i]]) < max(unlist(lapply(SCC, length)))) {
      ndis <- ndis + 1
      dis_nodes[[ndis]] <- names(SCC[[i]])
    }
  }

  res <- list(
    "nsub" = nsub,
    "subnetworks" = subnetworks,
    "disnet" = dis_nodes
  )

  res
}
