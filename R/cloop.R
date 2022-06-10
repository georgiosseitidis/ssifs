cloop <- function(data) {

  # Network's geometry
  data_2 <- unique(data.frame("treat1" = data$treat1, "treat2" = data$treat2))
  tt <- igraph::graph.data.frame(data_2, directed = F)

  ##
  # Find all closed loops of the network
  ##

  # A list with the closed loops of the network
  S <- list()

  for (i in 1:length(igraph::V(tt))) {
    SP <- SP2 <- NULL

    # All possible paths that start and end in the node V(tt)[i]
    SP <- igraph::all_simple_paths(tt, from = igraph::V(tt)[i], to = igraph::neighbors(tt, v = igraph::V(tt)[i]))

    # length SP>1 because if a node is only compared with one other node it stops
    if (!is.null(unlist(SP)) & length(SP) > 1) {

      # Exclude the paths that have length <=2 since they are not closed loops.
      SP2 <- SP[sapply(SP, function(p) length(p) > 2)]

      # Make the igraph object to character
      if (length(SP2) != 0) { # star networks

        S[[i]] <- lapply(SP2, names)
      }
    }
  }

  if (length(S) == 0) {
    stop("Lu & Ades model cannot applied in star networks")
  }

  ##
  # Exclude duplicate loops.
  ##

  loops_dub <- duplicateloops(S)

  UL <- loops_dub$Unique_loops # Unique loops
  min_loop <- loops_dub$min # length of the smallest loop
  max_loop <- loops_dub$max # length of the largest loop

  ##
  # Categorize the loops based on their length
  ##

  loops <- list()
  for (i in min_loop:dim(UL)[2]) {
    loops[[paste0(i, "loop")]] <- UL[which(apply(UL, MARGIN = 1, function(x) sum(is.na(x))) == c(max_loop - i)), ]
    # Omit the NA columns
    loops[[paste0(i, "loop")]] <- loops[[paste0(i, "loop")]][colSums(!is.na(loops[[paste0(i, "loop")]])) > 0]
  }

  loops
}
