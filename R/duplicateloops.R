duplicateloops <- function(S) {
  S <- lapply(rapply(S, enquote, how = "unlist"), eval)
  ##
  n_nodes <- unlist(lapply(S, length)) # Number of nodes in the loops
  max_loop <- max(n_nodes) # Maximum number of nodes in a loop
  min_loop <- min(n_nodes) # Minimum number of nodes in a loop

  ##
  # Find duplicate loops
  ##

  dub <- rep(NA, length(S))

  while (sum(!is.na(dub)) < length(S)) {

    # loop i
    i <- which(is.na(dub))[1]
    loop_i <- S[[i]]
    dub[i] <- FALSE
    n_loop_i <- length(loop_i)

    # Find the loops with the same nodes
    same_nodes <- sapply(S, FUN = function(x) {
      d <- sum(x %in% loop_i)
    })

    # Duplicate loops are the loops that have the same length with the loop i and the same nodes
    dub[which(is.na(dub) & same_nodes == n_loop_i & n_nodes == n_loop_i)] <- TRUE
  }

  # Exclude duplicate loops
  S <- S[which(dub == FALSE)]

  # Store the unique loops to a dataframe
  rows <- length(S) # number of rows
  cols <- max_loop # number of columns
  UL <- data.frame(matrix(rep(NA, rows * cols), nrow = rows, ncol = cols))

  for (i in 1:rows) {
    for (j in 1:length(S[[i]])) {
      UL[i, j] <- unlist(S[i])[j]
    }
  }

  res <- list(
    "Unique_loops" = UL,
    "min" = min_loop,
    "max" = max_loop
  )

  res
}
