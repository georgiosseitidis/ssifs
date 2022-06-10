NMAcomparisons <- function(data, ref) {
  nodes <- unique(c(data$treat1, data$treat2)) # Nodes of the network

  # Find all observed comparisons
  obs_comparisons <- data[, grep(paste(c("treat1", "treat2"), collapse = "|"), names(data))]
  obs_comparisons <- unique(obs_comparisons)

  ##
  # Exclude duplicate comparisons
  ##

  rownames(obs_comparisons) <- 1:dim(obs_comparisons)[1]
  colnames(obs_comparisons) <- c("V1", "V2")
  obs_comparisons$AvsB <- paste(obs_comparisons$V1, obs_comparisons$V2, sep = ",")
  obs_comparisons$BvsA <- paste(obs_comparisons$V2, obs_comparisons$V1, sep = ",")
  obs_comparisons$it <- 0

  i <- 1
  while (!sum(obs_comparisons$it) == dim(obs_comparisons)[1]) {
    omit <- which(obs_comparisons$BvsA == obs_comparisons$AvsB[i])

    if (is.na(omit[1])) {
      obs_comparisons[i, "it"] <- 1
      i <- i + 1
    } else {
      obs_comparisons[c(i, omit), "it"] <- 1
      obs_comparisons <- obs_comparisons[-omit, ]
      rownames(obs_comparisons) <- 1:dim(obs_comparisons)[1]
      i <- 0
    }
  }

  # Remove the last auxiliary column it
  obs_comparisons <- obs_comparisons[, -dim(obs_comparisons)[2]]

  ##
  # Change the order for the comparisons that have the reference treatment at the 2nd group
  ##

  obs_comparisons <- obs_comparisons

  # Comparisons that have the ref category in the second group
  zz <- obs_comparisons[which(obs_comparisons$V2 == ref), ]

  # Change the order
  if (dim(zz)[1] > 0) {
    for (i in 1:dim(zz)[1]) {
      z1 <- zz[i, 1]
      z2 <- zz[i, 2]
      z.i <- which(obs_comparisons$V1 == z1 & obs_comparisons$V2 == z2)
      obs_comparisons$V1[z.i] <- z2
      obs_comparisons$V2[z.i] <- z1
    }
  }

  bas_parameters <- obs_comparisons[which(obs_comparisons$V1 == ref | obs_comparisons$V2 == ref), ] # Basic parameters
  fun_parameters <- obs_comparisons[which(!obs_comparisons$V1 == ref & !obs_comparisons$V2 == ref), ] # Functional parameter

  res <- list(
    "comparisons" = obs_comparisons,
    "basic" = bas_parameters,
    "functional" = fun_parameters,
    "nodes" = nodes
  )

  res
}
