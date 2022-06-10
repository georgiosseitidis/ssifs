indloop <- function(IF, Close_loops, nodes, obs_comparisons) {
  maxloop <- as.numeric(max(unlist(strsplit(names(Close_loops), split = "loop")))) # max length of a loop
  ##
  nloop <- unlist(lapply(Close_loops, dim))
  nloop <- sum(nloop[seq(1, length(nloop), 2)]) # Number of close loops

  ##
  # Construct a table with all the close loops
  ##

  loop_table <- as.data.frame(matrix(nrow = nloop, ncol = maxloop + 1))
  colnames(loop_table) <- c(paste("X", 1:maxloop, sep = ""), "incl")
  m <- 0
  for (i in names(Close_loops)) {
    for (j in 1:dim(Close_loops[[i]])[1]) {
      m <- m + 1
      for (k in 1:dim(Close_loops[[i]])[2]) {
        loop_table[m, k] <- Close_loops[[i]][j, k]
      }
    }
  }


  ##
  # Find the independent loops
  ##

  # A loop is independent if it contains at least one edge which is not a part of any other independent loop

  loop_table$ind <- 0
  loop_table$ind[1] <- 1 # first closed loop is independent


  if (dim(loop_table)[1] > 1) {
    for (i in 2:dim(loop_table)[1]) {

      # A vector with the loop
      if (length(which(is.na(loop_table[i, 1:maxloop]))) > 0) {
        v <- loop_table[i, 1:maxloop][-which(is.na(loop_table[i, 1:maxloop]))]
      } else {
        v <- loop_table[i, 1:maxloop]
      }

      v <- plyr::mapvalues(unlist(v), from = nodes, to = 1:length(nodes), warn_missing = FALSE)

      loop_comp <- apply(gtools::permutations(n = length(v), r = 2, v = as.numeric(v)), 1, paste, collapse = ",") # loop comparisons
      # Exclude the unobserved comparisons from the loop comparison
      loop_comp <- loop_comp[which(loop_comp %in% obs_comparisons$AvsB)]

      # Find the comparisons of the previous independent loops
      loop_previous <- NULL
      for (j in which(loop_table$ind == 1)) {

        # A vector for the previous loop
        if (length(which(is.na(loop_table[j, 1:maxloop]))) > 0) {
          v_pr <- loop_table[j, 1:maxloop][-which(is.na(loop_table[j, 1:maxloop]))]
        } else {
          v_pr <- loop_table[j, 1:maxloop]
        }

        v_pr <- plyr::mapvalues(unlist(v_pr), from = nodes, to = 1:length(nodes), warn_missing = FALSE)

        loop_comp_pr <- apply(gtools::permutations(n = length(v_pr), r = 2, v = as.numeric(v_pr)), 1, paste, collapse = ",") # loop comparisons

        # Exclude the unobserved comparisons from the loop comparison
        loop_comp_pr <- loop_comp_pr[which(loop_comp_pr %in% obs_comparisons$AvsB)]
        loop_previous <- c(loop_previous, loop_comp_pr)
      }

      # Check if any of the above comparisons is not included in the previous independent loops
      if (sum(!loop_comp %in% loop_previous) > 0) {
        loop_table$ind[i] <- 1
      }
    }
  }

  # Exclude the non independent comparisons
  loop_table <- loop_table[which(loop_table$ind == 1), ]

  ##
  # Keep the comparisons that referred to a independent loops
  ##

  IF_final <- NULL
  IF_index <- data.frame(IF = IF, index = 0, stringsAsFactors = F)
  for (i in 1:dim(loop_table)[1]) {

    # A vector with the loop
    if (length(which(is.na(loop_table[i, 1:maxloop]))) > 0) {
      v <- loop_table[i, 1:maxloop][-which(is.na(loop_table[i, 1:maxloop]))]
    } else {
      v <- loop_table[i, 1:maxloop]
    }

    for (j in 1:dim(IF_index)[1]) {
      if (sum(strsplit(IF_index$IF[j], split = ",")[[1]] %in% v) == 2 & is.na(loop_table[i, maxloop + 1]) & IF_index$index[j] == 0) {
        loop_table[i, maxloop + 1] <- 1
        IF_index$index[j] <- 1
        IF_final <- c(IF_final, IF_index$IF[j])
      }
    }
  }

  if (is.null(IF_final)) {
    stop("Lu and Ades model cannot be applied in this network", call. = FALSE)
  }

  IF_final
}
