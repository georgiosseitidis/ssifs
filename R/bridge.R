bridge <- function(IF, Close_loops) {
  bridges_comp <- NULL

  for (i in 1:length(IF)) {
    bridge <- 1

    ##
    # Specification of bridges
    ##

    # If the element i is included in any closed loop, it is not a bridge

    for (j in 1:length(Close_loops)) {
      zz <- Close_loops[[j]]
      for (k in 1:dim(zz)[1]) {
        if (sum(strsplit(IF[i], ",")[[1]] %in% zz[k, ]) == 2) { # equal to 2 since both nodes must be included in the loop
          bridge <- 0
        }
      }
    }

    if (bridge == 1) {
      bridges_comp <- c(bridges_comp, i)
    }
  }

  # Exclude the bridges of the network from the IF vector
  if (length(bridges_comp) > 0) {
    IF <- IF[-bridges_comp]
  }

  if (length(IF) == 0) {
    stop("Lu and Ades model cannot be applied in this network", call. = FALSE)
  }

  IF
}
