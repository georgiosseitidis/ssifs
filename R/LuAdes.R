LuAdes <- function(data, ref, Multiarm_studies) {


  ##
  # Find all the closed loops of the network
  ##

  cl <- cloop(data)

  ##
  # Find network characteristics
  ##

  NMAchar <- NMAcomparisons(data, ref)
  ##
  basic <- NMAchar$basic # Basic parameters
  functional <- NMAchar$functional # Functional parameters

  if (dim(functional)[1] == 0) {
    stop("Lu & Ades model since the network does not contain functional parameters", call. = FALSE)
  }

  nodes <- NMAchar$nodes # Nodes of the network
  observed_comp <- NMAchar$comparisons # Observed comparisons

  ##
  # Define the comparisons where inconsistency factors should be added
  ##

  # Inconsistency Factors are added to all functional parameters
  if (length(Multiarm_studies) == 0) {
    IF <- paste(functional$V1, functional$V2, sep = ",")
    one_source_excluded <- NULL
  } else {

    # Exclude the comparisons with one source of evidence
    one_source_comp <- onesource(data, Multiarm_studies, functional)
    one_source_excluded <- one_source_comp$mult_comp
    IF <- one_source_comp$IFcomp
  }

  # Exclude the bridges of the network
  IF_bridge <- bridge(IF, cl)

  # Keep the inconsistency factors that referred to a independent loop
  IFs <- indloop(IF_bridge, cl, nodes, observed_comp)

  res <- list(
    "comparison" = IFs,
    "observed" = observed_comp[, 1:2],
    "basic" = basic[, 1:2],
    "functional" = functional[, 1:2],
    "one_source" = one_source_excluded
  )

  res
}
