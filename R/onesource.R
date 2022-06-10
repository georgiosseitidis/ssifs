onesource <- function(data, Multiarm_studies, funpar) {

  # Comparisons created from multi-arms studies
  comparisons <- data[which(data$studlab %in% Multiarm_studies), c(
    grep("treat1", names(data)),
    grep("treat2", names(data)),
    grep("studlab", names(data))
  )]
  comparisons$treat1 <- as.character(comparisons$treat1)
  comparisons$treat2 <- as.character(comparisons$treat2)

  ##
  # Identification of comparisons with one source of evidence
  ##

  comparisons$NotOnlyMultyArm <- FALSE
  i <- 1
  while (sum(is.na(comparisons$OnlyMultyArm)) <= dim(comparisons)[1] & i <= dim(comparisons)[1]) {
    zz <- which(data$treat1 %in% c(comparisons$treat1[i], comparisons$treat2[i]) &
      data$treat2 %in% c(comparisons$treat1[i], comparisons$treat2[i]) &
      data$studlab != comparisons$studlab[i])
    if (!is.na(zz[1])) {
      comparisons$NotOnlyMultyArm[i] <- TRUE
    }
    i <- i + 1
  }

  # Comparisons that are created from multi-arm studies with NotOnlyMultyArm==FALSE, have only one source of evidence
  mult_comp <- NULL
  if (sum(comparisons$NotOnlyMultyArm) != dim(comparisons)[1]) {
    zz <- which(comparisons$NotOnlyMultyArm == FALSE)
    for (i in 1:length(zz)) {
      mult_comp[i] <- paste(comparisons$treat1[zz[i]], comparisons$treat2[zz[i]], sep = ",")
    }
  }

  ##
  # Specify the comparisons in which inconsistency factor could be added
  ##

  IFcomp <- paste(funpar$V1, funpar$V2, sep = ",")
  if (length(mult_comp) > 0 & length(which(IFcomp %in% mult_comp)) > 0) { # comparisons with one source may be basic parameters
    IFcomp <- IFcomp[-which(IFcomp %in% mult_comp)]
    if (length(IFcomp) == 0) {
      stop("Lu and Ades model cannot be applied in this network", call. = FALSE)
    }
  }

  res <- list(
    "IFcomp" = IFcomp,
    "mult_comp" = mult_comp
  )

  res
}
