correct.data <- function(data) {
  study_freq <- as.data.frame(table(data$studlab))
  colnames(study_freq)[1] <- "study"

  # Studies with freq 1 are ok

  # Studies with freq 2 are wrong
  two <- which(study_freq$Freq == 2)

  if (length(two) == 1) {
    stop(paste0("Please check study ", as.character(study_freq$study[two]), ". Comparisons are missing!"), call. = FALSE)
  } else if (length(two) > 1) {
    stop(paste0(
      "Please check the following studies: ",
      paste(as.character(study_freq$study[two]), collapse = ", "), ". Comparisons are missing!"
    ),
    call. = FALSE
    )
  }

  # More than two comparisons
  more_than_two <- which(study_freq$Freq > 2)

  for (i in more_than_two) {

    # Check the number of treatments for each study
    study_i <- data[which(data$studlab == study_freq$study[i]), ]

    # Number of treatments
    n_trt <- length(unique(c(study_i$treat1, study_i$treat2)))

    # Check if the comparisons are unique
    comp_i_t1t2 <- unique(paste(study_i$treat1, "vs", study_i$treat2))
    comp_i_t2t1 <- unique(paste(study_i$treat2, "vs", study_i$treat1))

    if (length(comp_i_t1t2) != nrow(study_i) | length(comp_i_t2t1) != nrow(study_i)) {
      stop(paste("Study", study_i$studlab[1], "has duplicate comparisons"), call. = FALSE)
    }

    if (sum(comp_i_t1t2 %in% comp_i_t2t1) > 0) {
      stop(paste("Study", study_i$studlab[1], "has duplicate comparisons"), call. = FALSE)
    }

    # Check if all comparisons are present
    if (dim(utils::combn(n_trt, 2))[2] != nrow(study_i)) {
      stop(paste0("Please check study ", as.character(study_freq$study[i]), ". Comparisons are missing!"), call. = FALSE)
    }
  }
}
