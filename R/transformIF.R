transformIF <- function(IF, trt_trans, method) {
  trt_trans$from <- as.character(trt_trans$from)
  trt_trans$to <- as.character(trt_trans$to)

  if (method == "LuAdes") {
    IF <- unlist(strsplit(IF, split = ","))
    IF <- plyr::mapvalues(IF, from = trt_trans$to, to = trt_trans$from, warn_missing = FALSE)
    ##
    z_1 <- IF[seq(1, length(IF), 2)]
    z_2 <- IF[seq(2, length(IF), 2)]
    ##
    IF_trans <- paste(z_1, z_2, sep = ",")
  } else {
    IF <- strsplit(IF, split = " ; ")

    IF_trans <- NULL

    for (i in 1:length(IF)) {
      if (sum(IF[[i]] %in% trt_trans$to) == 2) {
        z <- paste(trt_trans$from[which(trt_trans$to == IF[[i]][1])],
          trt_trans$from[which(trt_trans$to == IF[[i]][2])],
          sep = " ; "
        )

        IF_trans <- c(IF_trans, z)
      } else {
        z <- strsplit(IF[[i]][2], split = "_")[[1]]
        z_1 <- trt_trans$from[which(trt_trans$to == z[1])]
        ##
        z_2 <- strsplit(z[2], split = "|")[[1]]
        z_2 <- z_2[-which(z_2 == "|")]
        z_2 <- plyr::mapvalues(z_2, from = trt_trans$to, to = trt_trans$from, warn_missing = FALSE)
        ##
        z_12 <- paste(z_1, paste(z_2, collapse = ""), sep = "_")
        z_f <- paste(trt_trans$from[which(trt_trans$to == IF[[i]][1])], z_12, sep = " ; ")
        ##
        IF_trans <- c(IF_trans, z_f)
      }
    }

    IF_trans
  }
}
