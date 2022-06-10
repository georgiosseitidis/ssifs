jags.summary <- function(model, trt_trans, Data_MCMC, IF, digits, method) {
  covariates <- model$parameters.to.save

  # jags output
  model_effects <- model$BUGSoutput
  model_effects_summary <- as.data.frame(model_effects$summary)

  # Table with the MCMC NMA results
  result <- data.frame(matrix(nrow = dim(model_effects_summary)[1], ncol = 4))
  colnames(result) <- c("Parameter", "Estimate", "lb", "ub")
  ##
  result$Parameter <- row.names(model_effects_summary)
  result$Estimate <- model_effects_summary$`50%`
  result$lb <- model_effects_summary$`2.5%`
  result$ub <- model_effects_summary$`97.5%`
  ##
  result[-grep(paste(covariates[-which(covariates == "d")], collapse = "|"), result$Parameter), "Parameter"] <- paste(trt_trans$from, "vs", trt_trans$from[as.numeric(Data_MCMC$ref)])
  ##
  result[which(result$Parameter == "precision"), 2:4] <- 1 / result[which(result$Parameter == "precision"), 2:4]
  result[which(result$Parameter == "precision"), "Parameter"] <- "tau2"
  result[which(result$Parameter == "sigma"), "Parameter"] <- "sigma2"

  ##
  # Posterior inclusion probabilities
  ##

  beta <- result[grep(paste("beta", sep = "|", collapse = "|"), result$Parameter), ]
  beta <- beta[(dim(beta)[1] - Data_MCMC$p + 1):dim(beta)[1], ]
  ##
  gamma <- result[grep(paste("gamma", sep = "|", collapse = "|"), result$Parameter), ]
  gamma$Estimate <- model$BUGSoutput$mean$gamma
  gamma <- gamma[(dim(gamma)[1] - Data_MCMC$p + 1):dim(gamma)[1], ]


  IF_PIP <- data.frame(matrix(ncol = 6, nrow = dim(gamma)[1]))
  colnames(IF_PIP) <- c("Comparison", "Design", "PIP", "b", "b.lb", "b.ub")

  if (method == "LuAdes") {
    IF_PIP$Comparison <- IF
    IF_PIP$Design <- NA
  } else {
    z <- strsplit(IF, split = "_")
    z_len <- sapply(z, length)

    IF_PIP$Comparison <- IF_PIP$Design <- sapply(z, FUN = function(x) {
      x[1]
    })

    if (sum(z_len > 1) > 0) {
      IF_PIP$Design[which(z_len > 1)] <- sapply(z[which(z_len > 1)], FUN = function(x) {
        x[2]
      })
    }

    if (sum(z_len == 1) > 0) {
      IF_PIP$Design[which(z_len == 1)] <- sapply(z[which(z_len == 1)], FUN = function(x) {
        paste0(trimws(unlist(strsplit(x, split = ";"))), collapse = "")
      })
    }
  }
  IF_PIP$PIP <- round(gamma$Estimate, digits = digits)
  IF_PIP$b <- round(beta$Estimate, digits = digits)
  IF_PIP$b.lb <- round(beta$lb, digits = digits)
  IF_PIP$b.ub <- round(beta$ub, digits = digits)



  ##
  # Posterior Odds
  ##

  drawns <- as.matrix(model$BUGSoutput$sims.matrix)
  drawns <- as.data.frame(drawns[, grep("gamma", colnames(drawns))])
  drawns$model <- NA
  colnames(drawns)[1:length(IF)] <- IF

  if (dim(drawns)[2] == 2) {
    drawns$model <- sapply(drawns[, -dim(drawns)[2]], FUN = function(x) {
      paste(colnames(drawns)[which(x == 1)], collapse = ",", sep = " ")
    })
  } else {
    drawns$model <- apply(drawns[, -dim(drawns)[2]], 1, function(x) {
      paste(colnames(drawns)[which(x == 1)], collapse = ",", sep = " ")
    })
  }

  models <- as.data.frame(table(drawns$model), stringsAsFactors = F)
  models <- models[order(models$Freq, decreasing = T), ]
  models$Percent <- models$Freq / dim(drawns)[1]
  models$Var1[which(models$Var1 == "")] <- "No IFs"
  colnames(models)[c(1, 3)] <- c("IFs", "f(m|y)")
  models$PO_IFCONS <- models$`f(m|y)`[which(models$IFs == "No IFs")] / models$`f(m|y)`
  rownames(models) <- 1:dim(models)[1]

  ##
  # Bayes Factor of the consistent NMA model versus the inconsistent NMA model
  ##

  fm1y <- models$`f(m|y)`[which(models$IFs == "No IFs")]
  fm2y <- 1 - fm1y
  fm1 <- result$Estimate[which(result$Parameter == "p.cons")]
  fm2 <- 1 - fm1
  ##
  BF_m1m2 <- round((fm1y * fm2) / (fm2y * fm1), digits = digits)

  # Round results
  for (i in 2:4) {
    result[, i] <- round(result[, i], digits = digits)
    models[, i] <- round(models[, i], digits = digits)
  }

  # Change the names according to inconsistency factors names
  result[grep("beta", result$Parameter), "Parameter"] <- paste("Effect of", IF, sep = " ")
  result[grep("gamma", result$Parameter), "Parameter"] <- paste("P.I.P of", IF, sep = " ")

  res <- list(
    "Summary" = result,
    "PIP" = IF_PIP,
    "PO" = models,
    "BF" = BF_m1m2
  )
}
