plotden <- function(betas, gammas, oneIF, IF_names) {
  p <- list()
  x <- y <- NULL

  if (oneIF) {
    plot_data <- data.frame("x" = betas, "y" = gammas)
    p_i <- ggplot(data = plot_data) +
      geom_density(aes(x = x, fill = y), alpha = 0.5) +
      xlab("Effect") +
      ylab("Density") +
      ggtitle(paste0("Spike and Slab for ", IF_names)) +
      guides(fill = guide_legend(title = "Inconsistency Factor")) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      )

    p[[IF_names]] <- p_i
  } else {
    colnames(betas) <- paste0("beta_", 1:dim(betas)[2])
    colnames(gammas) <- paste0("gamma_", 1:dim(betas)[2])

    for (i in 1:dim(betas)[2]) {
      plot_data <- data.frame("x" = betas[, i], "y" = gammas[, i])
      p_i <- ggplot(data = plot_data) +
        geom_density(aes(x = x, fill = y), alpha = 0.5) +
        xlab("Effect") +
        ylab("Density") +
        ggtitle(paste0("Spike and Slab for ", IF_names[i])) +
        guides(fill = guide_legend(title = "Inconsistency Factor")) +
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)
        )

      p[[IF_names[i]]] <- p_i
    }
  }
  p
}
