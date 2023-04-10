#' Inconsistency Factors' Spike and Slab
#'
#' @description
#' The function visualizes the inconsistency factor's effect when the inconsistency
#' factor is included in the Network Meta-Analysis (NMA) model and when is not.
#'
#' @details
#' The function creates two density plots for each inconsistency factor based on the inconsistency
#' factors' effects, which are obtained from the \code{ssifs} model. The former visualizes the effect when the
#' inconsistency factor is included in the NMA model (spike), while the latter when
#' is not (slab). A good mixing of the SSIFS model indicates that the spike has high density for values
#' close to zero whereas the slab is flatter.
#'
#' @param x An object of class \code{ssifs}.
#'
#' @return An object of class \code{ggplot}.
#' @importFrom ggplot2 ggplot aes %+% geom_density xlab ylab ggtitle guides guide_legend theme element_text
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(Alcohol)
#'
#' TE <- Alcohol$TE
#' seTE <- Alcohol$seTE
#' studlab <- Alcohol$studyid
#' treat1 <- Alcohol$treat2
#' treat2 <- Alcohol$treat1
#'
#' # Stochastic Search Inconsistency Factor Selection using intervention AO-CT as reference.
#' m <- ssifs(TE, seTE, treat1, treat2, studlab, ref = "AO-CT")
#' spike.slab(m)
#' }
#'
spike.slab <- function(x) {
  if (inherits(x, "ssifs") == FALSE) {
    stop("The class of x is not of ssifs", call. = FALSE)
  }

  # Get the inconsistency factors effects
  MCMC_draws <- x$MCMC_run$BUGSoutput$sims.matrix

  betas <- MCMC_draws[, grep("beta", colnames(MCMC_draws))] # Inconsistency factors effects
  gammas <- MCMC_draws[, grep("gamma", colnames(MCMC_draws))] # Posterior inclusion probabilities
  gammas <- plyr::mapvalues(gammas, from = c(0, 1), to = c("Not included", "Included"))

  # Inconsistency Factors names
  IF_names <- x$Posterior_inclusion_probabilities$Comparison
  if (sum(is.na(x$Posterior_inclusion_probabilities$Design)) == 0) {
    IF_names <- paste(x$Posterior_inclusion_probabilities$Comparison, x$Posterior_inclusion_probabilities$Design, sep = "_")
  }

  if (!is.null(dim(betas)[2])) {
    p <- plotden(betas, gammas, FALSE, IF_names)
  } else {
    p <- plotden(betas, gammas, TRUE, IF_names)
  }

  p
}
