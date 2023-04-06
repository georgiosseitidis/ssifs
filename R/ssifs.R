#' Stochastic Search Inconsistency Factor Selection
#'
#' @description
#' Stochastic Search Inconsistency Factor Selection evaluates the consistency assumption of Network Meta-Analysis in the Bayesian framework,
#' by treating the inconsistency detection as a variable selection problem. The consistency assumption is evaluated locally, but also globally.
#'
#' @param TE Estimate of treatment effect (e.g. log odds ratio, mean difference, or log hazard ratio).
#' @param seTE Standard error of the treatment estimate.
#' @param treat1 Label/Number of the first treatment.
#' @param treat2 Label/Number of the second treatment.
#' @param studlab Study labels.
#' @param ref Reference treatment.
#' @param method Method used for the specification of the inconsistency factors; Possible choices are: \itemize{
#' \item \code{"LuAdes"} for the Lu & Ades model (Lu & Ades, 2006)
#' \item \code{"DBT"} for the design-by-treatment method (Higgins et al., 2012)
#' \item \code{"Jackson"} for the random-effects implementation of the design-by-treatment model (Jackson et al., 2014)}
#' @param rpcons \code{"logical"}. If \code{TRUE}, an informative beta distribution Beta(157, 44) is used for the probability to have a consistent network.
#' @param pcons Probability to have a consistent network.
#' @param zellner \code{"logical"}. If \code{TRUE}, Zellner g-prior is used for the dependency between the inconsistency factors. If \code{FALSE}, inconsistency factors assumed independent.
#' @param c Tuning parameter.
#' @param psi Tuning parameter.
#' @param digits Digits of the exported results.
#' @param M Number of NMA MCMC iterations.
#' @param B Burn-in period of the NMA MCMC run.
#' @param n_thin Thinning interval of the NMA MCMC run.
#' @param n_chains Number of parallel chains for the NMA MCMC run.
#' @param M_pilot Number of pilot MCMC iterations.
#' @param B_pilot Burn-in period of the pilot MCMC run.
#' @param n_thin_pilot Thinning interval of the pilot MCMC run.
#' @param n_chains_pilot Number of parallel chains for the pilot MCMC run.
#'
#'
#' @details
#' Stochastic Search Inconsistency Factor Selection (SSIFS) is the extension of Stochastic Search Variable Selection (SSVS)
#' (George & McCulloch, 1993) for identifying inconsistencies in Network Meta-Analysis (NMA).
#'
#' SSIFS (Seitidis et al., 2022), is a two-step method in which the inconsistency factors are specified in the first step, and in the second step,
#' SSVS is performed on the inconsistency factors. The method used for the specification of the inconsistency factors, is
#' controlled by the argument \code{method}. Among the choices that may be considered are the Lu and Ades model (Lu & Ades, 2006),
#' the design-by-treatment model (Higgins et al., 2012), and the random-effects implementation of the design-by-treatment model (Jackson et al., 2014).
#'
#' After specifying the inconsistency factors, the random-effects NMA model is implemented in the Bayesian framework using the \code{R2jags} package.
#' An uninformative normal is assumed for the prior distribution of the treatment effects, while for the heterogeneity parameter tau, an uninformative half-normal is assumed.
#' The function provides the MCMC run of the NMA model (item \code{MCMC_run}), whereby the user can check the convergence of the MCMC run.
#'
#' SSIFS by default assumes that inconsistency factors are dependent by using a Zellner g-prior to describe this dependency (\code{zellner = TRUE}).
#' Parameter g in the Zellner g-prior is specified using the unit information criterion (Kass & Wasserman, 1995), which is translated in
#' SSIFS to the total number of observed comparisons that are included in the network. By setting the argument \code{zellner = FALSE}, inconsistency factors are assumed independent.
#' Regarding the inclusion probabilities, the function by default assumes an informative Beta distribution (Beta(157, 44)) for the probability to have a
#' consistent network (\code{rpcons = TRUE}). In the case where \code{rpcons = FALSE}, this probability is assumed fixed and equal to 0.5 (\code{pcons = 0.5}).
#' The user can modify this probability through the argument \code{pcons}.
#'
#' Tuning parameters in SSIFS are specified by the arguments \code{c} and \code{psi}. They should be specified in a way that, when
#' an inconsistency factor is included in the NMA model, the corresponding coefficient lies in an area close to zero,
#' and far away from this area when it is not included in the NMA model. Regarding the argument \code{c}, values between 10 and 100
#' usually perform well in most cases. Argument \code{psi} can be obtained either from a pilot MCMC run of the NMA model, as the
#' standard deviation of the inconsistency factors (\code{psi = NULL}), or can be set fixed a-priory by the analyst.
#'
#' In order to evaluate the consistency assumption, we can examine \itemize{
#' \item {the posterior inclusion probabilities
#' of the inconsistency factors (item \code{Posterior_inclusion_probabilities})}
#' \item{the posterior model probabilities (item \code{Posterior_Odds})}
#' \item{the posterior model odds (item \code{Posterior_Odds})}
#' \item{the Bayes factor of the consistent NMA model over the inconsistent NMA models (item \code{Bayes_Factor})} }
#' A posterior inclusion probability above 0.5 indicates inconsistency.
#' Also, an inconsistent NMA model with large posterior model probability suggests the presence of inconsistency. Item \code{Bayes_Factor}
#' provides a global test for testing the consistency assumption, by calculating the Bayes factor of the consistent NMA model (model without inconsistency factors) over the
#' rest inconsistent NMA models that  were observed in the MCMC run. An estimate above 1 favors the consistent NMA model.
#' For the calculation of the Bayes factor, the inconsistent NMA models are treated as a single model and the corresponding
#' posterior model probabilities are summed.
#'
#' @note
#' The function uses the random effects inverse-variance NMA model, and assumes common heterogeneity between different treatment comparisons
#' and no correlation between different studies. Also note that the function keeps only those studies that belong to
#' the largest sub-network, if the network is disconnected, in order to maintain one connected network.
#'
#' In a multi-arm study with T comparisons, T-1 are required for the NMA model since the rest are obtained as a linear combination.
#' The function automatically excludes the unnecessary comparisons, while maintaining the basic comparisons (if possible)
#' and comparisons in which an inconsistency factor has been added (if possible). Therefore, all possible comparisons
#' of multi-arm studies must be provided by the user.
#'
#' For the names of the inconsistency factor, treatments are separated by the symbol \code{" ; "}. For example, if an inconsistency
#' factor is added in the comparison between treatments A and B, the inconsistency factor name will be \code{A ; B}.
#' In the case where \code{method = "DBT"} or \code{method = "Jackson"}, inconsistency factors' names for multi-arm designs are
#' denoted as \emph{treatment.comparison}_\emph{design}. Thus, if an inconsistency factor is added in the comparison
#' between treatments A and B of the ABC design, the inconsistency factor name will be \code{A ; B_ABC}.
#'
#' In extremely large networks potentially there are exponentially many paths between two nodes, and you may run out of
#' memory when using the Lu and Ades model for the specification of the inconsistency factors, if your network is lattice-like.
#'
#' @importFrom Rdpack reprompt
#' @importFrom utils combn
#' @importFrom plyr mapvalues
#' @importFrom R2jags jags
#'
#' @return
#' A list containing the following components:
#' \item{MCMC_run}{An object of class \code{rjags} containing the MCMC run of the NMA model}
#' \item{Bayes_Factor}{Bayes factor of the consistent NMA model over the inconsistent NMA model}
#' \item{Posterior_inclusion_probabilities}{A \code{data.frame} containing the posterior inclusion
#' probabilities of the inconsistency factors. Columns \code{Comparison} and \code{Design} denote in which comparisons
#' inconsistency factors are added. When argument \code{method = "LuAdes"}, column \code{Design} is \code{NA},
#' because only loop inconsistencies are accounted.
#' \code{PIP} is the estimated posterior inclusion probability, \code{b} is the estimated median effect
#'  of the inconsistency factors, \code{b.lb} and \code{b.ub} the lower and upper bound of inconsistency factors' effect estimates,
#'  respectively}
#' \item{Posterior_Odds}{A \code{data.frame} containing the model posterior odds. Column \code{IFs} denotes
#'  in which comparisons inconsistency factors are added, \code{Freq} the number of times the model is observed in
#'  the MCMC run, \code{f(m|y)} the posterior model probability and \code{PO_IFCONS} the posterior model odds
#'  of the consistent NMA model (denoted as \code{NO IFs}) over the corresponding inconsistent NMA model.}
#' \item{Summary}{A \code{data.frame} containing the summary estimates of the MCMC run of the NMA model}
#' \item{Z_matrix}{A \code{data.frame} containing in the first 3 columns the treatment comparisons used for
#'  the Z matrix and the Z matrix in the rest columns}
#' \item{disconnected_studies}{A vector with the studies that were excluded, in order to have
#'  a connected network}
#' \item{n_subnetworks}{Number of sub-networks}
#' \item{subnetworks}{A list with the sub-networks of the original NMA data}
#'
#' @references
#' George, E. I., & McCulloch, R. E. (1993):
#' Variable selection via Gibbs sampling.
#' \emph{Journal of the American Statistical Association},
#' \bold{88}(423), 881-889.
#'
#' Seitidis, G., Nikolakopoulos, S., Ntzoufras, I., & Mavridis, D. (2022):
#' Inconsistency identification in network meta-analysis via stochastic search variable selection.
#' \emph{arXiv preprint},
#' \bold{arXiv:2211.07258}.
#'
#' Lu, G., & Ades, A. E. (2006):
#' Assessing evidence inconsistency in mixed treatment comparisons.
#' \emph{Journal of the American Statistical Association},
#' \bold{101}(474), 447-459.
#'
#' Higgins, J. P. T., Jackson, D., Barrett, J. K., Lu, G., Ades, A. E., & White, I. R. (2012):
#' Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies.
#' \emph{Research synthesis methods},
#' \bold{3}(2), 98-110.
#'
#' Jackson, D., Barrett, J. K., Rice, S., White, I. R., & Higgins, J. P. (2014):
#' A design-by-treatment interaction model for network meta-analysis with random inconsistency effects.
#' \emph{Statistics in medicine},
#' \bold{33}(21), 3639-3654.
#'
#' Kass, R. E., & Wasserman, L. (1995):
#' A Reference Bayesian Test for Nested Hypotheses and its Relationship to the Schwarz Criterion.
#' \emph{Journal of the American Statistical Association},
#' \bold{90}(431), 928–934.
#'
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
#' }
#'
ssifs <- function(TE, seTE, treat1, treat2, studlab, ref, method = "DBT", rpcons = TRUE, pcons = 0.5, zellner = TRUE, c = 3, psi = NULL, digits = 4,
                  M = 50000, B = 10000, n_thin = 1, n_chains = 2, M_pilot = 10000, B_pilot = 2000, n_thin_pilot = 1, n_chains_pilot = 1) {

  ##
  # Check arguments
  ##

  if (inherits(TE, "numeric") == FALSE) {
    stop("The class of TE must be numeric", call. = FALSE)
  } else if (inherits(seTE, "numeric") == FALSE) {
    stop("The class of seTE must be numeric", call. = FALSE)
  } else if (inherits(treat1, c("character", "numeric", "integer")) == FALSE) {
    stop("The class of treat1 must be either character, or numeric, or integer", call. = FALSE)
  } else if (inherits(treat2, c("character", "numeric", "integer")) == FALSE) {
    stop("The class of treat2 must be either character, or numeric, or integer", call. = FALSE)
  } else if (inherits(studlab, c("character", "numeric", "integer")) == FALSE) {
    stop("The class of studlab must be either character, or numeric, or integer", call. = FALSE)
  } else if (sum(is.na(studlab)) > 0) {
    stop("studlab must not contain missing values", call. = FALSE)
  } else if (!(length(TE) == length(seTE) & length(TE) == length(treat1) &
    length(TE) == length(treat2) & length(TE) == length(studlab))) {
    stop("TE, seTE, treat1, treat2 and studlab must have the same length")
  } else if (inherits(ref, c("character", "numeric", "integer")) == FALSE) {
    stop("The class of ref must be either character, or numeric, or integer", call. = FALSE)
  } else if (!ref %in% c(treat1, treat2)) {
    stop("The reference category must be a node of the network", call. = FALSE)
  } else if (!method %in% c("DBT", "LuAdes", "Jackson")) {
    stop("method must be either DBT, or LuAdes, or Jackson", call. = FALSE)
  } else if (length(method) > 1) {
    stop("The length of method must be one", call. = FALSE)
  } else if (inherits(rpcons, "logical") == FALSE) {
    stop("The class of rpcons is not logical", call. = FALSE)
  } else if (length(rpcons) > 1) {
    stop("The length of rpcons must be one", call. = FALSE)
  } else if (rpcons == FALSE) {
    if (inherits(pcons, c("numeric", "integer")) == FALSE) {
      stop("The class of pcons must be either numeric, or integer", call. = FALSE)
    } else if (length(pcons) > 1) {
      stop("The length of pcons must be one", call. = FALSE)
    } else if (pcons > 1 | pcons < 0) {
      stop("pcons must be between 0-1", call. = FALSE)
    }
  } else if (inherits(zellner, "logical") == FALSE) {
    stop("The class of zellner is not logical", call. = FALSE)
  } else if (length(zellner) > 1) {
    stop("The length of zellner must be one", call. = FALSE)
  } else if (inherits(c, c("numeric", "integer")) == FALSE) {
    stop("The class of c must be either numeric, or integer", call. = FALSE)
  } else if (length(c) > 1) {
    stop("The length of c must be one", call. = FALSE)
  } else if (c <= 0) {
    stop("c must be a positive number", call. = FALSE)
  } else if (is.null(psi) == FALSE) {
    if (inherits(psi, c("numeric", "integer")) == FALSE) {
      stop("The class of psi must be either numeric, or integer", call. = FALSE)
    } else if (length(psi) > 1) {
      stop("The length of psi must be one", call. = FALSE)
    } else if (psi <= 0) {
      stop("psi must be a positive number", call. = FALSE)
    }
  } else if (inherits(digits, c("numeric", "integer")) == FALSE) {
    stop("The class of digits must be either numeric, or integer", call. = FALSE)
  } else if (length(digits) > 1) {
    stop("The length of digits must be one", call. = FALSE)
  } else if (digits < 0 | digits %% 1 != 0) {
    stop("digits must be a positive integer number", call. = FALSE)
  } else if (inherits(M, c("numeric", "integer")) == FALSE) {
    stop("The class of M must be either numeric, or integer", call. = FALSE)
  } else if (length(M) > 1) {
    stop("The length of M must be one", call. = FALSE)
  } else if (M <= 0 | M %% 1 != 0) {
    stop("M must be a positive integer number", call. = FALSE)
  } else if (inherits(B, c("numeric", "integer")) == FALSE) {
    stop("The class of B must be either numeric, or integer", call. = FALSE)
  } else if (length(B) > 1) {
    stop("The length of B must be one", call. = FALSE)
  } else if (B <= 0 | B %% 1 != 0) {
    stop("B must be a positive integer number", call. = FALSE)
  } else if (B >= M) {
    stop("M must be larger than B", call. = FALSE)
  } else if (inherits(n_thin, c("numeric", "integer")) == FALSE) {
    stop("The class of n_thin must be either numeric, or integer", call. = FALSE)
  } else if (length(n_thin) > 1) {
    stop("The length of n_thin must be one", call. = FALSE)
  } else if (n_thin <= 0 | n_thin %% 1 != 0) {
    stop("n_thin must be a positive integer number", call. = FALSE)
  } else if (n_thin >= B) {
    stop("n_thin must be smaller than B", call. = FALSE)
  } else if (inherits(n_chains, c("numeric", "integer")) == FALSE) {
    stop("The class of n_chains must be either numeric, or integer", call. = FALSE)
  } else if (length(n_chains) > 1) {
    stop("The length of n_chains must be one", call. = FALSE)
  } else if (n_chains < 1 | n_chains %% 1 != 0) {
    stop("n_chains must be a positive integer number larger than one", call. = FALSE)
  } else if (inherits(M_pilot, c("numeric", "integer")) == FALSE) {
    stop("The class of M_pilot must be either numeric, or integer", call. = FALSE)
  } else if (length(M_pilot) > 1) {
    stop("The length of M_pilot must be one", call. = FALSE)
  } else if (M_pilot <= 0 | M_pilot %% 1 != 0) {
    stop("M_pilot must be a positive integer number", call. = FALSE)
  } else if (inherits(B_pilot, c("numeric", "integer")) == FALSE) {
    stop("The class of B_pilot must be either numeric, or integer", call. = FALSE)
  } else if (length(B_pilot) > 1) {
    stop("The length of B_pilot must be one", call. = FALSE)
  } else if (B_pilot <= 0 | B_pilot %% 1 != 0) {
    stop("B_pilot must be a positive integer number", call. = FALSE)
  } else if (B_pilot >= M_pilot) {
    stop("M_pilot must be larger than B_pilot", call. = FALSE)
  } else if (inherits(n_thin_pilot, c("numeric", "integer")) == FALSE) {
    stop("The class of n_thin_pilot must be either numeric, or integer", call. = FALSE)
  } else if (length(n_thin_pilot) > 1) {
    stop("The length of n_thin_pilot must be one", call. = FALSE)
  } else if (n_thin_pilot <= 0 | n_thin_pilot %% 1 != 0) {
    stop("n_thin_pilot must be a positive integer number", call. = FALSE)
  } else if (n_thin_pilot >= B_pilot) {
    stop("n_thin_pilot must be smaller than B_pilot", call. = FALSE)
  } else if (inherits(n_chains_pilot, c("numeric", "integer")) == FALSE) {
    stop("The class of n_chains_pilot must be either numeric, or integer", call. = FALSE)
  } else if (length(n_chains_pilot) > 1) {
    stop("The length of n_chains_pilot must be one", call. = FALSE)
  } else if (n_chains_pilot < 1 | n_chains_pilot %% 1 != 0) {
    stop("n_chains_pilot must be a positive integer number larger than one", call. = FALSE)
  }


  ##
  # Prepare network's data
  ##

  treat1 <- as.character(treat1)
  treat2 <- as.character(treat2)
  ref <- as.character(ref)
  studlab <- as.character(studlab)
  data <- data.frame(TE = TE, seTE = seTE, studlab = studlab, treat1 = treat1, treat2 = treat2, stringsAsFactors = FALSE)

  # Exclude studies with NAs
  NAs <- apply(data, 1, FUN = function(x) {
    sum(is.na(x))
  })
  omit_stud <- data$studlab[NAs > 0]
  if (length(omit_stud) > 0) {
    data <- data[-which(data$studlab %in% omit_stud), ]
  }

  # Keep the studies that belong to the largest sub-network
  cnet <- connet(data)
  if (!ref %in% unlist(cnet$data[, c("treat1", "treat2")])) {
    stop("The reference category must a node of the largest sub-network", call. = FALSE)
  }
  data <- cnet$data

  # Check if the data are correctly specified
  correct.data(data)

  # Make all treatments names numeric values
  trt_orig <- unique(c(data$treat1, data$treat2))
  trt_trans <- data.frame("from" = trt_orig, "to" = 1:length(trt_orig))
  data$treat1 <- plyr::mapvalues(data$treat1, from = trt_trans$from, to = trt_trans$to, warn_missing = FALSE)
  data$treat2 <- plyr::mapvalues(data$treat2, from = trt_trans$from, to = trt_trans$to, warn_missing = FALSE)
  ref <- plyr::mapvalues(ref, from = trt_trans$from, to = trt_trans$to, warn_missing = FALSE)

  # NMA MCMC model requires multi-arm studies to be at the end of the dataframe
  Multiarm_studies <- labels(which(table(data$studlab) >= 3))
  if (length(Multiarm_studies) > 0) {
    data_multi <- data[which(data$studlab %in% Multiarm_studies), ]
    data <- rbind(data[-which(data$studlab %in% Multiarm_studies), ], data_multi[order(data_multi$studlab), ])
  }

  ##
  # Generate NMA MCMC data
  ##

  if (method == "LuAdes") {
    NMA_data <- LuAdes.data(data, ref, Multiarm_studies)
  } else {
    NMA_data <- DBT.data(data, ref, Multiarm_studies, method)
  }

  # Data used for NMA MCMC run
  Data_MCMC <- list(
    N = NMA_data$N,
    NHtH = NMA_data$NHtH,
    p = NMA_data$p,
    NT = length(trt_orig),
    ref = ref,
    ##
    y = NMA_data$y,
    w = 1 / (data$seTE[1:NMA_data$NHtH]^2),
    t = NMA_data$Z_MCMC$treat1,
    b = NMA_data$Z_MCMC$treat2,
    ##
    x = NMA_data$Z,
    H = NMA_data$H,
    PREC = NMA_data$PREC,
    ZTZ = NMA_data$ZTZ
  )

  ##
  # Additional arguments of the NMA MCMC run
  ##

  if (length(Multiarm_studies) > 0) {
    multi <- TRUE
  } else {
    multi <- FALSE
    Data_MCMC <- Data_MCMC[-which(names(Data_MCMC) %in% c("H", "PREC"))]
  }

  ##
  par_save <- c("precision", "beta", "gamma", "sigma", "p.cons", "sd", "d") # parameteres saved in the MCMC

  if (NMA_data$p == 1) {
    one_IF <- TRUE
    Data_MCMC$x <- as.vector(Data_MCMC$x)
    par_save <- par_save[-which(par_save == "sigma")]

    if (zellner) {
      Data_MCMC <- Data_MCMC[-which(names(Data_MCMC) == "ZTZ")]
    }
  } else {
    one_IF <- FALSE
  }

  if (!zellner) {
    Data_MCMC <- Data_MCMC[-which(names(Data_MCMC) == "ZTZ")]

    if ("sigma" %in% par_save) {
      par_save <- par_save[-which(par_save == "sigma")]
    }
  }
  ##
  # Tunning
  ##

  if (is.null(psi) == TRUE) {

    # Pilot NMA MCMC run
    Data_MCMC_pilot <- Data_MCMC
    if (one_IF) {
      Data_MCMC_pilot <- Data_MCMC_pilot[-which(names(Data_MCMC_pilot) == "p")]
    }

    model <- R2jags::jags(
      model.file = textConnection(NMApilot(multi, one_IF, zellner)),
      n.iter = M_pilot, n.burnin = B_pilot, n.thin = n_thin_pilot, n.chains = n_chains_pilot,
      parameters.to.save = c("beta"), data = Data_MCMC_pilot
    )
    ##
    Data_MCMC[["psi"]] <- model$BUGSoutput$summary[-which(rownames(model$BUGSoutput$summary) == "deviance"), 2]
  } else {
    Data_MCMC[["psi"]] <- rep(psi, Data_MCMC$p)
  }
  ##
  Data_MCMC[["c"]] <- c

  ##
  # NMA MCMC run
  ##

  m_NMA <- NMAmodel(multi, one_IF, rpcons, pcons, zellner)

  ##
  model <- R2jags::jags(
    model.file = textConnection(m_NMA), n.iter = M, n.burnin = B, n.thin = n_thin, n.chains = n_chains,
    parameters.to.save = par_save, DIC = FALSE, data = Data_MCMC
  )

  ##
  # Export results
  ##

  # Transform numerical comparisons to the original comparison's names
  IF <- transformIF(NMA_data$IF, trt_trans, method)
  ##
  Z_matrix <- as.data.frame(NMA_data$Z)
  Z_t1 <- plyr::mapvalues(NMA_data$Z_MCMC$treat1, from = trt_trans$to, to = trt_trans$from, warn_missing = FALSE)
  Z_t2 <- plyr::mapvalues(NMA_data$Z_MCMC$treat2, from = trt_trans$to, to = trt_trans$from, warn_missing = FALSE)
  Z_matrix <- data.frame("treat1" = Z_t1, "treat2" = Z_t2, "studlab" = NMA_data$Z_MCMC$studlab, Z_matrix)
  colnames(Z_matrix)[4:dim(Z_matrix)[2]] <- IF
  ##
  Results <- jags.summary(model, trt_trans, Data_MCMC, IF, digits, method)


  res <- list(
    "MCMC_run" = model,
    "Bayes_Factor" = Results$BF,
    "Posterior_inclusion_probabilities" = Results$PIP,
    "Posterior_Odds" = Results$PO,
    "Summary" = Results$Summary,
    "Z_matrix" = Z_matrix,
    "disconnected_studies" = cnet$excluded_studies,
    "n_subnetworks" = cnet$n_sub,
    "subnetworks" = cnet$subnet
  )

  class(res) <- "ssifs"

  res
}
