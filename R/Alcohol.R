#' Stochastic Search Inconsistency Factor Selection of brief alcohol interventions.
#'
#' @description
#' Stochastic Search Inconsistency Factor Selection for the evaluation of the consistency
#' assumption for the network meta-analysis model.
#'
#' These data are used as an example in Seitidis et al. (2021).
#'
#' @name Alcohol
#'
#' @docType data
#'
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{studyid}}\tab study id \cr
#' \bold{\emph{treat1}}\tab treatment 1 \cr
#' \bold{\emph{treat2}}\tab treatment 2 \cr
#' \bold{\emph{m1}}\tab mean value of brief alcohol intervention in arm 1 \cr
#' \bold{\emph{m2}}\tab mean value of brief alcohol intervention in arm 2 \cr
#' \bold{\emph{n1}}\tab number of individuals in arm 1 \cr
#' \bold{\emph{n2}}\tab number of individuals in arm 2 \cr
#' \bold{\emph{sd1}}\tab standard deviation of brief alcohol intervention in arm 1 \cr
#' \bold{\emph{sd2}}\tab standard deviation of brief alcohol intervention in arm 2 \cr
#' \bold{\emph{TE}}\tab standardized mean difference of treat1 versus treat2 \cr
#' \bold{\emph{seTE}}\tab standard error of standardized mean difference \cr
#' }
#'
#' @source
#' Seitidis G, Nikolakopoulos S, Hennessy EA, Tanner-Smith EE, Mavridis D (2021):
#' Network Meta-Analysis Techniques for Synthesizing Prevention Science Evidence
#' \emph{Prevention Science},
#' 1-10
#'
#'
#' @examples
#' data(Alcohol)
#'
#' TE <- Alcohol$TE
#' seTE <- Alcohol$seTE
#' studlab <- Alcohol$studyid
#' treat1 <- Alcohol$treat2
#' treat2 <- Alcohol$treat1
#'
#' # Stochastic Search Inconsistency Factor Selection using as reference treatment AO-CT and the
#' # design-by-treatment method for the specification of the Z matrix.
#'
#' m <- ssifs(TE, seTE, treat1, treat2, studlab, ref = "AO-CT",
#' M = 1000, B = 100, M_pilot = 1000, B_pilot = 100)
#'
NULL
