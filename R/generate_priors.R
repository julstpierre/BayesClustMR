#' Generate a uniform eta prior weight matrix
#'
#' Returns a J × M matrix where every entry is 1/M, giving each instrument
#' equal prior weight across all M cluster columns. This is the only prior
#' appropriate for real (non-simulation) data; informative priors for
#' simulation studies are constructed directly in the simulation scripts.
#'
#' Column layout (M = K + num_null + num_junk):
#'   [null (if any)] | [K causal clusters] | [junk (if any)]
#'
#' @param K        Number of causal clusters.
#' @param J        Total number of instruments.
#' @param num_null 1 if a null cluster column is included, 0 otherwise.
#' @param num_junk 1 if a junk cluster column is included, 0 otherwise.
#'
#' @return A J × M numeric matrix with every entry equal to 1/M.
generate_priors <- function(K, J, num_null = 0, num_junk = 0) {
  M <- K + num_null + num_junk
  matrix(1 / M, nrow = J, ncol = M)
}
