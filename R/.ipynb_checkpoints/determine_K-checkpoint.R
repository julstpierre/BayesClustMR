#' Select the number of clusters K by WAIC
#'
#' Fits the model for each K in 2:K_max on already-prepared data (ratio
#' estimates and their SEs) and returns the K with the lowest WAIC.
#'
#' @param theta    Vector of ratio estimates (output of `.prepare_data()`)
#' @param sigma    Vector of ratio SEs
#' @param num_null 1 if a null cluster is present, 0 otherwise
#' @param num_junk 1 if a junk cluster is present, 0 otherwise
#' @param K_max    Maximum number of causal clusters to consider
#'
#' @return A list with `best_K` (integer), `waic_all` (named numeric vector),
#'   and `best_fit` (the `.fit_model` result for the selected K)
determine_K <- function(theta, sigma, num_null, num_junk, K_max = 5) {
  J <- length(theta)

  fits <- lapply(1:K_max, function(k) {
    etas <- generate_priors(
      K        = k,
      J        = J,
      num_null = num_null,
      num_junk = num_junk
    )
    .fit_model(theta, sigma, k, etas, num_null, num_junk)
  })

  waic_vals        <- vapply(fits, `[[`, numeric(1), "waic")
  names(waic_vals) <- paste0("K=", 1:K_max)

  best_idx <- which.min(waic_vals)

  list(
    best_K   = best_idx + 1,
    waic_all = waic_vals,
    best_fit = fits[[best_idx]]
  )
}
