#' Initialise cluster centers using precision-weighted quantiles
#'
#' Places K initial cluster centers at evenly-spaced quantiles of the ratio
#' estimate distribution, weighted by precision (1 / sigma^2). The centers
#' are nudged apart by a small epsilon if any two coincide, so that the
#' ordered-centers constraint in the Stan model is satisfied at initialisation.
#'
#' Called exclusively from `.fit_model()` in `fit_BayesClustMR.R`.
#'
#' @param theta  Numeric vector of ratio estimates (length J).
#' @param sigma  Numeric vector of ratio standard errors (length J).
#' @param K      Number of cluster centers to place.
#' @param jitter_sd  Not currently used; kept for API stability.
#'
#' @return Numeric vector of length K with ordered initial center values.
initialize_centers_weighted <- function(theta, sigma, K, jitter_sd = 1e-3) {
  w <- 1 / (sigma^2 + 1e-12)
  probs <- seq(0, 1, length.out = K + 2)[2:(K + 1)]
  centers <- as.numeric(weighted_quantile(theta, w, probs))
  
  # Align with ordered centers constraint in the Stan model.
  # The loop is only meaningful for K >= 2.
  if (K >= 2) {
    for (k in 2:K) {
      if (centers[k] <= centers[k - 1]) {
        centers[k] <- centers[k - 1] + 1e-4
      }
    }
  }
  
  centers
}

#' Precision-weighted empirical quantiles (internal helper)
#'
#' Computes empirical quantiles of `x` using normalized weights `w`. Used only
#' by `initialize_centers_weighted()`.
#'
#' @param x     Numeric vector of values.
#' @param w     Numeric vector of non-negative weights (need not sum to 1).
#' @param probs Numeric vector of probabilities in [0, 1].
#'
#' @return Numeric vector of the same length as `probs`.
weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  })
}