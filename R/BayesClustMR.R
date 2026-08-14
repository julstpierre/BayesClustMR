# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------
#' Fit the Stan model for a fixed number of clusters K
#'
#' @param theta  Vector of ratio estimates (length J)
#' @param sigma  Vector of ratio SEs (length J)
#' @param K      Number of causal clusters
#' @param etas   J x (K + num_null + num_junk) prior weight matrix
#' @param num_null,num_junk  Cluster-type flags
#' @param df     Degrees of freedom (passed to Stan data)
#'
#' @return A list with `df` (one row per SNP, columns: snp, cluster, post_prob,
#'   cluster_size, cluster_mean, cluster_se, cluster_lcl, cluster_ucl, theta_mean, theta_se)
#'   and `waic`
.fit_model <- function(theta, sigma, K, etas, num_null, num_junk, df = 4) {
  J <- length(theta)

  data_list <- list(
    J              = J,
    K              = K,
    rho            = rep(0, K),
    phi            = rep(1, K),
    theta_estimates = theta,
    sigma_estimates = sigma,
    num_null       = num_null,
    num_junk       = num_junk,
    deg_freedom    = df,
    etas           = etas
  )

  stan_file <- here::here("inst", "stan", "BayesClustMR.stan")
  model     <- cmdstanr::cmdstan_model(stan_file)

  base_init <- initialize_centers_weighted(theta, sigma, K)
  init_fun  <- function() list(cluster_center = base_init + rnorm(K, 0, 1e-3))

  fit <- model$sample(
    init          = init_fun,
    data          = data_list,
    iter_sampling = 1000,
    iter_warmup   = 1000,
    chains        = 4
  )

  log_lik     <- fit$draws("log_lik", format = "matrix")
  waic_result <- loo::waic(log_lik)
  waic_val    <- waic_result$estimates["waic", "Estimate"]

  # Cluster center summaries for the K causal clusters
  center_summary <- fit$summary(
    variables = paste0("cluster_center[", 1:K, "]"),
    "mean", "sd",
    ~quantile(.x, probs = c(0.025, 0.975))
  )

  # Build a lookup table for cluster center stats indexed by Stan cluster number.
  # Stan ordering: [null (opt)] [causal 1..K] [junk (opt)]
  # causal cluster k -> Stan index (num_null + k)
  num_clusters  <- K + num_null + num_junk
  cluster_mean  <- numeric(num_clusters)   # 0 for null/junk (fixed or undefined)
  cluster_se    <- rep(NA_real_, num_clusters)
  cluster_lcl   <- rep(NA_real_, num_clusters)
  cluster_ucl   <- rep(NA_real_, num_clusters)

  causal_idx <- seq_len(K) + num_null      # Stan cluster indices for causal clusters
  cluster_mean[causal_idx] <- center_summary$mean
  cluster_se[causal_idx]   <- center_summary$sd
  cluster_lcl[causal_idx]  <- center_summary$`2.5%`
  cluster_ucl[causal_idx]  <- center_summary$`97.5%`

  # Posterior cluster probabilities (J x num_clusters), averaged over draws
  pp_vars   <- outer(1:J, 1:num_clusters, function(j, k) sprintf("post_prob[%d,%d]", j, k))
  pp_draws  <- fit$draws(as.vector(pp_vars), format = "matrix")
  post_prob <- matrix(colMeans(pp_draws), nrow = J, ncol = num_clusters)

  # MAP cluster allocation and its posterior probability
  cluster_alloc <- apply(post_prob, 1, which.max)
  cluster_pp    <- post_prob[cbind(seq_len(J), cluster_alloc)]

  # Cluster sizes based on MAP allocation
  size_tbl      <- tabulate(cluster_alloc, nbins = num_clusters)

  # Cluster class labels: Stan ordering is [null (opt)] [causal 1..K] [junk (opt)]
  cluster_class_lut <- c(
    if (num_null) "Null",
    seq_len(K),
    if (num_junk) "Junk"
  )

  df_out <- data.frame(
    snp           = seq_len(J),
    cluster       = cluster_alloc,
    post_prob     = cluster_pp,
    cluster_size  = size_tbl[cluster_alloc],
    cluster_mean  = cluster_mean[cluster_alloc],
    cluster_se    = cluster_se[cluster_alloc],
    cluster_lcl   = cluster_lcl[cluster_alloc],
    cluster_ucl   = cluster_ucl[cluster_alloc],
    theta         = theta,
    theta_se      = sigma,
    cluster_class = cluster_class_lut[cluster_alloc]
  )

  list(results = df_out, waic = waic_val)
}


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------
#' Fit BayesClustMR
#'
#' @param bx,by,byse  Genetic summary statistics
#' @param etas        Optional J x M prior weight matrix. If NULL, the number
#'                    of clusters is selected automatically via `determine_K`.
#' @param num_null    1 if a null cluster is present, 0 otherwise
#' @param num_junk    1 if a junk cluster is present, 0 otherwise
#' @param df          Degrees of freedom for the junk distribution
#' @param K_max       Maximum number of causal clusters to consider
#'
#' @return A list with `center_estimates`, `center_ci` (K x 2 matrix of 95% credible
#'   intervals), and `waic`
BayesClustMR <- function(
    bx, by, byse,
    etas    = NULL,
    num_null = 0,
    num_junk = 0,
    df = 4,
    K_max = 10
) {
  # Data preparation happens exactly once
  dat <- list(theta = by / bx, sigma = byse / abs(bx), J = length(by))

  if (is.null(etas)) {
    K_result <- determine_K(
      dat$theta, dat$sigma,
      num_null = num_null,
      num_junk = num_junk,
      K_max = K_max
    )
    # determine_K already fitted the model for every K and kept the best result;
    return(K_result$best_fit)
  }

  K <- ncol(etas) - num_null - num_junk
  if (nrow(etas) != dat$J) stop("etas must have J rows")

  .fit_model(dat$theta, dat$sigma, K, etas, num_null, num_junk, df)
}
