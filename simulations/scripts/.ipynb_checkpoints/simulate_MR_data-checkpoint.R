simulate_MR_data <- function(
    cluster_sizes = c(10, 20, 40),  # variants per causal cluster
    theta         = c(0.4, -0.4, 0.8),  # causal effect per cluster
    N             = 10000,              # sample size
    tau           = 1,                  # overdispersion parameter
    h2            = 0.5,                # heritability of the exposure
    null_clust    = FALSE,              # include a null cluster (theta = 0)?
    null_size     = 10,                 # number of null instruments
    junk_clust    = FALSE,              # include a junk cluster (pleiotropic, heavy-tailed)?
    junk_size     = 10,                 # number of junk instruments
    df            = 4,                  # degrees of freedom for the junk t-distribution
    seed          = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  K        <- length(cluster_sizes)
  J_causal <- sum(cluster_sizes)
  h2_j   <- h2 / (J_causal + null_clust * null_size + junk_clust * junk_size)
  
  stopifnot(length(theta) == K)

  # Simulate bx, by, byse for J instruments with per-instrument true ratio theta_j
  .sim_block <- function(J, theta_j, h2_j) {
    MAF  <- runif(J, 0.05, 0.5)
    gamma <- rnorm(J, mean = 0, sd = sqrt(h2_j / (2 * MAF * (1 - MAF))))
    bxse <- sqrt(1 / (2 * N * MAF * (1 - MAF)))
    byse <- sqrt(tau / (2 * N * MAF * (1 - MAF)))
    bx   <- rnorm(J, mean = gamma, sd = bxse)
    by   <- rnorm(J, mean = theta_j * gamma,   sd = byse)
    list(bx = bx, by = by, bxse = bxse, byse = byse, MAF = MAF)
  }

  # --- Causal instruments ---
  labels_causal  <- rep(seq_len(K), times = cluster_sizes)
  causal         <- .sim_block(J_causal, theta[labels_causal], h2_j = h2_j)

  # --- Null instruments (true causal effect = 0) ---
  # by ~ N(0, byse): the instrument affects the exposure but not the outcome.
  null_block <- if (null_clust) .sim_block(null_size, rep(0, null_size), h2_j = h2_j) else NULL

  # --- Junk instruments (direct pleiotropic effect ~ t(df)) ---
  # The ratio estimate is heavy-tailed: theta_hat ~ N(alpha_j, sigma_j)
  # where alpha_j ~ t(df), mimicking the Stan junk-cluster likelihood.
  junk_block <- if (junk_clust) {
    alpha_j <- rt(junk_size, df = df)
    .sim_block(junk_size, alpha_j, h2_j = h2_j)
  } else NULL

  # --- Combine in Stan column order: [null | causal | junk] ---
  bx   <- c(null_block$bx,   causal$bx,   junk_block$bx)
  by   <- c(null_block$by,   causal$by,   junk_block$by)
  bxse <- c(null_block$bxse, causal$bxse, junk_block$bxse)
  byse <- c(null_block$byse, causal$byse, junk_block$byse)
  MAF  <- c(null_block$MAF,  causal$MAF,  junk_block$MAF)

  # Cluster labels: 0 = null, 1:K = causal, K+1 = junk
  labels <- c(
    if (null_clust)  rep(0L,     null_size),
    labels_causal,
    if (junk_clust)  rep(K + 1L, junk_size)
  )

  list(
    bx       = bx,
    by       = by,
    bxse     = bxse,
    byse     = byse,
    cluster  = labels,
    MAF      = MAF,
    theta    = theta,
    cluster_sizes = cluster_sizes,
    K        = K,
    J        = length(bx),
    N = N,
    num_null = as.integer(null_clust),
    num_junk = as.integer(junk_clust),
    null_size = if (null_clust) null_size else 0L,
    junk_size = if (junk_clust) junk_size else 0L
  )
}
