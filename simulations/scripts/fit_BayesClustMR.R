fit_BayesClustMR <- function(sim, prior) {
  
  num_null <- sim$num_null
  num_junk <- sim$num_junk
  
  if (prior == "none"){
    etas <- NULL
  } else {  
    # Generate priors
    M        <- sim$K + num_null + num_junk
    
    if (prior == "uniform") {
        true_cluster_weight <- 1 / M
    } else if (prior == "right_diffuse") {
        true_cluster_weight <- 0.5
    } else if (prior == "right_confident") {
        true_cluster_weight <- 0.8
    } else if (prior == "wrong_diffuse") {
        true_cluster_weight <- 0.5
    } else if (prior == "wrong_confident") {
        true_cluster_weight <- 0.8
    }
    
    cluster_rank <- rank(sim$theta)
    false_weight <- if (M > 1) (1 - true_cluster_weight) / (M - 1) else 0
    etas         <- matrix(false_weight, nrow = sim$J, ncol = M)
    
    for (j in seq_len(sim$J)) {
      if ( prior  %in% c("uniform", "right_diffuse", "right_confident") ) {
         label <- sim$cluster[j]
      } else if ( prior %in% c("wrong_diffuse", "wrong_confident") ) {
         label <- sample(setdiff(unique(sim$cluster), sim$cluster[j]), 1)
      }
      col <- if (label == 0L) {
        1L                                    # null column (always first)
      } else if (label == sim$K + 1L) {
        M                                     # junk column (always last)
      } else {
        cluster_rank[label] + num_null        # causal cluster, sorted position
      }
      etas[j, col] <- true_cluster_weight
    }
  }
  
  # Fit BayesClustMR
  res_bayes <- BayesClustMR(
      sim$bx, sim$by, sim$byse,
      etas     = etas,
      num_null = num_null,
      num_junk = num_junk,
      df       = 4,
      K_max = 10
    )
  
  return(
    tibble(
      theta_hat = list(sort(unique(res_bayes$res$cluster_mean[res_bayes$res$cluster_class != "Junk"]))),
      adj.rand = mclust::adjustedRandIndex(res_bayes$res$cluster, sim$cluster),
      rand = fossil::rand.index(res_bayes$res$cluster, sim$cluster),
      waic = res_bayes$waic,
      mse = mean((res_bayes$results$cluster_mean - c(0, sim$theta)[sim$cluster + 1])^2, na.rm = TRUE)
    )
  )
  
}
