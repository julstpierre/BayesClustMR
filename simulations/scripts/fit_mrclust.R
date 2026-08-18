#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
fit_mr_clust <- function(sim, prior){
  
  k <- if (prior == "none") NULL else sim$K
  k_max <- if (prior == "none") 10 else NULL
  
  ratios <- as.vector(sim$by / sim$bx)
  ratios.se <- sim$byse / abs(sim$bx)

  res_em <- mr_clust_em(theta = ratios,
                        theta_se = ratios.se,
                        bx = sim$bx,
                        by = sim$by,
                        bxse = sim$bxse,
                        byse = sim$byse,
                        k = k, k_max = k_max
  )
  
  res_em_best <- res_em$results$best
  
  return(
    tibble(
      theta_hat = list(sort(unique(res_em_best$cluster_mean[res_em_best$cluster_class != "Junk"]))),
      adj.rand = mclust::adjustedRandIndex(res_em_best$cluster, sim$cluster),
      rand = fossil::rand.index(res_em_best$cluster, sim$cluster),
      mse = mean((res_em_best$cluster_mean - c(0, sim$theta)[sim$cluster + 1])^2, na.rm = TRUE)
    )
  )
}
