#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
fit_mr_clust <- function(sim){
  ratios <- as.vector(sim$by / sim$bx)
  ratios.se <- sim$byse / abs(sim$bx)

  res_em <- mr_clust_em(theta = ratios,
                        theta_se = ratios.se,
                        bx = sim$bx,
                        by = sim$by,
                        bxse = sim$bxse,  # not returned by simulate_MR_data; unused by EM
                        byse = sim$byse,
                        k = sim$K
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
