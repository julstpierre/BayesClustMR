# introduce "scale" variable to shrink denisty to fit in the vertical plot range
gen_t_scale <- function(theta, df = 4, mu, sig, log = TRUE, scale = 1) {
  tmp <- ((gamma((df + 1) / 2) / gamma(df / 2) / sqrt(pi * df) / sig) *
            (1 + (theta - mu)^2 / sig^2 / df) ^ (- (df + 1) / 2))
  if (log) {
    tmp <- log(tmp)
  }
  return(scale * tmp)
}

junk_clust_plot <- function(res, xrange = NULL, sig = NULL,
                            mu = 0, mu_null = 0,
                            junk_mixture = TRUE, null_mixture = TRUE) {
  cbpalette <- c("#000000", "#999999", "#0072B2",
                 "#D55E00", "#CC79A7", "#009E73",
                 "#56B4E9", "#E69F00", "#F0E442")
  n_palette <- length(unique(res$cluster)) - length(cbpalette)
  if (n_palette >= 0) {
    cbpalette <- c(cbpalette,
                   RColorBrewer::brewer.pal(max(3, n_palette + 1),
                                            "Spectral")[1:(n_palette + 1)])
  }
  
  # redefine null cluster as cluster zero
  res$cluster[res$cluster_class == "Null" | res$cluster_class == "Junk"] <- 0
  res$clusters <- as.factor(res$cluster)
  
  junk_obs <- res$cluster_class == "Junk"
  null_obs <- res$cluster_class == "Null"
  clust_obs <- which(!res$cluster_class == "Null" &
                       !res$cluster_class == "Junk")
  theta <- res$theta
  theta_se <- res$theta_se
  
  if (is.null(sig)) {
    rng_thet <- range(theta)
    max.disp <- which.max(abs(theta) + 2 * theta_se)
    sig <- (rng_thet[2] - rng_thet[1] + theta_se[max.disp])
  }
  if (is.null(mu)) {
    mu <- 0
  }
  if (sum(null_obs) > 0) {
    sig_null <- max(theta_se[null_obs])
    mu_null <- 0
  } else {
    sig_null <- 1
    mu_null <- 0
  }
  
  jit_junk <- gen_t(theta[junk_obs], df = 4, mu = mu, sig = sig, log = F)
  jit_null <- stats::dnorm(theta[null_obs], mean = mu_null, sd = sig_null,
                           log = F)
  jit_clusts <- rep(0, length(theta))
  jit_clusts[junk_obs] <- jit_junk * stats::runif(sum(junk_obs), 0, 0.5)
  jit_clusts[null_obs] <- jit_null * stats::runif(sum(null_obs), 0, 0.5)
  if (sum(junk_obs) > 0 | sum(null_obs) > 0) {
    mny <- max(c(jit_clusts[junk_obs], jit_clusts[null_obs]))
  } else {
    mny <- 1
  }
  jit_clusts[clust_obs] <- -stats::runif(length(clust_obs), 0, mny)
  
  res$jit_clusts <- jit_clusts
  maxx <- max(c(theta + 2.5 * theta_se, theta - 2.5 * theta_se))
  minx <- min(c(theta + 2.5 * theta_se, theta - 2.5 * theta_se))
  
  maxy <- if (sum(junk_obs | null_obs) > 0) {
    max(jit_clusts[junk_obs | null_obs])
  } else if (sum(junk_obs | null_obs) == 0 & junk_mixture & !null_mixture) {
    max(gen_t(c(minx, maxx), mu = mu, sig = sig, log = FALSE))
  } else if (sum(junk_obs | null_obs) == 0 & !junk_mixture & null_mixture) {
    if (sum(sig_null) > 0) {
      stats::dnorm(0, 0, sig_null)
    } else {
      0
    }
  } else if (sum(junk_obs | null_obs) == 0 & junk_mixture & null_mixture) {
    max(max(gen_t(c(minx, maxx), mu = mu, sig = sig, log = FALSE)),
        if (sum(sig_null) > 0) {
          stats::dnorm(0, 0, sig_null)
        } else {
          0
        }
    )
  }
  
  miny <- if (length(clust_obs) > 0) {
    min(jit_clusts[clust_obs])
  } else {
    0
  }
  mx_dens <- if (junk_mixture & !null_mixture) {
    gen_t(mu, df = 4, mu = mu, sig = sig, log = FALSE)
  } else if (!junk_mixture & null_mixture) {
    maxy
  } else if (junk_mixture & null_mixture) {
    max(gen_t(mu, mu = mu, sig = sig, log = FALSE),
        stats::dnorm(0, 0, sig_null))
  } else {
    0
  }
  
  if (maxy < mx_dens) {
    shrink <- maxy / mx_dens
    res$jit_clusts[junk_obs | null_obs] <- (shrink *
                                              jit_clusts[junk_obs | null_obs])
    minx <- 1.25 * minx
  } else {
    shrink <- 1
  }
  
  if (!is.null(xrange)) {
    maxx <- xrange[2]
    minx <- xrange[1]
  }
  
  annotations <- data.frame(
    xpos = c(minx, minx),
    ypos = c(-miny, miny / 1.25),
    annotateText = c("Junk/null\nestimates", "Clustered\nestimates"),
    hjustvar = c(0, 0),
    vjustvar = c(1, 0)
  )
  
  
  p <- ggplot2::ggplot(data = res, ggplot2::aes(theta, jit_clusts))
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = cluster_class), size = 1.55)
  p <- p + ggplot2::scale_color_manual(values = cbpalette)
  p <- p + ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = theta - 1.96 * theta_se,
                 xmax = theta + 1.96 * theta_se,
                 color = cluster_class), linetype = "solid")
  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                               color = cbpalette[1], lwd = 0.5) 
  if (junk_mixture) {
    p <- p + ggplot2::geom_line(stat = "function", fun = gen_t_scale,
                                args = with(res, c(df = 4, mu = mu, sig = sig, log = FALSE,
                                                   scale = shrink)),
                                ggplot2::aes(colour = "Junk density")
    )
  }
  p <- p + ggplot2::ylim(miny, maxy)
  clusts <- unique(res$cluster_class)
  clusts <- clusts[order(clusts)]
  ind <- which(!clusts %in% "Null" & !clusts %in% "Junk")
  if (sum(ind) > 0) {
    clust_label <- c("Junk cluster")
  } else {
    clust_label <- NULL
  }
  if (sum(ind) > 0) {
    for (i in seq_len(length(ind))) {
      mn <- res$cluster_mean[res$cluster == as.numeric(clusts[ind[i]])][1]
      clust_label <- c(clust_label, paste0("cluster_", ind[i]))
      col <- cbpalette[i]
      p <- p + ggplot2::geom_segment(x = mn, y = miny, xend = mn, yend = 0,
                                     color = col, linetype = "dotted",
                                     lwd = 0.2)
    }
  }
  p <- p + ggplot2::geom_text(data = annotations,
                              ggplot2::aes(x = xpos, y = ypos, hjust = hjustvar,
                                           vjust = vjustvar,
                                           label = annotateText), size = 3) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Cluster",
                                                  labels = clust_label)) +
    ggplot2::xlim(minx, maxx)
  p <- p +
    ggplot2::ylab("") + ggplot2::xlab("Two-stage ratio estimate") +
    ggplot2::ggtitle("Null and/or junk observations (top)\nand cluster
                     partitioned ratio estimates (bottom)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(), ## <- this line
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(colour = "white"),
      plot.background = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  
  if (null_mixture & sum(null_obs) > 0) {
    p <- p + ggplot2::geom_line(stat = "function",
                                fun = scaled_dnorm,
                                args = with(res, c(mean = mu_null, sd = sig_null, log = FALSE,
                                                   scale = shrink)),
                                ggplot2::aes(colour = "Null density")) +
      ggplot2::ylab("") + ggplot2::xlab("Two-stage ratio estimate") +
      ggplot2::ggtitle("Null and/or junk observations (top)\nand cluster
                       partitioned ratio estimates (bottom)") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(), ## <- this line
        axis.text.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(colour = "white"),
        plot.background = ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }
  
  return(p)
}

clust_inc_list <- function(dta, by_prob = 0.2, bound = 0) {
  if (!((1 - bound) / by_prob) %% 1 == 0) {
    stop("(1-bound)/by_prob needs to be an integer")
  }
  unique_clust <- unique(dta$cluster)
  res <- vector("list", length(unique_clust))
  clust_nam <- vector("character", length(unique_clust))
  nms_tmp <- unique(dta$cluster_class)
  for (i in unique_clust) {
    if (bound == 0) {
      sub_grps <- (1 / by_prob)
    } else {
      sub_grps <- ceiling((1 - bound) / by_prob)
    }
    tmp <- vector("list", sub_grps)
    tmp_nam <- vector("character", sub_grps)
    for (j in 1:sub_grps) {
      p_mx <- 1 - (j - 1) * by_prob
      p_mn <- 1 - j * by_prob
      if (p_mx <= bound) {
        break
      }
      if ((!bound == 0) & (1 - j * by_prob < bound)) {
        p_mn <- bound
      }
      tmp[[j]] <- (dta$observation[(dta$probability <= p_mx)
                                   & (dta$probability > p_mn)
                                   & (dta$cluster == i)])
      tmp_nam[j] <- paste0("pr_range_", p_mx, "_to_", p_mn)
    }
    names(tmp) <- tmp_nam
    clust_nam[i] <- paste0("cluster_", nms_tmp[i])
    res[[i]] <- tmp
  }
  names(res) <- clust_nam
  return(res)
}

best_clust <- function(dta) {
  res <- NULL
  unique.obs <- unique(dta$observation)
  for (i in unique.obs) {
    tmp <- dta$observation == i
    mx <- which.max(dta$probability[tmp])
    res <- rbind(res, dta[tmp, ][mx, ])
  }
  return(res)
}

# prevents R check throwing errors about unassigned global variables
utils::globalVariables(c("annotateText", "cluster_class", "hjustvar",
                         "probability", "runif", "vjustvar", "xpos", "ypos",
                         "num_c"))

clust_prob <- function(bic_clust, clust_mn, clust_pi, theta, theta_sd, df,
                       junk_mixture = FALSE, junk_mean = NULL, junk_sd = NULL,
                       null_mixture = FALSE, null_mean = NULL, null_sd = NULL,
                       obs_names, clust_size_prior = FALSE, prior = NULL) {
  if (!clust_size_prior) {
    mx <- which.min(bic_clust)
  } else {
    m <- length(theta)
    pr <- stats::dbinom(1:m, m, prior)
    pr <- pr[seq_len(bic_clust)] / sum(pr[seq_len(bic_clust)])
    mx <- which.max(exp(-bic_clust / 2) * pr)
  }
  theta_mn <- clust_mn[[mx]]
  clust_probs <- clust_pi[[mx]]
  n_clust <- length(theta_mn)
  m <- length(theta)
  res1 <- round(as.numeric(t(rij(1:m, 1:n_clust, clust_probs,
                                 theta, theta_sd, theta_mn,
                                 junk_mixture, df, junk_mean, junk_sd,
                                 null_mixture, null_mean, null_sd))), 3)
  res2 <- round(rep(theta_mn, m), 3)
  res <- data.table::data.table(rep(obs_names, each = n_clust),
                                rep(1:n_clust, m),
                                round(clust_probs, 3),
                                res1,
                                res2,
                                rep(theta, each = n_clust),
                                rep(theta_sd, each = n_clust))
  names(res) <- c("observation", "cluster", "pi_hat","probability",
                  "cluster_mean", "theta", "theta_se")
  factor_clust <- res$cluster
  if (null_mixture & !junk_mixture) {
    tmp <- res$cluster == n_clust
    factor_clust[tmp] <- "Null"
  } else if (null_mixture & junk_mixture) {
    tmp <- res$cluster == n_clust - 1
    factor_clust[tmp] <- "Null"
  }
  if (junk_mixture) {
    tmp <- res$cluster == n_clust
    factor_clust[tmp] <- "Junk"
  }
  res$cluster_class <- as.character(factor_clust)
  return(list(res = res, r_ij = rij(1:m, 1:n_clust, clust_probs,
                                    theta, theta_sd, theta_mn,
                                    junk_mixture, df, junk_mean, junk_sd,
                                    null_mixture, null_mean, null_sd)))
}

pi_updt <- function(j, m, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
                    df, mu, sig, junk_prob, fix_junk_prob, null_mixture,
                    mu_null, sig_null, null_prob, fix_null_prob) {
  tmp <- rij(1:m, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
             df, mu, sig, null_mixture, mu_null, sig_null)
  num <- if (length(j) > 1) {
    colSums(tmp)
  } else {
    sum(tmp)
  }
  jk_null <- sum(junk_mixture) + sum(null_mixture)
  n_prob <- j_prob <- FALSE
  if (null_mixture & !is.null(null_prob)) {
    tmp_lgt <- length(num) - sum(junk_mixture)
    if (num[tmp_lgt] / m < null_prob | fix_null_prob == TRUE) {
      num[tmp_lgt] <- null_prob * m
      n_prob <- TRUE
    }
  }
  if (junk_mixture & !is.null(junk_prob)) {
    tmp_lgt <- length(num)
    if (num[length(num)] / m < junk_prob | fix_junk_prob == TRUE) {
      num[tmp_lgt] <- junk_prob * m
      j_prob <- TRUE
    }
  }
  if (j_prob | n_prob) {
    num[1:(length(num) - jk_null)] <- (m *
                                         (1 - sum(num[(length(num) - jk_null + 1):length(num)]) / m) *
                                         num[1:(length(num) - jk_null)] / (sum(num[1:(length(num) - jk_null)])))
  }
  return(num / m)
}

rij <- function(i, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
                df, mu, sig, null_mixture, mu_null, sig_null) {
  tmp_lgt <- length(i)
  if (junk_mixture | null_mixture) {
    tmp_lgt2 <- length(j) - sum(junk_mixture) - sum(null_mixture)
    tmp_mn <- theta_clust[1:tmp_lgt2]
    j <- 1:(length(j) - sum(junk_mixture) - sum(null_mixture))
  } else {
    tmp_lgt2 <- length(j)
    tmp_mn <- theta_clust
  }
  if (junk_mixture | null_mixture) {
    tmp_lg_num <- log_norm(i, j, tmp_lgt, tmp_lgt2, theta, theta_sd, tmp_mn)
    if (null_mixture) {
      tmp_lg_num <- c(tmp_lg_num, null_den(x = theta, mu = mu_null,
                                           sig = sig_null, log = TRUE))
    }
    if (junk_mixture) {
      tmp_lg_num <- c(tmp_lg_num,
                      gen_t(theta, df = df, mu = mu, sig = sig, log = TRUE))
    }
    num <- exp(matrix(rep(log(pi_clust), each = tmp_lgt) + tmp_lg_num,
                      nrow = tmp_lgt,
                      ncol = tmp_lgt2 + sum(junk_mixture) + sum(null_mixture)))
    num[num < 1e-300] <- 1e-300
    den <- rowSums(num)
  } else {
    tmp_lg_num <- log_norm(i, j = seq_len(length(theta_clust)), tmp_lgt,
                           tmp_lgt2 = length(theta_clust),
                           theta, theta_sd, theta_clust)
    num <- exp(matrix(rep(log(pi_clust[j]), each = tmp_lgt) + tmp_lg_num,
                      nrow = tmp_lgt, ncol = tmp_lgt2))
    num[num < 1e-300] <- 1e-300
    den <- rowSums(num)
  }
  res <- num / den
  return(res)
}

theta_updt <- function(j, m, theta, theta_sd, theta_clust, pi_clust,
                       junk_mixture, df, mu, sig, null_mixture, mu_null,
                       sig_null) {
  tmp <- rij(1:m, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture, df,
             mu, sig, null_mixture, mu_null, sig_null) / theta_sd[1:m]^2
  den <- if (length(j) > 1) {
    colSums(tmp)
  } else {
    sum(tmp)
  }
  num <- if (length(j) > 1) {
    colSums(tmp * theta)
  } else {
    sum(tmp * theta)
  }
  new_theta <- num / den
  null_jk <- sum(null_mixture) + sum(junk_mixture)
  if (null_mixture & null_jk == 2) {
    new_theta[length(j) - 1] <- mu_null
  } else if (null_mixture & null_jk == 1) {
    new_theta[length(j)] <- mu_null
  }
  if (junk_mixture) {
    new_theta[length(j)] <- mu
  }
  return(num / den)
}

gen_t <- function(x, df = 4, mu, sig, log = TRUE) {
  tmp <- ((gamma((df + 1) / 2) / gamma(df / 2) / sqrt(pi * df) / sig)
          * (1 + (x - mu)^2 / sig^2 / df) ^ (- (df + 1) / 2))
  if (log) {
    tmp <- log(tmp)
  }
  return(tmp)
}

null_den <- function(x, mu, sig, log = TRUE) {
  tmp <- -0.5 * log(2 * pi) - log(sig) - 0.5 * (x - mu)^2 / sig^2
  return(tmp)
}

log_norm <- function(i, j, tmp_lgt, tmp_lgt2, theta, theta_sd, theta_clust) {
  tmp_theta <- rep(theta[i], tmp_lgt2)
  tmp_theta_sd <- rep(theta_sd[i], tmp_lgt2)
  tmp_clust <- rep(theta_clust[j], each = tmp_lgt)
  res <- (-0.5 * log(2 * pi) - log(tmp_theta_sd)
          - 0.5 * (tmp_theta - tmp_clust)^2 / tmp_theta_sd^2)
  return(res)
}

loglik <- function(m, pi_clust, theta, theta_sd, theta_clust,
                   junk_mixture, df, mu, sig,
                   null_mixture, mu_null, sig_null) {
  i <- 1:m
  if (junk_mixture | null_mixture) {
    k <- length(theta_clust) - sum(junk_mixture) - sum(null_mixture)
  } else {
    k <- length(theta_clust)
  }
  tmp_lg_den <- log_norm(i, j = 1:k, m, tmp_lgt2 = k, theta, theta_sd,
                         theta_clust)
  if (null_mixture) {
    tmp_lg_den <- c(tmp_lg_den, null_den(x = theta, mu = mu_null,
                                         sig = sig_null, log = TRUE))
  }
  if (junk_mixture) {
    tmp_lg_den <- c(tmp_lg_den, gen_t(theta, df = df, mu = 0, sig = sig,
                                      log = TRUE))
  }
  
  tmp <- exp(matrix(rep(log(pi_clust), each = m) + tmp_lg_den, nrow = m,
                    ncol = length(pi_clust)))
  tmp_rowsum <- rowSums(tmp)
  res <- sum(log(tmp_rowsum))
  return(res)
}

mr_clust_em <- function (theta, theta_se, bx, by, bxse, byse, obs_names = NULL,
                         max_iter = 5000, tol = 1e-05, junk_sd = NULL, junk_mean = 0,
                         stop_bic_iter = 5, min_clust_search = 10, results_list = list("all",
                                                                                       "best"), cluster_membership = list(by_prob = 0.1, bound = 0),
                         plot_results = list("best", min_pr = 0.5), trait_search = FALSE,
                         trait_pvalue = 1e-05, proxy_r2 = 0.8, catalogue = "GWAS",
                         proxies = "None", build = 37, cluster_sizes = 0:length(theta),
                         k = NULL, k_max = NULL)
{
  if (!is.null(k)) cluster_sizes <- k
  if (!is.null(k_max)) cluster_sizes <- 0:k_max
  
  init_clust_means <- NULL
  init_clust_probs <- NULL
  junk_mixture <- TRUE
  df <- 4
  junk_prob <- NULL
  junk_null_aware_bic <- TRUE
  fix_junk_prob <- FALSE
  null_mixture <- TRUE
  null_sd <- NULL
  null_prob <- NULL
  fix_null_prob <- FALSE
  scale_grid_search <- FALSE
  grid_increment <- 0.1
  grid_max <- 2
  clust_size_prior <- FALSE
  bic_prior <- NULL
  rand_num <- 5
  rand_sample <- seq(0.05, 0.4, by = 0.05)
  m <- length(theta)

  if (is.null(obs_names)) {
    obs_names <- paste0("snp_", 1:m)
  }
  num_clust <- length(cluster_sizes)
  num_clust <- num_clust * rand_num
  bic_clust <- bic_clust_mx <- vector()
  clust_mn <- vector("list", num_clust)
  clust_pi <- vector("list", num_clust)
  log_like <- vector("list", num_clust)
  count0 <- count_k <- 1
  for (i in cluster_sizes) {
    for (itr in 1:(rand_num + 1)) {
      if (is.null(init_clust_means) | is.null(init_clust_probs)) {
        if (i > 0 & i != m) {
          init_conds <- stats::kmeans(x = theta, centers = i, 
                                      iter.max = 5000)
          clust_means <- as.numeric(init_conds$centers)
          clust_probs <- table(init_conds$cluster)/m
          init_conds
          clust_means
          clust_probs
        }
        else if (i == m) {
          clust_means <- theta
          clust_probs <- rep(1, m)/m
        }
        else if (i == 0) {
          clust_means <- clust_probs <- NULL
        }
      }
      if (junk_mixture & is.null(junk_sd)) {
        rng_thet <- range(theta)
        max_disp <- which.max(abs(theta) + 2 * theta_se)
        sig <- (rng_thet[2] - rng_thet[1] + theta_se[max_disp])
      }
      else if (junk_mixture & !is.null(junk_sd)) {
        sig <- junk_sd
      }
      else {
        sig <- NULL
      }
      if (junk_mixture & is.null(junk_mean)) {
        mu <- 0
      }
      else if (junk_mixture & !is.null(junk_mean)) {
        mu <- junk_mean
      }
      else {
        mu <- NULL
      }
      null_obs <- abs(theta/theta_se) < 1.96
      if (null_mixture & is.null(null_sd)) {
        sig_null <- theta_se
        mu_null <- 0
      }
      else if (null_mixture & !is.null(null_sd)) {
        sig_null <- null_sd
        mu_null <- 0
      }
      else {
        sig_null <- mu_null <- NULL
      }
      if (itr == 1) {
        if (junk_mixture | null_mixture) {
          sum_pr <- 0
          junk_obs <- sum(2 * (1 - stats::pnorm(theta - 
                                                  stats::median(theta), 0, theta_se)) < 0.05)
          if (junk_mixture & !fix_junk_prob) {
            if (sum(junk_obs) == 0) {
              jk_pr <- 1/m
              sum_pr <- jk_pr
            }
            else {
              jk_pr <- sum(junk_obs)/m
              sum_pr <- jk_pr
            }
          }
          else if (junk_mixture & fix_junk_prob) {
            jk_pr <- junk_prob
            sum_pr <- jk_pr
          }
          else {
            jk_pr <- NULL
          }
          if (null_mixture & !fix_null_prob) {
            if (is.null(jk_pr)) {
              tmp_jk_pr <- 0
            }
            else {
              tmp_jk_pr <- jk_pr
            }
            if (sum(null_obs) == 0) {
              null_pr <- 1/m
              sum_pr <- sum_pr + null_pr
            }
            else {
              null_pr <- sum(null_obs)/m
              if (null_pr + tmp_jk_pr > (1 - 1/m)) {
                tmp_sum <- null_pr + tmp_jk_pr
                if (is.null(jk_pr)) {
                  null_pr <- null_pr/tmp_sum - 1/m
                }
                else {
                  null_pr <- null_pr/tmp_sum - 1/(2 * 
                                                    m)
                  tmp_jk_pr <- tmp_jk_pr/tmp_sum - 1/(2 * 
                                                        m)
                  jk_pr <- tmp_jk_pr
                  sum_pr <- tmp_jk_pr
                }
              }
              sum_pr <- sum_pr + null_pr
            }
          }
          else if (null_mixture & fix_null_prob) {
            null_pr <- null_prob
            sum_pr <- sum_pr + null_pr
          }
          else {
            null_pr <- NULL
          }
          clust_means <- c(clust_means, mu_null, mu)
          clust_probs <- c(clust_probs * (1 - sum_pr), 
                           null_pr, jk_pr)
          k_clust <- cluster_sizes[count_k] + sum(junk_mixture) + 
            sum(null_mixture)
        }
        else {
          k_clust <- cluster_sizes[count_k]
          sig <- mu <- NULL
          sig_null <- mu_null <- NULL
          junk_null_aware_bic <- FALSE
        }
      }
      else {
        if (null_mixture) {
          null_pr <- sample(rand_sample, 1)
          sum_pr <- null_pr
        }
        else {
          null_pr <- NULL
          sum_pr <- 0
        }
        if (junk_mixture) {
          jk_pr <- sample(rand_sample, 1)
          sum_pr <- sum_pr + jk_pr
        }
        else {
          jk_pr <- NULL
        }
        clust_means <- c(clust_means, mu_null, mu)
        clust_probs <- c(clust_probs * (1 - sum_pr), 
                         null_pr, jk_pr)
        k_clust <- cluster_sizes[count_k] + sum(junk_mixture) + 
          sum(null_mixture)
      }
      if (scale_grid_search & junk_mixture) {
        sigs <- sig * seq(1, grid_max, by = grid_increment)
        pi_tmp <- vector("list", length(sigs))
        theta_tmp <- vector("list", length(sigs))
        loglik_tmp <- vector("numeric", length(sigs))
        count2 <- 1
        for (inc in sigs) {
          sig_tmp <- inc
          pi_clust <- clust_probs
          theta_clust <- clust_means
          count <- 1
          loglik_diff <- 1
          loglik_iter <- NULL
          while (loglik_diff > tol & count < max_iter) {
            tmp_loglik0 <- loglik(m, pi_clust, theta, 
                                  theta_se, theta_clust, junk_mixture, df, 
                                  mu, sig_tmp, null_mixture, mu_null, sig_null)
            tmp_clust <- theta_updt(1:k_clust, m, theta, 
                                    theta_se, theta_clust, pi_clust, junk_mixture, 
                                    df, mu, sig_tmp, null_mixture, mu_null, 
                                    sig_null)
            tmp_pi <- pi_updt(1:k_clust, m, pi_clust, 
                              theta, theta_se, theta_clust, junk_mixture, 
                              df, mu, sig_tmp, junk_prob, fix_junk_prob, 
                              null_mixture, mu_null, sig_null, null_prob, 
                              fix_null_prob)
            theta_clust <- tmp_clust
            pi_clust <- tmp_pi
            tmp_loglik <- loglik(m, pi_clust, theta, 
                                 theta_se, theta_clust, junk_mixture, df, 
                                 mu, sig_tmp, null_mixture, mu_null, sig_null)
            loglik_diff <- abs(tmp_loglik - tmp_loglik0)
            loglik_iter <- c(loglik_iter, tmp_loglik0)
            count <- count + 1
          }
          pi_tmp[[count2]] <- pi_clust
          theta_tmp[[count2]] <- theta_clust
          loglik_tmp[count2] <- tmp_loglik
          count2 <- count2 + 1
        }
        mx <- which.max(loglik_tmp)
        theta_clust <- theta_tmp[[mx]]
        pi_clust <- pi_tmp[[mx]]
        tmp_loglik <- loglik_tmp[mx]
        sig <- sigs[mx]
      }
      else {
        pi_clust <- clust_probs
        theta_clust <- clust_means
        count <- 1
        loglik_diff <- 1
        loglik_iter <- NULL
        while (loglik_diff > tol & count < max_iter) {
          tmp_loglik0 <- loglik(m, pi_clust, theta, theta_se, 
                                theta_clust, junk_mixture, df, mu, sig, null_mixture, 
                                mu_null, sig_null)
          tmp_clust <- theta_updt(1:k_clust, m, theta, 
                                  theta_se, theta_clust, pi_clust, junk_mixture, 
                                  df, mu, sig, null_mixture, mu_null, sig_null)
          tmp_pi <- pi_updt(1:k_clust, m, pi_clust, theta, 
                            theta_se, theta_clust, junk_mixture, df, 
                            mu, sig, junk_prob, fix_junk_prob, null_mixture, 
                            mu_null, sig_null, null_prob, fix_null_prob)
          theta_clust <- tmp_clust
          pi_clust <- tmp_pi
          tmp_loglik <- loglik(m, pi_clust, theta, theta_se, 
                               theta_clust, junk_mixture, df, mu, sig, null_mixture, 
                               mu_null, sig_null)
          loglik_diff <- abs(tmp_loglik - tmp_loglik0)
          loglik_iter <- c(loglik_iter, tmp_loglik0)
          count <- count + 1
        }
      }
      clust_mn[[count0]] <- theta_clust
      clust_pi[[count0]] <- pi_clust
      log_like[[count0]] <- loglik_iter
      bic_clust[count0] <- if (!junk_null_aware_bic) {
        -2 * tmp_loglik + (2 * k_clust - 1) * log(m)
      }
      else {
        -2 * tmp_loglik + (2 * k_clust - 1 - sum(junk_mixture) - 
                             sum(null_mixture)) * log(m)
      }
      count0 <- count0 + 1
    }
    bic_clust_mx[count_k] <- max(bic_clust[(count0 - (rand_num + 
                                                        2)):(count0 - 1)])
    if ((count_k > stop_bic_iter) & (i >= min_clust_search)) {
      tmp_cond <- TRUE
      for (j in 1:stop_bic_iter) {
        tmp_cond <- (tmp_cond & (bic_clust_mx[count_k - 
                                                j + 1] > bic_clust_mx[count_k - j]))
      }
      if (tmp_cond) {
        break
      }
    }
    count_k <- count_k + 1
  }
  results <- clust_prob(bic_clust, clust_mn, clust_pi, theta = theta, 
                        theta_sd = theta_se, junk_mixture = junk_mixture, df = df, 
                        junk_mean = mu, junk_sd = sig, null_mixture = null_mixture, 
                        null_mean = mu_null, null_sd = sig_null, obs_names = obs_names, 
                        clust_size_prior = clust_size_prior, prior = bic_prior)
  results_all <- results$res[order(results$res$cluster), ]
  results_best <- best_clust(results$res)
  if (!is.null(cluster_membership)) {
    variant_clusters <- clust_inc_list(results$res, by_prob = cluster_membership$by_prob, 
                                       bound = cluster_membership$bound)
  }
  # if (!is.null(plot_results)) {
  #   if (plot_results[[1]] == "best") {
  #     plot_1 <- junk_clust_plot(results_best, sig = sig, 
  #                               mu = mu, mu_null = mu_null)
  #     plot_2 <- two_stage_plot(results_best, bx, by, bxse, 
  #                              byse, obs_names)
  #   }
  #   else {
  #     tmp <- plot_results[[2]]
  #     tmp_res <- pr_clust(results, prob = tmp)
  #     plot_1 <- junk_clust_plot(tmp_res, sig = sig, mu = mu, 
  #                               mu_null = mu_null, junk_mixture = junk_mixture, 
  #                               null_mixture = null_mixture)
  #     plot_2 <- two_stage_plot(tmp_res, bx, by, bxse, byse, 
  #                              obs_names)
  #   }
  #   res <- list(results = list(all = results_all, best = results_best), 
  #               cluster_membership = variant_clusters, plots = list(two_stage = plot_2, 
  #                                                                   split_plot = plot_1), log_likelihood = log_like, 
  #               bic = bic_clust)
  # }
  # else {
    res <- list(results = list(all = results_all, best = results_best), 
                cluster_membership = variant_clusters, log_likelihood = log_like, 
                bic = bic_clust, r_ij = results$r_ij)
  # }
  if (trait_search) {
    trt_search <- pheno_search(variant_clusters = variant_clusters, 
                               p_value = trait_pvalue, r2 = proxy_r2, catalogue = catalogue, 
                               proxies = proxies, build = build)
    res <- list(results = list(all = results_all, best = results_best), 
                cluster_membership = variant_clusters, trait_search = trt_search, 
                log_likelihood = log_like, bic = bic_clust)
  }
  return(res)
}
