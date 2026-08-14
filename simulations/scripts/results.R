#!/usr/bin/env Rscript
rm(list=ls())

# -----------------------------
# Load packages
# -----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(ggplot2)
library(clue)

# -----------------------------
# Parameters
# -----------------------------
results_dir <- paste0("../results/rds/")
plot_dir <- "../results/plots/"
scenario_prefix <- "D_S"

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# Scenario metadata
# -----------------------------
scenario_metadata <- expand.grid(
  N = c(25000, 50000, 100000),
  null_clust = TRUE,
  junk_clust = TRUE
) %>% mutate(scenario = paste0("S", row_number()))

group_vars <- c("N", "method", "prior")

# -----------------------------
# Load and process all results
# -----------------------------
all_files <- list.files(
  results_dir,
  pattern = paste0("^", scenario_prefix, ".*\\.rds$"),
  full.names = TRUE,
  recursive = TRUE
)

read_and_expand <- function(path) {
  df <- readRDS(path)
  df$source_file <- path
  df$prior <- basename(dirname(path))
  
  return(df)
}

all_results <- map_dfr(all_files, read_and_expand) %>%
  mutate(scenario = str_extract(source_file, "S\\d+(_null)?")) %>%
  mutate(prior = factor(prior, levels = c("none", "right_confident", "right_diffuse", "uniform", 
                                          "wrong_diffuse", "wrong_confident"))) %>% 
  left_join(scenario_metadata, by = "scenario") %>%
  select(-source_file)

# -----------------------------
# Summary statistics
# -----------------------------
summary_stats <- all_results %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    nclust = mean(lengths(theta_hat)),
    adj.rand = mean(adj.rand),
    rand = mean(rand),
    mse = mean(mse),
    waic = mean(waic)
  )

# Print table
print(summary_stats, n = Inf)

# Plot number of clusters
ggplot(filter(all_results, prior == "none"), aes(x = method, y = lengths(theta_hat), fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_grid(. ~ N) +
  # coord_cartesian(ylim = c(0, 0.4)) + 
  labs(
    x = "Method",
    y = "Number of clusters",
    title = "Number of clusters by method and sample size (N)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

# Plot MSE values
ggplot(all_results, aes(x = method, y = mse, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_grid(N ~ prior) +
  coord_cartesian(ylim = c(0, 0.4)) + 
  labs(
    x = "Method",
    y = "MSE",
    title = "MSE by method, sample size (N), and prior"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

# Plot WAIC values
ggplot(filter(all_results, method == "BayesClustMR"), aes(x = method, y = waic, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_grid(N ~ prior) +
  labs(
    x = "",
    y = "WAIC",
    title = "WAIC by sample size (N), and prior"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

# -----------------------------
# Plots
# -----------------------------

## 0. Create color palette
method_levels <- c(
  "BayesClustMR", "mr_clust"
)

tol12 <- c(
  "#1F78B4", "grey"
)

method_colors <- setNames(tol12, method_levels)

## CLuster means density plot
selected_prior <- c("none", "right_confident", "right_diffuse", "uniform", 
           "wrong_diffuse", "wrong_confident")[1]

df_long <- all_results %>%
  filter(prior == selected_prior) %>%
  select(theta_hat, method, scenario, N, prior) %>%
  unnest_longer(theta_hat)

group_by(df_long, method) %>%
  summarise(
    max = max(theta_hat),
    min = min(theta_hat),
    sd = sd(theta_hat)
  )

true_thetas <- c(-0.4, 0, 0.4, 0.8)

p <- ggplot(filter(df_long, abs(theta_hat) < 3), 
            aes(x = theta_hat, fill = method, color = method)) +
  geom_density(alpha = 0.4, adjust = 0.5) +   # smaller `adjust` = less smoothing, sharper peaks
  geom_vline(xintercept = true_thetas, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_grid(N ~ ., scales = "free_y",
             labeller = labeller(N = function(x) paste0("N = ", format(as.numeric(x), scientific = FALSE, big.mark = " ")))) + 
  coord_cartesian(xlim = c(-0.8, 1.2)) +
  theme_minimal() +
  labs(x = expression(hat(theta)), y = "Density",
       title = bquote("Density of " ~ hat(theta) ~ " by method," ~ .(selected_prior) ~ "prior"))

print(p)
ggsave(paste0(plot_dir, "theta_hat_density.pdf"), plot = p, width = 8, height = 6, units = "in")
