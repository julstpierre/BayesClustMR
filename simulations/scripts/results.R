#!/usr/bin/env Rscript

# -----------------------------
# Load packages
# -----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(ggplot2)

# -----------------------------
# Parameters
# -----------------------------
results_dir <- "../results/rds/"
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

group_vars <- c("N", "method")

# -----------------------------
# Load and process all results
# -----------------------------
all_files <- list.files(results_dir, pattern = paste0("^", scenario_prefix, ".*\\.rds$"), full.names = TRUE)

read_and_expand <- function(path) {
  df <- readRDS(path)
  df$source_file <- path
  
  # df <- df %>%
  #   mutate(
  #     across(c(Pvalue, Estimate, StdError, LCL, UCL),
  #            ~ suppressWarnings(as.numeric(.)))
  #   )
  
  return(df)
}

all_results <- map_dfr(all_files, read_and_expand) %>%
  mutate(scenario = str_extract(source_file, "S\\d+(_null)?")) %>%
  left_join(scenario_metadata, by = "scenario") %>%
  select(-source_file)

# -----------------------------
# Summary statistics
# -----------------------------
summary_stats <- all_results %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    adj.rand = mean(adj.rand),
    rand = mean(rand)
  )

summary_stats

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
df_long <- all_results %>%
  select(theta_hat, method, scenario, N) %>%
  unnest_longer(theta_hat) %>%
  filter(theta_hat > -3, theta_hat < 3)

group_by(df_long, method) %>%
  summarise(
    max = max(theta_hat),
    min = min(theta_hat),
    sd = sd(theta_hat)
  )

true_thetas <- c(-0.4, 0, 0.4, 0.8)

p <- ggplot(df_long, aes(x = theta_hat, fill = method, color = method)) +
  geom_density(alpha = 0.4, adjust = 0.5) +   # smaller `adjust` = less smoothing, sharper peaks
  geom_vline(xintercept = true_thetas, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_grid(N ~ ., scales = "free_y",
             labeller = labeller(N = function(x) paste0("N = ", format(as.numeric(x), scientific = FALSE, big.mark = " ")))) + 
  coord_cartesian(xlim = c(-0.8, 1.2)) +
  theme_minimal() +
  labs(x = expression(hat(theta)), y = "Density",
       title = expression("Density of " ~ hat(theta) ~ " by method"))

print(p)
ggsave(paste0(plot_dir, "theta_hat_density.pdf"), plot = p, width = 8, height = 6, units = "in")

# -----------------------------
# Done
# -----------------------------
message("✅ Summary and plots saved in: ", results_dir, " and ", plot_dir)
