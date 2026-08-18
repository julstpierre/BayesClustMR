#!/usr/bin/env Rscript

# Load libraries
library(cmdstanr)
library(loo)
library(dplyr)

# Source functions
source(here::here("R", "generate_priors.R"))
source(here::here("R", "determine_K.R"))
source(here::here("R", "initialization.R"))
source(here::here("R", "BayesClustMR.R"))
source(here::here("simulations", "scripts", "fit_BayesClustMR.R"))
source(here::here("simulations", "scripts", "mrclust.R"))
source(here::here("simulations", "scripts", "fit_mrclust.R"))

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

# args <- c("D_S1", "1", "2", "none")
scenario <- args[1]               # e.g., "D_S1"
start_idx     <- as.numeric(args[2])  # e.g., 1
end_idx       <- as.numeric(args[3])  # e.g., 200
prior <- args[4] # e.g., "none", "uniform", "right_diffuse", "right_confident", "wrong_diffuse", "wrong_confident"

# -----------------------------
# Load simulation data
# -----------------------------
load(paste0(here::here("simulations", "data", scenario), ".RData"))
D <- get(scenario)                # Extract the dataset corresponding to scenario_name

# Safety check for indexing
end_idx <- min(end_idx, length(D))
chunk <- D[start_idx:end_idx]

# -----------------------------
# Run methods
# -----------------------------
tag_method <- function(res_list, method) {
  bind_rows(res_list) %>% 
    mutate(method = method)
}

res_all <- bind_rows(
  tag_method(lapply(chunk, fit_BayesClustMR, prior = prior), "BayesClustMR"),
  tag_method(lapply(chunk, fit_mr_clust, prior = prior), "mr_clust")
)

# -----------------------------
# Save results
# -----------------------------
out_dir <- here::here("simulations", "results", "rds", prior)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
out_file <- paste0(scenario, "_chunk_", start_idx, "_to_", end_idx, ".rds")
out_path <- file.path(out_dir, out_file)
saveRDS(res_all, out_path)

message("✅ Done: ", scenario, " | chunk ", start_idx, "–", end_idx, " | saved to ", out_path)
