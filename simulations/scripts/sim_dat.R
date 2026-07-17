rm(list = ls())
library(parallel)
library(dplyr)
source("simulate_MR_data.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 4
num_iter <- 1000
params <- expand.grid(
  N = c(25000, 50000, 100000),
  null_clust = TRUE,
  junk_clust = TRUE
)

#----------------------------------------------------------------------
# Simulate data
#----------------------------------------------------------------------

simulate_MR_data_wrapper <- function(params) {
  # Get default arguments of the target function
  defaults <- formals(simulate_MR_data)
  
  # Override defaults with whatever is in params
  final_params <- modifyList(as.list(defaults), params)
  
  # Call the function with the merged argument list
  do.call(simulate_MR_data, final_params)
}

for (scenario_id in 1:nrow(params)){
  df_param <- params[scenario_id, , drop = FALSE] %>% list()
  
  set.seed(2026052900 + scenario_id)
  
  res <- replicate(num_iter, df_param) %>%
    mclapply(mc.cores = num_jobs, simulate_MR_data_wrapper, mc.set.seed = TRUE)
  
  # Save file
  D_name <- paste0("D_S", scenario_id)
  assign(paste0("D_S", scenario_id), res)
  save(list = D_name, file = sprintf("../data/D_S%d.RData", scenario_id))
}
