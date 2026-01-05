# ============================================================================
# Premise Validation: Basic Multivariate Normal (Y only, no X)
# ============================================================================
# Purpose:
#   1. Calibrate (n, rho) for good model selection entropy
#   2. Validate imputation bias ≈ tr(RIV) empirically
#   3. Test scaling with model complexity Q
#
# Date: 2025-12-29
# ============================================================================

# Load required packages
library(future)
library(future.apply)
library(progressr)  # For progress bars
library(lavaan)
library(condMVNorm)
library(MASS)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

# Source utility functions
source("simulations/utils/generate_covariance.R")
source("simulations/utils/generate_data.R")
source("simulations/utils/fit_models.R")
source("simulations/utils/impute_data.R")
source("simulations/utils/compute_metrics.R")
source("simulations/utils/parallel_utils.R")
source("simulations/utils/utils-00-premise-basic-Y.R")

# ============================================================================
# PARALLEL PROCESSING SETUP
# ============================================================================
# Configure parallel processing once at the start
# All functions will use this configuration

cat("=== Configuring Parallel Processing ===\n")
n_workers <- setup_parallel(
  strategy = "multisession",
  workers = 25,  # Will use availableCores() - 1
  verbose = TRUE
)
cat("\n")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# # ---- PART 1: Grid Search ----
# # Uncomment to run grid search
#
# grid_results <- calibrate_selection_entropy(
#   p = 10,
#   true_structure = "Toeplitz",
#   n_grid = c(50, 100, 200, 500, 1000),
#   rho_grid = c(0.3, 0.5, 0.7, 0.9),
#   candidate_models = c("CS", "Toeplitz", "Unstructured"),
#   n_reps = 1000,
#   base_seed = 12345,
#   use_parallel = TRUE
# )
#
# # Save results
# dir.create("simulations/00-premise-basic-Y", showWarnings = FALSE, recursive = TRUE)
# saveRDS(grid_results, "simulations/00-premise-basic-Y/grid_search_results.rds")
# visualize_grid_search(grid_results)


# # ---- PART 2: Bias Validation ----
# 
# # Load calibrated parameters from grid search
grid_results <- readr::read_rds(file = "simulations/00-premise-basic-Y/grid_search_results.rds") %>%
  dplyr::arrange(desc(entropy))

n_calibrated <- grid_results$n[1]
rho_calibrated <- grid_results$rho[1]
p_calibrated <- 10
M_calibrated <- 200

# cat(sprintf("Using calibrated parameters:\n"))
# cat(sprintf("  n = %d\n", n_calibrated))
# cat(sprintf("  rho = %.2f\n", rho_calibrated))
# cat(sprintf("  p = %d\n", p_calibrated))
# cat(sprintf("  M = %d\n\n", M_calibrated))
# 
# # Run validation for different structures
# validation_results <- list()
# 
# for (struct in c("Toeplitz")) {
#   cat(sprintf("\n\n=================================================\n"))
#   cat(sprintf("Running: structure=%s, p=%d, M=%d\n", struct, p_calibrated, M_calibrated))
#   cat(sprintf("=================================================\n"))
# 
#   result <- validate_imputation_bias(
#     n = n_calibrated,
#     p = p_calibrated,
#     true_structure = struct,
#     rho = rho_calibrated,
#     missing_rate = 0.6,
#     M = M_calibrated,
#     scenarios = c("True", "MLE"),
#     n_reps = 1000,
#     base_seed = 54321,
#     use_parallel = TRUE
#   )
# 
#   validation_results[[struct]] <- result
# }
# 
# # Combine and save
# validation_df <- do.call(rbind, validation_results)
# saveRDS(validation_df, "simulations/00-premise-basic-Y/validation_results.rds")
# visualize_bias_validation(validation_df)

# ---- PART 3: Term Decomposition Validation ----
# Empirically verify Terms 1-3 from Taylor expansion to identify
# source of the 2×tr(RIV) discrepancy

library(numDeriv)  # Required for numerical differentiation

cat("\n\n=================================================\n")
cat("PART 3: Term Decomposition Validation\n")
cat("=================================================\n")

term_decomp_results <- validate_term_decomposition(
  n = n_calibrated,
  p = p_calibrated,
  true_structure = "Toeplitz",
  rho = rho_calibrated,
  missing_rate = 0.6,
  M = M_calibrated,
  n_reps = 1000,
  base_seed = 99999,
  use_parallel = TRUE
)

# Save results
saveRDS(term_decomp_results, "simulations/00-premise-basic-Y/term_decomposition_results.rds")

cat("\n\n=== ALL SIMULATIONS COMPLETE ===\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================
