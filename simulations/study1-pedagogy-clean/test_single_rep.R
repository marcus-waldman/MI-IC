# ============================================================================
# Test Single Replication for All Three Models
# ============================================================================

rm(list = ls())

cat("=== Testing Single Replication: Models A, B, C ===\n\n")

# Load functions
source("R/00_seed_management.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/03_analytical_functions.R")
source("R/04_compute_metrics.R")
source("R/05_run_simulation.R")

# Configuration
true_params <- list(
  a = 0.5,
  b = 0.5,
  sigma2_M = 1,
  sigma2_Y = 1
)

n <- 500
master_seed <- 12345

# Generate one seed
cat("Generating seed...\n")
seeds <- generate_seeds(master_seed = master_seed, n_reps = 1)

cat("Seed generated:\n")
cat(sprintf("  Data seed: %d\n", seeds[[1]]$data))
cat(sprintf("  Ampute seed: %d\n\n", seeds[[1]]$ampute))

# Test Model A
cat("=== MODEL A (Congenial) ===\n")
cat("Running 1 replication...\n")
results_A <- run_simulation(
  seeds = seeds,
  n = n,
  mechanism = "MCAR",
  gamma = 0.5,
  cut1 = 0.5,
  cut2 = 0.75,
  true_params = true_params,
  imputation_model = "A",
  n_cores = 1
)

cat("\nResults:\n")
print(results_A)

# Test Model B
cat("\n\n=== MODEL B (Uncongenial) ===\n")
cat("Running 1 replication...\n")
results_B <- run_simulation(
  seeds = seeds,
  n = n,
  mechanism = "MCAR",
  gamma = 0.5,
  cut1 = 0.5,
  cut2 = 0.75,
  true_params = true_params,
  imputation_model = "B",
  n_cores = 1
)

cat("\nResults:\n")
print(results_B)

# Test Model C
cat("\n\n=== MODEL C (Saturated MVN) ===\n")
cat("Running 1 replication...\n")
results_C <- run_simulation(
  seeds = seeds,
  n = n,
  mechanism = "MCAR",
  gamma = 0.5,
  cut1 = 0.5,
  cut2 = 0.75,
  true_params = true_params,
  imputation_model = "C",
  n_cores = 1
)

cat("\nResults:\n")
print(results_C)

# Summary comparison
cat("\n\n=== COMPARISON ACROSS MODELS ===\n\n")

comparison <- data.frame(
  Model = c("A (Congenial)", "B (Uncongenial)", "C (Saturated)"),
  term1_improper = c(results_A$term1_improper, results_B$term1_improper, results_C$term1_improper),
  term2 = c(results_A$term2, results_B$term2, results_C$term2),
  total_improper = c(results_A$total_improper, results_B$total_improper, results_C$total_improper),
  hessian_correction = c(results_A$hessian_correction, results_B$hessian_correction, results_C$hessian_correction),
  total_proper = c(results_A$total_proper, results_B$total_proper, results_C$total_proper),
  tr_RIV = c(results_A$tr_RIV, results_B$tr_RIV, results_C$tr_RIV)
)

print(comparison, row.names = FALSE)

cat("\n\nExpected patterns:\n")
cat("  - Model A: term1_improper ≈ tr(RIV) (congenial)\n")
cat("  - Model B: term1_improper ≠ tr(RIV) (uncongenial)\n")
cat("  - Model C: term1_improper ≈ tr(RIV) (congenial asymptotically)\n")
cat("  - All models: term2 ≈ -0.5 * tr(RIV)\n")
cat("  - All models: hessian_correction < 0\n")

cat("\nSingle replication test complete!\n")
