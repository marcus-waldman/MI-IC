# ============================================================================
# Test Script: Analytical Approach (M = Inf)
# ============================================================================
# Quick validation that the analytical approach works correctly
# ============================================================================

rm(list = ls())

# Set working directory to script location
if (interactive()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

cat("=== Study 1: Analytical Approach Test ===\n\n")

# ============================================================================
# 1. Load Functions
# ============================================================================

cat("1. Loading functions...\n")

source("R/00_setup_seeds.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/06_compute_metrics.R")
source("R/07_run_replication.R")

# Only need base R for analytical approach (no brms/lavaan required)
cat("   Functions loaded (no Stan/brms required for analytical approach).\n\n")

# ============================================================================
# 2. Create Test Configuration with M = Inf
# ============================================================================

cat("2. Creating analytical test configuration (M = Inf)...\n")

test_config <- list(
  # Sample size
  n = 500,
  prop_missing = 0.3,
  n_reps = 100,  # 100 reps for quick validation

  # KEY: M = Inf triggers analytical approach
  M = Inf,

  # Seeds
  master_seed = 12345,
  seeds_file = "results/test_analytical_seeds.rds",

  # Population parameters
  beta_xm = 0.5,
  beta_my = 0.5,
  beta_xy = 0,
  sigma2_m = 1,
  sigma2_y = 1,
  var_x = 1,

  # Test conditions
  mechanisms = c("MCAR"),
  imputation_models = c("matched"),  # Only matched for now (analytical assumes congenial)

  # Seed offsets
  mechanism_offsets = c(MCAR = 0, MAR = 1000, MNAR = 2000),
  model_offsets = c(matched = 0, saturated = 100, uncongenial = 200)
)

cat(sprintf("   n = %d, prop_missing = %.0f%%, n_reps = %d\n",
            test_config$n, test_config$prop_missing * 100, test_config$n_reps))
cat("   M = Inf (analytical approach - no MCMC sampling)\n\n")

# ============================================================================
# 3. Setup Seeds
# ============================================================================

cat("3. Setting up seeds...\n")
seeds <- setup_seeds(test_config)
cat(sprintf("   Generated %d seeds.\n\n", length(seeds)))

# ============================================================================
# 4. Run Single Replication Test
# ============================================================================

cat("4. Testing single replication...\n")

single_result <- run_single_replication(
  rep_id = 1,
  mechanism = "MCAR",
  imputation_model = "matched",
  config = test_config,
  seeds = seeds,
  verbose = TRUE
)

cat("\n   Single replication results:\n")
cat(sprintf("   Q (analytical):      %.4f\n", single_result$Q_bar))
cat(sprintf("   ell_com:             %.4f\n", single_result$ell_com))
cat(sprintf("   Empirical bias:      %.4f\n", single_result$empirical_bias))
cat(sprintf("   Theoretical bias:    %.4f (= 0.5 * tr(RIV))\n", single_result$theoretical_bias))
cat(sprintf("   tr(RIV):             %.4f\n", single_result$tr_RIV))
cat(sprintf("   Bias ratio:          %.4f (should ≈ 1.0)\n", single_result$bias_ratio))
cat("\n")

# ============================================================================
# 5. Run Multiple Replications
# ============================================================================

cat("5. Running", test_config$n_reps, "replications...\n")

results <- vector("list", test_config$n_reps)
t_start <- Sys.time()

for (rep in seq_len(test_config$n_reps)) {
  if (rep %% 20 == 0) {
    cat(sprintf("   Progress: %d/%d\n", rep, test_config$n_reps))
  }

  results[[rep]] <- run_single_replication(
    rep_id = rep,
    mechanism = "MCAR",
    imputation_model = "matched",
    config = test_config,
    seeds = seeds,
    verbose = FALSE
  )
}

t_end <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))

cat(sprintf("   Completed in %.1f seconds (%.2f sec/rep)\n\n", elapsed, elapsed / test_config$n_reps))

# ============================================================================
# 6. Summarize Results
# ============================================================================

cat("6. Summarizing results...\n\n")

# Extract key metrics
empirical_bias <- sapply(results, `[[`, "empirical_bias")
theoretical_bias <- sapply(results, `[[`, "theoretical_bias")
tr_RIV <- sapply(results, `[[`, "tr_RIV")
term_A1 <- sapply(results, `[[`, "term_A1")
term_A2 <- sapply(results, `[[`, "term_A2")
bias_ratio <- sapply(results, `[[`, "bias_ratio")

cat("=== BIAS FORMULA VALIDATION ===\n\n")

cat("Mean Values:\n")
cat(sprintf("  Empirical bias:    %.4f (SD: %.4f)\n", mean(empirical_bias), sd(empirical_bias)))
cat(sprintf("  Theoretical bias:  %.4f (SD: %.4f)\n", mean(theoretical_bias), sd(theoretical_bias)))
cat(sprintf("  tr(RIV):           %.4f\n", mean(tr_RIV)))
cat("\n")

cat("Bias Decomposition:\n")
cat(sprintf("  Term 1 (Q - ell_com at θ_obs):    %.4f (expected ≈ tr(RIV) = %.4f)\n",
            mean(term_A1), mean(tr_RIV)))
cat(sprintf("  Term 2 (ell_com at θ_obs - θ_com): %.4f (expected ≈ -0.5*tr = %.4f)\n",
            mean(term_A2), -0.5 * mean(tr_RIV)))
cat(sprintf("  Total (= 0.5*tr(RIV)):            %.4f (expected = %.4f)\n",
            mean(empirical_bias), 0.5 * mean(tr_RIV)))
cat("\n")

cat("Ratios (all should ≈ 1.0):\n")
cat(sprintf("  Term 1 / tr(RIV):            %.4f\n", mean(term_A1) / mean(tr_RIV)))
cat(sprintf("  Term 2 / (-0.5*tr(RIV)):     %.4f\n", mean(term_A2) / (-0.5 * mean(tr_RIV))))
cat(sprintf("  Total / (0.5*tr(RIV)):       %.4f\n", mean(empirical_bias) / (0.5 * mean(tr_RIV))))
cat(sprintf("  Mean bias_ratio:             %.4f\n", mean(bias_ratio)))
cat("\n")

# Parameter recovery
a_bias <- sapply(results, `[[`, "a_bias")
b_bias <- sapply(results, `[[`, "b_bias")
a_covered <- sapply(results, `[[`, "a_covered")
b_covered <- sapply(results, `[[`, "b_covered")

cat("Parameter Recovery (FIML estimates):\n")
cat(sprintf("  a: mean bias = %+.4f, coverage = %.1f%%\n", mean(a_bias), 100 * mean(a_covered)))
cat(sprintf("  b: mean bias = %+.4f, coverage = %.1f%%\n", mean(b_bias), 100 * mean(b_covered)))
cat("\n")

# Overall assessment
ratio_check <- abs(mean(bias_ratio) - 1.0) < 0.10  # Within 10% of 1.0
cat("=== VALIDATION STATUS ===\n")
if (ratio_check) {
  cat("✓ PASS: Bias formula validated (ratio within 10% of 1.0)\n")
} else {
  cat("✗ FAIL: Bias ratio outside tolerance\n")
}

cat("\n=== Test Complete ===\n")
