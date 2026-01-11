# ============================================================================
# Test Clean Analytical Implementation with Proper vs Improper MI
# ============================================================================

rm(list = ls())

cat("=== Testing Analytical Implementation: Improper vs Proper MI ===\n\n")

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

n_reps <- 5000
master_seed <- 12345

# Generate seeds with SeedMaker
cat("Generating seeds with SeedMaker...\n")
seeds <- generate_seeds(master_seed = master_seed, n_reps = n_reps)
print_seed_summary(seeds, master_seed)

# Run simulation
cat(sprintf("\nRunning %d replications (n=500, MCAR monotone)...\n\n", n_reps))

results <- run_simulation(
  seeds = seeds,
  n = 500,
  mechanism = "MCAR",
  gamma = 0.5,
  cut1 = 0.5,
  cut2 = 0.75,
  true_params = true_params
)

cat("\n=== RESULTS ===\n\n")

# ---- Improper MI Validation ----
cat("=== IMPROPER MI (θ̃ fixed at θ̂_obs) ===\n\n")

cat("Bias Formula Validation:\n")
cat(sprintf("  Mean Term 1:        %.4f\n", mean(results$term1_improper)))
cat(sprintf("  Mean tr(RIV):       %.4f\n", mean(results$tr_RIV)))
cat(sprintf("  Ratio 1:            %.4f (should be ~1.0)\n\n", mean(results$term1_improper) / mean(results$tr_RIV)))

cat(sprintf("  Mean Term 2:        %.4f\n", mean(results$term2)))
cat(sprintf("  Expected (-0.5*tr): %.4f\n", -0.5 * mean(results$tr_RIV)))
cat(sprintf("  Ratio 2:            %.4f (should be ~1.0)\n\n", mean(results$term2) / (-0.5 * mean(results$tr_RIV))))

cat(sprintf("  Mean Total Bias:    %.4f\n", mean(results$total_improper)))
cat(sprintf("  Expected (0.5*tr):  %.4f\n", 0.5 * mean(results$tr_RIV)))
cat(sprintf("  Ratio Total:        %.4f (should be ~1.0)\n\n", mean(results$total_improper) / (0.5 * mean(results$tr_RIV))))

# ---- Proper MI Validation ----
cat("=== PROPER MI (E_θ̃[Q] via delta method) ===\n\n")

cat("Hessian Correction:\n")
cat(sprintf("  Mean correction:    %.4f (should be negative)\n", mean(results$hessian_correction)))
cat(sprintf("  SD correction:      %.4f\n\n", sd(results$hessian_correction)))

cat("Bias Formula Validation:\n")
cat(sprintf("  Mean Term 1:        %.4f\n", mean(results$term1_proper)))
cat(sprintf("  Mean tr(RIV):       %.4f\n", mean(results$tr_RIV)))
cat(sprintf("  Ratio 1:            %.4f (should be ~1.0)\n\n", mean(results$term1_proper) / mean(results$tr_RIV)))

cat(sprintf("  Mean Total Bias:    %.4f\n", mean(results$total_proper)))
cat(sprintf("  Expected (0.5*tr):  %.4f\n", 0.5 * mean(results$tr_RIV)))
cat(sprintf("  Ratio Total:        %.4f (should be ~1.0)\n\n", mean(results$total_proper) / (0.5 * mean(results$tr_RIV))))

# ---- Comparison ----
cat("=== IMPROPER vs PROPER COMPARISON ===\n\n")
cat(sprintf("  Q_improper - Q_proper:  %.4f (= -correction)\n",
            mean(results$Q_improper - results$Q_proper)))
cat(sprintf("  Term1_improper - Term1_proper: %.4f\n",
            mean(results$term1_improper - results$term1_proper)))

# ---- Parameter Recovery ----
cat("\n=== PARAMETER RECOVERY (FIML) ===\n\n")
cat(sprintf("  Mean bias in a:     %+.4f\n", mean(results$a_bias)))
cat(sprintf("  Mean bias in b:     %+.4f\n", mean(results$b_bias)))
cat("\n")

# ---- Validation Summary ----
ratio1_imp <- mean(results$term1_improper) / mean(results$tr_RIV)
ratio2 <- mean(results$term2) / (-0.5 * mean(results$tr_RIV))
ratio_total_imp <- mean(results$total_improper) / (0.5 * mean(results$tr_RIV))
ratio1_prop <- mean(results$term1_proper) / mean(results$tr_RIV)
ratio_total_prop <- mean(results$total_proper) / (0.5 * mean(results$tr_RIV))

cat("=== VALIDATION SUMMARY ===\n\n")

imp_pass <- abs(ratio1_imp - 1.0) < 0.10 && abs(ratio2 - 1.0) < 0.10 && abs(ratio_total_imp - 1.0) < 0.10
prop_pass <- abs(ratio1_prop - 1.0) < 0.10 && abs(ratio_total_prop - 1.0) < 0.10

if (imp_pass) {
  cat("IMPROPER MI: PASSED (all ratios within 10% of 1.0)\n")
} else {
  cat("IMPROPER MI: FAILED (at least one ratio outside tolerance)\n")
}

if (prop_pass) {
  cat("PROPER MI:   PASSED (all ratios within 10% of 1.0)\n")
} else {
  cat("PROPER MI:   FAILED (at least one ratio outside tolerance)\n")
}

# Verify correction is negative
if (mean(results$hessian_correction) < 0) {
  cat("HESSIAN CORRECTION: PASSED (negative as expected)\n")
} else {
  cat("HESSIAN CORRECTION: FAILED (should be negative)\n")
}
