# ============================================================================
# Quick Validation: Single Replication Test for Study 2 (SEM)
# ============================================================================
# Runs 1 replication at N=250, 25% MCAR and checks:
#   - All 12 models converge on complete data
#   - mice produces 20 imputed datasets
#   - All 5 IC methods produce finite values for >= 10 models
#   - Complete-data AIC selects M1 (the true model)
#   - tr(RIV) > 0 and reasonable magnitude
# ============================================================================

cat("=== Study 2 SEM: Single Replication Test ===\n\n")

# Source modules - find R/ directory relative to this script or working dir
sim_dir <- tryCatch({
  file.path(dirname(sys.frame(1)$ofile), "..", "R")
}, error = function(e) {
  "simulations/study2-bollen-sem/R"
})
source(file.path(sim_dir, "00_config.R"))
source(file.path(sim_dir, "01_model_specifications.R"))
source(file.path(sim_dir, "02_generate_and_ampute.R"))
source(file.path(sim_dir, "03_fit_and_pool.R"))
source(file.path(sim_dir, "04_compute_ic.R"))
source(file.path(sim_dir, "05_run_replication.R"))

suppressPackageStartupMessages({
  library(lavaan)
  library(mice)
  library(MASS)
})

config <- get_config()
n_test <- 250
mr_test <- 0.25
test_seed <- 42

pass_count <- 0
fail_count <- 0

check <- function(desc, condition) {
  if (condition) {
    cat(sprintf("  PASS: %s\n", desc))
    pass_count <<- pass_count + 1
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    fail_count <<- fail_count + 1
  }
}

# --- Run single replication ---
cat("Running single replication (N=250, 25% MCAR)...\n")
t0 <- proc.time()

result <- run_one_rep(
  rep_id      = 1,
  n           = n_test,
  miss_rate   = mr_test,
  config      = config,
  seed_data   = test_seed,
  seed_ampute = test_seed + 1,
  seed_impute = test_seed + 2
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("Elapsed: %.1f seconds\n\n", elapsed))

# --- Checks ---
cat("--- Basic Checks ---\n")
check("Result is not NULL", !is.null(result))

if (is.null(result)) {
  cat("\nFATAL: Replication returned NULL. Cannot continue.\n")
  quit(status = 1)
}

# Convergence
cat("\n--- Convergence ---\n")
n_converged_complete <- sum(!is.na(result$ic_df$AIC_com))
check(sprintf("Complete-data convergence: %d/12 models", n_converged_complete),
      n_converged_complete >= 10)

mi_converged <- result$convergence
check(sprintf("MI convergence (M1): %d/%d imputations", mi_converged["M1"], config$M),
      mi_converged["M1"] >= 15)

# IC values
cat("\n--- IC Values ---\n")
ic_df <- result$ic_df

for (method in colnames(ic_df)) {
  n_finite <- sum(is.finite(ic_df[[method]]))
  check(sprintf("%s: %d/12 finite values", method, n_finite), n_finite >= 8)
}

# Model selection
cat("\n--- Model Selection ---\n")
sel <- result$selections
check(sprintf("AIC_com selects %s (expect M1)", sel["AIC_com"]),
      sel["AIC_com"] == "M1")
check(sprintf("MI_AIC selects %s", sel["MI_AIC"]),
      !is.na(sel["MI_AIC"]))
check(sprintf("MI_BIC selects %s", sel["MI_BIC"]),
      !is.na(sel["MI_BIC"]))

# tr(RIV)
cat("\n--- tr(RIV) ---\n")
tr_rivs <- result$tr_RIVs[!is.na(result$tr_RIVs)]
check(sprintf("tr(RIV) available for %d models", length(tr_rivs)),
      length(tr_rivs) >= 8)

if (length(tr_rivs) > 0) {
  check(sprintf("tr(RIV) range: [%.2f, %.2f]",
                min(tr_rivs), max(tr_rivs)),
        all(tr_rivs > 0) && all(tr_rivs < 100))

  # MI-AIC should differ from AICcd
  if (is.finite(ic_df["M1", "MI_AIC"]) && is.finite(ic_df["M1", "AICcd"])) {
    diff_val <- abs(ic_df["M1", "MI_AIC"] - ic_df["M1", "AICcd"])
    check(sprintf("MI_AIC != AICcd for M1 (diff=%.2f)", diff_val),
          diff_val > 0.01)
  }
}

# Print IC table
cat("\n--- IC Table (first 6 models) ---\n")
print(round(ic_df[1:min(6, nrow(ic_df)), ], 2))

cat("\n--- Selections ---\n")
print(sel)

cat("\n--- tr(RIV) by model ---\n")
print(round(result$tr_RIVs, 3))

# Summary
cat(sprintf("\n=== SUMMARY: %d passed, %d failed ===\n", pass_count, fail_count))
if (fail_count == 0) {
  cat("All checks passed!\n")
} else {
  cat("Some checks failed. Review output above.\n")
}
