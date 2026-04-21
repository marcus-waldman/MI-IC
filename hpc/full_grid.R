# ============================================================================
# Full Study 2 Simulation: 2000 reps x 12 conditions
# ============================================================================
# N      x missing rate grid: {100, 250, 500, 1000} x {0.10, 0.25, 0.40}
# Imputations:  M = 100
# Replications: 2000 per condition
# Cores:        honors SLURM_CPUS_PER_TASK (100 under run_full_grid.sh)
# Expected wall time: ~18-20 hours on 100 cores (5x the M=20 run)
# ============================================================================

RESULTS_DIR <- "/biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full-M100"
M_IMPUTATIONS <- 100L

# 1. Install/refresh miicsem from GitHub (always; catches version bumps)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::install_github("marcus-waldman/MI-IC",
                         subdir  = "packages/miicsem",
                         upgrade = "never",
                         force   = TRUE,
                         quiet   = TRUE)

library(miicsem)

cat("=== miicsem full-grid simulation ===\n")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))
cat(sprintf("Results directory: %s\n", RESULTS_DIR))
cat(sprintf("SLURM cpus-per-task: %s\n",
            Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>")))
cat(sprintf("default_n_cores()  : %d\n\n", default_n_cores()))

# 2. Run the simulation (runtime deps auto-install if missing)
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

t0 <- proc.time()
res <- run_simulation(
  n_reps       = 2000,
  sample_sizes = c(100, 250, 500, 1000),
  miss_rates   = c(0.10, 0.25, 0.40),
  M            = M_IMPUTATIONS,
  results_dir  = RESULTS_DIR,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]

cat(sprintf("\nTotal wall time: %.1f seconds (%.2f hours)\n",
            elapsed, elapsed / 3600))

# 3. Quick summary
cat("\n=== Conditions completed ===\n")
for (cond_label in names(res)) {
  n_ok <- sum(!vapply(res[[cond_label]], is.null, logical(1)))
  cat(sprintf("  %-20s: %d/%d reps succeeded\n",
              cond_label, n_ok, length(res[[cond_label]])))
}

# 4. Chi-square means for M1 across each condition (sanity check)
cat("\n=== chi-square means for M1 (true model, df=22) ===\n")
cat(sprintf("  %-14s  %8s  %8s  %8s\n",
            "condition", "chi2_com", "chi2_MI", "chi2_D3"))
for (cond_label in names(res)) {
  chi <- do.call(rbind, lapply(res[[cond_label]], function(r) {
    if (is.null(r) || is.null(r$chi2_df)) return(NULL)
    data.frame(chi2_com = r$chi2_df["M1", "chi2_com"],
               chi2_MI  = r$chi2_df["M1", "chi2_MI"],
               chi2_D3  = r$chi2_df["M1", "chi2_D3"])
  }))
  if (!is.null(chi) && nrow(chi) > 0) {
    cat(sprintf("  %-14s  %8.2f  %8.2f  %8.2f\n",
                cond_label,
                mean(chi$chi2_com, na.rm = TRUE),
                mean(chi$chi2_MI,  na.rm = TRUE),
                mean(chi$chi2_D3,  na.rm = TRUE)))
  }
}

cat(sprintf("\nResults written to %s\n", RESULTS_DIR))
cat("  results_<n>_<mr>.rds  per-condition checkpoints\n")
cat("  results_combined.rds  rolled-up list\n")

cat("\n=== Full grid simulation complete ===\n")
