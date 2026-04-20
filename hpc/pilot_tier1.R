# ============================================================================
# Tier 1 Pilot: 10 reps at N=250, 25% MCAR
# ============================================================================
# Quick smoke test of the miicsem pipeline. Installs the package + all
# runtime dependencies if missing. Expected wall time ~100s on 32 cores.
# ============================================================================

RESULTS_DIR <- "/biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results/"

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

# 2. Run the pilot (runtime deps auto-install if missing)
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

cat("=== Tier 1 pilot: 10 reps, N=250, 25% MCAR ===\n")
cat(sprintf("Results directory: %s\n\n", RESULTS_DIR))

t0 <- proc.time()
pilot <- run_simulation(
  n_reps       = 10,
  sample_sizes = 250,
  miss_rates   = 0.25,
  results_dir  = RESULTS_DIR
)
elapsed <- (proc.time() - t0)[3]

cat(sprintf("\nTotal wall time: %.1fs\n", elapsed))

# 3. Quick inspection of the first rep
cat("\n=== First replication, dev_df ===\n")
print(round(pilot[[1]][[1]]$dev_df, 2))

cat("\n=== First replication, chi2_df ===\n")
print(round(pilot[[1]][[1]]$chi2_df, 2))

cat("\n=== Chi-squares for M1 (true model, df=22) across reps ===\n")
chi_M1 <- do.call(rbind, lapply(pilot[[1]], function(r) {
  if (is.null(r) || is.null(r$chi2_df)) return(NULL)
  data.frame(chi2_com = r$chi2_df["M1", "chi2_com"],
             chi2_MI  = r$chi2_df["M1", "chi2_MI"],
             chi2_D3  = r$chi2_df["M1", "chi2_D3"])
}))
cat(sprintf("  chi2_com: mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_com), sd(chi_M1$chi2_com)))
cat(sprintf("  chi2_MI : mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_MI), sd(chi_M1$chi2_MI)))
cat(sprintf("  chi2_D3 : mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_D3), sd(chi_M1$chi2_D3)))

cat(sprintf("\nFiles written to %s:\n", RESULTS_DIR))
print(list.files(RESULTS_DIR))

cat("\n=== Pilot complete ===\n")
