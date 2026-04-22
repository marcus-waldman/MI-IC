# ============================================================================
# Three-Term Bias Decomposition HPC Run (scaled-down Study 2 grid)
# ============================================================================
# Goal: identify which of Term 1, 2, or 3 has the wrong sign driving the
# empirical Q_bar < ell_com result (opposite to theory's +1/2 tr(RIV)).
#
#   Term 1 (Imputation) = Q_bar(theta_pooled) - ell_com(theta_pooled)
#   Term 2 (Pooling)    = ell_com(theta_pooled) - ell_com(theta_obs_fiml)
#   Term 3 (Estimation) = ell_com(theta_obs_fiml) - ell_com(theta_com)
#
# Scaled design (targeted, not full grid):
#   N            in {100, 250, 1000}   (drop 500 — between 250 and 1000)
#   miss_rate    = 0.40                 (max tr(RIV), max signal)
#   M            = 50
#   n_reps       = 500
#
# Output:
#   /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-decomp/
#     results_n=100_mr=0.40.rds
#     results_n=250_mr=0.40.rds
#     results_n=1000_mr=0.40.rds
#     results_combined.rds
#
# Expected wall time: ~1 hour on 100 cores (generous).
# ============================================================================

RESULTS_DIR   <- "/biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-decomp"
M_IMPUTATIONS <- 50L
N_REPS        <- 500L
SAMPLE_SIZES  <- c(100, 250, 1000)
MISS_RATES    <- 0.40

# 1. Install/refresh miicsem from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::install_github("marcus-waldman/MI-IC",
                         subdir  = "packages/miicsem",
                         upgrade = "never",
                         force   = TRUE,
                         quiet   = TRUE)

library(miicsem)

cat("=== miicsem three-term decomposition simulation ===\n")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))
cat(sprintf("Results directory: %s\n", RESULTS_DIR))
cat(sprintf("SLURM cpus-per-task: %s\n",
            Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>")))
cat(sprintf("default_n_cores()  : %d\n\n", default_n_cores()))

cat("Design:\n")
cat(sprintf("  N            = {%s}\n", paste(SAMPLE_SIZES, collapse = ", ")))
cat(sprintf("  miss_rate    = %s\n", paste(MISS_RATES, collapse = ", ")))
cat(sprintf("  M            = %d\n", M_IMPUTATIONS))
cat(sprintf("  n_reps       = %d per condition\n\n", N_REPS))

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps       = N_REPS,
  sample_sizes = SAMPLE_SIZES,
  miss_rates   = MISS_RATES,
  M            = M_IMPUTATIONS,
  results_dir  = RESULTS_DIR,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]

cat(sprintf("\nTotal wall time: %.1f seconds (%.2f hours)\n",
            elapsed, elapsed / 3600))

# 3. Summary per condition
cat("\n=== Conditions completed ===\n")
for (cond_label in names(res)) {
  n_ok <- sum(!vapply(res[[cond_label]], is.null, logical(1)))
  cat(sprintf("  %-20s: %d/%d reps succeeded\n",
              cond_label, n_ok, length(res[[cond_label]])))
}

# 4. Three-term decomposition means for M1 (true model, log-lik scale)
cat("\n=== Mean three-term decomposition for M1 (log-lik scale) ===\n")
cat("  (Term 1 = Imputation, Term 2 = Pooling, Term 3 = Estimation)\n")
cat("  (Theory: Term 1 = +tr(RIV), Term 2 ~ 0, Term 3 = -tr(RIV)/2)\n\n")
cat(sprintf("  %-20s  %7s  %7s  %7s  %7s  %7s\n",
            "condition", "tr_RIV", "Term 1", "Term 2", "Term 3", "Total"))
for (cond_label in names(res)) {
  rows <- lapply(res[[cond_label]], function(r) {
    if (is.null(r) || is.null(r$dev_df)) return(NULL)
    d <- r$dev_df["M1", ]
    if (any(is.na(c(d$loglik_com_at_com, d$loglik_com_at_obs,
                    d$loglik_com_at_pooled, d$mean_loglik_at_pooled)))) {
      return(NULL)
    }
    data.frame(
      tr_RIV = d$tr_RIV,
      term1  = d$mean_loglik_at_pooled - d$loglik_com_at_pooled,
      term2  = d$loglik_com_at_pooled  - d$loglik_com_at_obs,
      term3  = d$loglik_com_at_obs     - d$loglik_com_at_com,
      total  = d$mean_loglik_at_pooled - d$loglik_com_at_com
    )
  })
  df <- do.call(rbind, rows)
  if (!is.null(df) && nrow(df) > 0) {
    cat(sprintf("  %-20s  %7.3f  %+7.3f  %+7.3f  %+7.3f  %+7.3f\n",
                cond_label,
                mean(df$tr_RIV),
                mean(df$term1), mean(df$term2),
                mean(df$term3), mean(df$total)))
  }
}

cat(sprintf("\nResults written to %s\n", RESULTS_DIR))
cat("  results_<n>_<mr>.rds  per-condition checkpoints\n")
cat("  results_combined.rds  rolled-up list\n")

cat("\n=== Decomposition simulation complete ===\n")
