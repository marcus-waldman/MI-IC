# ============================================================================
# Local Congeniality Test: pmm vs norm vs congenial MVN-from-saturated
# ============================================================================
# Runs the Study 2 decomposition pipeline three times with different
# imputation backends to isolate whether uncongeniality causes the
# empirical Term 1 < 0 result.
#
# Design:
#   N             = 250  (single condition — matches decomp baseline)
#   miss_rate     = 0.40 (max signal)
#   M             = 50
#   reps          = 2000 per method   (Study 2 convention)
#   cores         = 16   (local)
#
# Three methods run sequentially:
#   1. "pmm"      — current Study 2 default (predictive mean matching)
#   2. "norm"     — Bayesian Gaussian regression per variable
#   3. "mvn_msat" — congenial imputation from saturated MVN (FIML)
#
# Results saved to hpc/results-congeniality-<method>/ with per-condition
# and combined .rds files, identical layout to the main decomp run.
#
# Expected wall time on 16 cores: ~10-15 hours total (overnight).
# Per-method resume support: re-run the script and any completed
# per-condition .rds is loaded cached; only unfinished methods/conditions
# are computed.
# ============================================================================

library(miicsem)

stopifnot(as.character(packageVersion("miicsem")) >= "0.4.1")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

N_REPS       <- 2000L
N_CORES      <- 16L
SAMPLE_SIZES <- 250
MISS_RATES   <- 0.40
M_IMPUTATIONS <- 50L

METHODS <- c("pmm", "norm", "mvn_msat")

overall_t0 <- proc.time()

for (method in METHODS) {
  results_dir <- sprintf("hpc/results-congeniality-%s", method)
  cat("\n======================================================\n")
  cat(sprintf("  Imputation method: %s\n", method))
  cat(sprintf("  Results dir     : %s\n", results_dir))
  cat("======================================================\n")

  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  t0 <- proc.time()
  res <- run_simulation(
    n_reps       = N_REPS,
    n_cores      = N_CORES,
    sample_sizes = SAMPLE_SIZES,
    miss_rates   = MISS_RATES,
    M            = M_IMPUTATIONS,
    mice_method  = method,
    results_dir  = results_dir,
    verbose      = TRUE
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("\nMethod '%s' wall time: %.1f seconds (%.2f hours)\n",
              method, elapsed, elapsed / 3600))

  # Quick per-method M1 decomposition summary
  rows <- lapply(res[[1]], function(r) {
    if (is.null(r) || is.null(r$dev_df)) return(NULL)
    d <- r$dev_df["M1", ]
    need <- c("loglik_com_at_com", "loglik_com_at_obs",
              "loglik_com_at_pooled", "mean_loglik_at_pooled", "tr_RIV")
    if (any(is.na(d[need]))) return(NULL)
    data.frame(
      tr_RIV = d$tr_RIV,
      term1  = d$mean_loglik_at_pooled - d$loglik_com_at_pooled,
      term2  = d$loglik_com_at_pooled  - d$loglik_com_at_obs,
      term3  = d$loglik_com_at_obs     - d$loglik_com_at_com,
      total  = d$mean_loglik_at_pooled - d$loglik_com_at_com
    )
  })
  df <- do.call(rbind, rows)
  cat(sprintf("  Reps ok: %d / %d\n", nrow(df), N_REPS))
  cat(sprintf("  mean tr(RIV): %.3f\n", mean(df$tr_RIV)))
  cat(sprintf("  Term 1       : %+.3f  (theory congenial: +%.3f)\n",
              mean(df$term1), mean(df$tr_RIV)))
  cat(sprintf("  Term 2       : %+.3f  (theory: ~0)\n", mean(df$term2)))
  cat(sprintf("  Term 3       : %+.3f  (theory: %+.3f)\n",
              mean(df$term3), -mean(df$tr_RIV) / 2))
  cat(sprintf("  Total        : %+.3f  (theory: %+.3f)\n",
              mean(df$total), mean(df$tr_RIV) / 2))
}

total_elapsed <- (proc.time() - overall_t0)[3]
cat(sprintf("\n\n====================================================\n"))
cat(sprintf("  All 3 methods complete.\n"))
cat(sprintf("  Total wall time: %.1f seconds (%.2f hours)\n",
            total_elapsed, total_elapsed / 3600))
cat(sprintf("====================================================\n"))
cat("Analyze with: Rscript hpc/analyze_congeniality.R\n")
