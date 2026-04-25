# ============================================================================
# Full Study 2 Simulation, Amelia Backend (mirror of full_grid.R)
# ============================================================================
# N x miss-rate grid: {100, 250, 500, 1000} x {0.10, 0.25, 0.40}
# Imputations:  M = 100
# Replications: 2000 per condition
# Imputation:   amelia (joint MVN, EMB, empri = 0.01 * N for stability)
# Cores:        honors SLURM_CPUS_PER_TASK (100 under run_full_grid_amelia.sh)
# Expected wall time: ~8-12 hours on 100 cores
#
# Design rationale: amelia is congenial with the analysis CFA family
# (joint MVN nests every candidate model).  This run pairs with the
# existing PMM grid (results-full-M100/) for the manuscript's
# congenial-vs-uncongenial comparison.
# ============================================================================

RESULTS_DIR   <- "/biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full-M100-amelia"
M_IMPUTATIONS <- 100L

# 1. Install/refresh miicsem from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("Amelia", quietly = TRUE)) {
  install.packages("Amelia", repos = "https://cloud.r-project.org")
}
devtools::install_github("marcus-waldman/MI-IC",
                         subdir  = "packages/miicsem",
                         upgrade = "never",
                         force   = TRUE,
                         quiet   = TRUE)

library(miicsem)

cat("=== miicsem full-grid simulation (AMELIA backend) ===\n")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))
cat(sprintf("Amelia version : %s\n", as.character(packageVersion("Amelia"))))
cat(sprintf("Results directory: %s\n", RESULTS_DIR))
cat(sprintf("SLURM cpus-per-task: %s\n",
            Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>")))
cat(sprintf("default_n_cores()  : %d\n\n", default_n_cores()))

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps       = 2000,
  sample_sizes = c(100, 250, 500, 1000),
  miss_rates   = c(0.10, 0.25, 0.40),
  M            = M_IMPUTATIONS,
  mice_method  = "amelia",
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

# 4. Decomposition T1/T3 means for M1 (sanity check)
cat("\n=== Three-term decomp for M1 across conditions (log-lik scale) ===\n")
cat(sprintf("  %-14s  %8s  %8s  %8s  %8s\n",
            "condition", "tr(RIV)", "T1", "T3", "T1/tr(RIV)"))
for (cond_label in names(res)) {
  rows <- lapply(res[[cond_label]], function(r) {
    if (is.null(r) || is.null(r$dev_df)) return(NULL)
    d <- r$dev_df["M1", ]
    if (any(is.na(c(d$loglik_com_at_com, d$loglik_com_at_pooled,
                    d$mean_loglik_at_pooled)))) return(NULL)
    data.frame(
      tr_RIV = d$tr_RIV,
      t1     = d$mean_loglik_at_pooled - d$loglik_com_at_pooled,
      t3     = d$loglik_com_at_obs    - d$loglik_com_at_com
    )
  })
  df <- do.call(rbind, rows)
  if (!is.null(df) && nrow(df) > 0) {
    r1 <- mean(df$t1) / mean(df$tr_RIV)
    cat(sprintf("  %-14s  %8.3f  %+8.3f  %+8.3f  %+8.3f\n",
                cond_label, mean(df$tr_RIV), mean(df$t1), mean(df$t3), r1))
  }
}

cat(sprintf("\nResults written to %s\n", RESULTS_DIR))
cat("\n=== Full grid (amelia) simulation complete ===\n")
