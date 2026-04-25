# ============================================================================
# Study 2 Simulation: Amelia × empri Sweep (Congeniality Dose-Response)
# ============================================================================
# Design (12 conditions per empri level x 2 empri levels = 24 cells):
#   N           in {100, 500, 1000}
#   miss_rate   in {0.20, 0.40}
#   empri_frac  in {0.1, 100}     (multiplied by N to give Amelia's empri)
# Imputations:  M = 100
# Replications: 2000 per condition
# Imputer:      Amelia (joint MVN, EMB)
# Cores:        honors SLURM_CPUS_PER_TASK (100 under run_full_grid_amelia.sh)
# Expected wall time: ~10-12 hours on 100 cores
#
# Why empri varies as a simulation factor:
#   In Amelia, empri is the prior dof of an inverse-Wishart prior on
#   Sigma centered at a DIAGONAL matrix.  Small empri (e.g. 0.1*N) merely
#   stabilizes EM at small N / high missingness.  Large empri (e.g.
#   100*N) drowns out the data and shrinks Sigma* toward independence,
#   intentionally introducing UNCONGENIALITY with a structured analysis
#   model like the 3-factor CFA.  This produces a continuous handle on
#   congeniality within a single imputation framework — a much cleaner
#   contrast than amelia-vs-PMM (which differ in many ways at once).
#
# Predicted dose-response on Term 1 (Imputation Bias):
#   empri = 0.1*N  (mostly congenial)   ->  r_term1 ~ +1
#   empri = 100*N  (severely uncongenial) ->  r_term1 << 0  (sign-flipped)
#
# Output layout:
#   results-amelia-empri0.1/<cond>.rds       12 conds x 2000 reps each
#   results-amelia-empri100/<cond>.rds       12 conds x 2000 reps each
# ============================================================================

EMPRI_LEVELS  <- c(0.1, 100)
M_IMPUTATIONS <- 100L
N_REPS        <- 2000L
SAMPLE_SIZES  <- c(100, 500, 1000)
MISS_RATES    <- c(0.20, 0.40)
RESULTS_BASE  <- "/biostats_share/waldmanm/simulation-studies/MI-IC/SeM"

# 1. Install/refresh miicsem from GitHub.  As of miicsem 0.5.1, Amelia
# is in Imports, so install_github will pull it automatically.
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::install_github("marcus-waldman/MI-IC",
                         subdir       = "packages/miicsem",
                         upgrade      = "never",
                         force        = TRUE,
                         dependencies = TRUE,
                         quiet        = TRUE)

library(miicsem)

cat("=== miicsem amelia x empri sweep ===\n")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))
cat(sprintf("Amelia version : %s\n", as.character(packageVersion("Amelia"))))
cat(sprintf("SLURM cpus-per-task: %s\n",
            Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>")))
cat(sprintf("default_n_cores()  : %d\n\n", default_n_cores()))

overall_t0 <- proc.time()

for (efrac in EMPRI_LEVELS) {
  results_dir <- file.path(
    RESULTS_BASE,
    sprintf("results-amelia-empri%s", format(efrac, nsmall = 0, trim = TRUE))
  )
  cat("\n##########################################################\n")
  cat(sprintf("###  empri = %g * N    -> %s\n", efrac, results_dir))
  cat("##########################################################\n")
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  t0 <- proc.time()
  res <- run_simulation(
    n_reps            = N_REPS,
    sample_sizes      = SAMPLE_SIZES,
    miss_rates        = MISS_RATES,
    M                 = M_IMPUTATIONS,
    mice_method       = "amelia",
    amelia_empri_frac = efrac,
    results_dir       = results_dir,
    verbose           = TRUE
  )
  el <- (proc.time() - t0)[3]
  cat(sprintf("\n[empri=%g*N] wall time: %.1fs (%.2fh)\n",
              efrac, el, el / 3600))

  # Per-empri Term 1 summary across the 6 cells
  cat(sprintf("\n=== Term 1 summary by condition (empri = %g * N) ===\n", efrac))
  cat(sprintf("  %-14s  %8s  %8s  %8s  %10s\n",
              "condition", "n_reps", "tr(RIV)", "T1", "r_term1"))
  for (cond_label in names(res)) {
    rows <- lapply(res[[cond_label]], function(r) {
      if (is.null(r) || is.null(r$dev_df)) return(NULL)
      d <- r$dev_df["M1", ]
      if (any(is.na(c(d$loglik_com_at_pooled,
                      d$mean_loglik_at_pooled, d$tr_RIV)))) return(NULL)
      data.frame(
        tr_RIV = d$tr_RIV,
        t1     = d$mean_loglik_at_pooled - d$loglik_com_at_pooled
      )
    })
    df <- do.call(rbind, rows)
    if (!is.null(df) && nrow(df) > 0) {
      cat(sprintf("  %-14s  %8d  %8.3f  %+8.3f  %+10.3f\n",
                  cond_label, nrow(df),
                  mean(df$tr_RIV), mean(df$t1),
                  mean(df$t1) / mean(df$tr_RIV)))
    }
  }
}

total_el <- (proc.time() - overall_t0)[3]
cat(sprintf("\n\n======================================================\n"))
cat(sprintf("  All empri levels complete.\n"))
cat(sprintf("  Total wall time: %.1fs (%.2fh)\n",
            total_el, total_el / 3600))
cat(sprintf("======================================================\n"))
