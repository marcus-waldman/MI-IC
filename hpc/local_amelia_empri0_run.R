# ============================================================================
# Local verify: Amelia with empri = 0 (no prior).
#
# Question: does the formula  Δ ≈ 2 · [tr(RIV)_Msat − tr(RIV)_M1]  hold for
# Amelia once the empri prior is removed entirely?
#
# Why this matters: with empri = 0.01·N (mild prior), Δ_emp / Δ_pred ≈ 0.74
# across mr values — the formula systematically over-predicts. Hypothesis:
# the prior shrinks imputation noise in Φ\Θ directions by a fixed
# multiplicative factor. With empri = 0, the prior is gone and the formula
# should match within Monte Carlo error (the same way it matches PMM).
#
# Design (single N, varied mr):
#   N             = 500
#   miss_rate     ∈ {0.10, 0.25, 0.40}
#   M             = 50
#   reps          = 250 per cell
#   cores         = 16
#   empri_frac    = 0       (no prior)
#
# Output:
#   hpc/results-amelia-empri0/results_n=500_mr=0.10.rds
#   hpc/results-amelia-empri0/results_n=500_mr=0.25.rds
#   hpc/results-amelia-empri0/results_n=500_mr=0.40.rds
#   hpc/results-amelia-empri0/results_combined.rds
# ============================================================================

library(miicsem)
stopifnot(as.character(packageVersion("miicsem")) >= "0.4.2")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

N            <- 500L
MR_VALUES    <- c(0.10, 0.25, 0.40)
M_IMP        <- 50L
N_REPS       <- 250L
N_CORES      <- 16L
BASE_SEED    <- 32897891L
EMPRI_FRAC   <- 0           # the whole point: no prior

results_dir  <- "hpc/results-amelia-empri0"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps            = N_REPS,
  n_cores           = N_CORES,
  seed              = BASE_SEED,
  sample_sizes      = N,
  miss_rates        = MR_VALUES,
  M                 = M_IMP,
  mice_method       = "amelia",
  amelia_empri_frac = EMPRI_FRAC,
  results_dir       = results_dir,
  verbose           = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nWall time: %.1fs (%.2fh)\n", elapsed, elapsed / 3600))

# --- quick Δ summary --------------------------------------------------------
cat("\n=== Empirical Δ vs predicted 2·Δtr(RIV) ===\n")
cat(sprintf("%-6s  %5s  %8s  %10s  %8s  %10s  %10s  %8s  %8s\n",
            "mr", "n_rep", "chi2_com", "chi2_adhoc", "chi2_MI",
            "tr_RIV_M1", "tr_RIV_Msat", "Δ_emp", "Δ_pred"))

for (mr in MR_VALUES) {
  cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
  reps <- res[[cell_label]]
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    chi <- rr$chi2_df; dev <- rr$dev_df
    if (is.null(chi) || is.null(dev) || !"M1" %in% rownames(chi)) next
    if (!"Msat" %in% rownames(dev)) next
    sat_dev_adhoc <- dev["Msat", "DEV_adhoc"]
    if (is.na(sat_dev_adhoc)) next
    r <- chi["M1", ]; d <- dev["M1", ]; s <- dev["Msat", ]
    chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
    if (any(is.na(c(d$tr_RIV, s$tr_RIV, r$chi2_com, r$chi2_MI)))) next
    k <- k + 1
    rows[[k]] <- data.frame(
      chi2_com = r$chi2_com, chi2_adhoc = chi2_adhoc, chi2_MI = r$chi2_MI,
      tr_M1 = d$tr_RIV, tr_Msat = s$tr_RIV
    )
  }
  if (length(rows) == 0) {
    cat(sprintf("%-6.2f  %5s  (no usable reps)\n", mr, "0"))
    next
  }
  df <- do.call(rbind, rows)
  e_com   <- mean(df$chi2_com)
  e_adhoc <- mean(df$chi2_adhoc)
  e_MI    <- mean(df$chi2_MI)
  e_trM1  <- mean(df$tr_M1)
  e_trMs  <- mean(df$tr_Msat)
  delta_emp  <- e_adhoc - e_com
  delta_pred <- 2 * (e_trMs - e_trM1)
  cat(sprintf("%-6.2f  %5d  %8.2f  %10.2f  %8.2f  %10.3f  %10.3f  %8.2f  %8.2f\n",
              mr, nrow(df), e_com, e_adhoc, e_MI, e_trM1, e_trMs,
              delta_emp, delta_pred))
}

cat(sprintf("\nDone. Results saved to %s\n", results_dir))
