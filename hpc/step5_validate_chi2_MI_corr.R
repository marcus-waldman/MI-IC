# ============================================================================
# Step 5 (cozy-imagining-waterfall.md): end-to-end Type I validation of
# chi2_MI_corr (v4.5 §13 scaled-shifted correction).
#
# Why: miicsem 0.6.0 added chi2_MI_corr per rep. Existing results pre-date
# 0.6.0 and don't have the new column. To validate Type I error of
# chi2_MI_corr at nominal alpha=0.05, we need to re-run the high-stakes cells.
#
# Cells (validation grid):
#   PMM,    N in {250, 500}, mr in {0.25, 0.40}, M=50, 250 reps each
#   Amelia, N in {250, 500}, mr in {0.25, 0.40}, M=50, 250 reps each
# Total: 8 cells x 250 reps = 2000 reps
# Wall time: ~2-4 hours on 16 cores per local_amelia_empri0_run.R benchmark.
#
# Predicted Type I:
#   chi2_adhoc:    18-21% (4x nominal -- standard practice)
#   chi2_MI:       9.8-11% (v4 first-moment correction)
#   chi2_MI_corr:  ~5% (v4.5 finite-M correction -- this is what we test)
#   chi2_com:      6-7% (oracle, finite-N Bartlett residual)
#   chi2_D3:       14-15%
#
# Acceptance: chi2_MI_corr within +/- 1pp of nominal 5% across all cells.
#
# Output: hpc/results-step5/results_*.rds plus combined summary
# ============================================================================

library(miicsem)
stopifnot(as.character(packageVersion("miicsem")) >= "0.6.0")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

N_VALUES   <- c(250L, 500L)
MR_VALUES  <- c(0.25, 0.40)
M_IMP      <- 50L
N_REPS     <- 250L
N_CORES    <- 16L
BASE_SEED  <- 32897891L

out_root <- "hpc/results-step5"
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# ---- PMM cells -------------------------------------------------------------
pmm_dir <- file.path(out_root, "pmm")
if (!dir.exists(pmm_dir)) dir.create(pmm_dir, recursive = TRUE)

cat("\n========================================\n")
cat("  Cells 1-4: PMM, N in {250,500}, mr in {0.25,0.40}, M=50, 250 reps\n")
cat("========================================\n")
t0 <- proc.time()
res_pmm <- run_simulation(
  n_reps            = N_REPS,
  n_cores           = N_CORES,
  seed              = BASE_SEED,
  sample_sizes      = N_VALUES,
  miss_rates        = MR_VALUES,
  M                 = M_IMP,
  mice_method       = "pmm",
  results_dir       = pmm_dir,
  verbose           = TRUE
)
elapsed_pmm <- (proc.time() - t0)[3]
cat(sprintf("\nPMM cells wall time: %.1fs (%.2fh)\n",
            elapsed_pmm, elapsed_pmm / 3600))

# ---- Amelia empri=0 cells --------------------------------------------------
amelia_dir <- file.path(out_root, "amelia_empri0")
if (!dir.exists(amelia_dir)) dir.create(amelia_dir, recursive = TRUE)

cat("\n========================================\n")
cat("  Cells 5-8: Amelia empri=0, N in {250,500}, mr in {0.25,0.40}, M=50, 250 reps\n")
cat("========================================\n")
t0 <- proc.time()
res_amelia <- run_simulation(
  n_reps            = N_REPS,
  n_cores           = N_CORES,
  seed              = BASE_SEED,
  sample_sizes      = N_VALUES,
  miss_rates        = MR_VALUES,
  M                 = M_IMP,
  mice_method       = "amelia",
  amelia_empri_frac = 0,
  results_dir       = amelia_dir,
  verbose           = TRUE
)
elapsed_amelia <- (proc.time() - t0)[3]
cat(sprintf("\nAmelia cells wall time: %.1fs (%.2fh)\n",
            elapsed_amelia, elapsed_amelia / 3600))

# ---- Type I summary --------------------------------------------------------
cat("\n========================================\n")
cat("  Type I error at alpha = 0.05 (true model M1)\n")
cat("========================================\n\n")

extract_M1 <- function(reps) {
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    chi <- rr$chi2_df; dev <- rr$dev_df
    if (is.null(chi) || !"M1" %in% rownames(chi)) next
    if (!"Msat" %in% rownames(dev)) next
    sat_dev_adhoc <- dev["Msat", "DEV_adhoc"]
    r <- chi["M1", ]
    chi2_adhoc <- dev["M1", "DEV_adhoc"] - sat_dev_adhoc
    if (any(is.na(c(r$chi2_com, r$chi2_MI, r$chi2_MI_corr, chi2_adhoc, r$df)))) next
    k <- k + 1
    rows[[k]] <- data.frame(
      chi2_com     = r$chi2_com,
      chi2_adhoc   = chi2_adhoc,
      chi2_MI      = r$chi2_MI,
      chi2_MI_corr = r$chi2_MI_corr,
      chi2_D3      = r$chi2_D3,
      df           = r$df
    )
  }
  if (length(rows) == 0) NULL else do.call(rbind, rows)
}

ALPHA <- 0.05
cat(sprintf("%-12s %4s %5s %5s | %7s %7s %7s %7s %7s\n",
            "imputer", "N", "mr", "n_rep",
            "com", "adhoc", "MI", "MI_corr", "D3"))
cat(strrep("-", 80), "\n")

for (config_block in list(
       list(label = "PMM",           res = res_pmm),
       list(label = "Amelia_empri0", res = res_amelia))) {
  for (Ni in N_VALUES) for (mr in MR_VALUES) {
    cell_lbl <- sprintf("n=%d_mr=%.2f", Ni, mr)
    df_M1 <- extract_M1(config_block$res[[cell_lbl]])
    if (is.null(df_M1)) next
    thr <- stats::qchisq(1 - ALPHA, df = df_M1$df[1])
    rates <- sapply(c("chi2_com","chi2_adhoc","chi2_MI","chi2_MI_corr","chi2_D3"),
                    function(stat) mean(df_M1[[stat]] > thr) * 100)
    cat(sprintf("%-12s %4d %5.2f %5d | %6.1f%% %6.1f%% %6.1f%% %6.1f%% %6.1f%%\n",
                config_block$label, Ni, mr, nrow(df_M1),
                rates["chi2_com"], rates["chi2_adhoc"],
                rates["chi2_MI"], rates["chi2_MI_corr"], rates["chi2_D3"]))
  }
}

cat("\nAcceptance: chi2_MI_corr within +/- 1pp of nominal 5%.\n")
cat(sprintf("\nResults saved under: %s\n", out_root))
