# ============================================================================
# Local verification of miicsem 0.5.3's new fit-index pipeline.
#
# Runs PMM, N=500, mr ∈ {0.10, 0.25, 0.40}, M=50, 200 reps, 16 cores.
# Same seed scheme as hpc/verify_mi_cfi.R so the per-rep CFI_MI values
# computed here should match within Monte Carlo noise.
#
# Validation checks:
#   1. fit_indices_df is populated for every successful rep, with rows for
#      M1..M12 and columns cfi_*, tli_*, rmsea_*.
#   2. mean(cfi_MI) ≈ mean(cfi_MI from hpc/verify_mi_cfi.R) within 0.001.
#   3. mean(tli_MI) ≈ mean(cfi_MI) (similar magnitude — both incremental).
#   4. mean(rmsea_MI) is sensible (≈ 0 for correctly specified M1).
# ============================================================================

suppressPackageStartupMessages({
  library(miicsem); library(lavaan); library(mice)
})
stopifnot(as.character(packageVersion("miicsem")) >= "0.5.3")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

N            <- 500L
MR_VALUES    <- c(0.10, 0.25, 0.40)
M_IMP        <- 50L
N_REPS       <- 200L
N_CORES      <- 16L
BASE_SEED    <- 32897891L

results_dir <- "hpc/results-fitindices-N500"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps       = N_REPS,
  n_cores      = N_CORES,
  seed         = BASE_SEED,
  sample_sizes = N,
  miss_rates   = MR_VALUES,
  M            = M_IMP,
  mice_method  = "pmm",
  results_dir  = results_dir,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nWall time: %.1fs (%.2fh)\n", elapsed, elapsed / 3600))

# --- Quick validation summary ----------------------------------------------
cat("\n=== fit_indices summary (M1, true model) ===\n")
cat(sprintf("%-6s %5s | %8s %10s %8s | %8s %10s %8s | %10s %12s %10s\n",
            "mr", "n_rep",
            "cfi_com", "cfi_adhoc", "cfi_MI",
            "tli_com", "tli_adhoc", "tli_MI",
            "rmsea_com", "rmsea_adhoc", "rmsea_MI"))
for (mr in MR_VALUES) {
  cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
  reps <- res[[cell_label]]
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    fi <- rr$fit_indices_df
    if (is.null(fi) || !"M1" %in% rownames(fi)) next
    k <- k + 1
    rows[[k]] <- fi["M1", ]
  }
  if (length(rows) == 0) { cat(sprintf("%-6.2f  no usable reps\n", mr)); next }
  d <- do.call(rbind, rows)
  cat(sprintf("%-6.2f %5d | %8.4f %10.4f %8.4f | %8.4f %10.4f %8.4f | %10.4f %12.4f %10.4f\n",
              mr, nrow(d),
              mean(d$cfi_com),    mean(d$cfi_adhoc),    mean(d$cfi_MI),
              mean(d$tli_com),    mean(d$tli_adhoc),    mean(d$tli_MI),
              mean(d$rmsea_com),  mean(d$rmsea_adhoc),  mean(d$rmsea_MI)))
}

# --- Cross-check: cfi_MI here should match hpc/verify_mi_cfi.R ------------
cat("\n=== Cross-check vs hpc/verify_mi_cfi.R (same N, mr, seed scheme) ===\n")
v0 <- tryCatch(readRDS("hpc/results-mi-cfi/cfi_summary.rds"),
               error = function(e) NULL)
if (!is.null(v0)) {
  cat(sprintf("%-6s | %15s %15s %12s\n",
              "mr", "cfi_MI (0.5.3)", "cfi_MI (verify)", "diff"))
  for (mr in MR_VALUES) {
    cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
    reps <- res[[cell_label]]
    cfi_m1_new <- vapply(reps, function(r) {
      if (is.null(r) || isTRUE(r$failed) || is.null(r$fit_indices_df) ||
          !"M1" %in% rownames(r$fit_indices_df)) return(NA_real_)
      r$fit_indices_df["M1", "cfi_MI"]
    }, numeric(1))
    cfi_m1_new <- cfi_m1_new[!is.na(cfi_m1_new)]
    cfi_m1_v0  <- v0[[cell_label]]$cfi_MI
    diff <- mean(cfi_m1_new) - mean(cfi_m1_v0)
    cat(sprintf("%-6.2f | %15.5f %15.5f %12.5f\n",
                mr, mean(cfi_m1_new), mean(cfi_m1_v0), diff))
  }
} else {
  cat("  (verify_mi_cfi summary not found; skipping cross-check)\n")
}

cat("\nDone.\n")
