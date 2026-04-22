# ============================================================================
# tr(RIV) Monte Carlo Precision: M=20 vs M=100
# ============================================================================
# Compares the variability of the MI-based tr(RIV) estimate across
# replications. SD across reps is the direct measure of Monte Carlo
# sampling noise — a smaller SD at M=100 vs M=20 confirms that increasing
# M stabilises the RIV correction.
#
# Note: tr_RIV_fiml is only available in the M=100 results (the FIML
# reference feature was added after the M=20 run). The M=100 FIML target
# is therefore used as the reference for bias; for M=20 we report the
# MI-based distribution only.
# ============================================================================
library(miicsem)

RESULTS_M100 <- "hpc/results-full-M100"
RESULTS_M20  <- "hpc/results-full"
FIGURES_DIR  <- "hpc/figures-M100"

# ---------- load selection-level data (small, has tr_RIV_M1) ---------------
cat("Loading selection results...\n")
res100 <- load_results(RESULTS_M100)
res20  <- load_results(RESULTS_M20)
cat(sprintf("  M=100: %d reps  |  M=20: %d reps\n", nrow(res100), nrow(res20)))

# ---------- tr(RIV) summary per condition for M=100 and M=20 ---------------
summarise_triv <- function(res_df, M_label) {
  conditions <- unique(res_df[, c("n", "miss_rate")])
  conditions <- conditions[order(conditions$n, conditions$miss_rate), ]
  rows <- lapply(seq_len(nrow(conditions)), function(i) {
    cond <- conditions[i, ]
    vals <- res_df$tr_RIV_M1[res_df$n == cond$n &
                              res_df$miss_rate == cond$miss_rate]
    vals <- vals[!is.na(vals)]
    data.frame(
      M         = M_label,
      N         = cond$n,
      miss_rate = cond$miss_rate,
      n_reps    = length(vals),
      mean      = mean(vals),
      sd        = stats::sd(vals),
      cv        = stats::sd(vals) / mean(vals),
      q05       = stats::quantile(vals, 0.05),
      q95       = stats::quantile(vals, 0.95),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cat("\n--- tr(RIV) for M1: M=20 distribution ---\n")
s20  <- summarise_triv(res20,  M_label = 20L)
print(s20, digits = 3, row.names = FALSE)

cat("\n--- tr(RIV) for M1: M=100 distribution ---\n")
s100 <- summarise_triv(res100, M_label = 100L)
print(s100, digits = 3, row.names = FALSE)

# ---------- side-by-side SD comparison -------------------------------------
cat("\n--- SD reduction: M=100 vs M=20 (lower SD = more precise estimate) ---\n")
comp <- merge(
  s20[ , c("N", "miss_rate", "mean", "sd", "cv")],
  s100[, c("N", "miss_rate", "mean", "sd", "cv")],
  by = c("N", "miss_rate"),
  suffixes = c("_M20", "_M100")
)
comp$sd_ratio   <- round(comp$sd_M100  / comp$sd_M20,  3)  # should be ~sqrt(20/100) ≈ 0.447
comp$mean_delta <- round(comp$mean_M100 - comp$mean_M20, 4)
comp <- comp[order(comp$N, comp$miss_rate), ]
print(comp, digits = 3, row.names = FALSE)

# ---------- theoretical prediction: SD ratio = sqrt(20/100) ≈ 0.447 --------
expected_ratio <- sqrt(20 / 100)
cat(sprintf("\nExpected SD ratio (sqrt(20/100)): %.3f\n", expected_ratio))
cat(sprintf("Observed SD ratio: mean = %.3f, range [%.3f, %.3f]\n",
            mean(comp$sd_ratio), min(comp$sd_ratio), max(comp$sd_ratio)))

# ---------- save -----------------------------------------------------------
all_summ <- rbind(s20, s100)
write.csv(all_summ, file.path(FIGURES_DIR, "trRIV_by_M_summary.csv"),
          row.names = FALSE)
write.csv(comp, file.path(FIGURES_DIR, "trRIV_precision_M20_vs_M100.csv"),
          row.names = FALSE)
cat(sprintf("\n[saved] %s/trRIV_by_M_summary.csv\n", FIGURES_DIR))
cat(sprintf("[saved] %s/trRIV_precision_M20_vs_M100.csv\n", FIGURES_DIR))

# ---------- tr(RIV) FIML comparison (M=100 only) ---------------------------
cat("\n--- tr(RIV) MI vs FIML for M=100 (M1 only) ---\n")
cat("Loading M=100 deviances (this may take a moment)...\n")
dev100 <- load_deviances(RESULTS_M100)
triv100_comp <- tr_riv_comparison_table(dev100)
rm(dev100); gc()

triv100_M1 <- triv100_comp[triv100_comp$model == "M1", ]
triv100_M1$bias_pct <- round(triv100_M1$bias / triv100_M1$mean_tr_fiml * 100, 1)

cat("\n(Note: large bias at N=100 is expected — FIML vcov is unstable at small N)\n")
print(
  triv100_M1[, c("n", "miss_rate", "n_reps",
                  "mean_tr_mi", "mean_tr_fiml",
                  "bias", "bias_pct", "rmse", "cor")],
  digits = 3, row.names = FALSE
)

write.csv(triv100_M1, file.path(FIGURES_DIR, "trRIV_MI_vs_FIML_M100_M1only.csv"),
          row.names = FALSE)
cat(sprintf("[saved] %s/trRIV_MI_vs_FIML_M100_M1only.csv\n", FIGURES_DIR))
cat("\nDone.\n")
