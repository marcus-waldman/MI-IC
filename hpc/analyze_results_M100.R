# ============================================================================
# Study 2 SEM Results Analysis — M = 100
# ============================================================================
# Loads results from results-full-M100/ and produces:
#   1. Selection accuracy table  (all IC methods x 12 conditions)
#   2. Full selection frequency tables (per method)
#   3. Chi-square sanity check   (mean chi2 vs df = 22 for M1)
#   4. tr(RIV) summary           (MI vs FIML target, by condition & model)
#   5. Accuracy plots            (PDF to figures/)
#   6. Side-by-side M20 vs M100  (selection accuracy comparison)
#
# Run from repo root:
#   Rscript hpc/analyze_results_M100.R
# ============================================================================

library(miicsem)

RESULTS_M100 <- "hpc/results-full-M100"
RESULTS_M20  <- "hpc/results-full"         # M=20 baseline (for comparison)
FIGURES_DIR  <- "hpc/figures-M100"

if (!dir.exists(FIGURES_DIR)) dir.create(FIGURES_DIR, recursive = TRUE)

cat("==========================================================\n")
cat("  Study 2 SEM — M=100 Results Analysis\n")
cat(sprintf("  miicsem version: %s\n", as.character(packageVersion("miicsem"))))
cat(sprintf("  Results dir   : %s\n", RESULTS_M100))
cat("==========================================================\n\n")

# ============================================================
# 1. Load data
# ============================================================
cat("Loading M=100 results...\n")
res_df   <- load_results(RESULTS_M100)
dev_df   <- load_deviances(RESULTS_M100)
chi2_df  <- load_chi_squares(RESULTS_M100)

cat(sprintf("  Reps loaded (selection df) : %d\n", nrow(res_df)))
cat(sprintf("  Rows in deviance df        : %d\n", nrow(dev_df)))
cat(sprintf("  Rows in chi-square df      : %d\n", nrow(chi2_df)))

# ============================================================
# 2. Condition-level rep counts / failure rate
# ============================================================
cat("\n--- Rep counts & failure rate by condition ---\n")
cond_counts <- aggregate(rep_id ~ n + miss_rate, data = res_df, FUN = length)
colnames(cond_counts)[3] <- "n_succeeded"
cond_counts$pct_ok <- round(cond_counts$n_succeeded / 2000 * 100, 1)
print(cond_counts, row.names = FALSE)

# ============================================================
# 3. Selection accuracy (all methods)
# ============================================================
cat("\n--- Selection accuracy: % selecting true model M1 ---\n")
acc <- selection_accuracy_table(res_df)
print(acc, digits = 1, row.names = FALSE)

# Save as CSV
write.csv(acc, file.path(FIGURES_DIR, "selection_accuracy_M100.csv"), row.names = FALSE)
cat(sprintf("  [saved] %s/selection_accuracy_M100.csv\n", FIGURES_DIR))

# ============================================================
# 4. Full selection frequency tables (MI_AIC and MI_BIC)
# ============================================================
cat("\n--- Selection frequencies: MI_AIC ---\n")
freq_mi_aic <- selection_frequency_table(res_df, ic_method = "MI_AIC")
print(freq_mi_aic, row.names = FALSE)

cat("\n--- Selection frequencies: MI_BIC ---\n")
freq_mi_bic <- selection_frequency_table(res_df, ic_method = "MI_BIC")
print(freq_mi_bic, row.names = FALSE)

cat("\n--- Selection frequencies: AICcd ---\n")
freq_aiccd <- selection_frequency_table(res_df, ic_method = "AICcd")
print(freq_aiccd, row.names = FALSE)

# ============================================================
# 5. Chi-square sanity check (mean chi2 vs df = 22 for M1)
# ============================================================
cat("\n--- Chi-square means for M1 (true model, df = 22) ---\n")
chi2_M1 <- chi2_df[chi2_df$model == "M1", ]
chi2_M1_agg <- aggregate(
  cbind(chi2_com, chi2_MI, chi2_D3) ~ n + miss_rate,
  data = chi2_M1,
  FUN  = function(x) round(mean(x, na.rm = TRUE), 2)
)
chi2_M1_agg <- chi2_M1_agg[order(chi2_M1_agg$n, chi2_M1_agg$miss_rate), ]
print(chi2_M1_agg, row.names = FALSE)

# Full chi2 summary (all models)
chi2_summ <- chi2_summary_table(chi2_df)
write.csv(chi2_summ, file.path(FIGURES_DIR, "chi2_summary_M100.csv"), row.names = FALSE)
cat(sprintf("  [saved] %s/chi2_summary_M100.csv\n", FIGURES_DIR))

# ============================================================
# 6. tr(RIV) summary (MI-based estimate)
# ============================================================
cat("\n--- tr(RIV) summary for M1 (from selection results) ---\n")
triv_summ <- tr_riv_summary(res_df)
print(triv_summ, digits = 3, row.names = FALSE)

# ============================================================
# 7. tr(RIV) MI vs FIML comparison (deviance-level, all models)
# ============================================================
cat("\n--- tr(RIV) MI vs FIML comparison (M1 only, across conditions) ---\n")
triv_comp <- tr_riv_comparison_table(dev_df)
triv_M1   <- triv_comp[triv_comp$model == "M1", ]
triv_M1$bias_pct <- round(triv_M1$bias / triv_M1$mean_tr_fiml * 100, 2)
print(triv_M1[, c("n", "miss_rate", "n_reps",
                   "mean_tr_mi", "mean_tr_fiml", "bias",
                   "bias_pct", "rmse", "cor")],
      digits = 3, row.names = FALSE)

# Save full table
write.csv(triv_comp, file.path(FIGURES_DIR, "trRIV_MI_vs_FIML_M100.csv"), row.names = FALSE)
cat(sprintf("  [saved] %s/trRIV_MI_vs_FIML_M100.csv\n", FIGURES_DIR))

# ============================================================
# 8. Accuracy plots
# ============================================================
cat("\n--- Generating accuracy plots ---\n")

# AIC-class: AIC_com, AIC_adhoc, AICcd, MI_AIC
plot_accuracy_by_n(
  acc_table   = acc,
  output_file = file.path(FIGURES_DIR, "accuracy_AIC_M100.pdf")
)
cat(sprintf("  [saved] %s/accuracy_AIC_M100.pdf\n", FIGURES_DIR))

# BIC-class: BIC_com, BIC_adhoc, MI_BIC
pdf(file.path(FIGURES_DIR, "accuracy_BIC_M100.pdf"), width = 10, height = 6)
ic_bic    <- c("BIC_com", "BIC_adhoc", "MI_BIC")
colors_b  <- c("black", "blue", "darkgreen")
ltys_b    <- c(1, 2, 1)
miss_rates <- sort(unique(acc$miss_rate))
graphics::par(mfrow = c(1, length(miss_rates)), mar = c(4, 4, 3, 1))
for (mr in miss_rates) {
  sub <- acc[acc$miss_rate == mr, ]
  graphics::plot(NULL, xlim = range(sub$N), ylim = c(0, 100),
       xlab = "N", ylab = "% Selecting M1",
       main = sprintf("Miss rate = %.0f%%", mr * 100), log = "x")
  for (j in seq_along(ic_bic)) {
    m <- ic_bic[j]
    if (m %in% colnames(sub)) {
      graphics::lines(sub$N, sub[[m]], col = colors_b[j], lty = ltys_b[j], lwd = 2)
      graphics::points(sub$N, sub[[m]], col = colors_b[j], pch = 16)
    }
  }
  if (mr == miss_rates[1]) {
    graphics::legend("bottomright", legend = ic_bic, col = colors_b,
           lty = ltys_b, lwd = 2, cex = 0.8)
  }
}
grDevices::dev.off()
cat(sprintf("  [saved] %s/accuracy_BIC_M100.pdf\n", FIGURES_DIR))

# ============================================================
# 9. M=20 vs M=100 comparison (if M=20 results available)
# ============================================================
if (dir.exists(RESULTS_M20) &&
    file.exists(file.path(RESULTS_M20, "results_combined.rds"))) {

  cat("\n--- M=20 vs M=100 selection accuracy comparison ---\n")

  res20 <- load_results(RESULTS_M20)
  acc20 <- selection_accuracy_table(res20)
  acc20$M <- 20L
  acc$M   <- 100L

  shared_cols <- intersect(colnames(acc20), colnames(acc))
  shared_ic   <- setdiff(shared_cols, c("N", "miss_rate", "n_reps", "M"))

  comp_rows <- list()
  for (i in seq_len(nrow(acc20))) {
    row20 <- acc20[i, ]
    row100_idx <- which(acc$N == row20$N & acc$miss_rate == row20$miss_rate)
    if (length(row100_idx) == 0) next
    row100 <- acc[row100_idx, ]
    for (m in shared_ic) {
      comp_rows[[length(comp_rows) + 1]] <- data.frame(
        N         = row20$N,
        miss_rate = row20$miss_rate,
        method    = m,
        acc_M20   = row20[[m]],
        acc_M100  = row100[[m]],
        delta     = row100[[m]] - row20[[m]],
        stringsAsFactors = FALSE
      )
    }
  }
  comp_df <- do.call(rbind, comp_rows)
  comp_df <- comp_df[order(comp_df$method, comp_df$N, comp_df$miss_rate), ]

  cat("  Δ accuracy (M100 − M20) for MI_AIC and MI_BIC:\n")
  mi_comp <- comp_df[comp_df$method %in% c("MI_AIC", "MI_BIC"), ]
  print(mi_comp, digits = 1, row.names = FALSE)

  write.csv(comp_df, file.path(FIGURES_DIR, "M20_vs_M100_accuracy.csv"), row.names = FALSE)
  cat(sprintf("  [saved] %s/M20_vs_M100_accuracy.csv\n", FIGURES_DIR))
} else {
  cat("\n  (M=20 results not found — skipping M=20 vs M=100 comparison)\n")
}

# ============================================================
# 10. tr(RIV) precision comparison: M=20 vs M=100
# ============================================================
# Handled by the standalone script hpc/trRIV_precision_compare.R
# (M=20 results lack tr_RIV_fiml; comparison uses SD across reps.)
cat("\n  [Note] tr(RIV) precision comparison (M=20 vs M=100):\n")
cat("         run hpc/trRIV_precision_compare.R separately.\n")
cat("         Outputs land in hpc/figures-M100/trRIV_*.csv\n")

cat("\n==========================================================\n")
cat("  Analysis complete.\n")
cat(sprintf("  All outputs in: %s/\n", FIGURES_DIR))
cat("==========================================================\n")
