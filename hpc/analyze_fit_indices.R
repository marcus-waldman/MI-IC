# ============================================================================
# Analyze fit_indices_df from miicsem 0.5.3+ result files.
#
# For each (cell, candidate model), reports mean and SD of the three
# variants (com / adhoc / MI) of CFI, TLI, and RMSEA. Also reports the
# bias of the adhoc variant relative to oracle and the truncation-corrected
# analytic prediction from v4.3 §10.
#
# Default cells: hpc/results-fitindices-N500/. Pass --cells <dir1,dir2> to
# add others.
# ============================================================================

ALPHA <- 0.05

extract_fit_indices_M1 <- function(rds_path, label) {
  if (!file.exists(rds_path)) return(NULL)
  rs <- readRDS(rds_path)
  is_nested <- !is.null(rs[[1]]) && is.list(rs[[1]]) &&
               !("rep_id" %in% names(rs[[1]]))
  reps <- if (is_nested) unlist(rs, recursive = FALSE) else rs
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    fi <- rr$fit_indices_df
    chi <- rr$chi2_df
    if (is.null(fi) || !"M1" %in% rownames(fi)) next
    if (is.null(chi) || !"M1" %in% rownames(chi)) next
    if (!"Mnull" %in% rownames(chi)) next
    k <- k + 1
    rows[[k]] <- data.frame(
      cell_label = label, n = rr$n, mr = rr$miss_rate,
      cfi_com    = fi["M1", "cfi_com"],
      cfi_adhoc  = fi["M1", "cfi_adhoc"],
      cfi_MI     = fi["M1", "cfi_MI"],
      tli_com    = fi["M1", "tli_com"],
      tli_adhoc  = fi["M1", "tli_adhoc"],
      tli_MI     = fi["M1", "tli_MI"],
      rmsea_com  = fi["M1", "rmsea_com"],
      rmsea_adhoc = fi["M1", "rmsea_adhoc"],
      rmsea_MI   = fi["M1", "rmsea_MI"],
      chi2_com_M1   = chi["M1", "chi2_com"],
      chi2_adhoc_M1 = chi["M1", "chi2_adhoc"],
      chi2_MI_M1    = chi["M1", "chi2_MI"],
      chi2_com_M0   = chi["Mnull", "chi2_com"],
      df_M1         = chi["M1", "df"],
      df_M0         = chi["Mnull", "df"]
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# Default: results-fitindices-N500 cells
cells <- list(
  list(label = "PMM_e0_N500_mr0.10", path = "hpc/results-fitindices-N500/results_n=500_mr=0.10.rds"),
  list(label = "PMM_e0_N500_mr0.25", path = "hpc/results-fitindices-N500/results_n=500_mr=0.25.rds"),
  list(label = "PMM_e0_N500_mr0.40", path = "hpc/results-fitindices-N500/results_n=500_mr=0.40.rds")
)

cat("=== Fit indices summary (M1, true model) ===\n")
cat(sprintf("%-22s %5s | %8s %10s %8s | %8s %10s %8s | %10s %12s %10s\n",
            "cell", "n_rep",
            "cfi_com", "cfi_adhoc", "cfi_MI",
            "tli_com", "tli_adhoc", "tli_MI",
            "rmsea_com", "rmsea_adhoc", "rmsea_MI"))

results <- list(); ki <- 0
for (cell in cells) {
  d <- extract_fit_indices_M1(cell$path, cell$label)
  if (is.null(d)) {
    cat(sprintf("%-22s  (no data)\n", cell$label)); next
  }
  ki <- ki + 1
  cat(sprintf("%-22s %5d | %8.4f %10.4f %8.4f | %8.4f %10.4f %8.4f | %10.4f %12.4f %10.4f\n",
              cell$label, nrow(d),
              mean(d$cfi_com),    mean(d$cfi_adhoc),    mean(d$cfi_MI),
              mean(d$tli_com),    mean(d$tli_adhoc),    mean(d$tli_MI),
              mean(d$rmsea_com),  mean(d$rmsea_adhoc),  mean(d$rmsea_MI)))
  results[[ki]] <- d
}

cat("\n=== Bias relative to oracle ===\n")
cat(sprintf("%-22s | %10s %10s | %10s %10s | %10s %10s\n",
            "cell",
            "cfi_adhoc", "cfi_MI",
            "tli_adhoc", "tli_MI",
            "rmsea_adhoc", "rmsea_MI"))
for (d in results) {
  if (is.null(d) || nrow(d) == 0) next
  cat(sprintf("%-22s | %10.5f %10.5f | %10.5f %10.5f | %10.5f %10.5f\n",
              d$cell_label[1],
              mean(d$cfi_com)   - mean(d$cfi_adhoc),
              mean(d$cfi_com)   - mean(d$cfi_MI),
              mean(d$tli_com)   - mean(d$tli_adhoc),
              mean(d$tli_com)   - mean(d$tli_MI),
              mean(d$rmsea_com) - mean(d$rmsea_adhoc),
              mean(d$rmsea_com) - mean(d$rmsea_MI)))
}

cat("\n=== Truncation-corrected v4.3 §10 prediction for CFI_adhoc bias ===\n")
cat("bias ≈ {E[max(a + Δ, 0)] − E[max(a, 0)]} / [chi2_com(M0) − df_0]\n")
cat("       where a ≡ chi2_com(M1) − df_M1, Δ = chi2_adhoc(M1) − chi2_com(M1)\n\n")
cat(sprintf("%-22s | %12s %12s %12s\n",
            "cell", "bias_emp", "bias_pred", "ratio"))
for (d in results) {
  if (is.null(d) || nrow(d) == 0) next
  bias_emp <- mean(d$cfi_com) - mean(d$cfi_adhoc)
  delta_M1 <- mean(d$chi2_adhoc_M1) - mean(d$chi2_com_M1)
  null_disc <- mean(d$chi2_com_M0) - mean(d$df_M0)
  # truncation-corrected prediction:
  #   E[max(a + Δ, 0)] − E[max(a, 0)] for a ~ centered with var ≈ 2·df_M1
  sigma <- sqrt(2 * d$df_M1[1])
  trunc_diff <- delta_M1 * pnorm(delta_M1 / sigma) +
                sigma * (dnorm(delta_M1 / sigma) - dnorm(0))
  bias_pred <- trunc_diff / null_disc
  ratio <- if (abs(bias_pred) > 1e-8) bias_emp / bias_pred else NA
  cat(sprintf("%-22s | %12.5f %12.5f %12.3f\n",
              d$cell_label[1], bias_emp, bias_pred, ratio))
}
cat("\nDone.\n")
