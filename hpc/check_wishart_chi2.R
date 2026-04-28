# ============================================================================
# Apply Wishart correction post-hoc to chi2_MI and check Type I error
# ============================================================================
# Bias factor for inverse-Hessian-based covariance:
#   factor(n, p) = n / (n - p - 1)
# So a bias-CORRECTED tr(RIV) for a model with p_M params is:
#   tr_RIV_corr_M = tr_RIV_raw_M * (n - p_M - 1) / n
#
# For chi-square: chi2_MI uses tr_RIV(M1) - tr_RIV(Msat).  Since M1 and
# Msat have DIFFERENT p (32 vs 45), the Wishart factors differ, so
# correction changes the difference (not just a uniform shift).
# ============================================================================

ALPHA <- 0.05

extract_chi2 <- function(rds_path, label) {
  rs <- readRDS(rds_path)
  is_nested <- !is.null(rs[[1]]) && is.list(rs[[1]]) &&
               !("rep_id" %in% names(rs[[1]]))
  reps <- if (is_nested) unlist(rs, recursive = FALSE) else rs
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    chi <- rr$chi2_df
    dev <- rr$dev_df
    if (is.null(chi) || is.null(dev) || !"M1" %in% rownames(chi)) next
    if (!"Msat" %in% rownames(dev)) next
    sat_dev_adhoc <- dev["Msat", "DEV_adhoc"]
    if (is.na(sat_dev_adhoc)) next
    r <- chi["M1", ]; d <- dev["M1", ]; s <- dev["Msat", ]
    chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
    if (any(is.na(c(d$tr_RIV, s$tr_RIV, d$npar, s$npar, rr$n)))) next
    k <- k + 1
    rows[[k]] <- data.frame(
      method_label   = label,
      n              = rr$n,
      df             = r$df,
      chi2_com       = r$chi2_com,
      chi2_adhoc     = chi2_adhoc,
      chi2_MI        = r$chi2_MI,
      tr_RIV_M1_raw  = d$tr_RIV,
      tr_RIV_Msat_raw = s$tr_RIV,
      p_M1           = d$npar,
      p_Msat         = s$npar
    )
  }
  do.call(rbind, rows)
}

cases <- list(
  amelia = "hpc/results-congeniality-amelia/results_combined.rds",
  PMM    = "hpc/results-decomp/results_n=250_mr=0.40.rds",
  mvn_M1 = "hpc/results-mvn_M1-N250/results_combined.rds"
)

cat("Wishart correction post-hoc: chi-square Type I error at alpha=0.05\n")
cat("N=250, mr=0.40, M=50, true model M1 (df=22)\n\n")
cat(sprintf("%-12s  %6s  %8s  %10s  %8s  %12s\n",
            "imputer", "n_reps", "chi2_com", "chi2_adhoc",
            "chi2_MI", "chi2_MI_wish"))

for (lab in names(cases)) {
  df <- extract_chi2(cases[[lab]], lab)
  thresh <- stats::qchisq(1 - ALPHA, df = df$df[1])

  # Wishart-corrected per-rep tr(RIV)
  fac_M1   <- (df$n - df$p_M1   - 1) / df$n   # all should be 0.868 at N=250, p=32
  fac_Msat <- (df$n - df$p_Msat - 1) / df$n   # all should be 0.816 at N=250, p=45
  tr_M1_w   <- df$tr_RIV_M1_raw   * fac_M1
  tr_Msat_w <- df$tr_RIV_Msat_raw * fac_Msat

  # Reconstruct chi2_MI as: chi2_adhoc + (tr_RIV_M1 - tr_RIV_Msat)
  # Apply Wishart-corrected version analogously
  chi2_MI_wish <- df$chi2_adhoc + (tr_M1_w - tr_Msat_w)

  cat(sprintf("%-12s  %6d  %7.1f%%  %9.1f%%  %7.1f%%  %11.1f%%\n",
              lab, nrow(df),
              mean(df$chi2_com   > thresh) * 100,
              mean(df$chi2_adhoc > thresh) * 100,
              mean(df$chi2_MI    > thresh) * 100,
              mean(chi2_MI_wish  > thresh) * 100))
}

cat("\nMean chi2 (target df=22)\n\n")
cat(sprintf("%-12s  %8s  %10s  %8s  %12s\n",
            "imputer", "chi2_com", "chi2_adhoc", "chi2_MI", "chi2_MI_wish"))
for (lab in names(cases)) {
  df <- extract_chi2(cases[[lab]], lab)
  fac_M1   <- (df$n - df$p_M1   - 1) / df$n
  fac_Msat <- (df$n - df$p_Msat - 1) / df$n
  chi2_MI_wish <- df$chi2_adhoc + df$tr_RIV_M1_raw * fac_M1 -
                                  df$tr_RIV_Msat_raw * fac_Msat
  cat(sprintf("%-12s  %8.2f  %10.2f  %8.2f  %12.2f\n",
              lab,
              mean(df$chi2_com), mean(df$chi2_adhoc),
              mean(df$chi2_MI),  mean(chi2_MI_wish)))
}

cat("\nWishart factor diagnostics at N=250:\n")
cat(sprintf("  M1 (p=32)  : (n-p-1)/n = %.3f\n",   (250-32-1)/250))
cat(sprintf("  Msat (p=45): (n-p-1)/n = %.3f\n",   (250-45-1)/250))
cat(sprintf("  Difference factor:    %.3f\n",
            (250-45-1)/250 - (250-32-1)/250))
