# ============================================================================
# Map empirical flexibility-gap Δ across all available results-* cells.
#
# Empirical Δ ≡ E[chi2_adhoc(M1)] - E[chi2_com(M1)]   (raw inflation)
#
# This is what the flexibility-gap derivation must reproduce (modulo the
# Bartlett/Wishart pieces, which are subtracted off below).
#
# For each results-* cell:
#   - n_reps        : usable replications
#   - df            : nominal df = npar(Msat) - npar(M1)         (lavaan: 45-32=13?)
#                     [empirically df=22 reported in briefing — check this]
#   - p_M1, p_Msat  : analysis-model and saturated npar (lavaan)
#   - mean(chi2_com)         : oracle mean
#   - mean(chi2_adhoc)       : ad-hoc pooled (per-imputation own MLE)
#   - mean(chi2_MI)          : -2[Q_bar(M1)-Q_bar(Msat)] + tr_RIV_M1 - tr_RIV_Msat
#   - mean(tr_RIV_M1), mean(tr_RIV_Msat)
#   - bartlett_emp     : df / mean(chi2_com)         (oracle Bartlett scaling)
#   - delta_raw        : mean(chi2_adhoc) - df       (total inflation over df)
#   - delta_overOracle : mean(chi2_adhoc) - mean(chi2_com)  (over-and-above oracle)
#   - delta_W          : amount left over after Wishart corr applied to MI gap
# ============================================================================

ALPHA <- 0.05

extract_chi2 <- function(rds_path, label) {
  if (!file.exists(rds_path)) return(NULL)
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
      cell_label  = label,
      n           = rr$n,
      mr          = if (!is.null(rr$miss_rate)) rr$miss_rate else NA_real_,
      df          = r$df,
      chi2_com    = r$chi2_com,
      chi2_adhoc  = chi2_adhoc,
      chi2_MI     = r$chi2_MI,
      tr_RIV_M1   = d$tr_RIV,
      tr_RIV_Msat = s$tr_RIV,
      p_M1        = d$npar,
      p_Msat      = s$npar
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

cells <- list(
  list(label = "amelia_e0.1_n100_mr0.20",  path = "hpc/results-amelia-empri0.1/results_n=100_mr=0.20.rds"),
  list(label = "amelia_e0.1_n100_mr0.40",  path = "hpc/results-amelia-empri0.1/results_n=100_mr=0.40.rds"),
  list(label = "amelia_e0.1_n500_mr0.20",  path = "hpc/results-amelia-empri0.1/results_n=500_mr=0.20.rds"),
  list(label = "amelia_e0.1_n500_mr0.40",  path = "hpc/results-amelia-empri0.1/results_n=500_mr=0.40.rds"),
  list(label = "amelia_e0.1_n1000_mr0.20", path = "hpc/results-amelia-empri0.1/results_n=1000_mr=0.20.rds"),
  list(label = "amelia_e0.1_n1000_mr0.40", path = "hpc/results-amelia-empri0.1/results_n=1000_mr=0.40.rds"),
  list(label = "amelia_e100_n100_mr0.20",  path = "hpc/results-amelia-empri100/results_n=100_mr=0.20.rds"),
  list(label = "amelia_e100_n100_mr0.40",  path = "hpc/results-amelia-empri100/results_n=100_mr=0.40.rds"),
  list(label = "amelia_e100_n500_mr0.20",  path = "hpc/results-amelia-empri100/results_n=500_mr=0.20.rds"),
  list(label = "amelia_e100_n500_mr0.40",  path = "hpc/results-amelia-empri100/results_n=500_mr=0.40.rds"),
  list(label = "amelia_e100_n1000_mr0.20", path = "hpc/results-amelia-empri100/results_n=1000_mr=0.20.rds"),
  list(label = "amelia_e100_n1000_mr0.40", path = "hpc/results-amelia-empri100/results_n=1000_mr=0.40.rds"),
  list(label = "PMM_n100_mr0.40",          path = "hpc/results-decomp/results_n=100_mr=0.40.rds"),
  list(label = "PMM_n250_mr0.40",          path = "hpc/results-decomp/results_n=250_mr=0.40.rds"),
  list(label = "PMM_n1000_mr0.40",         path = "hpc/results-decomp/results_n=1000_mr=0.40.rds"),
  list(label = "PMM_M100_n100_mr0.10",     path = "hpc/results-full-M100/results_n=100_mr=0.10.rds"),
  list(label = "PMM_M100_n100_mr0.25",     path = "hpc/results-full-M100/results_n=100_mr=0.25.rds"),
  list(label = "PMM_M100_n100_mr0.40",     path = "hpc/results-full-M100/results_n=100_mr=0.40.rds"),
  list(label = "PMM_M100_n250_mr0.10",     path = "hpc/results-full-M100/results_n=250_mr=0.10.rds"),
  list(label = "PMM_M100_n250_mr0.25",     path = "hpc/results-full-M100/results_n=250_mr=0.25.rds"),
  list(label = "PMM_M100_n250_mr0.40",     path = "hpc/results-full-M100/results_n=250_mr=0.40.rds"),
  list(label = "PMM_M100_n500_mr0.10",     path = "hpc/results-full-M100/results_n=500_mr=0.10.rds"),
  list(label = "PMM_M100_n500_mr0.25",     path = "hpc/results-full-M100/results_n=500_mr=0.25.rds"),
  list(label = "PMM_M100_n500_mr0.40",     path = "hpc/results-full-M100/results_n=500_mr=0.40.rds"),
  list(label = "PMM_M100_n1000_mr0.10",    path = "hpc/results-full-M100/results_n=1000_mr=0.10.rds"),
  list(label = "PMM_M100_n1000_mr0.25",    path = "hpc/results-full-M100/results_n=1000_mr=0.25.rds"),
  list(label = "PMM_M100_n1000_mr0.40",    path = "hpc/results-full-M100/results_n=1000_mr=0.40.rds"),
  list(label = "mvn_M1_n250_mr0.40",       path = "hpc/results-mvn_M1-N250/results_n=250_mr=0.40.rds"),
  list(label = "mvn_M1_n1000_mr0.40",      path = "hpc/results-mvn_M1-N1000/results_combined.rds")
)

cat(sprintf("%-26s %5s %5s %4s %5s %6s %8s %9s %8s %10s %10s %9s %9s\n",
            "cell", "n_rep", "n", "df", "p_M1", "p_Msat",
            "chi_com", "chi_adhoc", "chi_MI",
            "tr_M1", "tr_Msat", "delta_oR", "delta_W"))

results_summary <- list(); ki <- 0
for (cell in cells) {
  df <- extract_chi2(cell$path, cell$label)
  if (is.null(df)) next
  ki <- ki + 1
  n_val <- df$n[1]; mr_val <- df$mr[1]; df_val <- df$df[1]
  p_M1   <- df$p_M1[1]; p_Msat <- df$p_Msat[1]
  e_com   <- mean(df$chi2_com)
  e_adhoc <- mean(df$chi2_adhoc)
  e_MI    <- mean(df$chi2_MI)
  e_trM1  <- mean(df$tr_RIV_M1)
  e_trMsat<- mean(df$tr_RIV_Msat)
  fac_M1   <- (n_val - p_M1   - 1) / n_val
  fac_Msat <- (n_val - p_Msat - 1) / n_val
  e_chi2_MI_W <- e_adhoc + e_trM1 * fac_M1 - e_trMsat * fac_Msat
  delta_oR    <- e_adhoc - e_com   # raw over oracle
  delta_W     <- e_chi2_MI_W - e_com  # leftover after Wishart applied
  cat(sprintf("%-26s %5d %5d %4d %5d %6d %8.2f %9.2f %8.2f %10.3f %10.3f %9.2f %9.2f\n",
              cell$label, nrow(df), n_val, df_val, p_M1, p_Msat,
              e_com, e_adhoc, e_MI, e_trM1, e_trMsat, delta_oR, delta_W))
  results_summary[[ki]] <- data.frame(
    cell = cell$label, n = n_val, mr = mr_val, df = df_val,
    n_rep = nrow(df), p_M1 = p_M1, p_Msat = p_Msat,
    e_com = e_com, e_adhoc = e_adhoc, e_MI = e_MI,
    e_trM1 = e_trM1, e_trMsat = e_trMsat,
    e_chi2_MI_W = e_chi2_MI_W,
    delta_oR = delta_oR, delta_W = delta_W
  )
}

summary_df <- do.call(rbind, results_summary)
saveRDS(summary_df, "hpc/empirical_delta_summary.rds")
cat("\nSaved hpc/empirical_delta_summary.rds (", nrow(summary_df), "cells)\n")
