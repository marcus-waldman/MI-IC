# ============================================================================
# Chi-square test analysis: does MI correction beat ad-hoc for SEM fit tests?
# ============================================================================
# For each rep x candidate model M, we have four chi-squares:
#   chi2_com     = oracle (target value, what we'd get on complete data)
#   chi2_adhoc   = mean per-imputation chi^2 (derivable from DEV_adhoc)
#   chi2_MI      = MI bias-corrected: -2[Q_bar(M) - Q_bar(sat)]
#                                     + (tr_RIV_M - tr_RIV_sat)
#   chi2_D3      = Meng-Rubin D_3 statistic
#
# Diagnostics:
#   1. Bias and RMSE vs chi2_com, by candidate model
#   2. Type I error rate at alpha = 0.05 for the TRUE model M1
#      (should be ~5% if calibrated)
#   3. Power for misspecified models (M2..M12)
# ============================================================================

CANDS <- paste0("M", 1:12)
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
    if (is.null(chi) || is.null(dev)) next
    if (!"Msat" %in% rownames(dev)) next
    sat_dev_adhoc <- dev["Msat", "DEV_adhoc"]
    if (is.na(sat_dev_adhoc)) next

    for (m in CANDS) {
      if (!(m %in% rownames(chi))) next
      r <- chi[m, ]
      d <- dev[m, ]
      if (any(is.na(c(r$chi2_com, r$chi2_MI, r$chi2_D3, d$DEV_adhoc)))) next
      chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
      k <- k + 1
      rows[[k]] <- data.frame(
        method_label = label,
        n            = rr$n,
        miss_rate    = rr$miss_rate,
        rep_id       = rr$rep_id,
        model        = m,
        df           = r$df,
        chi2_com     = r$chi2_com,
        chi2_adhoc   = chi2_adhoc,
        chi2_MI      = r$chi2_MI,
        chi2_D3      = r$chi2_D3,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cat("Loading chi-square data...\n")
chi_a <- extract_chi2("hpc/results-congeniality-amelia/results_combined.rds",
                       "amelia")
chi_p <- extract_chi2("hpc/results-decomp/results_n=250_mr=0.40.rds",
                       "PMM")
cat(sprintf("  amelia: %d (rep, model) rows\n", nrow(chi_a)))
cat(sprintf("  PMM   : %d (rep, model) rows\n", nrow(chi_p)))

# ---- Bias and RMSE vs chi2_com per (method, model) -----------------------
diag_table <- function(df, methods = c("chi2_adhoc", "chi2_MI", "chi2_D3")) {
  out <- list()
  for (m in CANDS) {
    sub <- df[df$model == m, ]
    if (nrow(sub) == 0) next
    row <- data.frame(model = m, df = sub$df[1], n_reps = nrow(sub))
    for (method in methods) {
      diff <- sub[[method]] - sub$chi2_com
      row[[paste0(method, ".bias")]] <- mean(diff)
      row[[paste0(method, ".rmse")]] <- sqrt(mean(diff^2))
    }
    out[[m]] <- row
  }
  do.call(rbind, out)
}

cat("\n========================================================\n")
cat("  Bias and RMSE of chi2 methods vs chi2_com (oracle)\n")
cat("  N=250, mr=0.40, M=50 — amelia (top) vs PMM (bottom)\n")
cat("========================================================\n")

ame_diag <- diag_table(chi_a)
pmm_diag <- diag_table(chi_p)

round_df <- function(df, digits = 2) {
  for (n in names(df)) if (is.numeric(df[[n]])) df[[n]] <- round(df[[n]], digits)
  df
}

cat("\n-- amelia (congenial) --\n")
print(round_df(ame_diag), row.names = FALSE)
cat("\n-- PMM (uncongenial) --\n")
print(round_df(pmm_diag), row.names = FALSE)

# ---- Type I error for true model M1 -------------------------------------
type1 <- function(df) {
  m1 <- df[df$model == "M1", ]
  if (nrow(m1) == 0) return(NULL)
  thresh <- stats::qchisq(1 - ALPHA, df = m1$df[1])
  out <- list(
    threshold = thresh,
    df        = m1$df[1],
    n_reps    = nrow(m1),
    rate_com   = mean(m1$chi2_com   > thresh) * 100,
    rate_adhoc = mean(m1$chi2_adhoc > thresh) * 100,
    rate_MI    = mean(m1$chi2_MI    > thresh) * 100,
    rate_D3    = mean(m1$chi2_D3    > thresh) * 100
  )
  data.frame(out)
}

cat("\n========================================================\n")
cat(sprintf("  Type I error at alpha = %.2f for TRUE model M1\n", ALPHA))
cat("  (should be 5%% if calibrated)\n")
cat("========================================================\n\n")
t1 <- rbind(
  cbind(method = "amelia", round_df(type1(chi_a), 2)),
  cbind(method = "PMM",    round_df(type1(chi_p), 2))
)
print(t1, row.names = FALSE)

# ---- Power for misspecified models --------------------------------------
power_for <- function(df, model, alpha = ALPHA) {
  sub <- df[df$model == model, ]
  if (nrow(sub) == 0) return(NULL)
  thresh <- stats::qchisq(1 - alpha, df = sub$df[1])
  data.frame(
    model     = model,
    df        = sub$df[1],
    n_reps    = nrow(sub),
    com       = mean(sub$chi2_com   > thresh) * 100,
    adhoc     = mean(sub$chi2_adhoc > thresh) * 100,
    MI        = mean(sub$chi2_MI    > thresh) * 100,
    D3        = mean(sub$chi2_D3    > thresh) * 100
  )
}

cat("\n========================================================\n")
cat(sprintf("  Power at alpha = %.2f for misspecified models\n", ALPHA))
cat("========================================================\n")
cat("\n-- amelia (congenial) --\n")
ame_power <- do.call(rbind, lapply(CANDS[-1], function(m) power_for(chi_a, m)))
print(round_df(ame_power), row.names = FALSE)
cat("\n-- PMM (uncongenial) --\n")
pmm_power <- do.call(rbind, lapply(CANDS[-1], function(m) power_for(chi_p, m)))
print(round_df(pmm_power), row.names = FALSE)

# ---- save ----------------------------------------------------------------
fdir <- "hpc/figures-selection"
if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE)
write.csv(round_df(ame_diag),  file.path(fdir, "chi2_bias_rmse_amelia.csv"), row.names = FALSE)
write.csv(round_df(pmm_diag),  file.path(fdir, "chi2_bias_rmse_PMM.csv"),    row.names = FALSE)
write.csv(t1,                  file.path(fdir, "chi2_type1_M1.csv"),         row.names = FALSE)
write.csv(round_df(ame_power), file.path(fdir, "chi2_power_amelia.csv"),     row.names = FALSE)
write.csv(round_df(pmm_power), file.path(fdir, "chi2_power_PMM.csv"),        row.names = FALSE)
cat(sprintf("\n[saved] %s/chi2_*.csv\n", fdir))
