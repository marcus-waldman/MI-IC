# Test Bartlett correction on chi2_com and combined Bartlett+Wishart on chi2_MI
ALPHA <- 0.05

extract_chi2 <- function(rds_path, label) {
  rs <- readRDS(rds_path)
  is_nested <- !is.null(rs[[1]]) && is.list(rs[[1]]) &&
               !("rep_id" %in% names(rs[[1]]))
  reps <- if (is_nested) unlist(rs, recursive = FALSE) else rs
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
    k <- k + 1
    rows[[k]] <- data.frame(
      method_label = label, n = rr$n, df = r$df,
      chi2_com   = r$chi2_com,
      chi2_adhoc = chi2_adhoc,
      chi2_MI    = r$chi2_MI,
      tr_RIV_M1  = d$tr_RIV,  tr_RIV_Msat = s$tr_RIV,
      p_M1       = d$npar,    p_Msat      = s$npar
    )
  }
  do.call(rbind, rows)
}

cases <- list(
  amelia = "hpc/results-congeniality-amelia/results_combined.rds",
  PMM    = "hpc/results-decomp/results_n=250_mr=0.40.rds",
  mvn_M1 = "hpc/results-mvn_M1-N250/results_combined.rds"
)

# --- Estimate Bartlett scaling empirically from chi2_com ---
cat("Estimating Bartlett scaling factor from chi2_com:\n")
cat("(target: scale chi2 such that mean = df)\n\n")

for (lab in names(cases)) {
  df <- extract_chi2(cases[[lab]], lab)
  bart <- df$df[1] / mean(df$chi2_com)
  cat(sprintf("  %-8s : empirical Bartlett factor = %.4f  (1 - %.3f/N where N=%d)\n",
              lab, bart, (1 - bart) * df$n[1], df$n[1]))
}

# Use a single estimated Bartlett factor (since data conditions identical)
df_a <- extract_chi2(cases[[1]], names(cases)[1])
bart_factor <- df_a$df[1] / mean(df_a$chi2_com)
cat(sprintf("\nUsing Bartlett factor = %.4f\n\n", bart_factor))

cat("Type I error at alpha=0.05 with various corrections (M1, df=22)\n\n")
cat(sprintf("%-12s  %8s  %10s  %10s  %8s  %12s  %12s\n",
            "imputer", "chi2_com", "chi2_com_B", "chi2_adhoc",
            "chi2_MI", "chi2_MI_W", "chi2_MI_BW"))
for (lab in names(cases)) {
  df <- extract_chi2(cases[[lab]], lab)
  thresh <- stats::qchisq(1 - ALPHA, df = df$df[1])
  fac_M1   <- (df$n - df$p_M1   - 1) / df$n
  fac_Msat <- (df$n - df$p_Msat - 1) / df$n
  chi2_MI_W  <- df$chi2_adhoc + df$tr_RIV_M1 * fac_M1 - df$tr_RIV_Msat * fac_Msat
  chi2_MI_BW <- chi2_MI_W * bart_factor
  cat(sprintf("%-12s  %7.1f%%  %9.1f%%  %9.1f%%  %7.1f%%  %11.1f%%  %11.1f%%\n",
              lab,
              mean(df$chi2_com   > thresh) * 100,
              mean(df$chi2_com * bart_factor > thresh) * 100,
              mean(df$chi2_adhoc > thresh) * 100,
              mean(df$chi2_MI    > thresh) * 100,
              mean(chi2_MI_W     > thresh) * 100,
              mean(chi2_MI_BW    > thresh) * 100))
}

cat("\n  chi2_com_B  = chi2_com * Bartlett        (oracle with Bartlett correction)\n")
cat("  chi2_MI_W   = chi2_adhoc + Wishart-corrected tr(RIV) gap\n")
cat("  chi2_MI_BW  = chi2_MI_W * Bartlett        (combined correction)\n")
