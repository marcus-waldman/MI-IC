# Test the user's hypothesis: imputing from M1 (the model being tested under H_0)
# gives chi2 Type I error closer to nominal 5% than imputing from Msat (amelia)
# or from chained PMM.
CANDS <- paste0("M", 1:12)
ALPHA <- 0.05

extract_chi2_M1 <- function(rds_path, label) {
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
    r <- chi["M1", ]; d <- dev["M1", ]
    chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
    k <- k + 1
    rows[[k]] <- data.frame(
      method_label = label, df = r$df,
      chi2_com   = r$chi2_com,
      chi2_adhoc = chi2_adhoc,
      chi2_MI    = r$chi2_MI,
      chi2_D3    = r$chi2_D3
    )
  }
  do.call(rbind, rows)
}

cases <- list(
  amelia = "hpc/results-congeniality-amelia/results_combined.rds",
  PMM    = "hpc/results-decomp/results_n=250_mr=0.40.rds",
  mvn_M1 = "hpc/results-mvn_M1-N250/results_combined.rds"
)

cat("Type I error at alpha=0.05 for true model M1 (df=22)\n")
cat("N=250, mr=0.40, M=50\n\n")
cat(sprintf("%-12s  %6s  %8s  %10s  %8s  %8s\n",
            "imputer", "n_reps", "chi2_com", "chi2_adhoc", "chi2_MI", "chi2_D3"))

for (lab in names(cases)) {
  df <- extract_chi2_M1(cases[[lab]], lab)
  thresh <- stats::qchisq(1 - ALPHA, df = df$df[1])
  cat(sprintf("%-12s  %6d  %7.1f%%  %9.1f%%  %7.1f%%  %7.1f%%\n",
              lab, nrow(df),
              mean(df$chi2_com   > thresh) * 100,
              mean(df$chi2_adhoc > thresh) * 100,
              mean(df$chi2_MI    > thresh) * 100,
              mean(df$chi2_D3    > thresh) * 100))
}

cat("\nMean chi-square (target: df = 22)\n\n")
cat(sprintf("%-12s  %8s  %10s  %8s  %8s\n",
            "imputer", "chi2_com", "chi2_adhoc", "chi2_MI", "chi2_D3"))
for (lab in names(cases)) {
  df <- extract_chi2_M1(cases[[lab]], lab)
  cat(sprintf("%-12s  %8.2f  %10.2f  %8.2f  %8.2f\n",
              lab, mean(df$chi2_com), mean(df$chi2_adhoc),
              mean(df$chi2_MI), mean(df$chi2_D3)))
}

cat("\nVar chi-square (target 2*df = 44)\n\n")
cat(sprintf("%-12s  %8s  %10s  %8s  %8s\n",
            "imputer", "chi2_com", "chi2_adhoc", "chi2_MI", "chi2_D3"))
for (lab in names(cases)) {
  df <- extract_chi2_M1(cases[[lab]], lab)
  cat(sprintf("%-12s  %8.2f  %10.2f  %8.2f  %8.2f\n",
              lab, stats::var(df$chi2_com), stats::var(df$chi2_adhoc),
              stats::var(df$chi2_MI), stats::var(df$chi2_D3)))
}
