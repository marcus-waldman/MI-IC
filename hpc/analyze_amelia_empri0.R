# ============================================================================
# Compare three Amelia variants and PMM at N=500, varying mr:
#   - amelia, empri = 0       (no prior; this run)
#   - amelia, empri = 0.01·N  (mild prior; existing hpc/results-amelia-empri0.1/)
#   - amelia, empri = 100·N   (severe prior; existing hpc/results-amelia-empri100/)
#   - PMM,    M=100           (baseline; existing hpc/results-full-M100/)
#
# For each cell, report:
#   chi2_com / chi2_adhoc / chi2_MI / tr_RIV_M1 / tr_RIV_Msat / Δ_emp / Δ_pred
# and the ratio Δ_emp / Δ_pred — the v0 derivation predicts this ratio = 1
# under saturated proper imputation with no prior.
# ============================================================================

extract_chi2 <- function(rds_path, label, n_target = NULL, mr_target = NULL) {
  if (!file.exists(rds_path)) return(NULL)
  rs <- readRDS(rds_path)
  is_nested <- !is.null(rs[[1]]) && is.list(rs[[1]]) &&
               !("rep_id" %in% names(rs[[1]]))
  reps <- if (is_nested) unlist(rs, recursive = FALSE) else rs
  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed)) next
    if (!is.null(n_target)  && !is.null(rr$n)         && rr$n         != n_target)  next
    if (!is.null(mr_target) && !is.null(rr$miss_rate) && abs(rr$miss_rate - mr_target) > 1e-6) next
    chi <- rr$chi2_df; dev <- rr$dev_df
    if (is.null(chi) || is.null(dev) || !"M1" %in% rownames(chi)) next
    if (!"Msat" %in% rownames(dev)) next
    sat_dev_adhoc <- dev["Msat", "DEV_adhoc"]
    if (is.na(sat_dev_adhoc)) next
    r <- chi["M1", ]; d <- dev["M1", ]; s <- dev["Msat", ]
    chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
    if (any(is.na(c(d$tr_RIV, s$tr_RIV, r$chi2_com, r$chi2_MI)))) next
    k <- k + 1
    rows[[k]] <- data.frame(
      label = label, n = rr$n, mr = rr$miss_rate, df = r$df,
      chi2_com = r$chi2_com, chi2_adhoc = chi2_adhoc, chi2_MI = r$chi2_MI,
      tr_M1 = d$tr_RIV, tr_Msat = s$tr_RIV
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

cells <- list(
  list(label = "amelia_e0_N500_mr0.10",    path = "hpc/results-amelia-empri0/results_n=500_mr=0.10.rds",   n=500, mr=0.10),
  list(label = "amelia_e0_N500_mr0.25",    path = "hpc/results-amelia-empri0/results_n=500_mr=0.25.rds",   n=500, mr=0.25),
  list(label = "amelia_e0_N500_mr0.40",    path = "hpc/results-amelia-empri0/results_n=500_mr=0.40.rds",   n=500, mr=0.40),
  list(label = "amelia_e0.1_N500_mr0.20",  path = "hpc/results-amelia-empri0.1/results_n=500_mr=0.20.rds", n=500, mr=0.20),
  list(label = "amelia_e0.1_N500_mr0.40",  path = "hpc/results-amelia-empri0.1/results_n=500_mr=0.40.rds", n=500, mr=0.40),
  list(label = "amelia_e100_N500_mr0.20",  path = "hpc/results-amelia-empri100/results_n=500_mr=0.20.rds", n=500, mr=0.20),
  list(label = "amelia_e100_N500_mr0.40",  path = "hpc/results-amelia-empri100/results_n=500_mr=0.40.rds", n=500, mr=0.40),
  list(label = "PMM_M100_N500_mr0.10",     path = "hpc/results-full-M100/results_n=500_mr=0.10.rds",        n=500, mr=0.10),
  list(label = "PMM_M100_N500_mr0.25",     path = "hpc/results-full-M100/results_n=500_mr=0.25.rds",        n=500, mr=0.25),
  list(label = "PMM_M100_N500_mr0.40",     path = "hpc/results-full-M100/results_n=500_mr=0.40.rds",        n=500, mr=0.40)
)

cat("Comparing imputation methods at N=500 across mr values\n")
cat("Δ_pred = 2 · [tr(RIV)_Msat − tr(RIV)_M1]\n")
cat("ratio  = Δ_emp / Δ_pred   (v0 predicts ≈ 1 for saturated proper MI without prior)\n\n")
cat(sprintf("%-25s %5s  %4s  %8s  %10s  %8s  %8s  %10s  %8s  %8s  %5s\n",
            "cell", "n_rep", "df", "chi2_com", "chi2_adhoc", "chi2_MI",
            "tr_M1", "tr_Msat", "Δ_emp", "Δ_pred", "ratio"))

for (cell in cells) {
  df <- extract_chi2(cell$path, cell$label, n_target = cell$n, mr_target = cell$mr)
  if (is.null(df)) {
    cat(sprintf("%-25s  (no data)\n", cell$label))
    next
  }
  e_com   <- mean(df$chi2_com)
  e_adhoc <- mean(df$chi2_adhoc)
  e_MI    <- mean(df$chi2_MI)
  e_trM1  <- mean(df$tr_M1)
  e_trMs  <- mean(df$tr_Msat)
  delta_emp  <- e_adhoc - e_com
  delta_pred <- 2 * (e_trMs - e_trM1)
  ratio <- if (abs(delta_pred) > 0.01) delta_emp / delta_pred else NA
  cat(sprintf("%-25s %5d  %4.0f  %8.2f  %10.2f  %8.2f  %8.3f  %10.3f  %8.2f  %8.2f  %5.2f\n",
              cell$label, nrow(df), df$df[1],
              e_com, e_adhoc, e_MI, e_trM1, e_trMs,
              delta_emp, delta_pred, ratio))
}

cat("\nInterpretation:\n")
cat("  ratio ≈ 1.00  -> formula matches; v0 derivation correct for this cell\n")
cat("  ratio < 1.00  -> imputation noise damped (e.g., by prior shrinkage)\n")
cat("  ratio NA     -> Δ_pred ≈ 0 (no flexibility gap to test)\n")
