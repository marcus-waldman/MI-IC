# ============================================================================
# Aggregate Var[chi2_MI] across all existing simulation cells
# ============================================================================
# For each (imputer, N, mr, M, candidate model), compute:
#   var(chi2_MI), var(chi2_com), var(chi2_adhoc), var(chi2_D3)
#   mean of each, plus mean(tr_RIV_M), mean(tr_RIV_Msat), p_M, p_Msat, df, n_reps
#
# This is the empirical input to the finite-M variance correction
# (Step 1 of plan cozy-imagining-waterfall.md). The output saved as
# hpc/var_chi2_mi_summary.rds is then consumed by Step 3 fit_var_correction.R.
#
# Sanity check: M=50, N=250, mr=0.40, PMM, M1 should give Var[chi2_MI] ≈ 64.1
# (per claude/notes/2026-04-26-chisquare-recalibration.md).
# ============================================================================

CANDS <- paste0("M", 1:12)

# Per-directory metadata. Each entry: list(dir, imputer, M, files = NULL or vector)
SOURCES <- list(
  list(dir = "hpc/results-full-M100",          imputer = "PMM",            M = 100),
  list(dir = "hpc/results-decomp",             imputer = "PMM",            M = 50),
  list(dir = "hpc/results-amelia-empri0",      imputer = "amelia_empri0",  M = 50),
  list(dir = "hpc/results-congeniality-amelia",imputer = "amelia_default", M = 50),
  list(dir = "hpc/results-congeniality-mvn_msat", imputer = "mvn_Msat",    M = 50),
  list(dir = "hpc/results-mvn_M1-N250",        imputer = "mvn_M1",         M = 50),
  list(dir = "hpc/results-mvn_M1-N1000",       imputer = "mvn_M1",         M = 50)
)

# ---- per-rep extraction --------------------------------------------------
extract_rep_rows <- function(rds_path, imputer, M_imp) {
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
    sat_tr_riv    <- dev["Msat", "tr_RIV"]
    sat_npar      <- dev["Msat", "npar"]
    if (is.na(sat_dev_adhoc)) next

    for (m in CANDS) {
      if (!(m %in% rownames(chi))) next
      r <- chi[m, ]
      d <- dev[m, ]
      if (any(is.na(c(r$chi2_com, r$chi2_MI, r$chi2_D3, d$DEV_adhoc)))) next

      chi2_adhoc <- d$DEV_adhoc - sat_dev_adhoc
      k <- k + 1
      rows[[k]] <- data.frame(
        imputer    = imputer,
        M_imp      = M_imp,
        n          = rr$n,
        miss_rate  = rr$miss_rate,
        rep_id     = rr$rep_id,
        model      = m,
        df         = r$df,
        p_M        = d$npar,
        p_Msat     = sat_npar,
        chi2_com   = r$chi2_com,
        chi2_adhoc = chi2_adhoc,
        chi2_MI    = r$chi2_MI,
        chi2_D3    = r$chi2_D3,
        tr_RIV_M   = d$tr_RIV,
        tr_RIV_Msat = sat_tr_riv,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# ---- enumerate per-cell rds files ----------------------------------------
enumerate_files <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  # Prefer per-cell files when present (they have the cell suffix in the name);
  # if only results_combined.rds is available, use that.
  per_cell <- list.files(dir, pattern = "^results_n=.*\\.rds$",
                         full.names = TRUE)
  if (length(per_cell) > 0) return(per_cell)
  combined <- file.path(dir, "results_combined.rds")
  if (file.exists(combined)) return(combined)
  character(0)
}

# ---- driver --------------------------------------------------------------
all_rows <- list()
for (src in SOURCES) {
  files <- enumerate_files(src$dir)
  if (length(files) == 0) {
    cat(sprintf("[skip] no rds in %s\n", src$dir)); next
  }
  for (f in files) {
    cat(sprintf("[load] %s\n", f))
    df <- tryCatch(extract_rep_rows(f, src$imputer, src$M),
                   error = function(e) { cat("  ERROR:", conditionMessage(e),"\n"); NULL })
    if (!is.null(df)) all_rows[[length(all_rows) + 1]] <- df
  }
}
rep_df <- do.call(rbind, all_rows)
cat(sprintf("\nTotal (rep, model) rows: %d across %d source files\n",
            nrow(rep_df), length(all_rows)))

# ---- aggregate to (imputer, M, n, mr, model) ----------------------------
key_cols <- c("imputer", "M_imp", "n", "miss_rate", "model")
key <- do.call(paste, c(rep_df[key_cols], list(sep = "|")))

agg <- function(idx) {
  s <- rep_df[idx, ]
  data.frame(
    imputer    = s$imputer[1],
    M_imp      = s$M_imp[1],
    n          = s$n[1],
    miss_rate  = s$miss_rate[1],
    model      = s$model[1],
    df         = s$df[1],
    p_M        = s$p_M[1],
    p_Msat     = s$p_Msat[1],
    n_reps     = nrow(s),
    mean_chi2_com    = mean(s$chi2_com,    na.rm = TRUE),
    var_chi2_com     = stats::var(s$chi2_com,    na.rm = TRUE),
    mean_chi2_adhoc  = mean(s$chi2_adhoc,  na.rm = TRUE),
    var_chi2_adhoc   = stats::var(s$chi2_adhoc,  na.rm = TRUE),
    mean_chi2_MI     = mean(s$chi2_MI,     na.rm = TRUE),
    var_chi2_MI      = stats::var(s$chi2_MI,     na.rm = TRUE),
    mean_chi2_D3     = mean(s$chi2_D3,     na.rm = TRUE),
    var_chi2_D3      = stats::var(s$chi2_D3,     na.rm = TRUE),
    mean_tr_RIV_M    = mean(s$tr_RIV_M,    na.rm = TRUE),
    var_tr_RIV_M     = stats::var(s$tr_RIV_M,    na.rm = TRUE),
    mean_tr_RIV_Msat = mean(s$tr_RIV_Msat, na.rm = TRUE),
    var_tr_RIV_Msat  = stats::var(s$tr_RIV_Msat, na.rm = TRUE),
    cov_chi2_com_chi2_MI    = stats::cov(s$chi2_com,   s$chi2_MI, use = "complete.obs"),
    cov_chi2_com_chi2_adhoc = stats::cov(s$chi2_com,   s$chi2_adhoc, use = "complete.obs"),
    cov_chi2_MI_chi2_adhoc  = stats::cov(s$chi2_MI,    s$chi2_adhoc, use = "complete.obs"),
    cov_tr_RIV_M_Msat       = stats::cov(s$tr_RIV_M,   s$tr_RIV_Msat, use = "complete.obs"),
    cov_chi2_adhoc_dT       = stats::cov(s$chi2_adhoc, s$tr_RIV_M - s$tr_RIV_Msat, use = "complete.obs"),
    var_dT                  = stats::var(s$tr_RIV_M - s$tr_RIV_Msat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cells <- split(seq_len(nrow(rep_df)), key)
summary_df <- do.call(rbind, lapply(cells, agg))
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$imputer, summary_df$M_imp,
                                summary_df$n, summary_df$miss_rate,
                                summary_df$model), ]

cat(sprintf("\nAggregated summary: %d rows (cell x candidate model)\n", nrow(summary_df)))

# ---- sanity check vs note's reported 64.1 -------------------------------
sanity <- summary_df[summary_df$imputer == "PMM" &
                     summary_df$M_imp == 50 &
                     summary_df$n == 250 &
                     summary_df$miss_rate == 0.40 &
                     summary_df$model == "M1", ]
cat("\nSanity check (target: PMM M=50 N=250 mr=0.40 M1, expect var_chi2_MI ~ 64.1):\n")
if (nrow(sanity) > 0) {
  print(sanity[, c("imputer","M_imp","n","miss_rate","model",
                   "n_reps","mean_chi2_MI","var_chi2_MI","var_chi2_com")],
        row.names = FALSE)
} else {
  cat("  [not found]\n")
}

# ---- save ---------------------------------------------------------------
out_path <- "hpc/var_chi2_mi_summary.rds"
saveRDS(summary_df, out_path)
cat(sprintf("\n[saved] %s\n", out_path))

out_csv <- "hpc/var_chi2_mi_summary.csv"
utils::write.csv(summary_df, out_csv, row.names = FALSE)
cat(sprintf("[saved] %s\n", out_csv))
