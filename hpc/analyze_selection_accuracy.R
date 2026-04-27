# ============================================================================
# Selection-accuracy aggregation
# ============================================================================
# For each (rep, IC method), find the selected model (argmin of IC), then
# tabulate:
#   - % selecting the true model M1
#   - % matching the AIC_com oracle's selection (within-rep agreement)
#   - % matching the BIC_com oracle's selection
#
# Sources used:
#   amelia at N=250, mr=0.40, M=50  -> hpc/results-congeniality-amelia/
#   PMM    at N=250, mr=0.40, M=50  -> hpc/results-decomp/results_n=250_mr=0.40.rds
#   PMM    at the full M=100 grid   -> hpc/results-full-M100/   (for breadth)
# ============================================================================

CANDS <- paste0("M", 1:12)
IC_METHODS <- c("AIC_com", "BIC_com",
                "AIC_adhoc", "BIC_adhoc",
                "AICcd", "MI_AIC", "MI_BIC")

extract_selections <- function(rds_path, label) {
  rs <- readRDS(rds_path)
  # Two possible layouts:
  #   (a) nested:  rs[[cond_label]][[rep_idx]]   (results_combined.rds)
  #   (b) flat:    rs[[rep_idx]]                  (per-condition .rds)
  is_nested <- !is.null(rs[[1]]) && is.list(rs[[1]]) &&
               !("rep_id" %in% names(rs[[1]]))
  reps <- if (is_nested) unlist(rs, recursive = FALSE) else rs

  rows <- list(); k <- 0
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$ic_df)) next
    ic <- rr$ic_df
    if (!all(IC_METHODS %in% colnames(ic))) next
    cand <- if ("model" %in% colnames(ic)) ic$model else rownames(ic)
    sel <- vapply(IC_METHODS, function(m) {
      v <- ic[[m]]
      if (all(is.na(v))) return(NA_character_)
      cand[which.min(v)]
    }, character(1))
    k <- k + 1
    rows[[k]] <- data.frame(
      method_label = label,
      n            = rr$n,
      miss_rate    = rr$miss_rate,
      rep_id       = rr$rep_id,
      t(sel),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

cat("Loading selections...\n")
sel_amelia <- extract_selections(
  "hpc/results-congeniality-amelia/results_combined.rds", "amelia (M=50)"
)
sel_pmm <- extract_selections(
  "hpc/results-decomp/results_n=250_mr=0.40.rds", "PMM (M=50)"
)
cat(sprintf("  amelia: %d reps; PMM: %d reps\n", nrow(sel_amelia), nrow(sel_pmm)))

# ---- % selecting M1 (true-model accuracy) -------------------------------
pct_M1 <- function(df) {
  vapply(IC_METHODS, function(m) {
    sum(df[[m]] == "M1", na.rm = TRUE) / sum(!is.na(df[[m]])) * 100
  }, numeric(1))
}
acc_amelia <- pct_M1(sel_amelia)
acc_pmm    <- pct_M1(sel_pmm)

cat("\n=========================================================\n")
cat("  % selecting true model M1 — N=250, mr=0.40, M=50\n")
cat("=========================================================\n\n")
out_M1 <- data.frame(
  IC      = IC_METHODS,
  amelia  = round(acc_amelia, 1),
  PMM     = round(acc_pmm,    1),
  diff    = round(acc_amelia - acc_pmm, 1)
)
print(out_M1, row.names = FALSE)

# ---- Oracle agreement (within-rep match with AIC_com / BIC_com) --------
oracle_agree <- function(df, oracle_col, method_cols) {
  vapply(method_cols, function(m) {
    valid <- !is.na(df[[m]]) & !is.na(df[[oracle_col]])
    if (sum(valid) == 0) return(NA_real_)
    sum(df[[m]][valid] == df[[oracle_col]][valid]) / sum(valid) * 100
  }, numeric(1))
}
aic_partners <- c("AIC_adhoc", "AICcd", "MI_AIC")
bic_partners <- c("BIC_adhoc", "MI_BIC")

cat("\n=========================================================\n")
cat("  % within-rep agreement with AIC_com (oracle)\n")
cat("=========================================================\n\n")
agr_aic <- data.frame(
  IC      = aic_partners,
  amelia  = round(oracle_agree(sel_amelia, "AIC_com", aic_partners), 1),
  PMM     = round(oracle_agree(sel_pmm,    "AIC_com", aic_partners), 1)
)
print(agr_aic, row.names = FALSE)

cat("\n=========================================================\n")
cat("  % within-rep agreement with BIC_com (oracle)\n")
cat("=========================================================\n\n")
agr_bic <- data.frame(
  IC      = bic_partners,
  amelia  = round(oracle_agree(sel_amelia, "BIC_com", bic_partners), 1),
  PMM     = round(oracle_agree(sel_pmm,    "BIC_com", bic_partners), 1)
)
print(agr_bic, row.names = FALSE)

# ---- save ---------------------------------------------------------------
fdir <- "hpc/figures-selection"
if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE)
write.csv(out_M1, file.path(fdir, "pct_M1_amelia_vs_PMM.csv"), row.names = FALSE)
write.csv(agr_aic, file.path(fdir, "agreement_AIC_com.csv"),   row.names = FALSE)
write.csv(agr_bic, file.path(fdir, "agreement_BIC_com.csv"),   row.names = FALSE)
cat(sprintf("\n[saved] %s/{pct_M1, agreement_AIC_com, agreement_BIC_com}.csv\n", fdir))
