# ============================================================================
# Analyze the amelia x empri congeniality dose-response
# ============================================================================
# Loads results from results-amelia-empri{0.1, 100}/ and produces:
#   1. Per-(empri, N, mr) M1 decomposition table
#   2. r_term1 dose-response plot data
#   3. Comparison vs the prior congeniality run (N=250, M=50, empri=0.01*N)
# ============================================================================

library(miicsem)

EMPRI_LEVELS <- c("0.1", "100")
FIGURES <- "hpc/figures-empri"
if (!dir.exists(FIGURES)) dir.create(FIGURES, recursive = TRUE)

extract_M1 <- function(rds_path, empri_label) {
  all_res <- readRDS(rds_path)
  rows <- list()
  k <- 0
  for (cond_label in names(all_res)) {
    for (rr in all_res[[cond_label]]) {
      if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
      d <- rr$dev_df
      need <- c("loglik_com_at_com", "loglik_com_at_obs",
                "loglik_com_at_pooled", "mean_loglik_at_pooled", "tr_RIV")
      if (!all(need %in% colnames(d)) || !"M1" %in% rownames(d)) next
      m1 <- d["M1", ]
      if (any(is.na(m1[need]))) next
      k <- k + 1
      rows[[k]] <- data.frame(
        empri      = empri_label,
        n          = rr$n,
        miss_rate  = rr$miss_rate,
        rep_id     = rr$rep_id,
        tr_RIV     = m1$tr_RIV,
        term1      = m1$mean_loglik_at_pooled - m1$loglik_com_at_pooled,
        term2      = m1$loglik_com_at_pooled  - m1$loglik_com_at_obs,
        term3      = m1$loglik_com_at_obs     - m1$loglik_com_at_com,
        total      = m1$mean_loglik_at_pooled - m1$loglik_com_at_com,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cat("Loading empri-sweep results...\n")
dfs <- lapply(EMPRI_LEVELS, function(e) {
  rds <- sprintf("hpc/results-amelia-empri%s/results_combined.rds", e)
  cat(sprintf("  empri=%s*N from %s\n", e, rds))
  extract_M1(rds, e)
})
dlong <- do.call(rbind, dfs)
cat(sprintf("Total: %d rep rows across %d empri levels\n",
            nrow(dlong), length(EMPRI_LEVELS)))

# Also pull the prior congeniality run for comparison (empri = 0.01*N at N=250)
prior_path <- "hpc/results-congeniality-amelia/results_combined.rds"
if (file.exists(prior_path)) {
  prior_df <- extract_M1(prior_path, "0.01")
  cat(sprintf("Plus baseline (empri=0.01*N at N=250): %d rep rows\n",
              nrow(prior_df)))
  dlong <- rbind(dlong, prior_df)
}

# ---- summarize ------------------------------------------------------------
agg <- aggregate(
  cbind(tr_RIV, term1, term2, term3, total) ~ empri + n + miss_rate,
  data = dlong, FUN = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)))
)
# Flatten the nested matrices that aggregate() returns
flat <- data.frame(
  empri     = agg$empri,
  n         = agg$n,
  miss_rate = agg$miss_rate,
  tr_RIV    = agg$tr_RIV[, "mean"],
  T1        = agg$term1[, "mean"],
  se_T1     = agg$term1[, "se"],
  T2        = agg$term2[, "mean"],
  T3        = agg$term3[, "mean"],
  Total     = agg$total[, "mean"],
  stringsAsFactors = FALSE
)
flat$r_T1 <- flat$T1 / flat$tr_RIV
flat$r_T3 <- flat$T3 / flat$tr_RIV
n_reps <- aggregate(rep_id ~ empri + n + miss_rate, data = dlong, length)
flat$n_reps <- n_reps$rep_id[match(
  paste(flat$empri, flat$n, flat$miss_rate),
  paste(n_reps$empri, n_reps$n, n_reps$miss_rate)
)]
flat <- flat[order(flat$miss_rate, flat$n,
                   factor(flat$empri, levels = c("0.01", "0.1", "100"))), ]
flat[, sapply(flat, is.numeric)] <-
  lapply(flat[, sapply(flat, is.numeric)], function(x) round(x, 3))

cat("\n======================================================\n")
cat("  M1 three-term decomposition: amelia x empri sweep\n")
cat("  All M=100 except baseline (empri=0.01, M=50, N=250 only)\n")
cat("  Theory (congenial): r_T1 = +1, r_T3 = -0.5\n")
cat("======================================================\n\n")
print(flat[, c("empri", "n", "miss_rate", "n_reps",
                "tr_RIV", "T1", "se_T1", "r_T1", "r_T3")],
      row.names = FALSE)

# ---- pivot for dose-response visualization data --------------------------
dose <- flat[, c("empri", "n", "miss_rate", "tr_RIV", "T1", "r_T1")]
dose$empri_abs <- ifelse(dose$empri == "0.01", 0.01,
                         ifelse(dose$empri == "0.1", 0.1, 100)) * dose$n

cat("\n-- r_T1 vs absolute empri (T1 should slide negative as empri rises) --\n")
print(dose[, c("empri", "empri_abs", "n", "miss_rate", "tr_RIV", "r_T1")],
      row.names = FALSE)

write.csv(dlong, file.path(FIGURES, "empri_terms_long.csv"),
          row.names = FALSE)
write.csv(flat,  file.path(FIGURES, "empri_M1_summary.csv"),
          row.names = FALSE)
write.csv(dose,  file.path(FIGURES, "empri_dose_response.csv"),
          row.names = FALSE)
cat(sprintf("\n[saved] %s/empri_{terms_long,M1_summary,dose_response}.csv\n",
            FIGURES))
