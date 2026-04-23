# ============================================================================
# Analyze congeniality test: pmm vs norm vs mvn_msat
# ============================================================================
# Loads results-congeniality-<method>/ for each method and compares the
# three-term bias decomposition for M1.  Key prediction:
#
#   If congeniality is the cause of Term 1 < 0:
#     pmm       -> Term 1 ~ -1/2 tr(RIV)   (matches Study 2 decomp)
#     norm      -> Term 1 intermediate
#     mvn_msat  -> Term 1 ~ +tr(RIV)       (matches theory)
#
# All three use the same seeds, so any difference is purely due to the
# imputation method.
# ============================================================================

library(miicsem)

METHODS <- c("pmm", "norm", "mvn_msat")
FIGURES <- "hpc/figures-congeniality"
if (!dir.exists(FIGURES)) dir.create(FIGURES, recursive = TRUE)

extract_decomp <- function(results_dir, method_label) {
  combined <- file.path(results_dir, "results_combined.rds")
  if (!file.exists(combined)) {
    cat(sprintf("  [skip] %s not found\n", combined))
    return(NULL)
  }
  all_res <- readRDS(combined)
  rows <- list()
  k <- 0
  for (cond_label in names(all_res)) {
    for (rr in all_res[[cond_label]]) {
      if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
      d <- rr$dev_df
      need <- c("loglik_com_at_com", "loglik_com_at_obs",
                "loglik_com_at_pooled", "mean_loglik_at_pooled", "tr_RIV")
      if (!all(need %in% colnames(d))) next
      for (m in rownames(d)) {
        row <- d[m, ]
        if (any(is.na(row[need]))) next
        k <- k + 1
        rows[[k]] <- data.frame(
          method    = method_label,
          n         = rr$n,
          miss_rate = rr$miss_rate,
          rep_id    = rr$rep_id,
          model     = m,
          tr_RIV    = row$tr_RIV,
          term1     = row$mean_loglik_at_pooled - row$loglik_com_at_pooled,
          term2     = row$loglik_com_at_pooled  - row$loglik_com_at_obs,
          term3     = row$loglik_com_at_obs     - row$loglik_com_at_com,
          total     = row$mean_loglik_at_pooled - row$loglik_com_at_com,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

cat("Loading decomposition data for each method...\n")
dlong_list <- lapply(METHODS, function(m) {
  extract_decomp(sprintf("hpc/results-congeniality-%s", m), m)
})
names(dlong_list) <- METHODS
dlong <- do.call(rbind, dlong_list)
cat(sprintf("  %d rep x model rows across %d methods\n",
            nrow(dlong), length(METHODS)))

# ---- M1 comparison across methods ----------------------------------------
m1 <- dlong[dlong$model == "M1", ]

summarise <- function(sub) {
  se <- function(x) stats::sd(x) / sqrt(length(x))
  data.frame(
    n_reps      = nrow(sub),
    mean_tr_RIV = mean(sub$tr_RIV),
    term1       = mean(sub$term1), se1 = se(sub$term1),
    term2       = mean(sub$term2), se2 = se(sub$term2),
    term3       = mean(sub$term3), se3 = se(sub$term3),
    total       = mean(sub$total), se_total = se(sub$total),
    r_term1     = mean(sub$term1) / mean(sub$tr_RIV),
    r_term2     = mean(sub$term2) / mean(sub$tr_RIV),
    r_term3     = mean(sub$term3) / mean(sub$tr_RIV),
    r_total     = mean(sub$total) / mean(sub$tr_RIV)
  )
}

summ_rows <- lapply(METHODS, function(m) {
  cbind(method = m, summarise(m1[m1$method == m, ]))
})
summ <- do.call(rbind, summ_rows)

cat("\n==========================================================\n")
cat("  M1 (true model) three-term decomposition by method\n")
cat("  N = 250, miss_rate = 0.40, M = 50\n")
cat("  Theory congenial: Term 1 = +tr(RIV), Term 3 = -tr(RIV)/2\n")
cat("==========================================================\n\n")
print(round(summ[, c("method", "n_reps", "mean_tr_RIV",
                      "term1", "se1", "term2", "se2",
                      "term3", "se3", "total", "se_total")], 3),
      row.names = FALSE)

cat("\n-- Terms as multiples of tr(RIV) --\n")
print(round(summ[, c("method", "r_term1", "r_term2",
                      "r_term3", "r_total")], 3),
      row.names = FALSE)

cat("\nTheoretical prediction (congenial): r_term1=+1, r_term2=0, ",
    "r_term3=-0.5, r_total=+0.5\n", sep = "")

# ---- save --------------------------------------------------------------
write.csv(dlong, file.path(FIGURES, "congeniality_terms_long.csv"),
          row.names = FALSE)
write.csv(summ,  file.path(FIGURES, "congeniality_M1_summary.csv"),
          row.names = FALSE)
cat(sprintf("\n[saved] %s/congeniality_terms_long.csv\n", FIGURES))
cat(sprintf("[saved] %s/congeniality_M1_summary.csv\n", FIGURES))
