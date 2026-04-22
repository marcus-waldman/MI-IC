# ============================================================================
# Analyze the three-term bias decomposition results
# ============================================================================
# Loads results-decomp/results_combined.rds, extracts per-rep Term 1/2/3
# values for every (condition, model), and compares empirical means to
# theoretical predictions.
#
#   Term 1 (Imputation) = Q_bar(theta_pooled) - ell_com(theta_pooled)
#                          theory: +tr(RIV)
#   Term 2 (Pooling)    = ell_com(theta_pooled) - ell_com(theta_obs)
#                          theory: ~0
#   Term 3 (Estimation) = ell_com(theta_obs) - ell_com(theta_com)
#                          theory: -tr(RIV)/2
#   Total               = Q_bar(theta_pooled) - ell_com(theta_com)
#                          theory: +tr(RIV)/2
#
# All quantities are on the log-likelihood scale (NOT deviance).
# Run:   Rscript hpc/analyze_decomp.R
# ============================================================================
library(miicsem)

RESULTS <- "hpc/results-decomp"      # transfer target on local machine
FIGURES <- "hpc/figures-decomp"
if (!dir.exists(FIGURES)) dir.create(FIGURES, recursive = TRUE)

cat("Loading combined results...\n")
all_res <- readRDS(file.path(RESULTS, "results_combined.rds"))
cat(sprintf("  %d conditions\n", length(all_res)))

# ---- per-rep long-form table with all 4 logliks + 3 terms -----------------
cat("Building long-form terms table...\n")
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
        n               = rr$n,
        miss_rate       = rr$miss_rate,
        rep_id          = rr$rep_id,
        model           = m,
        tr_RIV          = row$tr_RIV,
        ll_com_at_com   = row$loglik_com_at_com,
        ll_com_at_obs   = row$loglik_com_at_obs,
        ll_com_at_pool  = row$loglik_com_at_pooled,
        Q_bar_at_pool   = row$mean_loglik_at_pooled,
        stringsAsFactors = FALSE
      )
    }
  }
}
dlong <- do.call(rbind, rows)
dlong$term1 <- dlong$Q_bar_at_pool  - dlong$ll_com_at_pool   # Imputation
dlong$term2 <- dlong$ll_com_at_pool - dlong$ll_com_at_obs    # Pooling
dlong$term3 <- dlong$ll_com_at_obs  - dlong$ll_com_at_com    # Estimation
dlong$total <- dlong$Q_bar_at_pool  - dlong$ll_com_at_com
cat(sprintf("  %d rep x model rows\n", nrow(dlong)))

# ---- M1 summary (true model) — the quantity theory is derived for -------
cat("\n=========================================================\n")
cat("  Three-term decomposition for M1 (true model, log-lik scale)\n")
cat("    Term 1 = Imputation        (theory: +tr(RIV))\n")
cat("    Term 2 = Pooling           (theory: ~0)\n")
cat("    Term 3 = Estimation        (theory: -tr(RIV)/2)\n")
cat("    Total  = Q_bar - ell_com   (theory: +tr(RIV)/2)\n")
cat("=========================================================\n")

summarise_M1 <- function(df) {
  m1 <- df[df$model == "M1", ]
  conds <- unique(m1[, c("n", "miss_rate")])
  conds <- conds[order(conds$n, conds$miss_rate), ]
  out <- lapply(seq_len(nrow(conds)), function(i) {
    cc  <- conds[i, ]
    sub <- m1[m1$n == cc$n & m1$miss_rate == cc$miss_rate, ]
    n_reps <- nrow(sub)
    se <- function(x) stats::sd(x) / sqrt(length(x))
    data.frame(
      N          = cc$n,
      miss_rate  = cc$miss_rate,
      n_reps     = n_reps,
      mean_trRIV = mean(sub$tr_RIV),
      term1      = mean(sub$term1),   se1 = se(sub$term1),
      term2      = mean(sub$term2),   se2 = se(sub$term2),
      term3      = mean(sub$term3),   se3 = se(sub$term3),
      total      = mean(sub$total),   se_total = se(sub$total),
      pred_t1    =  mean(sub$tr_RIV),
      pred_t2    = 0,
      pred_t3    = -mean(sub$tr_RIV) / 2,
      pred_tot   =  mean(sub$tr_RIV) / 2,
      diff_t1    = mean(sub$term1)  - mean(sub$tr_RIV),
      diff_t3    = mean(sub$term3)  + mean(sub$tr_RIV) / 2,
      diff_tot   = mean(sub$total)  - mean(sub$tr_RIV) / 2,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
m1_summ <- summarise_M1(dlong)

cat("\n-- Means (empirical) and paired SE --\n")
print(
  round(m1_summ[, c("N", "miss_rate", "n_reps", "mean_trRIV",
                     "term1", "se1", "term2", "se2",
                     "term3", "se3", "total", "se_total")], 3),
  row.names = FALSE
)

cat("\n-- Theory predictions (at empirical mean_trRIV) --\n")
print(
  round(m1_summ[, c("N", "miss_rate", "mean_trRIV",
                     "pred_t1", "pred_t2", "pred_t3", "pred_tot")], 3),
  row.names = FALSE
)

cat("\n-- Difference: empirical - theory (highlights which term is off) --\n")
print(
  round(m1_summ[, c("N", "miss_rate", "mean_trRIV",
                     "diff_t1", "diff_t3", "diff_tot")], 3),
  row.names = FALSE
)

# ---- Also report term-vs-tr(RIV) ratios to see scaling structure ---------
cat("\n-- Terms expressed as multiples of tr(RIV) --\n")
ratios <- data.frame(
  N          = m1_summ$N,
  miss_rate  = m1_summ$miss_rate,
  r_term1    = round(m1_summ$term1 / m1_summ$mean_trRIV, 3),
  r_term2    = round(m1_summ$term2 / m1_summ$mean_trRIV, 3),
  r_term3    = round(m1_summ$term3 / m1_summ$mean_trRIV, 3),
  r_total    = round(m1_summ$total / m1_summ$mean_trRIV, 3),
  pred_r1    = 1.0,
  pred_r2    = 0.0,
  pred_r3    = -0.5,
  pred_total = 0.5
)
print(ratios, row.names = FALSE)

# ---- Full per-model table (may be large) ----------------------------------
cat("\n-- All 13 models (M1..M12 + Msat): mean terms per condition --\n")
agg <- stats::aggregate(
  cbind(tr_RIV, term1, term2, term3, total) ~ n + miss_rate + model,
  data = dlong, FUN = function(x) mean(x, na.rm = TRUE)
)
agg <- agg[order(agg$n, agg$miss_rate, agg$model), ]
# Print just M1 and Msat as examples
cat("\n   M1 (true model):\n")
print(round(agg[agg$model == "M1", ], 3), row.names = FALSE)
cat("\n   Msat (saturated):\n")
print(round(agg[agg$model == "Msat", ], 3), row.names = FALSE)

# ---- save ----------------------------------------------------------------
write.csv(dlong,   file.path(FIGURES, "decomp_terms_long.csv"),     row.names = FALSE)
write.csv(m1_summ, file.path(FIGURES, "decomp_M1_summary.csv"),     row.names = FALSE)
write.csv(agg,     file.path(FIGURES, "decomp_all_models_mean.csv"), row.names = FALSE)
write.csv(ratios,  file.path(FIGURES, "decomp_M1_ratios.csv"),      row.names = FALSE)

cat(sprintf("\n[saved] %s/decomp_terms_long.csv\n", FIGURES))
cat(sprintf("[saved] %s/decomp_M1_summary.csv\n", FIGURES))
cat(sprintf("[saved] %s/decomp_all_models_mean.csv\n", FIGURES))
cat(sprintf("[saved] %s/decomp_M1_ratios.csv\n", FIGURES))
cat("\nDone.\n")
