# ============================================================================
# Empirical log-likelihood bias check (Study 2, M=100)
# ============================================================================
# Compares per-rep deviances for the TRUE model M1 to the oracle DEV_com:
#
#   bias_adhoc     = E[DEV_adhoc]     - E[DEV_com]
#   bias_raw_pool  = E[DEV_rawPool]   - E[DEV_com]     (DEV_rawPool = MI_DEV - tr_RIV)
#   bias_corrected = E[MI_DEVIANCE]   - E[DEV_com]
#   bias_MR        = E[MR_DEVIANCE]   - E[DEV_com]     (Meng-Rubin anchored)
#
# Theory for the corrected estimator:
#   E[-2*Q_bar(theta_pooled)] - E[-2*ell_com(theta_com)] = -tr(RIV) on the
#   deviance scale, so adding tr(RIV) should restore unbiasedness.
# Pairs within rep (reps see the same dataset), so the bias estimate is
# much more precise than the raw means would suggest.
# ============================================================================
library(miicsem)

RESULTS_M100 <- "hpc/results-full-M100"
FIGURES_DIR  <- "hpc/figures-M100"

cat("Loading deviances...\n")
dev_df <- load_deviances(RESULTS_M100)
cat(sprintf("  %d rows\n", nrow(dev_df)))

# Restrict to the true model M1, with all four deviance variants present
d <- dev_df[dev_df$model == "M1", ]
d <- d[complete.cases(d[, c("DEV_com", "DEV_adhoc",
                            "MI_DEVIANCE", "MR_DEVIANCE", "tr_RIV")]), ]
cat(sprintf("  %d reps with complete M1 deviances\n", nrow(d)))

# Raw pooled deviance (without the tr(RIV) correction)
d$DEV_rawPool <- d$MI_DEVIANCE - d$tr_RIV

# Paired differences from the complete-data oracle
d$diff_adhoc    <- d$DEV_adhoc    - d$DEV_com
d$diff_rawPool  <- d$DEV_rawPool  - d$DEV_com
d$diff_MI       <- d$MI_DEVIANCE  - d$DEV_com
d$diff_MR       <- d$MR_DEVIANCE  - d$DEV_com

summarise <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)))
}

conds <- unique(d[, c("n", "miss_rate")])
conds <- conds[order(conds$n, conds$miss_rate), ]

rows <- lapply(seq_len(nrow(conds)), function(i) {
  cc <- conds[i, ]
  sub <- d[d$n == cc$n & d$miss_rate == cc$miss_rate, ]

  mean_trRIV <- mean(sub$tr_RIV, na.rm = TRUE)
  # Theory: bias_rawPool ≈ -mean(tr_RIV)  (raw pooled dev too small)
  #         bias_MI      ≈ 0              (after correction)

  s_adhoc  <- summarise(sub$diff_adhoc)
  s_raw    <- summarise(sub$diff_rawPool)
  s_MI     <- summarise(sub$diff_MI)
  s_MR     <- summarise(sub$diff_MR)

  data.frame(
    N            = cc$n,
    miss_rate    = cc$miss_rate,
    n_reps       = nrow(sub),
    mean_tr_RIV  = round(mean_trRIV, 3),
    bias_adhoc   = round(s_adhoc["mean"], 3),
    bias_rawPool = round(s_raw["mean"],   3),
    bias_MI      = round(s_MI["mean"],    3),
    bias_MR      = round(s_MR["mean"],    3),
    se_MI        = round(s_MI["se"],      3),
    # How close is the MI correction to theory?  Should be ≈ 0:
    bias_MI_over_se = round(s_MI["mean"] / s_MI["se"], 2),
    # Raw pooled bias vs. -tr(RIV) prediction:
    pred_vs_obs  = round(s_raw["mean"] + mean_trRIV, 3),
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, rows)

cat("\n=========================================================\n")
cat("Deviance bias relative to oracle DEV_com (true model M1)\n")
cat("  bias = E[estimator] - E[DEV_com]  (paired within rep)\n")
cat("  Theory predicts: bias_rawPool ≈ -mean_tr_RIV,  bias_MI ≈ 0\n")
cat("=========================================================\n\n")
print(out, row.names = FALSE)

cat("\nColumn glossary:\n")
cat("  bias_adhoc    : E[-2*mean_m ell_m(theta_m)]      - E[DEV_com]\n")
cat("  bias_rawPool  : E[-2*mean_m ell_m(theta_pooled)] - E[DEV_com]\n")
cat("  bias_MI       : E[MI_DEVIANCE]                    - E[DEV_com]\n")
cat("                  (MI_DEVIANCE = rawPool + tr(RIV))\n")
cat("  bias_MR       : E[MR_DEVIANCE (Meng-Rubin)]       - E[DEV_com]\n")
cat("  pred_vs_obs   : bias_rawPool + mean_tr_RIV\n")
cat("                  (should be ~0 if theory holds)\n")
cat("  bias_MI_over_se: paired-bias t-stat for corrected estimator\n")

write.csv(out, file.path(FIGURES_DIR, "loglik_bias_M100.csv"),
          row.names = FALSE)
cat(sprintf("\n[saved] %s/loglik_bias_M100.csv\n", FIGURES_DIR))
