# ============================================================================
# Verify diagnostic claim: delta = tr(RIV)_MI / tr(RIV)_FIML tracks r_1
# ============================================================================
# Hypothesis (item 1 in the bad-news list):
#   delta near 1 -> imputation is well-congenial -> r_1 near +1
#   delta >> 1   -> uncongenial inflation of B -> r_1 << +1
# If this holds across the empri sweep + prior congeniality run, we have
# an actionable diagnostic.  If delta and r_1 don't co-move, the
# "diagnostic" claim collapses to wishful thinking.
# ============================================================================

library(miicsem)

extract_M1 <- function(rds_path, empri_label) {
  all_res <- readRDS(rds_path)
  rows <- list()
  k <- 0
  for (cond_label in names(all_res)) {
    for (rr in all_res[[cond_label]]) {
      if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
      d <- rr$dev_df
      need <- c("loglik_com_at_com", "loglik_com_at_pooled",
                "mean_loglik_at_pooled", "tr_RIV", "tr_RIV_fiml")
      if (!all(need %in% colnames(d)) || !"M1" %in% rownames(d)) next
      m1 <- d["M1", ]
      if (any(is.na(m1[need]))) next
      k <- k + 1
      rows[[k]] <- data.frame(
        empri      = empri_label,
        n          = rr$n,
        miss_rate  = rr$miss_rate,
        rep_id     = rr$rep_id,
        tr_RIV_mi  = m1$tr_RIV,
        tr_RIV_fiml = m1$tr_RIV_fiml,
        term1      = m1$mean_loglik_at_pooled - m1$loglik_com_at_pooled,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

paths <- list(
  list(rds = "hpc/results-amelia-empri0.1/results_combined.rds",     lab = "0.1"),
  list(rds = "hpc/results-amelia-empri100/results_combined.rds",     lab = "100"),
  list(rds = "hpc/results-congeniality-amelia/results_combined.rds", lab = "0.01")
)

cat("Loading per-rep tr(RIV) decompositions...\n")
dlong <- do.call(rbind, lapply(paths, function(p) extract_M1(p$rds, p$lab)))
cat(sprintf("  %d total rows\n", nrow(dlong)))

# Per-cell summary
agg <- aggregate(
  cbind(tr_RIV_mi, tr_RIV_fiml, term1) ~ empri + n + miss_rate,
  data = dlong, FUN = mean
)
agg$delta   <- agg$tr_RIV_mi / agg$tr_RIV_fiml
agg$r1_mi   <- agg$term1     / agg$tr_RIV_mi      # r_1 vs sample tr(RIV)
agg$r1_fiml <- agg$term1     / agg$tr_RIV_fiml    # r_1 vs theoretical tr(RIV)
agg$n_reps  <- aggregate(rep_id ~ empri + n + miss_rate,
                          data = dlong, length)$rep_id
agg <- agg[order(agg$miss_rate, agg$n,
                 factor(agg$empri, levels = c("0.01", "0.1", "100"))), ]

# Round numeric columns
agg[sapply(agg, is.numeric)] <-
  lapply(agg[sapply(agg, is.numeric)], function(x) round(x, 3))

cat("\n======================================================\n")
cat("  delta = tr(RIV)_MI / tr(RIV)_FIML  vs r_1 (M1, both definitions)\n")
cat("======================================================\n\n")
print(agg[, c("empri", "n", "miss_rate", "n_reps",
               "tr_RIV_mi", "tr_RIV_fiml", "delta",
               "term1", "r1_mi", "r1_fiml")],
      row.names = FALSE)

# Correlation diagnostic
cat("\n-- Correlation between delta and r_1 across the 13 cells --\n")
cat(sprintf("  cor(delta, r_1_mi)   = %+.3f\n",
            cor(agg$delta, agg$r1_mi)))
cat(sprintf("  cor(delta, r_1_fiml) = %+.3f\n",
            cor(agg$delta, agg$r1_fiml)))
cat(sprintf("  cor(log10(delta), r_1_fiml) = %+.3f\n",
            cor(log10(agg$delta), agg$r1_fiml)))

# Save
fdir <- "hpc/figures-empri"
if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE)
write.csv(agg, file.path(fdir, "delta_diagnostic.csv"), row.names = FALSE)
cat(sprintf("\n[saved] %s/delta_diagnostic.csv\n", fdir))
