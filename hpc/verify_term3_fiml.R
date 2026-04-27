# ============================================================================
# Verify item #2: T_3 = -1/2 * tr(RIV)_FIML across all cells
# ============================================================================
# Item #1 already showed tr(RIV)_FIML blows up at N=100 (140, 27, 16) due
# to FIML instability with 32 CFA params on ~80 effectively-complete
# cases.  Question now: does Term 3 still ride that inflation, or is it
# a stable -1.16-ish absolute value (suggesting it ISN'T governed by
# tr(RIV)_FIML, but by something else)?
#
# Compute per-cell:
#   T_3        = mean(loglik_com_at_obs - loglik_com_at_com)
#   tr_RIV_fiml
#   ratio      = -2 * T_3 / tr_RIV_fiml      (should be 1 if claim holds)
# ============================================================================

extract_T3 <- function(rds_path, empri_label) {
  all_res <- readRDS(rds_path)
  rows <- list()
  k <- 0
  for (cond_label in names(all_res)) {
    for (rr in all_res[[cond_label]]) {
      if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
      d <- rr$dev_df
      need <- c("loglik_com_at_com", "loglik_com_at_obs",
                "tr_RIV_fiml")
      if (!all(need %in% colnames(d)) || !"M1" %in% rownames(d)) next
      m1 <- d["M1", ]
      if (any(is.na(m1[need]))) next
      k <- k + 1
      rows[[k]] <- data.frame(
        empri       = empri_label,
        n           = rr$n,
        miss_rate   = rr$miss_rate,
        rep_id      = rr$rep_id,
        tr_RIV_fiml = m1$tr_RIV_fiml,
        term3       = m1$loglik_com_at_obs - m1$loglik_com_at_com,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

paths <- list(
  list(rds = "hpc/results-amelia-empri0.1/results_combined.rds",     lab = "0.1"),
  list(rds = "hpc/results-amelia-empri100/results_combined.rds",     lab = "100"),
  list(rds = "hpc/results-congeniality-amelia/results_combined.rds", lab = "0.01_baseline")
)
dlong <- do.call(rbind, lapply(paths, function(p) extract_T3(p$rds, p$lab)))
cat(sprintf("Loaded %d total rep rows\n", nrow(dlong)))

# Per-cell: mean T_3, mean tr_RIV_fiml, and the implied ratio
agg <- aggregate(
  cbind(tr_RIV_fiml, term3) ~ empri + n + miss_rate,
  data = dlong,
  FUN  = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)))
)
flat <- data.frame(
  empri        = agg$empri,
  n            = agg$n,
  miss_rate    = agg$miss_rate,
  trRIV_fiml   = agg$tr_RIV_fiml[, "mean"],
  T3           = agg$term3[, "mean"],
  se_T3        = agg$term3[, "se"],
  neg2T3       = -2 * agg$term3[, "mean"],
  ratio        = -2 * agg$term3[, "mean"] / agg$tr_RIV_fiml[, "mean"]
)
flat <- flat[order(flat$miss_rate, flat$n,
                   factor(flat$empri,
                          levels = c("0.01_baseline", "0.1", "100"))), ]
flat[sapply(flat, is.numeric)] <-
  lapply(flat[sapply(flat, is.numeric)], function(x) round(x, 3))

cat("\n========================================================\n")
cat("  T_3 vs tr(RIV)_FIML — does -2 T_3 / tr(RIV)_FIML = 1?\n")
cat("========================================================\n\n")
print(flat[, c("empri", "n", "miss_rate", "trRIV_fiml",
                "T3", "se_T3", "neg2T3", "ratio")],
      row.names = FALSE)

# Is T_3 approximately CONSTANT in absolute units when stratified by mr?
cat("\n-- T_3 absolute values, stratified by mr (looking for invariance) --\n")
abs_summ <- aggregate(T3 ~ miss_rate, data = flat,
                      FUN = function(x) c(min = min(x), max = max(x),
                                          ratio = max(x) / min(x)))
print(abs_summ)

cat("\n-- Same comparison restricted to N >= 500 (where FIML is stable) --\n")
big_n <- flat[flat$n >= 500, ]
print(big_n[, c("empri", "n", "miss_rate", "trRIV_fiml", "T3",
                 "neg2T3", "ratio")], row.names = FALSE)

write.csv(flat, "hpc/figures-empri/term3_vs_trRIV_fiml.csv", row.names = FALSE)
cat("\n[saved] hpc/figures-empri/term3_vs_trRIV_fiml.csv\n")
