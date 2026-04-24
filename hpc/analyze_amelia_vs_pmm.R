# ============================================================================
# Compare amelia (2000 reps, congenial joint MVN) vs pmm (500 reps from
# earlier decomp run at N=250, uncongenial chained PMM)
# ============================================================================
library(miicsem)

extract_decomp <- function(results_file, method_label) {
  all_res <- readRDS(results_file)
  rows <- list()
  k <- 0
  for (cond_label in names(all_res)) {
    for (rr in all_res[[cond_label]]) {
      if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
      if (rr$n != 250 || rr$miss_rate != 0.40) next   # restrict to matched cond
      d <- rr$dev_df
      need <- c("loglik_com_at_com", "loglik_com_at_obs",
                "loglik_com_at_pooled", "mean_loglik_at_pooled", "tr_RIV")
      if (!all(need %in% colnames(d))) next
      for (m in rownames(d)) {
        row <- d[m, ]
        if (any(is.na(row[need]))) next
        k <- k + 1
        rows[[k]] <- data.frame(
          method = method_label,
          model  = m,
          tr_RIV = row$tr_RIV,
          term1  = row$mean_loglik_at_pooled - row$loglik_com_at_pooled,
          term2  = row$loglik_com_at_pooled  - row$loglik_com_at_obs,
          term3  = row$loglik_com_at_obs     - row$loglik_com_at_com,
          total  = row$mean_loglik_at_pooled - row$loglik_com_at_com,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

cat("Loading amelia (2000 reps, congenial joint MVN)...\n")
amelia_df <- extract_decomp(
  "hpc/results-congeniality-amelia/results_combined.rds", "amelia")

cat("Loading pmm baseline (from decomp run, 500 reps at N=250)...\n")
pmm_df <- extract_decomp(
  "hpc/results-decomp/results_combined.rds", "pmm")

dlong <- rbind(amelia_df, pmm_df)
cat(sprintf("  amelia: %d rep x model rows\n", nrow(amelia_df)))
cat(sprintf("  pmm   : %d rep x model rows\n", nrow(pmm_df)))

# ---- M1 summary -----------------------------------------------------------
summarise_M1 <- function(sub) {
  se <- function(x) stats::sd(x) / sqrt(length(x))
  data.frame(
    n_reps      = nrow(sub),
    mean_tr_RIV = mean(sub$tr_RIV),
    term1       = mean(sub$term1),  se1 = se(sub$term1),
    term2       = mean(sub$term2),  se2 = se(sub$term2),
    term3       = mean(sub$term3),  se3 = se(sub$term3),
    total       = mean(sub$total),  se_total = se(sub$total),
    r_term1     = mean(sub$term1) / mean(sub$tr_RIV),
    r_term2     = mean(sub$term2) / mean(sub$tr_RIV),
    r_term3     = mean(sub$term3) / mean(sub$tr_RIV),
    r_total     = mean(sub$total) / mean(sub$tr_RIV)
  )
}

m1 <- dlong[dlong$model == "M1", ]
summ <- do.call(rbind, lapply(unique(m1$method), function(m) {
  cbind(method = m, summarise_M1(m1[m1$method == m, ]))
}))

cat("\n======================================================================\n")
cat("  M1 (true model) three-term decomposition: amelia vs pmm\n")
cat("  N = 250, miss_rate = 0.40, M = 50\n")
cat("  Theory (congenial): Term 1 = +tr(RIV), Term 3 = -tr(RIV)/2\n")
cat("======================================================================\n\n")

round_numeric <- function(df, digits = 3) {
  for (nm in names(df)) if (is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], digits)
  df
}

cat("-- Absolute values (log-lik scale) --\n")
print(round_numeric(summ[, c("method", "n_reps", "mean_tr_RIV",
                              "term1", "se1", "term2", "se2",
                              "term3", "se3", "total", "se_total")]),
      row.names = FALSE)

cat("\n-- As multiples of tr(RIV) --\n")
cat("   Theory: r_term1=+1, r_term2=0, r_term3=-0.5, r_total=+0.5\n\n")
print(round_numeric(summ[, c("method", "r_term1", "r_term2",
                              "r_term3", "r_total")]),
      row.names = FALSE)

# ---- Formal sign test on Term 1 under amelia ------------------------------
ame <- m1[m1$method == "amelia", ]
t_stat <- mean(ame$term1) / (stats::sd(ame$term1) / sqrt(nrow(ame)))
predicted <- mean(ame$tr_RIV)
diff_from_theory <- mean(ame$term1) - predicted
se_diff <- stats::sd(ame$term1) / sqrt(nrow(ame))

cat("\n-- Amelia formal tests --\n")
cat(sprintf("  H0: E[Term 1] = 0   |   t = %+.2f   (n=%d)\n", t_stat, nrow(ame)))
cat(sprintf("  H0: E[Term 1] = +tr(RIV) = +%.3f\n", predicted))
cat(sprintf("     observed E[Term 1] = %+.3f  (diff = %+.3f, SE = %.3f,  t = %+.2f)\n",
            mean(ame$term1), diff_from_theory, se_diff, diff_from_theory / se_diff))

# ---- Save ----------------------------------------------------------------
write.csv(dlong, "hpc/figures-M100/amelia_vs_pmm_terms_long.csv",
          row.names = FALSE)
write.csv(summ,  "hpc/figures-M100/amelia_vs_pmm_M1_summary.csv",
          row.names = FALSE)
cat("\n[saved] hpc/figures-M100/amelia_vs_pmm_{terms_long,M1_summary}.csv\n")
