ame <- readRDS('hpc/results-congeniality-amelia/results_combined.rds')
mvn <- readRDS('hpc/results-mvn_M1-N250/results_combined.rds')

extract <- function(rs) {
  rows <- list(); k <- 0
  for (cl in names(rs)) for (rr in rs[[cl]]) {
    if (is.null(rr) || is.null(rr$dev_df)) next
    d <- rr$dev_df['M1', ]
    if (any(is.na(c(d$loglik_com_at_pooled, d$mean_loglik_at_pooled,
                    d$tr_RIV, d$tr_RIV_fiml)))) next
    k <- k+1
    rows[[k]] <- data.frame(
      tr_RIV_mi = d$tr_RIV, tr_RIV_fiml = d$tr_RIV_fiml,
      term1 = d$mean_loglik_at_pooled - d$loglik_com_at_pooled
    )
  }
  do.call(rbind, rows)
}
a <- extract(ame); m <- extract(mvn)
se <- function(x) sd(x) / sqrt(length(x))

cat(sprintf("amelia : n=%d  T_1=%+.3f (SE %.3f)  tr_RIV_MI=%.3f  tr_RIV_FIML=%.3f  delta=%.3f\n",
            nrow(a), mean(a$term1), se(a$term1),
            mean(a$tr_RIV_mi), mean(a$tr_RIV_fiml),
            mean(a$tr_RIV_mi) / mean(a$tr_RIV_fiml)))
cat(sprintf("mvn_M1 : n=%d  T_1=%+.3f (SE %.3f)  tr_RIV_MI=%.3f  tr_RIV_FIML=%.3f  delta=%.3f\n",
            nrow(m), mean(m$term1), se(m$term1),
            mean(m$tr_RIV_mi), mean(m$tr_RIV_fiml),
            mean(m$tr_RIV_mi) / mean(m$tr_RIV_fiml)))

cat(sprintf("\nDifference (amelia - mvn_M1):\n"))
cat(sprintf("  Delta T_1     = %+.3f\n", mean(a$term1) - mean(m$term1)))
cat(sprintf("  Delta tr_RIV_MI = %+.3f\n",
            mean(a$tr_RIV_mi) - mean(m$tr_RIV_mi)))
cat(sprintf("  Delta tr_RIV_FIML = %+.3f\n",
            mean(a$tr_RIV_fiml) - mean(m$tr_RIV_fiml)))
