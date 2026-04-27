# ============================================================================
# Local 500-rep mvn_M1 run: does true-model imputation give r_1 = +1?
# ============================================================================
# N=250, mr=0.40, M=50, 500 reps on 16 cores (~36 min)
# Compares against existing amelia at empri=0.01*N (saturated, MI-proper)
# at the same N=250, mr=0.40, M=50, n=1989 reps (results-congeniality-amelia).
#
# Predictions:
#   - delta = tr(RIV)_MI / tr(RIV)_FIML -> 1 (no B inflation, no slack)
#   - r_1 -> +1 if the +0.85 amelia gap was saturated slack
#   - r_1 -> +0.85 if the gap is finite-N / finite-M / improper-MI
# ============================================================================

library(miicsem)
stopifnot(as.character(packageVersion("miicsem")) >= "0.5.2")

results_dir <- "hpc/results-mvn_M1-N250"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps       = 500L,
  n_cores      = 16L,
  seed         = 32897891L,
  sample_sizes = 250,
  miss_rates   = 0.40,
  M            = 50L,
  mice_method  = "mvn_M1",
  results_dir  = results_dir,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nWall time: %.1fs (%.2fh)\n", elapsed, elapsed / 3600))

# Quick summary
rows <- lapply(res[[1]], function(r) {
  if (is.null(r) || is.null(r$dev_df)) return(NULL)
  d <- r$dev_df["M1", ]
  need <- c("loglik_com_at_com", "loglik_com_at_pooled",
            "mean_loglik_at_pooled", "tr_RIV", "tr_RIV_fiml")
  if (any(is.na(d[need]))) return(NULL)
  data.frame(
    tr_RIV_mi   = d$tr_RIV,
    tr_RIV_fiml = d$tr_RIV_fiml,
    term1       = d$mean_loglik_at_pooled - d$loglik_com_at_pooled
  )
})
df <- do.call(rbind, rows)
se <- function(x) stats::sd(x) / sqrt(length(x))
cat(sprintf("\n=== mvn_M1 at N=250, mr=0.40, M=50 ===\n"))
cat(sprintf("  reps converged: %d / 500\n", nrow(df)))
cat(sprintf("  mean tr(RIV)_MI   = %.3f\n", mean(df$tr_RIV_mi)))
cat(sprintf("  mean tr(RIV)_FIML = %.3f\n", mean(df$tr_RIV_fiml)))
cat(sprintf("  mean delta        = %.3f\n", mean(df$tr_RIV_mi / df$tr_RIV_fiml)))
cat(sprintf("  mean Term 1       = %+.3f  (SE %.3f)\n",
            mean(df$term1), se(df$term1)))
cat(sprintf("  r_1 (vs MI tr_RIV)   = %+.3f\n",
            mean(df$term1) / mean(df$tr_RIV_mi)))
cat(sprintf("  r_1 (vs FIML tr_RIV) = %+.3f\n",
            mean(df$term1) / mean(df$tr_RIV_fiml)))
cat("\nCompare to amelia(empri=0.01*N) at same cell:  r_1 = +0.85\n")
cat("Theory prediction:                              r_1 = +1.00\n")
