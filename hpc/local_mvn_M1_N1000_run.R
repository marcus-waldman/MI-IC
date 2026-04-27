# ============================================================================
# Local mvn_M1 at N=1000: test the finite-N scaling hypothesis
# ============================================================================
# Question: does r_1 -> 1 as N -> infinity, and at what rate?
#
# Combined with the existing N=250 mvn_M1 result (r_1 = +0.77, n=499),
# this gives us TWO N values to:
#   (a) check whether the assumed form r_1 = 1 - c/N is consistent
#       across N, vs e.g. r_1 = 1 - c/sqrt(N) or some other scaling.
#   (b) estimate the coefficient c if (a) checks out, with proper SE.
#
# Design:
#   N = 1000, mr = 0.40, M = 50, 500 reps
#   mice_method = "mvn_M1" (true-model congenial imputation)
#   16 local cores
#
# Expected wall time ~4.2h (CFA fit cost ~4x N=250)
# ============================================================================

library(miicsem)
stopifnot(as.character(packageVersion("miicsem")) >= "0.5.2")

results_dir <- "hpc/results-mvn_M1-N1000"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

t0 <- proc.time()
res <- run_simulation(
  n_reps       = 500L,
  n_cores      = 16L,
  seed         = 32897891L,         # SAME master seed as N=250 run
  sample_sizes = 1000,
  miss_rates   = 0.40,
  M            = 50L,
  mice_method  = "mvn_M1",
  results_dir  = results_dir,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nWall time: %.1fs (%.2fh)\n", elapsed, elapsed / 3600))

# Pull both N=250 and N=1000 mvn_M1 results, do paired scaling test
extract <- function(rds_path, N_label) {
  rs <- readRDS(rds_path)
  rows <- list(); k <- 0
  for (cl in names(rs)) for (rr in rs[[cl]]) {
    if (is.null(rr) || is.null(rr$dev_df)) next
    d <- rr$dev_df['M1', ]
    if (any(is.na(c(d$loglik_com_at_pooled, d$mean_loglik_at_pooled,
                    d$tr_RIV, d$tr_RIV_fiml)))) next
    k <- k+1
    rows[[k]] <- data.frame(
      N           = N_label,
      tr_RIV_mi   = d$tr_RIV,
      tr_RIV_fiml = d$tr_RIV_fiml,
      term1       = d$mean_loglik_at_pooled - d$loglik_com_at_pooled,
      term3       = d$loglik_com_at_obs    - d$loglik_com_at_com
    )
  }
  do.call(rbind, rows)
}

n250  <- extract("hpc/results-mvn_M1-N250/results_combined.rds",  250)
n1000 <- extract("hpc/results-mvn_M1-N1000/results_combined.rds", 1000)
all   <- rbind(n250, n1000)

se <- function(x) stats::sd(x) / sqrt(length(x))

cat("\n=========================================================\n")
cat("  mvn_M1 (true-model congenial imputation) at N=250 vs N=1000\n")
cat("=========================================================\n")
for (NN in c(250, 1000)) {
  d <- all[all$N == NN, ]
  r1_mi   <- mean(d$term1) / mean(d$tr_RIV_mi)
  r1_fiml <- mean(d$term1) / mean(d$tr_RIV_fiml)
  c1_est  <- NN * (1 - r1_mi)
  cat(sprintf("\nN = %d   (n_reps = %d)\n", NN, nrow(d)))
  cat(sprintf("  T_1 = %+.3f (SE %.3f)\n",  mean(d$term1), se(d$term1)))
  cat(sprintf("  T_3 = %+.3f (SE %.3f)\n",  mean(d$term3), se(d$term3)))
  cat(sprintf("  tr(RIV)_MI   = %.3f\n", mean(d$tr_RIV_mi)))
  cat(sprintf("  tr(RIV)_FIML = %.3f\n", mean(d$tr_RIV_fiml)))
  cat(sprintf("  delta        = %.3f\n",
              mean(d$tr_RIV_mi) / mean(d$tr_RIV_fiml)))
  cat(sprintf("  r_1 (MI)     = %+.3f\n", r1_mi))
  cat(sprintf("  r_1 (FIML)   = %+.3f\n", r1_fiml))
  cat(sprintf("  Implied c1   = N*(1 - r_1_MI) = %.1f\n", c1_est))
}

cat("\n-- Scaling consistency check --\n")
r1_250  <- mean(all$term1[all$N == 250])  / mean(all$tr_RIV_mi[all$N == 250])
r1_1000 <- mean(all$term1[all$N == 1000]) / mean(all$tr_RIV_mi[all$N == 1000])
c1_250  <- 250  * (1 - r1_250)
c1_1000 <- 1000 * (1 - r1_1000)
cat(sprintf("  c1 estimate from N=250  : %.1f\n", c1_250))
cat(sprintf("  c1 estimate from N=1000 : %.1f\n", c1_1000))
cat(sprintf("  Ratio (consistency)     : %.2f\n", c1_1000 / c1_250))
cat("  (1.0 = perfect 1/N scaling; >1 = falls slower; <1 = falls faster)\n")
