# ============================================================================
# Local Congeniality Test: amelia (proper, joint-coherent) vs mvn_msat
#                          (improper MLE plug-in, joint-coherent)
# ============================================================================
# Tests whether Term 1 becomes +tr(RIV) under congenial joint-MVN
# imputation.  Both methods use a joint Gaussian model; they differ in
# proper (amelia, bootstrap-based parameter uncertainty) vs improper
# (mvn_msat, MLE plug-in) MI.
#
# Design:
#   N             = 250
#   miss_rate     = 0.40
#   M             = 50
#   reps          = 2000 per method, SPLIT INTO 4 BATCHES OF 500
#   cores         = 16
#
# Batching rationale: the previous 2000-rep single-batch run crashed at
# the per-condition save step (~2h14m in), losing all work.  By chunking
# into 4 x 500 reps with per-batch .rds files, at most one batch's work
# is lost on a crash, and incomplete methods can resume from where they
# left off.
#
# Output layout:
#   hpc/results-congeniality-<method>/batch<b>/
#     results_n=250_mr=0.40.rds     per-batch 500 reps
#   hpc/results-congeniality-<method>/combined_2000.rds
#     rolled up across all 4 batches after completion
# ============================================================================

library(miicsem)
stopifnot(as.character(packageVersion("miicsem")) >= "0.4.2")
cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

METHODS <- c("amelia", "mvn_msat")

N            <- 250
MR           <- 0.40
M_IMP        <- 50L
TOTAL_REPS   <- 2000L
BATCH_SIZE   <- 500L
N_BATCHES    <- TOTAL_REPS %/% BATCH_SIZE    # 4
N_CORES      <- 16L
BASE_SEED    <- 32897891L

overall_t0 <- proc.time()

for (method in METHODS) {
  cat("\n######################################################\n")
  cat(sprintf("###  METHOD: %s\n", method))
  cat("######################################################\n")

  method_dir <- sprintf("hpc/results-congeniality-%s", method)
  if (!dir.exists(method_dir)) dir.create(method_dir, recursive = TRUE)

  for (b in seq_len(N_BATCHES)) {
    batch_dir <- file.path(method_dir, sprintf("batch%d", b))
    batch_combined <- file.path(batch_dir, "results_combined.rds")

    if (file.exists(batch_combined)) {
      cat(sprintf("\n[batch %d/%d] already complete -> %s (skip)\n",
                  b, N_BATCHES, batch_combined))
      next
    }
    if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)

    cat(sprintf("\n------------------------------------------------\n"))
    cat(sprintf(" [%s] Batch %d / %d  (500 reps, seed = %d)\n",
                method, b, N_BATCHES, BASE_SEED + (b - 1) * BATCH_SIZE))
    cat(sprintf("------------------------------------------------\n"))

    t0 <- proc.time()
    res <- tryCatch({
      run_simulation(
        n_reps       = BATCH_SIZE,
        n_cores      = N_CORES,
        seed         = BASE_SEED + (b - 1) * BATCH_SIZE,
        sample_sizes = N,
        miss_rates   = MR,
        M            = M_IMP,
        mice_method  = method,
        results_dir  = batch_dir,
        verbose      = TRUE
      )
    }, error = function(e) {
      cat(sprintf("  [ERROR] %s\n", conditionMessage(e)))
      NULL
    })
    elapsed <- (proc.time() - t0)[3]
    cat(sprintf("  batch %d wall: %.1fs (%.2fh)\n",
                b, elapsed, elapsed / 3600))
    if (is.null(res)) break
  }

  # --- roll up completed batches for this method ---------------------------
  cat(sprintf("\n[%s] Rolling up batches...\n", method))
  batch_files <- file.path(method_dir, sprintf("batch%d", seq_len(N_BATCHES)),
                           "results_combined.rds")
  existing <- batch_files[file.exists(batch_files)]
  cat(sprintf("  %d/%d batches complete\n", length(existing), N_BATCHES))
  if (length(existing) == 0) next

  merged_list <- list()
  cond_label <- sprintf("n=%d_mr=%.2f", N, MR)
  merged_list[[cond_label]] <- list()
  rep_offset <- 0L
  for (f in existing) {
    chunk <- readRDS(f)
    reps <- chunk[[cond_label]]
    for (i in seq_along(reps)) {
      r <- reps[[i]]
      if (!is.null(r) && !is.null(r$rep_id)) {
        r$rep_id <- r$rep_id + rep_offset
      }
      merged_list[[cond_label]][[length(merged_list[[cond_label]]) + 1]] <- r
    }
    rep_offset <- rep_offset + BATCH_SIZE
  }
  out_file <- file.path(method_dir, "results_combined.rds")
  saveRDS(merged_list, out_file)
  cat(sprintf("  saved combined -> %s (%d reps)\n",
              out_file, length(merged_list[[cond_label]])))

  # --- quick M1 decomp summary --------------------------------------------
  rows <- lapply(merged_list[[cond_label]], function(r) {
    if (is.null(r) || is.null(r$dev_df)) return(NULL)
    d <- r$dev_df["M1", ]
    need <- c("loglik_com_at_com", "loglik_com_at_obs",
              "loglik_com_at_pooled", "mean_loglik_at_pooled", "tr_RIV")
    if (any(is.na(d[need]))) return(NULL)
    data.frame(
      tr_RIV = d$tr_RIV,
      term1  = d$mean_loglik_at_pooled - d$loglik_com_at_pooled,
      term2  = d$loglik_com_at_pooled  - d$loglik_com_at_obs,
      term3  = d$loglik_com_at_obs     - d$loglik_com_at_com,
      total  = d$mean_loglik_at_pooled - d$loglik_com_at_com
    )
  })
  df <- do.call(rbind, rows)
  if (!is.null(df) && nrow(df) > 0) {
    cat(sprintf("\n[%s] M1 decomposition over %d reps (mr=%.2f, N=%d, M=%d):\n",
                method, nrow(df), MR, N, M_IMP))
    cat(sprintf("  mean tr(RIV) = %.3f\n", mean(df$tr_RIV)))
    cat(sprintf("  Term 1 = %+.3f  (theory congenial: +%.3f)\n",
                mean(df$term1), mean(df$tr_RIV)))
    cat(sprintf("  Term 2 = %+.3f  (theory ~ 0)\n", mean(df$term2)))
    cat(sprintf("  Term 3 = %+.3f  (theory: %+.3f)\n",
                mean(df$term3), -mean(df$tr_RIV) / 2))
    cat(sprintf("  Total  = %+.3f  (theory: %+.3f)\n",
                mean(df$total), mean(df$tr_RIV) / 2))
    cat(sprintf("  r_term1 / r_term3 = %+.3f / %+.3f\n",
                mean(df$term1) / mean(df$tr_RIV),
                mean(df$term3) / mean(df$tr_RIV)))
  }
}

total_elapsed <- (proc.time() - overall_t0)[3]
cat(sprintf("\n\n======================================================\n"))
cat(sprintf("  All methods complete.\n"))
cat(sprintf("  Total wall time: %.1fs (%.2fh)\n",
            total_elapsed, total_elapsed / 3600))
cat(sprintf("======================================================\n"))
