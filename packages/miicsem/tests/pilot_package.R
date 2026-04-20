library(miicsem)

cat("=== Testing miicsem::run_simulation() with 10 reps ===\n\n")

t0 <- proc.time()
res <- run_simulation(
  n_reps       = 10,
  sample_sizes = 250,
  miss_rates   = 0.25,
  verbose      = TRUE
)
elapsed <- (proc.time() - t0)[3]

cat(sprintf("\nElapsed: %.1f sec (%.1f sec/rep)\n", elapsed, elapsed / 10))
cat(sprintf("Conditions: %d\n", length(res)))
cat(sprintf("Reps in first condition: %d (%d non-null)\n",
            length(res[[1]]),
            sum(!sapply(res[[1]], is.null))))

sel_mat <- do.call(rbind, lapply(res[[1]],
  function(r) if (!is.null(r)) r$selections))
cat("\n=== Selection frequency ===\n")
for (m in colnames(sel_mat)) {
  tab <- table(sel_mat[, m])
  cat(sprintf("  %-10s: %s\n", m,
              paste(sprintf("%s=%d", names(tab), as.integer(tab)),
                    collapse = ", ")))
}

cat("\n=== Package pilot: PASS ===\n")
