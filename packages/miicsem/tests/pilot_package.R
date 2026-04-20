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

# Chi-square summary for M1 (true model) across reps
chi_M1 <- do.call(rbind, lapply(res[[1]], function(r) {
  if (is.null(r) || is.null(r$chi2_df)) return(NULL)
  data.frame(chi2_com = r$chi2_df["M1", "chi2_com"],
             chi2_MI  = r$chi2_df["M1", "chi2_MI"],
             chi2_D3  = r$chi2_df["M1", "chi2_D3"],
             df       = r$chi2_df["M1", "df"])
}))
cat(sprintf("\n=== Chi-squares for M1 (df=%d) across %d reps ===\n",
            chi_M1$df[1], nrow(chi_M1)))
cat(sprintf("  chi2_com: mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_com), sd(chi_M1$chi2_com)))
cat(sprintf("  chi2_MI : mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_MI), sd(chi_M1$chi2_MI)))
cat(sprintf("  chi2_D3 : mean=%.2f, sd=%.2f\n",
            mean(chi_M1$chi2_D3), sd(chi_M1$chi2_D3)))

cat("\n=== Package pilot: PASS ===\n")
