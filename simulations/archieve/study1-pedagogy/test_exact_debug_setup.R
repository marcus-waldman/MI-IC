# Test with EXACT same setup as debug script
rm(list = ls())

source("R/06_compute_metrics.R")

# Run 50 reps with exact debug script setup
set.seed(12345)

n_reps <- 50
results <- matrix(NA, n_reps, 4)
colnames(results) <- c("term1", "term2", "total", "tr_RIV")

for (rep in 1:n_reps) {
  # Generate data (EXACTLY like debug script)
  a_true <- 0.5
  b_true <- 0.5
  sigma2_M_true <- 1
  sigma2_Y_true <- 1
  N <- 500
  prop_miss <- 0.3

  X <- rnorm(N)
  M <- a_true * X + rnorm(N, 0, sqrt(sigma2_M_true))
  Y <- b_true * M + rnorm(N, 0, sqrt(sigma2_Y_true))
  data_com <- data.frame(X = X, M = M, Y = Y)

  # Ampute (EXACTLY like debug script)
  data_mis <- data_com
  miss_M <- sample(1:N, round(N * prop_miss))
  miss_Y <- sample(1:N, round(N * prop_miss))
  data_mis$M[miss_M] <- NA
  data_mis$Y[miss_Y] <- NA

  # Compute using study1 functions
  riv_result <- compute_tr_RIV_analytical(data_com, data_mis)

  theta_obs <- riv_result$theta_obs
  theta_com <- riv_result$theta_com

  Q_at_obs <- compute_Q_analytical(data_mis, theta_obs, theta_obs)
  ell_com_at_com <- riv_result$ell_com
  ell_com_at_obs <- -negloglik_complete(c(theta_obs[1:2], log(theta_obs[3:4])), data_com)

  term1 <- Q_at_obs - ell_com_at_obs
  term2 <- ell_com_at_obs - ell_com_at_com
  total <- term1 + term2

  results[rep, ] <- c(term1, term2, total, riv_result$tr_RIV)

  if (rep %% 10 == 0) cat("Rep", rep, "\n")
}

cat("\n=== Results (Exact Debug Setup) ===\n")
cat("Mean Term 1:", round(mean(results[,"term1"]), 4), "\n")
cat("Mean tr(RIV):", round(mean(results[,"tr_RIV"]), 4), "\n")
cat("Ratio 1:", round(mean(results[,"term1"]) / mean(results[,"tr_RIV"]), 4), "\n\n")

cat("Mean Term 2:", round(mean(results[,"term2"]), 4), "\n")
cat("Expected:", round(-0.5 * mean(results[,"tr_RIV"]), 4), "\n")
cat("Ratio 2:", round(mean(results[,"term2"]) / (-0.5 * mean(results[,"tr_RIV"])), 4), "\n\n")

cat("Mean Total:", round(mean(results[,"total"]), 4), "\n")
cat("Expected:", round(0.5 * mean(results[,"tr_RIV"]), 4), "\n")
cat("Ratio Total:", round(mean(results[,"total"]) / (0.5 * mean(results[,"tr_RIV"])), 4), "\n")
