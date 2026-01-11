# Test with SIMPLE random amputation (like debug script)
rm(list = ls())

source("R/01_generate_data.R")
source("R/06_compute_metrics.R")

set.seed(810366)

# Generate data (same as debug script)
N <- 500
X <- rnorm(N)
M <- 0.5 * X + rnorm(N)
Y <- 0.5 * M + rnorm(N)
data_complete <- data.frame(X=X, M=M, Y=Y)

# Simple random amputation (NOT mice::ampute)
prop_miss <- 0.3
data_miss <- data_complete
miss_M <- sample(1:N, round(N * prop_miss))
miss_Y <- sample(1:N, round(N * prop_miss))
data_miss$M[miss_M] <- NA
data_miss$Y[miss_Y] <- NA

cat("Missing patterns:\n")
cat("  Complete:", sum(!is.na(data_miss$M) & !is.na(data_miss$Y)), "\n")
cat("  M miss:", sum(is.na(data_miss$M) & !is.na(data_miss$Y)), "\n")
cat("  Y miss:", sum(!is.na(data_miss$M) & is.na(data_miss$Y)), "\n")
cat("  Both miss:", sum(is.na(data_miss$M) & is.na(data_miss$Y)), "\n\n")

# Use analytical functions
riv_result <- compute_tr_RIV_analytical(data_complete, data_miss)
theta_obs <- riv_result$theta_obs
theta_com <- riv_result$theta_com

Q_at_obs <- compute_Q_analytical(data_miss, theta_obs, theta_obs)
ell_com_at_obs <- -negloglik_complete(
  c(theta_obs[1], theta_obs[2], log(theta_obs[3]), log(theta_obs[4])),
  data_complete
)
ell_com_at_com <- riv_result$ell_com

term1 <- Q_at_obs - ell_com_at_obs
term2 <- ell_com_at_obs - ell_com_at_com

cat("Results:\n")
cat(sprintf("tr(RIV): %.4f\n", riv_result$tr_RIV))
cat(sprintf("Term 1: %.4f (expected ≈ %.4f)\n", term1, riv_result$tr_RIV))
cat(sprintf("Term 2: %.4f (expected ≈ %.4f)\n", term2, -0.5 * riv_result$tr_RIV))
cat(sprintf("Ratio 1: %.4f (should ≈ 1.0)\n", term1 / riv_result$tr_RIV))
cat(sprintf("Ratio 2: %.4f (should ≈ 1.0)\n", term2 / (-0.5 * riv_result$tr_RIV)))
