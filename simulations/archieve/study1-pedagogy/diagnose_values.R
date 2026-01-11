# Diagnose what values we're getting
rm(list = ls())

source("R/01_generate_data.R")
source("R/06_compute_metrics.R")

set.seed(810366)

N <- 500
X <- rnorm(N)
M <- 0.5 * X + rnorm(N)
Y <- 0.5 * M + rnorm(N)
data_complete <- data.frame(X=X, M=M, Y=Y)

prop_miss <- 0.3
data_miss <- data_complete
miss_M <- sample(1:N, round(N * prop_miss))
miss_Y <- sample(1:N, round(N * prop_miss))
data_miss$M[miss_M] <- NA
data_miss$Y[miss_Y] <- NA

riv_result <- compute_tr_RIV_analytical(data_complete, data_miss)

cat("=== Parameter Estimates ===\n")
cat("theta_obs:", round(riv_result$theta_obs, 4), "\n")
cat("theta_com:", round(riv_result$theta_com, 4), "\n\n")

cat("=== Log-Likelihoods ===\n")
cat("ℓ_com(θ_com):", round(riv_result$ell_com, 4), "\n")
cat("ℓ_obs(θ_obs):", round(riv_result$ell_obs, 4), "\n\n")

# Compute Q
Q_at_obs <- compute_Q_analytical(data_miss, riv_result$theta_obs, riv_result$theta_obs)
cat("Q(θ_obs | θ_obs):", round(Q_at_obs, 4), "\n\n")

# Compute ℓ_com at θ_obs
theta_obs_log <- c(riv_result$theta_obs[1:2], log(riv_result$theta_obs[3:4]))
nll_com_at_obs <- negloglik_complete(theta_obs_log, data_complete)
ell_com_at_obs <- -nll_com_at_obs

cat("ℓ_com(θ_obs):", round(ell_com_at_obs, 4), "\n\n")

cat("=== Bias Decomposition ===\n")
term1 <- Q_at_obs - ell_com_at_obs
term2 <- ell_com_at_obs - riv_result$ell_com

cat("Term 1 [Q - ℓ_com(θ_obs)]:", round(term1, 4), "\n")
cat("  Expected (tr(RIV)):", round(riv_result$tr_RIV, 4), "\n")
cat("  Ratio:", round(term1 / riv_result$tr_RIV, 4), "\n\n")

cat("Term 2 [ℓ_com(θ_obs) - ℓ_com(θ_com)]:", round(term2, 4), "\n")
cat("  Expected (-0.5*tr):", round(-0.5 * riv_result$tr_RIV, 4), "\n")
cat("  Ratio:", round(term2 / (-0.5 * riv_result$tr_RIV), 4), "\n\n")

# Double-check: manually compute ℓ_com at θ_com
nll_com_at_com_manual <- negloglik_complete(
  c(riv_result$theta_com[1:2], log(riv_result$theta_com[3:4])),
  data_complete
)
ell_com_at_com_manual <- -nll_com_at_com_manual

cat("=== Verification ===\n")
cat("ℓ_com(θ_com) from fit:", round(riv_result$ell_com, 4), "\n")
cat("ℓ_com(θ_com) manual:", round(ell_com_at_com_manual, 4), "\n")
cat("Difference:", round(abs(riv_result$ell_com - ell_com_at_com_manual), 6), "\n")
