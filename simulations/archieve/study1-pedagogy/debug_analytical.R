# Direct comparison: debug script vs new analytical functions
rm(list = ls())

# Load new functions
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/06_compute_metrics.R")

set.seed(810366)  # Same seed as test_analytical first rep

# Parameters
a_true <- 0.5
b_true <- 0.5
sigma2_M_true <- 1
sigma2_Y_true <- 1
N <- 500
prop_miss <- 0.3

# Generate data
X <- rnorm(N)
M <- a_true * X + rnorm(N, 0, sqrt(sigma2_M_true))
Y <- b_true * M + rnorm(N, 0, sqrt(sigma2_Y_true))
data_complete <- data.frame(X = X, M = M, Y = Y)

# Ampute
data_miss <- data_complete
miss_M <- sample(1:N, round(N * prop_miss))
miss_Y <- sample(1:N, round(N * prop_miss))
data_miss$M[miss_M] <- NA
data_miss$Y[miss_Y] <- NA

cat("=== Using NEW analytical functions ===\n")

# Use new function
riv_result <- compute_tr_RIV_analytical(data_complete, data_miss)
cat("tr(RIV):", round(riv_result$tr_RIV, 4), "\n")

theta_obs <- riv_result$theta_obs
theta_com <- riv_result$theta_com

# Compute Q
Q_at_obs <- compute_Q_analytical(data_miss, theta_obs, theta_obs)
cat("Q(θ_obs|θ_obs):", round(Q_at_obs, 4), "\n")

# Compute ell_com at theta_obs
ell_com_at_obs <- -negloglik_complete(
  c(theta_obs[1], theta_obs[2], log(theta_obs[3]), log(theta_obs[4])),
  data_complete
)
cat("ℓ_com(θ_obs):", round(ell_com_at_obs, 4), "\n")

# Compute ell_com at theta_com
ell_com_at_com <- riv_result$ell_com
cat("ℓ_com(θ_com):", round(ell_com_at_com, 4), "\n\n")

term1 <- Q_at_obs - ell_com_at_obs
term2 <- ell_com_at_obs - ell_com_at_com

cat("Term 1:", round(term1, 4), "(expected ≈", round(riv_result$tr_RIV, 4), ")\n")
cat("Term 2:", round(term2, 4), "(expected ≈", round(-0.5 * riv_result$tr_RIV, 4), ")\n")
cat("Ratio 1:", round(term1 / riv_result$tr_RIV, 4), "\n")
cat("Ratio 2:", round(term2 / (-0.5 * riv_result$tr_RIV), 4), "\n")
