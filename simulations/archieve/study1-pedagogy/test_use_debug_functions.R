# Test using debug script functions directly
rm(list = ls())

# Load study1 data generation
source("R/01_generate_data.R")

# Load DEBUG script functions (overriding study1 versions)
source("C:/Users/marcu/git-repositories/MI-IC/claude/derivations/debug_custom_likelihood.R")

# Just keep the functions, remove the execution
# (the script runs code at the end, but we only want the functions)

set.seed(810366)

# Generate data using study1 method
N <- 500
X <- rnorm(N)
M <- 0.5 * X + rnorm(N)
Y <- 0.5 * M + rnorm(N)
data_complete <- data.frame(X=X, M=M, Y=Y)

# Simple amputation
prop_miss <- 0.3
data_miss <- data_complete
miss_M <- sample(1:N, round(N * prop_miss))
miss_Y <- sample(1:N, round(N * prop_miss))
data_miss$M[miss_M] <- NA
data_miss$Y[miss_Y] <- NA

cat("Using DEBUG script functions:\n")

# Fit models
theta_start <- c(0.5, 0.5, log(1), log(1))

fit_com <- optim(theta_start, negloglik_complete, data = data_complete,
                 method = "BFGS", hessian = TRUE)
fit_obs <- optim(theta_start, negloglik_observed, data = data_miss,
                 method = "BFGS", hessian = TRUE)

theta_com <- fit_com$par
theta_com[3:4] <- exp(theta_com[3:4])

theta_obs <- fit_obs$par
theta_obs[3:4] <- exp(theta_obs[3:4])

# Compute tr(RIV)
tr_RIV <- sum(diag(fit_com$hessian %*% solve(fit_obs$hessian))) - 4

# Compute Q and likelihoods
Q_at_obs <- compute_Q_analytical(data_miss, theta_obs, theta_obs)
ell_com_at_com <- -fit_com$value
ell_com_at_obs <- -negloglik_complete(c(theta_obs[1:2], log(theta_obs[3:4])), data_complete)

term1 <- Q_at_obs - ell_com_at_obs
term2 <- ell_com_at_obs - ell_com_at_com

cat(sprintf("tr(RIV): %.4f\n", tr_RIV))
cat(sprintf("Term 1: %.4f (expected ≈ %.4f)\n", term1, tr_RIV))
cat(sprintf("Term 2: %.4f (expected ≈ %.4f)\n", term2, -0.5 * tr_RIV))
cat(sprintf("Ratio 1: %.4f (should ≈ 1.0)\n", term1 / tr_RIV))
cat(sprintf("Ratio 2: %.4f (should ≈ 1.0)\n", term2 / (-0.5 * tr_RIV)))
