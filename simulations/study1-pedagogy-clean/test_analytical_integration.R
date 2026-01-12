# ============================================================================
# Validate Analytical Integration Against Monte Carlo Integration
# ============================================================================
# This script validates that our analytical delta method integration matches
# Monte Carlo integration from brms posterior samples.
#
# Approach:
# 1. Simulate one large dataset (n = 5000) with missingness
# 2. Fit imputation models in brms (A, B, C)
# 3. Extract M posterior samples of imputation parameters
# 4. Compute Q via Monte Carlo: Q_MC = (1/M) Σ_m Q(θ̂_obs | φ̃⁽ᵐ⁾)
# 5. Compute Q via analytical delta method: Q_analytical = Q(θ̂_obs | φ̂) + (1/2) tr(H_Q · V_φ)
# 6. Verify: |Q_MC - Q_analytical| < 2 × se(Q_MC)

library(brms)
library(posterior)

# Set working directory to script location
setwd("C:/Users/marcu/git-repositories/MI-IC/simulations/study1-pedagogy-clean")

# Source simulation functions
source("R/00_seed_management.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/03_analytical_functions.R")
source("R/04_compute_metrics.R")
source("R/05_run_simulation.R")

# ============================================================================
# Configuration
# ============================================================================

set.seed(20250111)  # Fixed seed for reproducibility

# True parameters
true_params <- list(a = 0.5, b = 0.5, sigma2_M = 1, sigma2_Y = 1)

# Sample size
n <- 200

# Moderate missingness (MCAR, ~25% missing)
mechanism <- "MCAR"
gamma <- 0.5
cut1 <- 0.5
cut2 <- 0.75

cat("=== Analytical Integration Validation ===\n\n")
cat("Configuration:\n")
cat(sprintf("  Sample size: %d\n", n))
cat(sprintf("  Mechanism: %s\n", mechanism))
cat(sprintf("  True a: %.2f, b: %.2f, σ²_M: %.2f, σ²_Y: %.2f\n\n",
            true_params$a, true_params$b, true_params$sigma2_M, true_params$sigma2_Y))

# ============================================================================
# Step 1: Generate Data with Missingness
# ============================================================================

cat("Step 1: Generating data...\n")

X <- rnorm(n)
M <- true_params$a * X + rnorm(n, 0, sqrt(true_params$sigma2_M))
Y <- true_params$b * M + rnorm(n, 0, sqrt(true_params$sigma2_Y))
data_com <- data.frame(X = X, M = M, Y = Y)

# Ampute data (monotone missingness)
data_mis <- ampute_data(data_com, mechanism = mechanism, gamma = gamma,
                        cut1 = cut1, cut2 = cut2, seed = 12345)

# Count missing patterns
n_complete <- sum(!is.na(data_mis$M) & !is.na(data_mis$Y))
n_Y_miss <- sum(!is.na(data_mis$M) & is.na(data_mis$Y))
n_MY_miss <- sum(is.na(data_mis$M) & is.na(data_mis$Y))

cat(sprintf("  Pattern 1 (complete): %d (%.1f%%)\n", n_complete, 100 * n_complete / n))
cat(sprintf("  Pattern 2 (Y missing): %d (%.1f%%)\n", n_Y_miss, 100 * n_Y_miss / n))
cat(sprintf("  Pattern 3 (M,Y missing): %d (%.1f%%)\n\n", n_MY_miss, 100 * n_MY_miss / n))

# ============================================================================
# Step 2: Fit Observed-Data Model (FIML)
# ============================================================================

cat("Step 2: Fitting observed-data model (FIML)...\n")

theta_start <- c(0.5, 0.5, log(1), log(1))
fit_o <- optim(theta_start, negloglik_observed, data = data_mis, method = "BFGS", hessian = TRUE)

theta_obs <- fit_o$par
theta_obs[3:4] <- exp(theta_obs[3:4])

cat(sprintf("  θ̂_obs = (a: %.4f, b: %.4f, σ²_M: %.4f, σ²_Y: %.4f)\n\n",
            theta_obs[1], theta_obs[2], theta_obs[3], theta_obs[4]))

# Compute posterior covariance
V_theta <- solve(fit_o$hessian)
J <- diag(c(1, 1, theta_obs[3], theta_obs[4]))
V_theta_orig <- J %*% V_theta %*% t(J)

# ============================================================================
# Validation Function
# ============================================================================

validate_model <- function(model_name, brm_formula, model_code) {

  cat(sprintf("\n=== MODEL %s ===\n", model_name))

  # --------------------------------------------------------------------------
  # Step 3: Fit Imputation Model in brms
  # --------------------------------------------------------------------------

  cat("  Fitting imputation model in brms...\n")

  fit_brms <- brm(
    brm_formula,
    data = data_mis,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 20250111,
    silent = 2,
    refresh = 0
  )

  cat("    Done.\n")

  # --------------------------------------------------------------------------
  # Step 4: Extract Posterior Samples
  # --------------------------------------------------------------------------

  cat("  Extracting posterior samples...\n")

  post_samples <- as_draws_df(fit_brms)
  M <- nrow(post_samples)

  cat(sprintf("    M = %d posterior samples\n", M))

  # Extract imputation parameters for each model
  if (model_code == "A") {
    # Model A: (a, b, sigma2_M, sigma2_Y)
    theta_impute_samples <- lapply(1:M, function(i) {
      list(
        a = post_samples$b_M_X[i],
        b = post_samples$b_Y_M[i],
        sigma2_M = post_samples$sigma_M[i]^2,
        sigma2_Y = post_samples$sigma_Y[i]^2
      )
    })

  } else if (model_code == "B") {
    # Model B: (a, sigma2_M, mu_Y, tau2_Y)
    theta_impute_samples <- lapply(1:M, function(i) {
      list(
        a = post_samples$b_M_X[i],
        sigma2_M = post_samples$sigma_M[i]^2,
        mu_Y = post_samples$b_Y_Intercept[i],
        tau2_Y = post_samples$sigma_Y[i]^2
      )
    })

  } else if (model_code == "C") {
    # Model C: Extract from multivariate brms model
    # Check what parameter names are available
    param_names <- names(post_samples)

    # Try to find correlation parameters with different naming conventions
    has_rescor_double <- any(grepl("rescor__", param_names))
    has_rescor_single <- any(grepl("rescor_[^_]", param_names))

    theta_impute_samples <- lapply(1:M, function(i) {
      # Extract means (intercepts in brms multivariate models are prefixed with b_)
      mu_X <- post_samples[[paste0("b_X_Intercept")]][i]
      mu_M <- post_samples[[paste0("b_M_Intercept")]][i]
      mu_Y <- post_samples[[paste0("b_Y_Intercept")]][i]

      # Extract standard deviations
      sigma_X <- post_samples$sigma_X[i]
      sigma_M <- post_samples$sigma_M[i]
      sigma_Y <- post_samples$sigma_Y[i]

      # Extract correlations (try different naming conventions)
      if (has_rescor_double) {
        # Double underscore format: rescor__X__M
        cor_XM <- post_samples[[paste0("rescor__X__M")]][i]
        cor_XY <- post_samples[[paste0("rescor__X__Y")]][i]
        cor_MY <- post_samples[[paste0("rescor__M__Y")]][i]
      } else if (has_rescor_single) {
        # Single underscore format: rescor_X_M
        cor_XM <- post_samples[[paste0("rescor_X_M")]][i]
        cor_XY <- post_samples[[paste0("rescor_X_Y")]][i]
        cor_MY <- post_samples[[paste0("rescor_M_Y")]][i]
      } else {
        # Fallback: compute from empirical data
        cor_XM <- cov(data_mis$X, data_mis$M, use = "pairwise.complete.obs") /
                  (sd(data_mis$X, na.rm = TRUE) * sd(data_mis$M, na.rm = TRUE))
        cor_XY <- cov(data_mis$X, data_mis$Y, use = "pairwise.complete.obs") /
                  (sd(data_mis$X, na.rm = TRUE) * sd(data_mis$Y, na.rm = TRUE))
        cor_MY <- cov(data_mis$M, data_mis$Y, use = "pairwise.complete.obs") /
                  (sd(data_mis$M, na.rm = TRUE) * sd(data_mis$Y, na.rm = TRUE))
      }

      # Construct covariance matrix
      Sigma <- matrix(0, 3, 3)
      Sigma[1, 1] <- sigma_X^2
      Sigma[2, 2] <- sigma_M^2
      Sigma[3, 3] <- sigma_Y^2
      Sigma[1, 2] <- Sigma[2, 1] <- cor_XM * sigma_X * sigma_M
      Sigma[1, 3] <- Sigma[3, 1] <- cor_XY * sigma_X * sigma_Y
      Sigma[2, 3] <- Sigma[3, 2] <- cor_MY * sigma_M * sigma_Y

      list(mu_X = mu_X, mu_M = mu_M, mu_Y = mu_Y, Sigma = Sigma)
    })
  }

  # --------------------------------------------------------------------------
  # Step 5: Monte Carlo Integration
  # --------------------------------------------------------------------------

  cat("  Computing Q via Monte Carlo integration...\n")

  Q_MC_samples <- numeric(M)

  for (i in 1:M) {
    Q_MC_samples[i] <- compute_Q_analytical(
      data_mis, theta_obs, theta_impute_samples[[i]], model = model_code
    )
  }

  Q_MC <- mean(Q_MC_samples)
  Q_MC_se <- sd(Q_MC_samples) / sqrt(M)

  cat(sprintf("    Q_MC = %.6f (SE = %.6f)\n", Q_MC, Q_MC_se))

  # --------------------------------------------------------------------------
  # Step 6: Analytical Integration (Delta Method)
  # --------------------------------------------------------------------------

  cat("  Computing Q via analytical delta method...\n")

  # Posterior mean of imputation parameters
  if (model_code == "A") {
    theta_impute_mean <- list(
      a = mean(sapply(theta_impute_samples, function(x) x$a)),
      b = mean(sapply(theta_impute_samples, function(x) x$b)),
      sigma2_M = mean(sapply(theta_impute_samples, function(x) x$sigma2_M)),
      sigma2_Y = mean(sapply(theta_impute_samples, function(x) x$sigma2_Y))
    )

  } else if (model_code == "B") {
    theta_impute_mean <- list(
      a = mean(sapply(theta_impute_samples, function(x) x$a)),
      sigma2_M = mean(sapply(theta_impute_samples, function(x) x$sigma2_M)),
      mu_Y = mean(sapply(theta_impute_samples, function(x) x$mu_Y)),
      tau2_Y = mean(sapply(theta_impute_samples, function(x) x$tau2_Y))
    )

  } else if (model_code == "C") {
    theta_impute_mean <- list(
      mu_X = mean(sapply(theta_impute_samples, function(x) x$mu_X)),
      mu_M = mean(sapply(theta_impute_samples, function(x) x$mu_M)),
      mu_Y = mean(sapply(theta_impute_samples, function(x) x$mu_Y)),
      Sigma = Reduce("+", lapply(theta_impute_samples, function(x) x$Sigma)) / M
    )
  }

  # Estimate V_theta from brms posterior
  V_theta_brms <- estimate_V_theta(data_mis, theta_impute_mean, V_theta_orig, model = model_code)

  # Compute analytical Q
  Q_results <- compute_Q_analytical_proper(
    data_mis, theta_obs, theta_impute_mean, V_theta_brms, model = model_code
  )

  Q_analytical <- Q_results$Q_proper

  cat(sprintf("    Q_analytical = %.6f (Hessian correction = %.6f)\n",
              Q_analytical, Q_results$hessian_correction))

  # --------------------------------------------------------------------------
  # Step 7: Compare
  # --------------------------------------------------------------------------

  diff <- abs(Q_MC - Q_analytical)
  z_score <- diff / Q_MC_se
  threshold <- 2 * Q_MC_se

  cat("\n  === VALIDATION RESULT ===\n")
  cat(sprintf("    |Q_MC - Q_analytical| = %.6f\n", diff))
  cat(sprintf("    Monte Carlo SE = %.6f\n", Q_MC_se))
  cat(sprintf("    z-score = %.2f\n", z_score))
  cat(sprintf("    Threshold (2×SE) = %.6f\n", threshold))

  if (diff < threshold) {
    cat(sprintf("    MODEL %s: PASS ✓\n", model_name))
    pass <- TRUE
  } else {
    cat(sprintf("    MODEL %s: FAIL ✗\n", model_name))
    pass <- FALSE
  }

  # Return results for summary table
  data.frame(
    Model = model_name,
    Q_MC = Q_MC,
    Q_analytical = Q_analytical,
    Difference = diff,
    MC_SE = Q_MC_se,
    z_score = z_score,
    Pass = pass
  )
}

# ============================================================================
# Model A: Congenial (Matches DGP)
# ============================================================================

results_A <- validate_model(
  model_name = "A (Congenial)",
  brm_formula = bf(M ~ X) + bf(Y ~ M) + set_rescor(FALSE),
  model_code = "A"
)

# ============================================================================
# Model B: Uncongenial (Ignores M→Y)
# ============================================================================

results_B <- validate_model(
  model_name = "B (Uncongenial)",
  brm_formula = bf(M ~ X) + bf(Y ~ 1) + set_rescor(FALSE),
  model_code = "B"
)

# ============================================================================
# Model C: Saturated MVN
# ============================================================================

results_C <- validate_model(
  model_name = "C (Saturated MVN)",
  brm_formula = mvbf(X ~ 1, M ~ 1, Y ~ 1),
  model_code = "C"
)

# ============================================================================
# Summary Table
# ============================================================================

cat("\n\n=== SUMMARY TABLE ===\n\n")

results_all <- rbind(results_A, results_B, results_C)
print(results_all, row.names = FALSE, digits = 6)

cat("\n")

# Overall pass/fail
if (all(results_all$Pass)) {
  cat("ALL MODELS PASSED ✓\n")
  cat("Analytical integration is validated against Monte Carlo integration.\n")
} else {
  cat("SOME MODELS FAILED ✗\n")
  cat("Check model implementation and Hessian corrections.\n")
}
