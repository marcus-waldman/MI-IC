# ============================================================================
# Test Error Handling in compute_rubin_variance
# ============================================================================
# Verify that the function handles NULL vcov matrices gracefully

library(future)
library(future.apply)
library(lavaan)
library(MASS)
library(mvtnorm)

# Source utilities
source("simulations/utils/generate_covariance.R")
source("simulations/utils/generate_data.R")
source("simulations/utils/fit_models.R")
source("simulations/utils/impute_data.R")
source("simulations/utils/compute_metrics.R")

# Generate test data
set.seed(12345)
n <- 100
p <- 10
M <- 10

# Create CS structure
Sigma_true <- generate_covariance(p = p, structure = "CS", sigma2 = 1, rho = 0.5)
mu_true <- rep(0, p)

# Generate complete data and impose missingness
data_complete <- generate_mvn_data(n, mu_true, Sigma_true, seed = 123)
miss_result <- impose_missingness(data_complete, missing_rate = 0.6,
                                   pattern = "monotone", prop_complete = 0.4,
                                   seed = 456)
data_miss <- miss_result$data_miss

cat("Testing error handling in compute_rubin_variance\n\n")

# Generate imputations
cat("Generating", M, "imputations...\n")
impute_result <- impute_wrapper(
  data_miss = data_miss,
  scenario = "MLE",
  structure = "CS",
  M = M,
  seed_base = 789
)

completed_datasets <- impute_result$completed_datasets

# Test 1: All imputations succeed
cat("\nTest 1: Normal case (all imputations should succeed)\n")
riv_result <- compute_rubin_variance(completed_datasets, structure = "CS")
cat(sprintf("  M = %d, M_successful = %d, tr(RIV) = %.3f\n",
            riv_result$M, riv_result$M_successful, riv_result$tr_RIV))
cat(sprintf("  Success rate: %.1f%%\n",
            riv_result$M_successful / riv_result$M * 100))

# Test 2: Simulate some failures by using difficult data
cat("\nTest 2: Challenging case (some imputations may fail)\n")
cat("Creating pathological dataset with high collinearity...\n")

# Create highly collinear data
Sigma_bad <- diag(p)
Sigma_bad[1, 2] <- Sigma_bad[2, 1] <- 0.999  # Nearly perfect correlation
data_bad <- generate_mvn_data(n = 50, mu = mu_true, Sigma = Sigma_bad, seed = 999)
miss_result_bad <- impose_missingness(data_bad, missing_rate = 0.7,
                                      pattern = "monotone", prop_complete = 0.3,
                                      seed = 888)

impute_result_bad <- impute_wrapper(
  data_miss = miss_result_bad$data_miss,
  scenario = "MLE",
  structure = "CS",
  M = M,
  seed_base = 777
)

riv_result_bad <- tryCatch(
  {
    compute_rubin_variance(impute_result_bad$completed_datasets, structure = "CS")
  },
  error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  }
)

if (!is.null(riv_result_bad)) {
  cat(sprintf("  M = %d, M_successful = %d, tr(RIV) = %.3f\n",
              riv_result_bad$M, riv_result_bad$M_successful, riv_result_bad$tr_RIV))
  cat(sprintf("  Success rate: %.1f%%\n",
              riv_result_bad$M_successful / riv_result_bad$M * 100))
}

cat("\n✓ Error handling tests complete\n")
