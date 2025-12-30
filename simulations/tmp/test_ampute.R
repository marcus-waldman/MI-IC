# ============================================================================
# Test mice::ampute Implementation for Monotone Missingness
# ============================================================================
# Quick test to verify the refactored impose_missingness() function

library(mice)
library(MASS)

# Source utilities
source("simulations/utils/generate_data.R")

# Generate test data
set.seed(12345)
n <- 100
p <- 10
mu <- rep(0, p)
Sigma <- diag(p)

data_complete <- generate_mvn_data(n, mu, Sigma, seed = 123)

# Test monotone missingness with mice::ampute
cat("Testing monotone missingness with mice::ampute\n")
cat("n =", n, ", p =", p, "\n\n")

result <- impose_missingness(
  data = data_complete,
  missing_rate = 0.6,
  pattern = "monotone",
  prop_complete = 0.4,
  seed = 456
)

# Print missingness summary
cat("Missingness by variable:\n")
miss_df <- data.frame(
  Variable = paste0("Y", 1:p),
  N_Missing = result$n_miss_by_var,
  Prop_Missing = round(result$prop_miss_by_var, 3),
  N_Observed = n - result$n_miss_by_var,
  Prop_Observed = round(1 - result$prop_miss_by_var, 3)
)
print(miss_df)

cat("\n\nExpected pattern:\n")
cat("- Y1: Always observed (100% observed)\n")
cat("- Y10: ~40-50% observed (much better than old ~40% × 1/9 = 4.4%)\n")
cat("- Overall: Monotone pattern with 40% complete observations\n")

# Verify monotonicity
cat("\n\nChecking monotone property:\n")
obs_pattern <- result$missing_pattern
monotone_violations <- 0

for (i in 1:n) {
  row_obs <- obs_pattern[i, ]
  # Find first missing value
  first_miss <- which(!row_obs)[1]
  if (!is.na(first_miss) && first_miss < p) {
    # All subsequent values should be missing
    if (any(row_obs[(first_miss+1):p])) {
      monotone_violations <- monotone_violations + 1
    }
  }
}

cat("Monotone violations:", monotone_violations, "out of", n, "observations\n")
if (monotone_violations == 0) {
  cat("✓ Perfect monotone pattern!\n")
} else {
  cat("✗ Violations detected\n")
}

# Check proportion of complete cases
n_complete <- sum(complete.cases(result$data_miss))
prop_complete_actual <- n_complete / n
cat("\n\nComplete cases:\n")
cat("Expected: ~40% (", 0.4 * n, "observations)\n")
cat("Actual:", n_complete, "observations (", round(prop_complete_actual, 3), ")\n")
