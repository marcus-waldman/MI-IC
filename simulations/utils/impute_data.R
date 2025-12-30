# ============================================================================
# Multiple Imputation Functions
# ============================================================================
# Functions for generating multiple imputations using conditional MVN
# Supports different scenarios:
#   (a) True parameters (θ̃ = θ*)
#   (b) MLE parameters (θ̃ = θ̂_MLE)
# ============================================================================

#' Impute Missing Values Using Conditional MVN
#'
#' Core imputation function using conditional multivariate normal distribution
#'
#' @param data_miss Numeric matrix with NA values (n x p)
#' @param mu Numeric vector. Mean parameter (length p)
#' @param Sigma Numeric matrix. Covariance parameter (p x p)
#' @param M Integer. Number of imputations to generate
#' @param seed_base Integer. Base seed for reproducibility (optional)
#' @return List of M completed datasets (each n x p matrix)
#' @details
#'   For each imputation m and each observation i with missing values:
#'   - Identify observed and missing indices
#'   - Draw Z_mis | Z_obs ~ N(μ_cond, Σ_cond) using condMVNorm
#'   Seeds: set.seed(seed_base + m) for each imputation m
#' @examples
#'   data_miss <- matrix(rnorm(100*5), 100, 5)
#'   data_miss[sample(500, 100)] <- NA
#'   mu <- rep(0, 5)
#'   Sigma <- diag(5)
#'   imputed_list <- impute_conditional_mvn(data_miss, mu, Sigma, M = 10, seed_base = 123)
impute_conditional_mvn <- function(data_miss, mu, Sigma, M, seed_base = NULL) {

  if (!requireNamespace("condMVNorm", quietly = TRUE)) {
    stop("Package 'condMVNorm' required for conditional MVN imputation")
  }

  # Input validation
  if (!is.matrix(data_miss)) {
    data_miss <- as.matrix(data_miss)
  }

  n <- nrow(data_miss)
  p <- ncol(data_miss)

  if (length(mu) != p || nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("Dimension mismatch between data_miss, mu, and Sigma")
  }

  if (M < 1) stop("M must be at least 1")

  # Identify observations with any missing values
  has_missing <- apply(data_miss, 1, function(row) any(is.na(row)))
  miss_indices <- which(has_missing)

  if (length(miss_indices) == 0) {
    # No missing data; return M copies of original data
    warning("No missing values in data_miss; returning M copies of complete data")
    return(replicate(M, data_miss, simplify = FALSE))
  }

  # Initialize list to store M completed datasets
  completed_datasets <- vector("list", M)

  # Generate M imputations
  for (m in 1:M) {

    # Set seed for this imputation
    if (!is.null(seed_base)) {
      set.seed(seed_base + m)
    }

    # Start with a copy of the incomplete data
    data_complete <- data_miss

    # Impute each observation with missing values
    for (i in miss_indices) {

      # Identify observed and missing indices for this observation
      obs_idx <- which(!is.na(data_miss[i, ]))
      miss_idx <- which(is.na(data_miss[i, ]))

      if (length(miss_idx) == 0) next  # Skip if no missing (shouldn't happen)

      # Extract observed values
      z_obs <- data_miss[i, obs_idx]

      # Compute conditional distribution parameters
      # Z_miss | Z_obs ~ N(μ_cond, Σ_cond)
      cond_params <- condMVNorm::condMVN(
        mean = mu,
        sigma = Sigma,
        dependent.ind = miss_idx,
        given.ind = obs_idx,
        X.given = z_obs
      )

      mu_cond <- cond_params$condMean
      Sigma_cond <- cond_params$condVar

      # Draw from conditional distribution
      # Handle univariate case
      if (length(miss_idx) == 1) {
        z_miss_imputed <- rnorm(1, mean = mu_cond, sd = sqrt(Sigma_cond))
      } else {
        # Multivariate case
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("Package 'MASS' required for multivariate imputation")
        }
        z_miss_imputed <- MASS::mvrnorm(n = 1, mu = mu_cond, Sigma = Sigma_cond)
      }

      # Fill in imputed values
      data_complete[i, miss_idx] <- z_miss_imputed
    }

    # Store completed dataset
    completed_datasets[[m]] <- data_complete
  }

  return(completed_datasets)
}


#' Multiple Imputation Using True Parameters (Scenario a)
#'
#' Imputes using known population parameters θ̃ = θ*
#' Expected imputation bias ≈ 0 (control scenario)
#'
#' @param data_miss Numeric matrix with NA values
#' @param mu_true Numeric vector. True mean parameter
#' @param Sigma_true Numeric matrix. True covariance parameter
#' @param M Integer. Number of imputations
#' @param seed_base Integer. Base seed for reproducibility
#' @return List of M completed datasets
#' @details
#'   This is the control scenario. By iterated expectations, if imputation
#'   uses the true parameters, there should be no imputation bias.
#' @examples
#'   result <- impute_scenario_true(data_miss, mu_true, Sigma_true, M = 100, seed_base = 1000)
impute_scenario_true <- function(data_miss, mu_true, Sigma_true, M, seed_base = NULL) {

  completed_datasets <- impute_conditional_mvn(
    data_miss = data_miss,
    mu = mu_true,
    Sigma = Sigma_true,
    M = M,
    seed_base = seed_base
  )

  return(completed_datasets)
}


#' Multiple Imputation Using MLE Parameters (Scenario b)
#'
#' Fits model to data_miss via lavaan, extracts θ̂_MLE, then imputes
#' using the same θ̂ for all M imputations (maximum correlation)
#' Expected imputation bias ≈ tr(RIV)
#'
#' @param data_miss Numeric matrix with NA values
#' @param structure Character. "CS", "Toeplitz", or "Unstructured"
#' @param M Integer. Number of imputations
#' @param seed_base Integer. Base seed for reproducibility
#' @return List with components:
#'   \item{completed_datasets}{List of M completed datasets}
#'   \item{fit_result}{Result from fit_lavaan_model()}
#'   \item{mu_hat}{Estimated mean vector}
#'   \item{Sigma_hat}{Estimated covariance matrix}
#' @details
#'   This scenario uses plug-in MLE: same θ̂ for all M imputations.
#'   This creates maximum correlation between imputation and estimation,
#'   leading to imputation bias.
#' @examples
#'   result <- impute_scenario_mle(data_miss, structure = "CS", M = 100, seed_base = 2000)
#'   completed_datasets <- result$completed_datasets
impute_scenario_mle <- function(data_miss, structure, M, seed_base = NULL) {

  # Source fit_models.R functions if not already loaded
  if (!exists("fit_lavaan_model", mode = "function")) {
    source("simulations/utils/fit_models.R")
  }

  # Fit model to data_miss to get MLE estimates
  fit_result <- fit_lavaan_model(data_miss, structure = structure, return_vcov = TRUE)

  if (!fit_result$converged) {
    warning("Model fitting did not converge; imputation may be unreliable")
  }

  mu_hat <- fit_result$mu_hat
  Sigma_hat <- fit_result$Sigma_hat

  if (is.null(mu_hat) || is.null(Sigma_hat)) {
    stop("Model fitting failed; cannot proceed with imputation")
  }

  # Impute using MLE estimates (same θ̂ for all M imputations)
  completed_datasets <- impute_conditional_mvn(
    data_miss = data_miss,
    mu = mu_hat,
    Sigma = Sigma_hat,
    M = M,
    seed_base = seed_base
  )

  return(list(
    completed_datasets = completed_datasets,
    fit_result = fit_result,
    mu_hat = mu_hat,
    Sigma_hat = Sigma_hat
  ))
}


#' Wrapper for Multiple Imputation (All Scenarios)
#'
#' Convenience function that handles both scenarios
#'
#' @param data_miss Numeric matrix with NA values
#' @param scenario Character. "True" or "MLE"
#' @param structure Character. "CS", "Toeplitz", or "Unstructured" (for MLE scenario)
#' @param mu_true Numeric vector. True mean (for True scenario)
#' @param Sigma_true Numeric matrix. True covariance (for True scenario)
#' @param M Integer. Number of imputations
#' @param seed_base Integer. Base seed
#' @return List with completed_datasets and scenario-specific components
#' @examples
#'   # Scenario (a): True parameters
#'   result_true <- impute_wrapper(
#'     data_miss, scenario = "True",
#'     mu_true = mu, Sigma_true = Sigma,
#'     M = 100, seed_base = 1000
#'   )
#'
#'   # Scenario (b): MLE
#'   result_mle <- impute_wrapper(
#'     data_miss, scenario = "MLE", structure = "CS",
#'     M = 100, seed_base = 2000
#'   )
impute_wrapper <- function(data_miss, scenario,
                           structure = NULL,
                           mu_true = NULL, Sigma_true = NULL,
                           M, seed_base = NULL) {

  scenario <- match.arg(scenario, c("True", "MLE"))

  if (scenario == "True") {
    if (is.null(mu_true) || is.null(Sigma_true)) {
      stop("Must provide mu_true and Sigma_true for True scenario")
    }

    result <- list(
      completed_datasets = impute_scenario_true(data_miss, mu_true, Sigma_true, M, seed_base),
      scenario = "True",
      mu_imputation = mu_true,
      Sigma_imputation = Sigma_true
    )

  } else if (scenario == "MLE") {
    if (is.null(structure)) {
      stop("Must provide structure for MLE scenario")
    }

    result <- impute_scenario_mle(data_miss, structure, M, seed_base)
    result$scenario <- "MLE"
    result$mu_imputation <- result$mu_hat
    result$Sigma_imputation <- result$Sigma_hat
  }

  return(result)
}


#' Check Imputation Quality
#'
#' Diagnostic function to assess imputation quality
#'
#' @param completed_datasets List of M completed datasets
#' @param data_complete Original complete data (before imposing missingness)
#' @param data_miss Data with missing values
#' @return List with diagnostic statistics
#' @details
#'   Computes:
#'   - Mean imputed values vs true values (for missing locations)
#'   - Variance of imputations across M datasets
#'   - RMSE of imputed vs true values
#' @examples
#'   diagnostics <- check_imputation_quality(completed_datasets, data_complete, data_miss)
check_imputation_quality <- function(completed_datasets, data_complete, data_miss) {

  M <- length(completed_datasets)
  n <- nrow(data_miss)
  p <- ncol(data_miss)

  # Identify missing locations
  is_missing <- is.na(data_miss)

  # Extract imputed values across M datasets
  # For each missing location, get M imputed values
  imputed_values_list <- lapply(1:M, function(m) {
    completed_datasets[[m]][is_missing]
  })

  # Stack into matrix: rows = missing locations, columns = M imputations
  imputed_matrix <- do.call(cbind, imputed_values_list)

  # True values at missing locations
  true_values <- data_complete[is_missing]

  # Mean imputed value per location
  mean_imputed <- rowMeans(imputed_matrix)

  # Bias: mean imputed - true
  bias_per_location <- mean_imputed - true_values
  mean_bias <- mean(bias_per_location)

  # RMSE
  rmse <- sqrt(mean((mean_imputed - true_values)^2))

  # Variance of imputations
  var_imputed <- apply(imputed_matrix, 1, var)
  mean_var <- mean(var_imputed)

  return(list(
    mean_bias = mean_bias,
    rmse = rmse,
    mean_variance_imputation = mean_var,
    n_missing_locations = sum(is_missing)
  ))
}
