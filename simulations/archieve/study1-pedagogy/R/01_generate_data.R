# ============================================================================
# Generate Mediation Data for Study 1
# ============================================================================
# Functions for generating X -> M -> Y mediation data with specified parameters
# DGM: Full mediation (no direct X -> Y effect by default)
# ============================================================================

#' Generate Mediation Data
#'
#' Generates complete data from a mediation model: X -> M -> Y
#'
#' @param n Integer. Sample size
#' @param beta_xm Numeric. Effect of X on M (default: 0.5)
#' @param beta_my Numeric. Effect of M on Y (default: 0.5)
#' @param beta_xy Numeric. Direct effect of X on Y (default: 0, full mediation)
#' @param sigma2_m Numeric. Residual variance for M (default: 1)
#' @param sigma2_y Numeric. Residual variance for Y (default: 1)
#' @param var_x Numeric. Variance of X (default: 1)
#' @param seed Integer. Random seed for reproducibility (optional)
#' @return data.frame with columns X, M, Y
#' @details
#'   Data generating model:
#'   - X ~ N(0, var_x)
#'   - M = beta_xm * X + epsilon_M, where epsilon_M ~ N(0, sigma2_m)
#'   - Y = beta_my * M + beta_xy * X + epsilon_Y, where epsilon_Y ~ N(0, sigma2_y)
#'
#' @examples
#'   data <- generate_mediation_data(n = 100, seed = 123)
#'   # Full mediation: indirect effect = 0.5 * 0.5 = 0.25
generate_mediation_data <- function(n,
                                     beta_xm = 0.5,
                                     beta_my = 0.5,
                                     beta_xy = 0,
                                     sigma2_m = 1,
                                     sigma2_y = 1,
                                     var_x = 1,
                                     seed = NULL) {

  # Input validation
 if (n < 1) stop("Sample size n must be at least 1")
  if (sigma2_m <= 0) stop("sigma2_m must be positive")
  if (sigma2_y <= 0) stop("sigma2_y must be positive")
  if (var_x <= 0) stop("var_x must be positive")

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Generate X
  X <- rnorm(n, mean = 0, sd = sqrt(var_x))

  # Generate M = beta_xm * X + epsilon_M
  epsilon_M <- rnorm(n, mean = 0, sd = sqrt(sigma2_m))
  M <- beta_xm * X + epsilon_M

  # Generate Y = beta_my * M + beta_xy * X + epsilon_Y
  epsilon_Y <- rnorm(n, mean = 0, sd = sqrt(sigma2_y))
  Y <- beta_my * M + beta_xy * X + epsilon_Y

  # Return as data.frame
  data <- data.frame(X = X, M = M, Y = Y)

  return(data)
}


#' Generate Mediation Data from Config
#'
#' Convenience wrapper that uses config parameters
#'
#' @param config List. Configuration from config.R
#' @param seed Integer. Random seed
#' @return data.frame with columns X, M, Y
generate_mediation_data_from_config <- function(config, seed = NULL) {
  generate_mediation_data(
    n = config$n,
    beta_xm = config$beta_xm,
    beta_my = config$beta_my,
    beta_xy = config$beta_xy,
    sigma2_m = config$sigma2_m,
    sigma2_y = config$sigma2_y,
    var_x = config$var_x,
    seed = seed
  )
}


#' Derive Population Covariance Matrix for Mediation Model
#'
#' Computes the theoretical covariance matrix Sigma for (X, M, Y)
#'
#' @param beta_xm Numeric. Effect of X on M
#' @param beta_my Numeric. Effect of M on Y
#' @param beta_xy Numeric. Direct effect of X on Y
#' @param sigma2_m Numeric. Residual variance for M
#' @param sigma2_y Numeric. Residual variance for Y
#' @param var_x Numeric. Variance of X (default: 1)
#' @return 3x3 covariance matrix with row/col names c("X", "M", "Y")
#' @details
#'   Derivation:
#'   - Var(X) = var_x
#'   - Var(M) = beta_xm^2 * var_x + sigma2_m
#'   - Var(Y) = beta_my^2 * Var(M) + beta_xy^2 * var_x + 2*beta_my*beta_xy*Cov(M,X) + sigma2_y
#'   - Cov(X, M) = beta_xm * var_x
#'   - Cov(X, Y) = beta_my * Cov(X, M) + beta_xy * var_x
#'   - Cov(M, Y) = beta_my * Var(M) + beta_xy * Cov(X, M)
derive_mediation_covariance <- function(beta_xm, beta_my, beta_xy,
                                         sigma2_m, sigma2_y, var_x = 1) {

  # Variance of X
  var_X <- var_x

  # Variance of M = beta_xm^2 * Var(X) + sigma2_m
  var_M <- beta_xm^2 * var_X + sigma2_m

  # Covariance of X and M
  cov_XM <- beta_xm * var_X

  # Covariance of X and Y = beta_my * Cov(X,M) + beta_xy * Var(X)
  cov_XY <- beta_my * cov_XM + beta_xy * var_X

  # Covariance of M and Y = beta_my * Var(M) + beta_xy * Cov(X,M)
  cov_MY <- beta_my * var_M + beta_xy * cov_XM

  # Variance of Y
  var_Y <- beta_my^2 * var_M + beta_xy^2 * var_X + 2 * beta_my * beta_xy * cov_XM + sigma2_y

  # Construct covariance matrix
  Sigma <- matrix(c(
    var_X,  cov_XM, cov_XY,
    cov_XM, var_M,  cov_MY,
    cov_XY, cov_MY, var_Y
  ), nrow = 3, ncol = 3, byrow = TRUE)

  rownames(Sigma) <- colnames(Sigma) <- c("X", "M", "Y")

  return(Sigma)
}


#' Derive Population Covariance from Config
#'
#' @param config List. Configuration from config.R
#' @return 3x3 covariance matrix
derive_mediation_covariance_from_config <- function(config) {
  derive_mediation_covariance(
    beta_xm = config$beta_xm,
    beta_my = config$beta_my,
    beta_xy = config$beta_xy,
    sigma2_m = config$sigma2_m,
    sigma2_y = config$sigma2_y,
    var_x = config$var_x
  )
}


#' Get Population Parameters as Named List
#'
#' Returns population parameters in a structured format for use in analysis
#'
#' @param config List. Configuration from config.R (or individual parameters)
#' @return List with population parameters and derived quantities
get_population_params <- function(config) {

  # Derived quantities
  indirect_effect <- config$beta_xm * config$beta_my
  direct_effect <- config$beta_xy
  total_effect <- indirect_effect + direct_effect

  # Covariance matrix
  Sigma <- derive_mediation_covariance_from_config(config)

  return(list(
    # Structural parameters
    beta_xm = config$beta_xm,
    beta_my = config$beta_my,
    beta_xy = config$beta_xy,
    sigma2_m = config$sigma2_m,
    sigma2_y = config$sigma2_y,
    var_x = config$var_x,

    # Derived effects
    indirect_effect = indirect_effect,
    direct_effect = direct_effect,
    total_effect = total_effect,

    # Covariance structure
    Sigma = Sigma,
    mu = c(X = 0, M = 0, Y = 0)  # Means are zero
  ))
}


#' Validate Generated Data Against Population Parameters
#'
#' Compares sample statistics to population values (for testing/diagnostics)
#'
#' @param data data.frame with X, M, Y columns
#' @param pop_params List from get_population_params()
#' @return List with comparison statistics
validate_generated_data <- function(data, pop_params) {

  # Sample covariance
  sample_cov <- cov(data)

  # Sample means
  sample_means <- colMeans(data)

  # Differences
  cov_diff <- sample_cov - pop_params$Sigma
  mean_diff <- sample_means - pop_params$mu

  # Max absolute difference
  max_cov_diff <- max(abs(cov_diff))
  max_mean_diff <- max(abs(mean_diff))

  return(list(
    sample_cov = sample_cov,
    population_cov = pop_params$Sigma,
    cov_difference = cov_diff,
    max_cov_diff = max_cov_diff,
    sample_means = sample_means,
    population_means = pop_params$mu,
    mean_difference = mean_diff,
    max_mean_diff = max_mean_diff
  ))
}
