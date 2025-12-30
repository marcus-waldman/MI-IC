# ============================================================================
# Generate Data for Simulation
# ============================================================================
# Functions for:
#   - Generating complete multivariate normal data
#   - Imposing missingness patterns (monotone MCAR)
# ============================================================================

#' Generate Multivariate Normal Data
#'
#' Generates complete data from multivariate normal distribution
#'
#' @param n Integer. Sample size
#' @param mu Numeric vector. Mean vector of length p
#' @param Sigma Numeric matrix. p x p covariance matrix
#' @param seed Integer. Random seed for reproducibility (optional)
#' @return n x p matrix of generated data
#' @details
#'   Uses MASS::mvrnorm for generation
#' @examples
#'   Sigma <- generate_sigma_cs(p = 5, sigma2 = 1, rho = 0.5)
#'   data <- generate_mvn_data(n = 100, mu = rep(0, 5), Sigma = Sigma, seed = 123)
generate_mvn_data <- function(n, mu, Sigma, seed = NULL) {

  # Load required package
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' required for mvrnorm()")
  }

  # Input validation
  p <- length(mu)
  if (!is.matrix(Sigma) || nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("Sigma must be a p x p matrix where p = length(mu)")
  }
  if (n < 1) stop("Sample size n must be at least 1")

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Generate data
  data <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

  # Ensure matrix format even if n=1
  if (n == 1) {
    data <- matrix(data, nrow = 1)
  }

  # Add column names
  colnames(data) <- paste0("Y", 1:p)

  return(data)
}


#' Impose Missingness Pattern on Complete Data
#'
#' Creates missing data by imposing a specific missingness pattern using mice::ampute.
#' Default is monotone MCAR for clean theoretical validation.
#'
#' @param data Numeric matrix. Complete n x p data matrix
#' @param missing_rate Numeric. Overall target proportion of missing values (0 to 1)
#' @param pattern Character. Missingness pattern: "monotone" (default) or "random"
#' @param prop_complete Numeric. Proportion of observations with complete data (default: 0.4)
#' @param seed Integer. Random seed for reproducibility (optional)
#' @return List with components:
#'   \item{data_miss}{n x p matrix with NA values}
#'   \item{missing_pattern}{n x p logical matrix (TRUE = observed, FALSE = missing)}
#'   \item{n_miss_by_var}{Vector of length p with count of missing values per variable}
#'   \item{prop_miss_by_var}{Vector of length p with proportion missing per variable}
#' @details
#'   Monotone pattern using mice::ampute:
#'   - prop_complete (default 40%) of observations are fully observed
#'   - Remaining observations have monotone missingness (Y1 always observed)
#'   - Incomplete observations distributed evenly across monotone patterns
#'
#'   Random pattern: Each value missing independently with probability missing_rate
#'
#' @examples
#'   data_complete <- generate_mvn_data(100, rep(0, 10), diag(10), seed = 123)
#'   result <- impose_missingness(data_complete, missing_rate = 0.6, pattern = "monotone")
#'   data_miss <- result$data_miss
impose_missingness <- function(data, missing_rate, pattern = "monotone",
                               prop_complete = 0.4, seed = NULL) {

  # Input validation
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  n <- nrow(data)
  p <- ncol(data)

  if (missing_rate < 0 || missing_rate > 1) {
    stop("missing_rate must be between 0 and 1")
  }

  if (prop_complete < 0 || prop_complete > 1) {
    stop("prop_complete must be between 0 and 1")
  }

  if (missing_rate == 0) {
    warning("missing_rate = 0; returning original data with no missingness")
    return(list(
      data_miss = data,
      missing_pattern = matrix(TRUE, nrow = n, ncol = p),
      n_miss_by_var = rep(0, p),
      prop_miss_by_var = rep(0, p)
    ))
  }

  pattern <- match.arg(pattern, c("monotone", "random"))

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Apply missingness pattern
  if (pattern == "monotone") {
    # Use mice::ampute for monotone MCAR pattern
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' required for ampute()")
    }

    # Create monotone patterns matrix (for incomplete observations only)
    # Pattern 1: Only Y_p missing (1 1 1 ... 1 0)
    # Pattern 2: Y_{p-1} and Y_p missing (1 1 1 ... 0 0)
    # ...
    # Pattern p-1: Only Y_1 observed (1 0 0 ... 0 0)

    patterns_matrix <- matrix(0, nrow = p - 1, ncol = p)
    for (k in 1:(p - 1)) {
      # Pattern k: first (p - k) variables observed, last k variables missing
      n_obs <- p - k
      patterns_matrix[k, 1:n_obs] <- 1
    }

    # Frequency: distribute incomplete observations evenly across patterns
    freq <- rep(1 / (p - 1), p - 1)

    # Use mice::ampute
    # prop: proportion of cases with ANY missing values
    # mech: "MCAR" for missing completely at random
    ampute_result <- mice::ampute(
      data = as.data.frame(data),
      prop = 1 - prop_complete,
      patterns = patterns_matrix,
      freq = freq,
      mech = "MCAR"
    )

    data_miss <- as.matrix(ampute_result$amp)

  } else if (pattern == "random") {
    # Random pattern: each value independently missing with probability missing_rate
    obs_pattern <- matrix(
      rbinom(n * p, size = 1, prob = 1 - missing_rate) == 1,
      nrow = n, ncol = p
    )

    # Ensure at least one variable is always observed for identifiability
    # (make Y1 fully observed)
    obs_pattern[, 1] <- TRUE

    data_miss <- data
    data_miss[!obs_pattern] <- NA
  }

  # Compute missingness summaries
  obs_pattern <- !is.na(data_miss)
  n_miss_by_var <- colSums(!obs_pattern)
  prop_miss_by_var <- n_miss_by_var / n

  # Add column names if present in original data
  if (!is.null(colnames(data))) {
    colnames(data_miss) <- colnames(data)
    colnames(obs_pattern) <- colnames(data)
    names(n_miss_by_var) <- colnames(data)
    names(prop_miss_by_var) <- colnames(data)
  }

  return(list(
    data_miss = data_miss,
    missing_pattern = obs_pattern,
    n_miss_by_var = n_miss_by_var,
    prop_miss_by_var = prop_miss_by_var
  ))
}


#' Generate Complete MVN Data with Missingness (Wrapper)
#'
#' Convenience function that generates complete data and imposes missingness in one step
#'
#' @param n Integer. Sample size
#' @param mu Numeric vector. Mean vector
#' @param Sigma Numeric matrix. Covariance matrix
#' @param missing_rate Numeric. Target proportion of missing values
#' @param pattern Character. "monotone" or "random"
#' @param prop_complete Numeric. Proportion of complete observations (default: 0.4)
#' @param seed_data Integer. Seed for data generation
#' @param seed_miss Integer. Seed for missingness generation
#' @return List with components:
#'   \item{data_complete}{Complete n x p data matrix}
#'   \item{data_miss}{Data matrix with NA values}
#'   \item{missing_pattern}{Logical matrix indicating observed values}
#'   \item{n_miss_by_var}{Count of missing per variable}
#'   \item{prop_miss_by_var}{Proportion missing per variable}
#' @examples
#'   Sigma <- generate_sigma_cs(p = 10, sigma2 = 1, rho = 0.5)
#'   result <- generate_data_with_missingness(
#'     n = 200, mu = rep(0, 10), Sigma = Sigma,
#'     missing_rate = 0.6, pattern = "monotone",
#'     prop_complete = 0.4,
#'     seed_data = 123, seed_miss = 456
#'   )
generate_data_with_missingness <- function(n, mu, Sigma,
                                            missing_rate,
                                            pattern = "monotone",
                                            prop_complete = 0.4,
                                            seed_data = NULL,
                                            seed_miss = NULL) {

  # Generate complete data
  data_complete <- generate_mvn_data(n, mu, Sigma, seed = seed_data)

  # Impose missingness
  miss_result <- impose_missingness(
    data_complete,
    missing_rate = missing_rate,
    pattern = pattern,
    prop_complete = prop_complete,
    seed = seed_miss
  )

  # Combine results
  result <- c(
    list(data_complete = data_complete),
    miss_result
  )

  return(result)
}


#' Get Missingness Summary Statistics
#'
#' Computes summary statistics for a data matrix with missing values
#'
#' @param data_miss Numeric matrix with NA values
#' @return List with summary statistics
#' @examples
#'   # After generating data with missingness
#'   summary_stats <- get_missingness_summary(data_miss)
get_missingness_summary <- function(data_miss) {

  if (!is.matrix(data_miss)) {
    data_miss <- as.matrix(data_miss)
  }

  n <- nrow(data_miss)
  p <- ncol(data_miss)

  # Overall missingness
  total_values <- n * p
  total_missing <- sum(is.na(data_miss))
  prop_missing_overall <- total_missing / total_values

  # Per-variable missingness
  n_miss_by_var <- colSums(is.na(data_miss))
  prop_miss_by_var <- n_miss_by_var / n

  # Per-observation missingness
  n_miss_by_obs <- rowSums(is.na(data_miss))
  prop_miss_by_obs <- n_miss_by_obs / p

  # Missingness patterns
  n_complete_cases <- sum(complete.cases(data_miss))
  prop_complete_cases <- n_complete_cases / n

  list(
    n = n,
    p = p,
    total_missing = total_missing,
    prop_missing_overall = prop_missing_overall,
    n_miss_by_var = n_miss_by_var,
    prop_miss_by_var = prop_miss_by_var,
    prop_complete_cases = prop_complete_cases,
    summary_table = data.frame(
      variable = if (!is.null(colnames(data_miss))) colnames(data_miss) else paste0("Y", 1:p),
      n_missing = n_miss_by_var,
      prop_missing = round(prop_miss_by_var, 3)
    )
  )
}
