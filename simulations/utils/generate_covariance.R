# ============================================================================
# Generate Covariance Matrices for Simulation
# ============================================================================
# Functions to create true covariance matrices for different structures:
#   - Compound Symmetry (CS)
#   - Toeplitz (general, with distance-dependent correlations)
#   - Unstructured (random positive definite)
# ============================================================================

#' Generate Compound Symmetry Covariance Matrix
#'
#' Creates a covariance matrix with equal variances and equal covariances.
#' All pairwise correlations are identical (exchangeable structure).
#'
#' @param p Integer. Dimension of the covariance matrix
#' @param sigma2 Numeric. Variance parameter (default = 1)
#' @param rho Numeric. Correlation parameter, must be in (-1/(p-1), 1)
#' @return p x p covariance matrix with compound symmetry structure
#' @details
#'   Sigma_ij = sigma2 * rho   if i != j
#'   Sigma_ii = sigma2          if i == j
#'   Number of parameters Q = p + 2 (means + variance + correlation)
#' @examples
#'   Sigma_cs <- generate_sigma_cs(p = 10, sigma2 = 1, rho = 0.5)
generate_sigma_cs <- function(p, sigma2 = 1, rho) {

  # Input validation
  if (p < 2) stop("Dimension p must be at least 2")
  if (sigma2 <= 0) stop("Variance sigma2 must be positive")
  if (rho <= -1/(p-1) || rho >= 1) {
    stop(sprintf("Correlation rho must be in (%.3f, 1) for positive definiteness", -1/(p-1)))
  }

  # Create correlation matrix
  R <- matrix(rho, nrow = p, ncol = p)
  diag(R) <- 1

  # Scale to covariance matrix
  Sigma <- sigma2 * R

  # Verify positive definiteness
  min_eig <- min(eigen(Sigma, only.values = TRUE)$values)
  if (min_eig < 1e-8) {
    stop(sprintf("Generated matrix not positive definite (min eigenvalue = %.2e)", min_eig))
  }

  return(Sigma)
}


#' Generate Toeplitz Covariance Matrix
#'
#' Creates a covariance matrix where correlations depend only on distance (lag).
#' Uses a decay function to generate lag-specific correlations.
#'
#' @param p Integer. Dimension of the covariance matrix
#' @param sigma2 Numeric. Variance parameter (default = 1)
#' @param rho Numeric. Baseline correlation parameter (used in default decay function)
#' @param rho_fn Function. Custom decay function rho_fn(k, p) returning correlation at lag k.
#'        Default is linear decay: rho * max(0, 1 - k/p)
#' @param seed Integer. Random seed for reproducibility if generating random decay (optional)
#' @return p x p covariance matrix with Toeplitz structure
#' @details
#'   Sigma_ij = sigma2 * rho_{|i-j|}
#'   Number of parameters Q = 2p (means + variance + (p-1) lag correlations)
#'   Default decay ensures positive definiteness through monotone decrease
#' @examples
#'   # Linear decay (default)
#'   Sigma_toep <- generate_sigma_toeplitz(p = 10, sigma2 = 1, rho = 0.7)
#'
#'   # Custom decay function
#'   Sigma_toep2 <- generate_sigma_toeplitz(
#'     p = 10, sigma2 = 1, rho = 0.7,
#'     rho_fn = function(k, p) 0.8 * exp(-0.3 * k)
#'   )
generate_sigma_toeplitz <- function(p, sigma2 = 1, rho, rho_fn = NULL, seed = NULL) {

  # Input validation
  if (p < 2) stop("Dimension p must be at least 2")
  if (sigma2 <= 0) stop("Variance sigma2 must be positive")
  if (rho < 0 || rho >= 1) stop("Baseline correlation rho must be in [0, 1)")

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Default decay function: linear decay
  if (is.null(rho_fn)) {
    rho_fn <- function(k, p) {
      rho * pmax(0, 1 - k/p)
    }
  }

  # Generate lag correlations
  rho_vec <- sapply(1:(p-1), rho_fn, p = p)

  # Check monotone decrease (sufficient for PD)
  if (any(diff(rho_vec) > 1e-10)) {
    warning("Correlations not monotonically decreasing; PD not guaranteed")
  }

  # Create correlation matrix using toeplitz()
  R <- toeplitz(c(1, rho_vec))

  # Scale to covariance matrix
  Sigma <- sigma2 * R

  # Verify positive definiteness
  min_eig <- min(eigen(Sigma, only.values = TRUE)$values)
  if (min_eig < 1e-8) {
    stop(sprintf(
      "Generated Toeplitz matrix not positive definite (min eigenvalue = %.2e).
      Try smaller rho or different decay function.",
      min_eig
    ))
  }

  return(Sigma)
}


#' Generate Unstructured Covariance Matrix
#'
#' Creates a random positive definite covariance matrix with approximate
#' average correlation level controlled by rho.
#'
#' @param p Integer. Dimension of the covariance matrix
#' @param sigma2 Numeric. Variance parameter (default = 1)
#' @param rho Numeric. Target average correlation level (0 to 1)
#' @param seed Integer. Random seed for reproducibility
#' @return p x p unstructured covariance matrix
#' @details
#'   Uses Wishart distribution to generate random PD matrix, then scales
#'   to achieve approximate target correlation.
#'   Number of parameters Q = p + p(p+1)/2 (means + full covariance matrix)
#' @examples
#'   Sigma_unstr <- generate_sigma_unstructured(p = 10, sigma2 = 1, rho = 0.5, seed = 123)
generate_sigma_unstructured <- function(p, sigma2 = 1, rho = 0.5, seed = NULL) {

  # Input validation
  if (p < 2) stop("Dimension p must be at least 2")
  if (sigma2 <= 0) stop("Variance sigma2 must be positive")
  if (rho < 0 || rho >= 1) stop("Target correlation rho must be in [0, 1)")

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Generate random PD matrix using Wishart distribution
  # Higher df gives more stable/less extreme correlations
  df <- p + 10
  S_base <- diag(p)
  Sigma_raw <- rWishart(1, df = df, Sigma = S_base)[,,1]

  # Convert to correlation matrix
  R_raw <- cov2cor(Sigma_raw)

  # Scale correlations toward target rho
  # New correlation = rho + (1-rho) * (original correlation - target)
  # This centers correlations around rho while preserving PD
  R_scaled <- R_raw
  diag(R_scaled) <- 1

  # Simple scaling approach: multiply off-diagonals by factor to achieve target mean
  mean_corr <- mean(R_raw[lower.tri(R_raw)])
  if (mean_corr > 1e-6) {
    scale_factor <- rho / mean_corr
    # Don't scale too aggressively to maintain PD
    scale_factor <- min(scale_factor, 2)
    R_scaled[lower.tri(R_scaled)] <- R_raw[lower.tri(R_raw)] * scale_factor
    R_scaled[upper.tri(R_scaled)] <- t(R_scaled)[upper.tri(R_scaled)]
  }

  # Ensure still positive definite after scaling
  min_eig <- min(eigen(R_scaled, only.values = TRUE)$values)
  if (min_eig < 1e-8) {
    # If not PD, use original unscaled version
    warning("Scaled matrix not PD; using unscaled Wishart-generated matrix")
    R_scaled <- R_raw
  }

  # Scale to covariance matrix
  Sigma <- sigma2 * R_scaled

  return(Sigma)
}


#' Get Number of Parameters for Covariance Structure
#'
#' Returns the number of parameters Q for a given covariance structure
#' (excluding means, which add p parameters)
#'
#' @param p Integer. Dimension
#' @param structure Character. One of "CS", "Toeplitz", "Unstructured"
#' @return Integer. Number of covariance parameters
#' @examples
#'   get_cov_nparams(10, "CS")          # Returns 2
#'   get_cov_nparams(10, "Toeplitz")    # Returns 10
#'   get_cov_nparams(10, "Unstructured") # Returns 55
get_cov_nparams <- function(p, structure) {
  switch(structure,
    "CS" = 2,                    # sigma2, rho
    "Toeplitz" = p,              # sigma2, rho_1, ..., rho_{p-1}
    "Unstructured" = p*(p+1)/2,  # full covariance matrix
    stop("Unknown structure. Must be 'CS', 'Toeplitz', or 'Unstructured'")
  )
}


#' Generate Covariance Matrix (Wrapper)
#'
#' Convenient wrapper function that calls the appropriate generator based on structure
#'
#' @param p Integer. Dimension
#' @param structure Character. One of "CS", "Toeplitz", "Unstructured"
#' @param sigma2 Numeric. Variance parameter
#' @param rho Numeric. Correlation parameter
#' @param seed Integer. Random seed (used for Toeplitz and Unstructured)
#' @param ... Additional arguments passed to specific generators
#' @return p x p covariance matrix
#' @examples
#'   Sigma1 <- generate_covariance(10, "CS", sigma2 = 1, rho = 0.5)
#'   Sigma2 <- generate_covariance(10, "Toeplitz", sigma2 = 1, rho = 0.7, seed = 123)
#'   Sigma3 <- generate_covariance(10, "Unstructured", sigma2 = 1, rho = 0.5, seed = 456)
generate_covariance <- function(p, structure, sigma2 = 1, rho, seed = NULL, ...) {

  structure <- match.arg(structure, c("CS", "Toeplitz", "Unstructured"))

  Sigma <- switch(structure,
    "CS" = generate_sigma_cs(p, sigma2, rho),
    "Toeplitz" = generate_sigma_toeplitz(p, sigma2, rho, seed = seed, ...),
    "Unstructured" = generate_sigma_unstructured(p, sigma2, rho, seed = seed)
  )

  return(Sigma)
}
