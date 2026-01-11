# ============================================================================
# Analytical Likelihood and Q-Functions for Mediation Model
# ============================================================================
# Assumes MONOTONE missingness: Y missing alone, or M+Y missing (never M alone)
# This gives 3 patterns where all θ̃-dependent terms are quadratic

#' Complete-Data Negative Log-Likelihood
#'
#' @param theta c(a, b, log_sigma2_M, log_sigma2_Y)
#' @param data data.frame with X, M, Y
#' @return Negative log-likelihood
negloglik_complete <- function(theta, data) {
  a <- theta[1]
  b <- theta[2]
  sigma2_M <- exp(theta[3])
  sigma2_Y <- exp(theta[4])

  X <- data$X
  M <- data$M
  Y <- data$Y

  ll_M <- sum(dnorm(M, mean = a * X, sd = sqrt(sigma2_M), log = TRUE))
  ll_Y <- sum(dnorm(Y, mean = b * M, sd = sqrt(sigma2_Y), log = TRUE))

  -(ll_M + ll_Y)
}

#' Observed-Data Negative Log-Likelihood (FIML) for Monotone Missingness
#'
#' @param theta c(a, b, log_sigma2_M, log_sigma2_Y)
#' @param data data.frame with X, M, Y (may contain NA)
#' @return Negative log-likelihood
negloglik_observed <- function(theta, data) {
  a <- theta[1]
  b <- theta[2]
  sigma2_M <- exp(theta[3])
  sigma2_Y <- exp(theta[4])

  X <- data$X
  M <- data$M
  Y <- data$Y
  n <- nrow(data)

  ll_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: Complete
      ll_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
              dnorm(Y[i], b * M[i], sqrt(sigma2_Y), log = TRUE)

    } else if (M_obs && !Y_obs) {
      # Pattern 2: Y missing only
      ll_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE)

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: Both M and Y missing
      ll_i <- 0

    } else {
      # M missing, Y observed - should not occur under monotone missingness
      stop("Non-monotone pattern detected: M missing but Y observed at row ", i)
    }

    ll_total <- ll_total + ll_i
  }

  -ll_total
}


#' Analytical Q-Function for Monotone Missingness (Improper MI)
#'
#' Computes Q(θ|θ̃) where θ̃ is fixed (improper MI).
#' Only handles the 3 monotone patterns.
#'
#' @param data_miss data.frame with missing values (monotone pattern)
#' @param theta_eval Parameters to evaluate at c(a, b, sigma2_M, sigma2_Y)
#' @param theta_impute Parameters for imputation c(a, b, sigma2_M, sigma2_Y)
#' @return Q-function value
compute_Q_analytical <- function(data_miss, theta_eval, theta_impute) {
  a <- theta_eval[1]
  b <- theta_eval[2]
  sigma2_M <- theta_eval[3]
  sigma2_Y <- theta_eval[4]

  a_tilde <- theta_impute[1]
  b_tilde <- theta_impute[2]
  sigma2_M_tilde <- theta_impute[3]
  sigma2_Y_tilde <- theta_impute[4]

  X <- data_miss$X
  M <- data_miss$M
  Y <- data_miss$Y
  n <- nrow(data_miss)

  Q_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: Complete - no dependence on θ̃
      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             dnorm(Y[i], b * M[i], sqrt(sigma2_Y), log = TRUE)

    } else if (M_obs && !Y_obs) {
      # Pattern 2: Y missing only
      # E[Y|M] = b̃*M, Var[Y|M] = σ²_Ỹ
      E_Y <- b_tilde * M[i]
      Var_Y <- sigma2_Y_tilde

      # E[(Y - bM)²|M] = Var[Y|M] + (E[Y|M] - bM)²
      E_Y_resid_sq <- Var_Y + (E_Y - b * M[i])^2

      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: Both M and Y missing
      # Posterior moments from imputation model
      E_M <- a_tilde * X[i]
      Var_M <- sigma2_M_tilde
      E_Y <- a_tilde * b_tilde * X[i]
      Var_Y <- b_tilde^2 * sigma2_M_tilde + sigma2_Y_tilde
      Cov_MY <- b_tilde * sigma2_M_tilde

      # E[(M - aX)²|X] = Var[M|X] + (E[M|X] - aX)²
      E_M_resid_sq <- Var_M + (E_M - a * X[i])^2

      # E[(Y - bM)²|X] using Var(Y - bM) = Var(Y) + b²Var(M) - 2b*Cov(M,Y)
      Var_Y_minus_bM <- Var_Y + b^2 * Var_M - 2 * b * Cov_MY
      E_Y_minus_bM <- E_Y - b * E_M
      E_Y_resid_sq <- Var_Y_minus_bM + E_Y_minus_bM^2

      Q_i <- -0.5 * log(2 * pi * sigma2_M) - E_M_resid_sq / (2 * sigma2_M) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)

    } else {
      stop("Non-monotone pattern detected: M missing but Y observed at row ", i)
    }

    Q_total <- Q_total + Q_i
  }

  Q_total
}


#' Analytical Q-Function with Proper Bayesian MI (Delta Method Correction)
#'
#' Computes E_θ̃[Q(θ|θ̃)] where θ̃ ~ N(θ̂, V_θ) using second-order Taylor expansion.
#' The correction is: E[Q] ≈ Q(θ|θ̂) + (1/2) tr(H_Q · V_θ)
#'
#' Under monotone missingness, all θ̃-dependent terms are quadratic, so the
#' Hessian H_Q is constant and the Taylor expansion is exact.
#'
#' @param data_miss data.frame with missing values (monotone pattern)
#' @param theta_eval Parameters to evaluate at c(a, b, sigma2_M, sigma2_Y)
#' @param theta_hat Posterior mean of θ̃ (= θ̂_obs)
#' @param V_theta Posterior covariance of θ̃ (= I_obs^{-1}), 4x4 matrix
#' @return List with Q_improper, hessian_correction, Q_proper
compute_Q_analytical_proper <- function(data_miss, theta_eval, theta_hat, V_theta) {

  # First compute improper Q at θ̃ = θ̂
  Q_improper <- compute_Q_analytical(data_miss, theta_eval, theta_hat)

  # Extract parameters for Hessian computation
  a <- theta_eval[1]
  b <- theta_eval[2]
  sigma2_M <- theta_eval[3]
  sigma2_Y <- theta_eval[4]

  a_hat <- theta_hat[1]
  sigma2_M_hat <- theta_hat[3]

  X <- data_miss$X
  M <- data_miss$M
  Y <- data_miss$Y
  n <- nrow(data_miss)

  # Accumulate Hessian correction by pattern
  # Correction = (1/2) Σ_i tr(H_i · V_θ)
  # For diagonal H_i, this is (1/2) Σ_i Σ_j H_i[j,j] * V_θ[j,j]

  correction_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: Complete - no dependence on θ̃, H_i = 0
      # No correction

    } else if (M_obs && !Y_obs) {
      # Pattern 2: Y missing only
      # Q_i contains: -E[(Y - bM)²]/(2σ²_Y) = -(σ²_Ỹ + (b̃-b)²M²)/(2σ²_Y)
      # ∂²/∂b̃² = -M²/σ²_Y
      # H_i = diag(0, -M_i²/σ²_Y, 0, 0)

      H_bb <- -M[i]^2 / sigma2_Y
      correction_i <- 0.5 * H_bb * V_theta[2, 2]
      correction_total <- correction_total + correction_i

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: Both M and Y missing
      # M-residual: -(σ²_M̃ + (ã-a)²X²)/(2σ²_M)
      #   ∂²/∂ã² = -X²/σ²_M
      #
      # Y-residual (simplified at θ̃ = θ̂ where cross-deriv vanishes):
      #   ∂²/∂b̃² = -(σ²_M̃ + ã²X²)/σ²_Y  evaluated at ã = â
      #
      # H_i = diag(-X²/σ²_M, -(σ²_M_hat + a_hat²X²)/σ²_Y, 0, 0)

      H_aa <- -X[i]^2 / sigma2_M
      H_bb <- -(sigma2_M_hat + a_hat^2 * X[i]^2) / sigma2_Y

      correction_i <- 0.5 * (H_aa * V_theta[1, 1] + H_bb * V_theta[2, 2])
      correction_total <- correction_total + correction_i

    } else {
      stop("Non-monotone pattern detected: M missing but Y observed at row ", i)
    }
  }

  Q_proper <- Q_improper + correction_total

  list(
    Q_improper = Q_improper,
    hessian_correction = correction_total,
    Q_proper = Q_proper
  )
}
