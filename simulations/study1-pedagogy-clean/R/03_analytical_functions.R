# ============================================================================
# Analytical Likelihood and Q-Functions for Mediation Model
# ============================================================================
# Assumes MONOTONE missingness: Y missing alone, or M+Y missing (never M alone)
# This gives 3 patterns where all θ̃-dependent terms are quadratic


#' Estimate Imputation Model Parameters
#'
#' @param data_miss Data frame with missing values (X always observed, M and Y may be missing)
#' @param theta_obs Observed-data MLE c(a, b, sigma2_M, sigma2_Y) on natural scale
#' @param model One of "A", "B", "C"
#' @return List of imputation model parameters
#' @details
#' Model A (Congenial): M|X, Y|M structure matching DGP
#'   Returns: list(a, b, sigma2_M, sigma2_Y)
#'
#' Model B (Uncongenial): M|X, Y marginal (ignores M->Y)
#'   Returns: list(a, sigma2_M, mu_Y, tau2_Y)
#'
#' Model C (Saturated MVN): Full (X,M,Y) covariance
#'   Returns: list(mu_X, mu_M, mu_Y, Sigma) where Sigma is 3x3 matrix
estimate_imputation_params <- function(data_miss, theta_obs, model = "A") {

  if (model == "A") {
    # Model A: Use observed-data MLE (already congenial)
    # theta_obs from FIML is the correct imputation parameter
    list(
      a = theta_obs[1],
      b = theta_obs[2],
      sigma2_M = theta_obs[3],
      sigma2_Y = theta_obs[4]
    )

  } else if (model == "B") {
    # Model B: Fit M ~ X, compute marginal Y statistics

    # Fit M ~ X on observed M only
    obs_M <- !is.na(data_miss$M)
    if (sum(obs_M) < 2) {
      stop("Model B requires at least 2 observed M values")
    }

    fit_M <- lm(M ~ X, data = data_miss[obs_M, ])
    a_tilde <- coef(fit_M)[2]
    sigma2_M_tilde <- sigma(fit_M)^2

    # Marginal Y statistics (observed Y only)
    obs_Y <- !is.na(data_miss$Y)
    if (sum(obs_Y) < 2) {
      stop("Model B requires at least 2 observed Y values")
    }

    mu_Y <- mean(data_miss$Y[obs_Y])
    tau2_Y <- var(data_miss$Y[obs_Y])

    list(
      a = a_tilde,
      sigma2_M = sigma2_M_tilde,
      mu_Y = mu_Y,
      tau2_Y = tau2_Y
    )

  } else if (model == "C") {
    # Model C: Estimate full (X, M, Y) covariance matrix

    # Means (pairwise deletion handles NA)
    mu_X <- mean(data_miss$X, na.rm = TRUE)
    mu_M <- mean(data_miss$M, na.rm = TRUE)
    mu_Y <- mean(data_miss$Y, na.rm = TRUE)

    # Covariance matrix (pairwise deletion)
    Sigma_hat <- matrix(0, 3, 3)

    # Variances
    Sigma_hat[1, 1] <- var(data_miss$X, na.rm = TRUE)
    Sigma_hat[2, 2] <- var(data_miss$M, na.rm = TRUE)
    Sigma_hat[3, 3] <- var(data_miss$Y, na.rm = TRUE)

    # Covariances (pairwise complete observations)
    Sigma_hat[1, 2] <- Sigma_hat[2, 1] <- cov(data_miss$X, data_miss$M, use = "pairwise.complete.obs")
    Sigma_hat[1, 3] <- Sigma_hat[3, 1] <- cov(data_miss$X, data_miss$Y, use = "pairwise.complete.obs")
    Sigma_hat[2, 3] <- Sigma_hat[3, 2] <- cov(data_miss$M, data_miss$Y, use = "pairwise.complete.obs")

    list(
      mu_X = mu_X,
      mu_M = mu_M,
      mu_Y = mu_Y,
      Sigma = Sigma_hat
    )

  } else {
    stop("model must be one of 'A', 'B', or 'C'")
  }
}


#' Estimate Posterior Covariance Matrix for Imputation Parameters
#'
#' @param data_miss Data frame with missing values
#' @param theta_impute Imputation parameters from estimate_imputation_params()
#' @param V_theta_obs Observed-data parameter covariance (4x4, from FIML Hessian)
#' @param model One of "A", "B", "C"
#' @return Covariance matrix for imputation parameters
#' @details
#' Model A: Returns V_theta_obs directly (4x4)
#' Model B: Block-diagonal (4x4) from lm() and sample variance
#' Model C: Approximate asymptotic covariance (9x9) for MVN MLE
estimate_V_theta <- function(data_miss, theta_impute, V_theta_obs, model = "A") {

  if (model == "A") {
    # Model A: Use observed-data FIML covariance
    return(V_theta_obs)

  } else if (model == "B") {
    # Model B: Estimate V_theta for (ã, σ̃²_M, μ̃_Y, τ̃²_Y)

    # Fit M ~ X to get covariance for (ã, σ̃²_M)
    obs_M <- !is.na(data_miss$M)
    fit_M <- lm(M ~ X, data = data_miss[obs_M, ])

    # vcov(fit_M) gives variance for (intercept, slope)
    # We need variance for (slope, sigma^2)
    vcov_M <- vcov(fit_M)
    V_a <- vcov_M[2, 2]  # Variance of slope (ã)

    # Variance of σ̃²_M: use asymptotic variance 2σ⁴/n
    n_M <- sum(obs_M)
    V_sigma2_M <- 2 * theta_impute$sigma2_M^2 / n_M

    # Variance for (μ̃_Y, τ̃²_Y)
    obs_Y <- !is.na(data_miss$Y)
    n_Y <- sum(obs_Y)
    V_mu_Y <- theta_impute$tau2_Y / n_Y  # Variance of sample mean
    V_tau2_Y <- 2 * theta_impute$tau2_Y^2 / n_Y  # Variance of sample variance

    # Construct block-diagonal 4x4 matrix
    # Order: (ã, σ̃²_M, μ̃_Y, τ̃²_Y)
    V_theta <- matrix(0, 4, 4)
    V_theta[1, 1] <- V_a
    V_theta[2, 2] <- V_sigma2_M
    V_theta[3, 3] <- V_mu_Y
    V_theta[4, 4] <- V_tau2_Y

    return(V_theta)

  } else if (model == "C") {
    # Model C: Asymptotic covariance for MVN MLE
    # Parameters: (μ_X, μ_M, μ_Y, Σ_11, Σ_22, Σ_33, Σ_12, Σ_13, Σ_23)

    # Get sample size (use complete X since it's always observed)
    n <- nrow(data_miss)

    # Asymptotic theory for MVN MLE:
    # Var(μ̂) = Σ/n
    # Var(vec(Σ̂)) ≈ (2/n) * (Σ ⊗ Σ) for unique elements

    Sigma <- theta_impute$Sigma

    # Simplified approximation: block-diagonal between means and covariances
    # V_theta is 9x9: [V_means (3x3), 0; 0, V_cov (6x6)]

    V_theta <- matrix(0, 9, 9)

    # Variance of means (diagonal of Σ/n)
    V_theta[1, 1] <- Sigma[1, 1] / n  # Var(μ̂_X)
    V_theta[2, 2] <- Sigma[2, 2] / n  # Var(μ̂_M)
    V_theta[3, 3] <- Sigma[3, 3] / n  # Var(μ̂_Y)

    # Variance of covariance elements (approximate)
    # For Σ_ij, Var(Σ̂_ij) ≈ (1/n)[Σ_ii*Σ_jj + Σ_ij²]

    # Σ_11 (variance of X)
    V_theta[4, 4] <- (2 / n) * Sigma[1, 1]^2

    # Σ_22 (variance of M)
    V_theta[5, 5] <- (2 / n) * Sigma[2, 2]^2

    # Σ_33 (variance of Y)
    V_theta[6, 6] <- (2 / n) * Sigma[3, 3]^2

    # Σ_12 (covariance X,M)
    V_theta[7, 7] <- (1 / n) * (Sigma[1, 1] * Sigma[2, 2] + Sigma[1, 2]^2)

    # Σ_13 (covariance X,Y)
    V_theta[8, 8] <- (1 / n) * (Sigma[1, 1] * Sigma[3, 3] + Sigma[1, 3]^2)

    # Σ_23 (covariance M,Y)
    V_theta[9, 9] <- (1 / n) * (Sigma[2, 2] * Sigma[3, 3] + Sigma[2, 3]^2)

    return(V_theta)

  } else {
    stop("model must be one of 'A', 'B', or 'C'")
  }
}


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
#' Supports three imputation models: A (congenial), B (uncongenial), C (saturated MVN)
#'
#' @param data_miss data.frame with missing values (monotone pattern)
#' @param theta_eval Parameters to evaluate at c(a, b, sigma2_M, sigma2_Y)
#' @param theta_impute List of imputation parameters (from estimate_imputation_params)
#' @param model One of "A", "B", "C"
#' @return Q-function value
compute_Q_analytical <- function(data_miss, theta_eval, theta_impute, model = "A") {
  # Analysis model parameters
  a <- theta_eval[1]
  b <- theta_eval[2]
  sigma2_M <- theta_eval[3]
  sigma2_Y <- theta_eval[4]

  X <- data_miss$X
  M <- data_miss$M
  Y <- data_miss$Y
  n <- nrow(data_miss)

  Q_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: Complete - no imputation needed (same for all models)
      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             dnorm(Y[i], b * M[i], sqrt(sigma2_Y), log = TRUE)

    } else if (M_obs && !Y_obs) {
      # Pattern 2: Y missing only - posterior moments depend on model

      if (model == "A") {
        # Model A: E[Y|M] = b̃*M, Var[Y|M] = σ̃²_Y
        E_Y <- theta_impute$b * M[i]
        Var_Y <- theta_impute$sigma2_Y

      } else if (model == "B") {
        # Model B: E[Y|M] = μ̃_Y (ignores M!), Var[Y|M] = τ̃²_Y
        E_Y <- theta_impute$mu_Y
        Var_Y <- theta_impute$tau2_Y

      } else if (model == "C") {
        # Model C: Conditional Y|X,M from MVN
        Sigma_YXM <- c(theta_impute$Sigma[3, 1], theta_impute$Sigma[3, 2])
        Sigma_XM <- theta_impute$Sigma[1:2, 1:2]

        resid <- c(X[i] - theta_impute$mu_X, M[i] - theta_impute$mu_M)
        E_Y <- as.numeric(theta_impute$mu_Y + Sigma_YXM %*% solve(Sigma_XM) %*% resid)
        Var_Y <- as.numeric(theta_impute$Sigma[3, 3] - Sigma_YXM %*% solve(Sigma_XM) %*% Sigma_YXM)

      } else {
        stop("model must be one of 'A', 'B', or 'C'")
      }

      # E[(Y - bM)²|M] = Var[Y|M] + (E[Y|M] - bM)²
      E_Y_resid_sq <- Var_Y + (E_Y - b * M[i])^2

      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: Both M and Y missing - posterior moments depend on model

      if (model == "A") {
        # Model A: Mediation structure
        E_M <- theta_impute$a * X[i]
        Var_M <- theta_impute$sigma2_M
        E_Y <- theta_impute$a * theta_impute$b * X[i]
        Var_Y <- theta_impute$b^2 * theta_impute$sigma2_M + theta_impute$sigma2_Y
        Cov_MY <- theta_impute$b * theta_impute$sigma2_M

      } else if (model == "B") {
        # Model B: M|X, Y marginal (Cov=0!)
        E_M <- theta_impute$a * X[i]
        Var_M <- theta_impute$sigma2_M
        E_Y <- theta_impute$mu_Y
        Var_Y <- theta_impute$tau2_Y
        Cov_MY <- 0  # Incorrectly assumes independence

      } else if (model == "C") {
        # Model C: Conditional (M,Y)|X from MVN
        Sigma_MYX <- theta_impute$Sigma[2:3, 1, drop = FALSE]
        Sigma_X <- theta_impute$Sigma[1, 1]

        resid_X <- X[i] - theta_impute$mu_X
        E_MY <- c(theta_impute$mu_M, theta_impute$mu_Y) +
                as.vector((Sigma_MYX / Sigma_X) * resid_X)

        E_M <- E_MY[1]
        E_Y <- E_MY[2]

        Sigma_MY_cond <- theta_impute$Sigma[2:3, 2:3] -
                         (Sigma_MYX %*% t(Sigma_MYX)) / Sigma_X

        Var_M <- Sigma_MY_cond[1, 1]
        Var_Y <- Sigma_MY_cond[2, 2]
        Cov_MY <- Sigma_MY_cond[1, 2]

      } else {
        stop("model must be one of 'A', 'B', or 'C'")
      }

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
#' Supports three imputation models: A (congenial), B (uncongenial), C (saturated MVN)
#' - Model A: Analytical Hessian (existing implementation)
#' - Model B: Analytical Hessian for (ã, σ̃²_M, μ̃_Y, τ̃²_Y)
#' - Model C: Numerical Hessian via finite differences
#'
#' @param data_miss data.frame with missing values (monotone pattern)
#' @param theta_eval Parameters to evaluate at c(a, b, sigma2_M, sigma2_Y)
#' @param theta_impute List of imputation parameters (from estimate_imputation_params)
#' @param V_theta Posterior covariance of θ̃ (model-specific dimensions)
#' @param model One of "A", "B", "C"
#' @return List with Q_improper, hessian_correction, Q_proper
compute_Q_analytical_proper <- function(data_miss, theta_eval, theta_impute,
                                         V_theta, model = "A") {

  # First compute improper Q at θ̃ = θ̂
  Q_improper <- compute_Q_analytical(data_miss, theta_eval, theta_impute, model = model)

  # Compute Hessian correction based on model
  if (model == "A") {
    correction_total <- compute_hessian_correction_A(data_miss, theta_eval,
                                                       theta_impute, V_theta)

  } else if (model == "B") {
    correction_total <- compute_hessian_correction_B(data_miss, theta_eval,
                                                       theta_impute, V_theta)

  } else if (model == "C") {
    correction_total <- compute_hessian_correction_C(data_miss, theta_eval,
                                                       theta_impute, V_theta)

  } else {
    stop("model must be one of 'A', 'B', or 'C'")
  }

  Q_proper <- Q_improper + correction_total

  list(
    Q_improper = Q_improper,
    hessian_correction = correction_total,
    Q_proper = Q_proper
  )
}


#' Hessian Correction for Model A (Congenial)
#'
#' @keywords internal
compute_hessian_correction_A <- function(data_miss, theta_eval, theta_impute, V_theta) {
  a <- theta_eval[1]
  b <- theta_eval[2]
  sigma2_M <- theta_eval[3]
  sigma2_Y <- theta_eval[4]

  a_hat <- theta_impute$a
  sigma2_M_hat <- theta_impute$sigma2_M

  X <- data_miss$X
  M <- data_miss$M
  n <- nrow(data_miss)

  correction_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(data_miss$Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: No correction

    } else if (M_obs && !Y_obs) {
      # Pattern 2: ∂²/∂b̃² = -M²/σ²_Y
      H_bb <- -M[i]^2 / sigma2_Y
      correction_total <- correction_total + 0.5 * H_bb * V_theta[2, 2]

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: ∂²/∂ã² and ∂²/∂b̃²
      H_aa <- -X[i]^2 / sigma2_M
      H_bb <- -(sigma2_M_hat + a_hat^2 * X[i]^2) / sigma2_Y

      correction_total <- correction_total + 0.5 * (H_aa * V_theta[1, 1] + H_bb * V_theta[2, 2])
    }
  }

  correction_total
}


#' Hessian Correction for Model B (Uncongenial)
#'
#' @keywords internal
compute_hessian_correction_B <- function(data_miss, theta_eval, theta_impute, V_theta) {
  a <- theta_eval[1]
  b <- theta_eval[2]
  sigma2_M <- theta_eval[3]
  sigma2_Y <- theta_eval[4]

  X <- data_miss$X
  M <- data_miss$M
  n <- nrow(data_miss)

  correction_total <- 0

  for (i in 1:n) {
    M_obs <- !is.na(M[i])
    Y_obs <- !is.na(data_miss$Y[i])

    if (M_obs && Y_obs) {
      # Pattern 1: No correction

    } else if (M_obs && !Y_obs) {
      # Pattern 2: Q depends on μ̃_Y (not b̃)
      # E[(Y - bM)²] = τ̃²_Y + (μ̃_Y - bM)²
      # ∂²/∂μ̃_Y² = -1/σ²_Y

      H_mu_Y <- -1 / sigma2_Y

      # V_theta for Model B has structure: (ã, σ̃²_M, μ̃_Y, τ̃²_Y)
      # Assuming block-diagonal: V_M for (ã, σ̃²_M), V_Y for (μ̃_Y, τ̃²_Y)
      # We need V_theta[3,3] for μ̃_Y

      correction_total <- correction_total + 0.5 * H_mu_Y * V_theta[3, 3]

    } else if (!M_obs && !Y_obs) {
      # Pattern 3: ∂²/∂ã² (same as Model A), ∂²/∂μ̃_Y²
      H_aa <- -X[i]^2 / sigma2_M
      H_mu_Y <- -1 / sigma2_Y

      correction_total <- correction_total + 0.5 * (H_aa * V_theta[1, 1] + H_mu_Y * V_theta[3, 3])
    }
  }

  correction_total
}


#' Hessian Correction for Model C (Saturated MVN) via Numerical Differentiation
#'
#' @keywords internal
compute_hessian_correction_C <- function(data_miss, theta_eval, theta_impute, V_theta) {

  # For Model C, use numerical Hessian via finite differences
  # Parameters: (μ_X, μ_M, μ_Y, Σ_11, Σ_22, Σ_33, Σ_12, Σ_13, Σ_23)
  # That's 9 parameters total

  # Extract parameter vector
  params <- c(
    theta_impute$mu_X,
    theta_impute$mu_M,
    theta_impute$mu_Y,
    theta_impute$Sigma[1, 1],  # σ²_X
    theta_impute$Sigma[2, 2],  # σ²_M
    theta_impute$Sigma[3, 3],  # σ²_Y
    theta_impute$Sigma[1, 2],  # σ_XM
    theta_impute$Sigma[1, 3],  # σ_XY
    theta_impute$Sigma[2, 3]   # σ_MY
  )

  # Function to convert parameter vector back to theta_impute list
  params_to_impute <- function(p) {
    Sigma_new <- matrix(0, 3, 3)
    Sigma_new[1, 1] <- p[4]
    Sigma_new[2, 2] <- p[5]
    Sigma_new[3, 3] <- p[6]
    Sigma_new[1, 2] <- Sigma_new[2, 1] <- p[7]
    Sigma_new[1, 3] <- Sigma_new[3, 1] <- p[8]
    Sigma_new[2, 3] <- Sigma_new[3, 2] <- p[9]

    list(mu_X = p[1], mu_M = p[2], mu_Y = p[3], Sigma = Sigma_new)
  }

  # Q-function wrapper for numerical differentiation
  Q_func <- function(p) {
    theta_imp_temp <- params_to_impute(p)
    compute_Q_analytical(data_miss, theta_eval, theta_imp_temp, model = "C")
  }

  # Compute numerical Hessian diagonal (second derivatives)
  h <- 1e-5  # Step size
  H_diag <- numeric(9)

  for (j in 1:9) {
    params_plus <- params_minus <- params
    step <- ifelse(abs(params[j]) > 1e-10, h * abs(params[j]), h)

    params_plus[j] <- params[j] + step
    params_minus[j] <- params[j] - step

    Q_plus <- Q_func(params_plus)
    Q_center <- Q_func(params)
    Q_minus <- Q_func(params_minus)

    H_diag[j] <- (Q_plus - 2 * Q_center + Q_minus) / (step^2)
  }

  # Hessian correction: (1/2) tr(H · V_theta) ≈ (1/2) Σ H_diag[i] * V_theta[i,i]
  correction_total <- 0.5 * sum(H_diag * diag(V_theta))

  correction_total
}
