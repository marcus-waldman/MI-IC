# ============================================================================
# Compute Metrics for Study 1
# ============================================================================
# Functions for computing:
# - Empirical bias (Q_bar - ell_com)
# - Theoretical bias (0.5 * tr(RIV))
# - Comparison metrics
#
# Supports two modes:
# - M = Inf: Analytical approach (exact, fast)
# - M = finite: Monte Carlo via brms (current approach)
# ============================================================================

# ============================================================================
# ANALYTICAL LIKELIHOOD FUNCTIONS (for M = Inf)
# ============================================================================

#' Complete-Data Negative Log-Likelihood for Mediation Model
#'
#' @param theta Numeric vector c(a, b, log_sigma2_M, log_sigma2_Y)
#' @param data data.frame with columns X, M, Y
#' @return Negative log-likelihood (scalar)
negloglik_complete <- function(theta, data) {
  a <- theta[1]
  b <- theta[2]
  sigma2_M <- exp(theta[3])  # Log transform for positivity
  sigma2_Y <- exp(theta[4])

  X <- data$X
  M <- data$M
  Y <- data$Y

  # M | X ~ N(a*X, sigma2_M)
  ll_M <- sum(dnorm(M, mean = a * X, sd = sqrt(sigma2_M), log = TRUE))

  # Y | M ~ N(b*M, sigma2_Y)
  ll_Y <- sum(dnorm(Y, mean = b * M, sd = sqrt(sigma2_Y), log = TRUE))

  -(ll_M + ll_Y)
}

#' Observed-Data Negative Log-Likelihood (FIML) for Mediation Model
#'
#' Computes likelihood by marginalizing over missing values
#'
#' @param theta Numeric vector c(a, b, log_sigma2_M, log_sigma2_Y)
#' @param data data.frame with columns X, M, Y (may contain NA)
#' @return Negative log-likelihood (scalar)
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
      # Pattern 1: Complete - full likelihood
      ll_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
              dnorm(Y[i], b * M[i], sqrt(sigma2_Y), log = TRUE)

    } else if (!M_obs && Y_obs) {
      # Pattern 2: M missing - marginalize over M
      # p(Y | X) = integral p(Y | M) p(M | X) dM
      # Y | X ~ N(a*b*X, b^2*sigma2_M + sigma2_Y)
      mean_Y_given_X <- a * b * X[i]
      var_Y_given_X <- b^2 * sigma2_M + sigma2_Y
      ll_i <- dnorm(Y[i], mean_Y_given_X, sqrt(var_Y_given_X), log = TRUE)

    } else if (M_obs && !Y_obs) {
      # Pattern 3: Y missing - just M likelihood
      ll_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE)

    } else {
      # Pattern 4: Both missing - no likelihood contribution
      ll_i <- 0
    }

    ll_total <- ll_total + ll_i
  }

  -ll_total
}

#' Analytical Q-Function for Mediation Model
#'
#' Computes E[ell_com(theta_eval) | Y_obs, theta_impute] analytically
#'
#' @param data_miss data.frame with missing values
#' @param theta_eval Numeric vector. Parameters to evaluate likelihood at
#' @param theta_impute Numeric vector. Parameters for imputation distribution
#' @return Q-function value (scalar)
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
      # Pattern 1: Complete - no expectation needed
      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             dnorm(Y[i], b * M[i], sqrt(sigma2_Y), log = TRUE)

    } else if (!M_obs && Y_obs) {
      # Pattern 2: M missing
      # Posterior: M | X, Y ~ N(mu_post, sigma2_post)
      sigma2_post <- (sigma2_M_tilde * sigma2_Y_tilde) /
                     (sigma2_Y_tilde + b_tilde^2 * sigma2_M_tilde)
      mu_post <- sigma2_post * (a_tilde * X[i] / sigma2_M_tilde +
                                b_tilde * Y[i] / sigma2_Y_tilde)

      # E[(M - a*X)^2] = Var(M) + (E[M] - a*X)^2
      E_M_resid_sq <- sigma2_post + (mu_post - a * X[i])^2
      # E[(Y - b*M)^2] = b^2*Var(M) + (Y - b*E[M])^2
      E_Y_resid_sq <- b^2 * sigma2_post + (Y[i] - b * mu_post)^2

      Q_i <- -0.5 * log(2 * pi * sigma2_M) - E_M_resid_sq / (2 * sigma2_M) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)

    } else if (M_obs && !Y_obs) {
      # Pattern 3: Y missing
      # Posterior: Y | M ~ N(b_tilde * M, sigma2_Y_tilde)
      E_Y <- b_tilde * M[i]
      Var_Y <- sigma2_Y_tilde
      E_Y_resid_sq <- Var_Y + (E_Y - b * M[i])^2

      Q_i <- dnorm(M[i], a * X[i], sqrt(sigma2_M), log = TRUE) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)

    } else {
      # Pattern 4: Both missing
      # Joint posterior: (M, Y) | X from imputation model
      E_M <- a_tilde * X[i]
      E_Y <- a_tilde * b_tilde * X[i]
      Var_M <- sigma2_M_tilde
      Var_Y <- b_tilde^2 * sigma2_M_tilde + sigma2_Y_tilde
      Cov_MY <- b_tilde * sigma2_M_tilde

      E_M_resid_sq <- Var_M + (E_M - a * X[i])^2
      # Var(Y - b*M) = Var(Y) + b^2*Var(M) - 2*b*Cov(M,Y)
      Var_Y_minus_bM <- Var_Y + b^2 * Var_M - 2 * b * Cov_MY
      E_Y_minus_bM <- E_Y - b * E_M
      E_Y_resid_sq <- Var_Y_minus_bM + E_Y_minus_bM^2

      Q_i <- -0.5 * log(2 * pi * sigma2_M) - E_M_resid_sq / (2 * sigma2_M) +
             -0.5 * log(2 * pi * sigma2_Y) - E_Y_resid_sq / (2 * sigma2_Y)
    }

    Q_total <- Q_total + Q_i
  }

  Q_total
}

#' Compute tr(RIV) Analytically Using FIML
#'
#' Uses Missing Information Principle: tr(RIV) = tr(I_com * I_obs^-1) - Q
#'
#' @param data_complete data.frame. Complete data
#' @param data_miss data.frame. Data with missing values
#' @return List with tr_RIV, theta_com, theta_obs, I_com, I_obs
compute_tr_RIV_analytical <- function(data_complete, data_miss) {

  theta_start <- c(0.5, 0.5, log(1), log(1))


  # Fit complete-data model
  fit_com <- optim(
    theta_start,
    negloglik_complete,
    data = data_complete,
    method = "BFGS",
    hessian = TRUE
  )

  # Fit observed-data model (FIML)
  fit_obs <- optim(
    theta_start,
    negloglik_observed,
    data = data_miss,
    method = "BFGS",
    hessian = TRUE
  )

  # Extract parameters (back-transform variances)
  theta_com <- fit_com$par
  theta_com[3:4] <- exp(theta_com[3:4])
  names(theta_com) <- c("a", "b", "sigma2_M", "sigma2_Y")

  theta_obs <- fit_obs$par
  theta_obs[3:4] <- exp(theta_obs[3:4])
  names(theta_obs) <- c("a", "b", "sigma2_M", "sigma2_Y")

  # Information matrices (Hessian at MLE)
  I_com <- fit_com$hessian
  I_obs <- fit_obs$hessian

  # tr(RIV) via Missing Information Principle
  Q_params <- 4
  ratio_matrix <- I_com %*% solve(I_obs)
  tr_RIV <- sum(diag(ratio_matrix)) - Q_params

  # Log-likelihoods
  ell_com <- -fit_com$value
  ell_obs <- -fit_obs$value

  list(
    tr_RIV = tr_RIV,
    theta_com = theta_com,
    theta_obs = theta_obs,
    I_com = I_com,
    I_obs = I_obs,
    ell_com = ell_com,
    ell_obs = ell_obs,
    fit_com = fit_com,
    fit_obs = fit_obs
  )
}

#' Compute Study 1 Metrics - Analytical Version (M = Inf)
#'
#' Computes all metrics analytically without Monte Carlo sampling
#'
#' @param data_complete data.frame. Complete data
#' @param data_miss data.frame. Data with missing values
#' @param pop_params List. Population parameters
#' @return List with all computed metrics
compute_study1_metrics_analytical <- function(data_complete, data_miss, pop_params) {

  # 1. Compute tr(RIV) and get MLEs
  riv_result <- compute_tr_RIV_analytical(data_complete, data_miss)

  theta_com <- riv_result$theta_com
  theta_obs <- riv_result$theta_obs

  # 2. Compute Q at theta_obs (with imputation from theta_obs)
  Q_at_obs <- compute_Q_analytical(data_miss, theta_obs, theta_obs)

  # 3. Complete-data log-likelihood at complete-data MLE
  ell_com_at_com <- riv_result$ell_com

  # 4. Complete-data log-likelihood at FIML MLE
  ell_com_at_obs <- -negloglik_complete(
    c(theta_obs[1], theta_obs[2], log(theta_obs[3]), log(theta_obs[4])),
    data_complete
  )

  # 5. Bias decomposition
  term_A1 <- Q_at_obs - ell_com_at_obs  # Imputation bias (should ≈ tr(RIV))
  term_A2 <- ell_com_at_obs - ell_com_at_com  # Estimation mismatch (should ≈ -0.5*tr(RIV))
  empirical_bias <- term_A1 + term_A2  # Total (should ≈ 0.5*tr(RIV))
  theoretical_bias <- 0.5 * riv_result$tr_RIV

  # 6. Parameter bias
  theta_true <- c(
    a = pop_params$beta_xm,
    b = pop_params$beta_my,
    sigma2_M = pop_params$sigma2_m,
    sigma2_Y = pop_params$sigma2_y
  )

  a_bias <- theta_obs["a"] - theta_true["a"]
  b_bias <- theta_obs["b"] - theta_true["b"]
  indirect_true <- theta_true["a"] * theta_true["b"]
  indirect_est <- theta_obs["a"] * theta_obs["b"]
  indirect_bias <- indirect_est - indirect_true

  # 7. Approximate SEs from inverse Hessian (for coverage)
  # Note: Hessian is at log scale for variances, so SEs for a,b are direct
  vcov_obs <- tryCatch(solve(riv_result$I_obs), error = function(e) NULL)

  if (!is.null(vcov_obs)) {
    se_a <- sqrt(vcov_obs[1, 1])
    se_b <- sqrt(vcov_obs[2, 2])

    # Coverage for a and b
    ci_a_lower <- theta_obs["a"] - 1.96 * se_a
    ci_a_upper <- theta_obs["a"] + 1.96 * se_a
    a_covered <- (theta_true["a"] >= ci_a_lower) & (theta_true["a"] <= ci_a_upper)

    ci_b_lower <- theta_obs["b"] - 1.96 * se_b
    ci_b_upper <- theta_obs["b"] + 1.96 * se_b
    b_covered <- (theta_true["b"] >= ci_b_lower) & (theta_true["b"] <= ci_b_upper)

    # Delta method for indirect effect SE
    grad_indirect <- c(theta_obs["b"], theta_obs["a"], 0, 0)
    se_indirect <- sqrt(t(grad_indirect) %*% vcov_obs %*% grad_indirect)
    ci_indirect_lower <- indirect_est - 1.96 * se_indirect
    ci_indirect_upper <- indirect_est + 1.96 * se_indirect
    indirect_covered <- (indirect_true >= ci_indirect_lower) & (indirect_true <= ci_indirect_upper)
  } else {
    se_a <- se_b <- se_indirect <- NA
    a_covered <- b_covered <- indirect_covered <- NA
  }

  list(
    # Bias metrics (primary outcomes)
    Q_bar = Q_at_obs,  # Named Q_bar for compatibility
    ell_com = ell_com_at_com,
    empirical_bias = empirical_bias,
    theoretical_bias = theoretical_bias,
    bias_difference = empirical_bias - theoretical_bias,
    bias_ratio = empirical_bias / theoretical_bias,

    # RIV components
    tr_RIV = riv_result$tr_RIV,
    param_RIV = NULL,  # Not computed in analytical version

    # Level 1 Diagnostics: Bias Decomposition
    ell_com_at_pooled = ell_com_at_obs,  # theta_obs plays role of theta_pooled
    term_A1 = term_A1,
    term_A2 = term_A2,
    term_A1_ratio = term_A1 / riv_result$tr_RIV,
    term_A2_ratio = term_A2 / (-0.5 * riv_result$tr_RIV),
    theta_pooled = theta_obs,
    theta_complete = theta_com,
    theta_true = theta_true,

    # Parameter recovery
    a_est = as.numeric(theta_obs["a"]),
    a_se = as.numeric(se_a),
    a_bias = as.numeric(a_bias),
    a_covered = as.logical(a_covered),

    b_est = as.numeric(theta_obs["b"]),
    b_se = as.numeric(se_b),
    b_bias = as.numeric(b_bias),
    b_covered = as.logical(b_covered),

    indirect_est = as.numeric(indirect_est),
    indirect_se = as.numeric(se_indirect),
    indirect_bias = as.numeric(indirect_bias),
    indirect_covered = as.logical(indirect_covered),

    # Metadata
    M = Inf,
    n = nrow(data_complete),

    # Analytical-specific
    method = "analytical"
  )
}

# ============================================================================
# ORIGINAL FUNCTIONS (for M = finite)
# ============================================================================

#' Compute Log-Likelihood at Specific Parameters
#'
#' Evaluates the log-likelihood of a dataset at a specific parameter vector.
#' Uses fixed-parameter fitting approach (setting lower=upper bounds).
#'
#' @param data data.frame. Dataset to evaluate (typically complete data)
#' @param params Named numeric vector. Parameter values to evaluate at
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @return Numeric. Log-likelihood value
compute_loglik_at_params <- function(data, params, model_type = "full_mediation") {

  # Get model syntax (labels are already included in the syntax)
  syntax <- get_mediation_syntax(model_type)

  # Create template fit to get parameter table
  fit_template <- lavaan::sem(
    syntax,
    data = data,
    do.fit = FALSE,
    meanstructure = TRUE
  )
  partable <- lavaan::parameterTable(fit_template)

  # Identify free parameters
  free_idx <- partable$free > 0

  # Construct parameter names matching lavaan's conventions
  # Use label if available, otherwise construct from lhs op rhs
  param_names <- ifelse(
    partable$label[free_idx] != "",
    partable$label[free_idx],
    paste0(partable$lhs[free_idx], partable$op[free_idx], partable$rhs[free_idx])
  )

  # Match provided parameter values to parameter table
  # Handle potential name mismatches
  param_values <- numeric(length(param_names))
  for (i in seq_along(param_names)) {
    if (param_names[i] %in% names(params)) {
      param_values[i] <- params[param_names[i]]
    } else {
      # Try alternative naming conventions
      alt_name <- paste0(partable$lhs[free_idx][i], partable$op[free_idx][i],
                         partable$rhs[free_idx][i])
      if (alt_name %in% names(params)) {
        param_values[i] <- params[alt_name]
      } else {
        stop(sprintf("Parameter '%s' not found in provided params vector", param_names[i]))
      }
    }
  }

  # Fix parameters by setting lower=upper bounds
  partable_fixed <- partable
  partable_fixed$ustart[free_idx] <- param_values
  partable_fixed$lower[free_idx] <- param_values
  partable_fixed$upper[free_idx] <- param_values

  # Fit with fixed parameters (no SE computation needed)
  fit_fixed <- lavaan::sem(
    model = partable_fixed,
    data = data,
    do.fit = TRUE,
    se = "none",
    test = "standard",
    verbose = FALSE
  )

  # Extract and return log-likelihood
  as.numeric(logLik(fit_fixed))
}

#' Compute Study 1 Metrics
#'
#' Computes all metrics for one replication
#'
#' @param data_complete data.frame. Complete data (before missingness)
#' @param imputed_datasets List. Imputed datasets from brms
#' @param fit_mi lavaan.mi fit object
#' @param riv_result List from compute_RIV_from_mi()
#' @param pop_params List. Population parameters
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @return List with all computed metrics
compute_study1_metrics <- function(data_complete,
                                    imputed_datasets,
                                    fit_mi,
                                    riv_result,
                                    pop_params,
                                    model_type = "full_mediation") {

  M <- length(imputed_datasets)

  # 1. Compute Q_bar: average log-likelihood across imputed datasets
  #    EVALUATED AT THE POOLED RUBIN'S RULES ESTIMATE
  #    Q̄ = (1/M)Σℓ(θ̄_pooled | Y^(m))
  logliks_imputed <- get_loglik_per_imputation(imputed_datasets, fit_mi, model_type)
  Q_bar <- mean(logliks_imputed, na.rm = TRUE)

  # 2. Compute ell_com: complete-data log-likelihood at complete-data MLE
  fit_complete <- fit_mediation_complete(data_complete, model_type)
  ell_com <- lavaan::fitMeasures(fit_complete, "logl")

  # 3. Empirical bias
  empirical_bias <- Q_bar - ell_com

  # 4. Theoretical bias
  theoretical_bias <- 0.5 * riv_result$tr_RIV

  # 5. Difference between empirical and theoretical
  bias_difference <- empirical_bias - theoretical_bias

  # 6. Get pooled parameter estimates and check against truth
  pooled_params <- get_pooled_estimates(fit_mi)
  pooled_params <- check_parameter_coverage(pooled_params, pop_params)

  # 7. Compute per-parameter metrics
  param_a <- pooled_params[pooled_params$label == "a", ]
  param_b <- pooled_params[pooled_params$label == "b", ]
  param_indirect <- pooled_params[pooled_params$label == "indirect", ]

  # ============================================================================
  # LEVEL 1 DIAGNOSTICS: Bias Decomposition (Approach A)
  # ============================================================================
  # Extract parameter vectors
  theta_pooled <- coef(fit_mi)
  theta_complete <- coef(fit_complete)
  theta_true <- c(
    a = pop_params$beta_xm,
    b = pop_params$beta_my,
    "M~~M" = pop_params$sigma2_m,
    "Y~~Y" = pop_params$sigma2_y,
    "M~1" = 0,   # Intercepts assumed 0 in data generation
    "Y~1" = 0
  )

  # Compute ell_com at pooled estimate (KEY missing quantity)
  # This enables decomposition: Total Bias = Term A1 + Term A2
  ell_com_at_pooled <- compute_loglik_at_params(
    data = data_complete,
    params = theta_pooled,
    model_type = model_type
  )

  # Term decomposition
  # Term A1 (Imputation Bias): Q̄(θ̄_pooled) - ℓ_com(θ̄_pooled) ≈ tr(RIV)
  # Term A2 (Estimation Mismatch): ℓ_com(θ̄_pooled) - ℓ_com(θ̂_com) ≈ -½tr(RIV)
  term_A1 <- Q_bar - ell_com_at_pooled
  term_A2 <- ell_com_at_pooled - ell_com

  # Ratios for diagnostics (should both ≈ 1.0 in expectation)
  term_A1_ratio <- term_A1 / riv_result$tr_RIV
  term_A2_ratio <- term_A2 / (-0.5 * riv_result$tr_RIV)

  # Verify arithmetic: term_A1 + term_A2 should equal empirical_bias
  term_sum_check <- abs((term_A1 + term_A2) - empirical_bias)
  if (term_sum_check > 1e-6) {
    warning(sprintf("Term decomposition check failed: |sum - total| = %.2e", term_sum_check))
  }

  # ============================================================================
  # LEVEL 2 DIAGNOSTICS: Taylor Approximation Check
  # ============================================================================
  # Check quality of quadratic approximation for Term A2
  # Term A2 ≈ -½(θ̄_pooled - θ̂_com)' × J_com × (θ̄_pooled - θ̂_com)
  # where J_com is the observed information matrix at θ̂_com

  # Get observed information matrix (inverse of covariance matrix)
  vcov_complete <- vcov(fit_complete)
  J_com <- tryCatch(
    solve(vcov_complete),
    error = function(e) {
      warning("Could not invert vcov matrix for J_com")
      return(NULL)
    }
  )

  if (!is.null(J_com)) {
    theta_diff <- theta_pooled - theta_complete
    term_A2_taylor <- -0.5 * as.numeric(t(theta_diff) %*% J_com %*% theta_diff)
  } else {
    term_A2_taylor <- NA
  }

  # ============================================================================
  # LEVEL 3 DIAGNOSTICS: Information Matrices and Bilinear Form
  # ============================================================================
  # Compute information matrices to verify Term 1 derivation
  # Term 1 ≈ (θ̃_obs - θ*)' × ℐ_mis|obs × (θ̂_obs - θ*)
  # where ℐ_mis|obs = ℐ_com - ℐ_obs

  # Information matrices (inverse of covariance matrices)
  I_com_matrix <- J_com  # Already computed above
  I_obs_matrix <- tryCatch({
    vcov_obs <- vcov(fit_mi)
    solve(vcov_obs)
  }, error = function(e) {
    warning("Could not compute I_obs from fit_mi vcov")
    return(NULL)
  })

  # Missing information matrix
  if (!is.null(I_com_matrix) && !is.null(I_obs_matrix)) {
    I_mis_obs_matrix <- I_com_matrix - I_obs_matrix

    # Step 7 bilinear form: (θ̄_pooled - θ*)' × ℐ_mis|obs × (θ̄_pooled - θ*)
    # Using θ̄_pooled as the observed-data estimate (from Rubin's rules)
    theta_pooled_deviation <- theta_pooled - theta_true
    term1_bilinear <- as.numeric(
      t(theta_pooled_deviation) %*% I_mis_obs_matrix %*% theta_pooled_deviation
    )
  } else {
    I_mis_obs_matrix <- NULL
    term1_bilinear <- NA
  }

  return(list(
    # Bias metrics (primary outcomes)
    Q_bar = Q_bar,
    ell_com = ell_com,
    empirical_bias = empirical_bias,
    theoretical_bias = theoretical_bias,
    bias_difference = bias_difference,
    bias_ratio = empirical_bias / theoretical_bias,

    # RIV components
    tr_RIV = riv_result$tr_RIV,
    param_RIV = riv_result$param_RIV,

    # Level 1 Diagnostics: Bias Decomposition
    ell_com_at_pooled = ell_com_at_pooled,
    term_A1 = term_A1,
    term_A2 = term_A2,
    term_A1_ratio = term_A1_ratio,
    term_A2_ratio = term_A2_ratio,
    theta_pooled = theta_pooled,
    theta_complete = theta_complete,
    theta_true = theta_true,

    # Level 2 Diagnostics: Taylor Approximation
    term_A2_taylor = term_A2_taylor,

    # Level 3 Diagnostics: Information Matrices and Bilinear Form
    term1_bilinear = term1_bilinear,
    tr_I_com = if (!is.null(I_com_matrix)) sum(diag(I_com_matrix)) else NA,
    tr_I_obs = if (!is.null(I_obs_matrix)) sum(diag(I_obs_matrix)) else NA,
    tr_I_mis_obs = if (!is.null(I_mis_obs_matrix)) sum(diag(I_mis_obs_matrix)) else NA,

    # Parameter recovery
    a_est = if (nrow(param_a) > 0) param_a$est else NA,
    a_se = if (nrow(param_a) > 0) param_a$se else NA,
    a_bias = if (nrow(param_a) > 0) param_a$bias else NA,
    a_covered = if (nrow(param_a) > 0) param_a$covered else NA,

    b_est = if (nrow(param_b) > 0) param_b$est else NA,
    b_se = if (nrow(param_b) > 0) param_b$se else NA,
    b_bias = if (nrow(param_b) > 0) param_b$bias else NA,
    b_covered = if (nrow(param_b) > 0) param_b$covered else NA,

    indirect_est = if (nrow(param_indirect) > 0) param_indirect$est else NA,
    indirect_se = if (nrow(param_indirect) > 0) param_indirect$se else NA,
    indirect_bias = if (nrow(param_indirect) > 0) param_indirect$bias else NA,
    indirect_covered = if (nrow(param_indirect) > 0) param_indirect$covered else NA,

    # Metadata
    M = M,
    n = nrow(data_complete)
  ))
}


#' Summarize Metrics Across Replications
#'
#' Aggregates metrics across all replications for a condition
#'
#' @param results_list List of results from compute_study1_metrics()
#' @return data.frame with summary statistics
summarize_condition_results <- function(results_list) {

  # Filter out NULL/failed results
  results_list <- results_list[!sapply(results_list, is.null)]
  n_reps <- length(results_list)

  if (n_reps == 0) {
    return(data.frame(n_reps = 0, note = "All replications failed"))
  }

  # Extract vectors
  empirical_bias <- sapply(results_list, `[[`, "empirical_bias")
  theoretical_bias <- sapply(results_list, `[[`, "theoretical_bias")
  bias_difference <- sapply(results_list, `[[`, "bias_difference")
  tr_RIV <- sapply(results_list, `[[`, "tr_RIV")

  a_bias <- sapply(results_list, `[[`, "a_bias")
  b_bias <- sapply(results_list, `[[`, "b_bias")
  indirect_bias <- sapply(results_list, `[[`, "indirect_bias")

  a_covered <- sapply(results_list, `[[`, "a_covered")
  b_covered <- sapply(results_list, `[[`, "b_covered")
  indirect_covered <- sapply(results_list, `[[`, "indirect_covered")

  summary_df <- data.frame(
    n_reps = n_reps,

    # Bias formula validation
    mean_empirical_bias = mean(empirical_bias, na.rm = TRUE),
    sd_empirical_bias = sd(empirical_bias, na.rm = TRUE),
    mean_theoretical_bias = mean(theoretical_bias, na.rm = TRUE),
    sd_theoretical_bias = sd(theoretical_bias, na.rm = TRUE),
    mean_bias_difference = mean(bias_difference, na.rm = TRUE),
    sd_bias_difference = sd(bias_difference, na.rm = TRUE),

    # Correlation between empirical and theoretical
    cor_emp_theo = cor(empirical_bias, theoretical_bias, use = "complete.obs"),

    # tr(RIV)
    mean_tr_RIV = mean(tr_RIV, na.rm = TRUE),
    sd_tr_RIV = sd(tr_RIV, na.rm = TRUE),

    # Parameter recovery
    mean_a_bias = mean(a_bias, na.rm = TRUE),
    mean_b_bias = mean(b_bias, na.rm = TRUE),
    mean_indirect_bias = mean(indirect_bias, na.rm = TRUE),

    # Coverage
    coverage_a = mean(a_covered, na.rm = TRUE),
    coverage_b = mean(b_covered, na.rm = TRUE),
    coverage_indirect = mean(indirect_covered, na.rm = TRUE)
  )

  return(summary_df)
}


#' Create Results Table for All Conditions
#'
#' Combines results across all conditions into a single table
#'
#' @param all_results List with structure: all_results[[mechanism]][[model]] = results_list
#' @param config Configuration list
#' @return data.frame with one row per condition
create_results_table <- function(all_results, config) {

  results_rows <- list()

  for (mechanism in config$mechanisms) {
    for (model in config$imputation_models) {
      if (!is.null(all_results[[mechanism]][[model]])) {
        summary <- summarize_condition_results(all_results[[mechanism]][[model]])
        summary$mechanism <- mechanism
        summary$imputation_model <- model
        results_rows[[paste(mechanism, model, sep = "_")]] <- summary
      }
    }
  }

  results_df <- do.call(rbind, results_rows)

  # Reorder columns
  col_order <- c("mechanism", "imputation_model", "n_reps",
                 setdiff(names(results_df), c("mechanism", "imputation_model", "n_reps")))
  results_df <- results_df[, col_order]

  rownames(results_df) <- NULL

  return(results_df)
}


#' Check Bias Formula Validity
#'
#' Evaluates whether empirical bias ≈ 0.5 * tr(RIV) for each condition
#'
#' @param results_df data.frame from create_results_table()
#' @param tolerance Numeric. Acceptable relative difference (default: 0.2 = 20%)
#' @return data.frame with validity assessment
check_bias_formula <- function(results_df, tolerance = 0.2) {

  results_df$rel_diff <- with(results_df,
    (mean_empirical_bias - mean_theoretical_bias) / abs(mean_theoretical_bias)
  )

  results_df$formula_holds <- abs(results_df$rel_diff) < tolerance

  # Expected pattern:
  # - matched, saturated: formula should hold

  # - uncongenial: formula should NOT hold
  results_df$expected_to_hold <- results_df$imputation_model != "uncongenial"
  results_df$matches_expectation <- results_df$formula_holds == results_df$expected_to_hold

  return(results_df)
}


#' Check Parameter Recovery Quality
#'
#' Evaluates whether pooled parameter estimates (from Rubin's rules via sem.mi)
#' show expected bias/unbiasedness patterns based on mechanism and imputation model.
#'
#' Expected patterns:
#'   - Unbiased: (MCAR or MAR) + (matched or saturated) - ignorable + congenial
#'   - Biased:
#'       - MNAR (any model): non-ignorable missingness
#'       - Uncongenial (any mechanism): imputation model omits M from Y
#'
#' @param results_df data.frame from create_results_table()
#' @param bias_tolerance Numeric. Max |mean bias| to consider "unbiased" (default: 0.05)
#' @param coverage_tolerance Numeric. Coverage range around 0.95 (default: 0.10, so 0.85-1.0)
#' @return data.frame with parameter recovery assessment
check_parameter_recovery <- function(results_df,
                                      bias_tolerance = 0.05,
                                      coverage_tolerance = 0.10) {


  # Determine expected unbiasedness for each condition
  # Unbiased if: (mechanism is MCAR or MAR) AND (model is matched or saturated)
  results_df$expect_unbiased <- with(results_df,
    (mechanism %in% c("MCAR", "MAR")) & (imputation_model %in% c("matched", "saturated"))
  )

  # Check actual bias for key parameters
  # Parameter 'a' (X -> M): Should be unbiased in all ignorable conditions
  #   - Only biased under MNAR if missingness in M depends on M itself
  results_df$a_unbiased <- abs(results_df$mean_a_bias) < bias_tolerance

  # Parameter 'b' (M -> Y): Most sensitive to congeniality
  #   - Uncongenial model omits M from Y imputation, so b estimates are biased
  results_df$b_unbiased <- abs(results_df$mean_b_bias) < bias_tolerance

  # Indirect effect (a*b): Combines both sources of bias

  results_df$indirect_unbiased <- abs(results_df$mean_indirect_bias) < bias_tolerance

  # Coverage checks (should be ~95% when unbiased)
  nominal_coverage <- 0.95
  results_df$a_coverage_ok <- abs(results_df$coverage_a - nominal_coverage) < coverage_tolerance
  results_df$b_coverage_ok <- abs(results_df$coverage_b - nominal_coverage) < coverage_tolerance
  results_df$indirect_coverage_ok <- abs(results_df$coverage_indirect - nominal_coverage) < coverage_tolerance

  # Overall parameter recovery assessment
  # For conditions expected to be unbiased: check that they ARE unbiased with good coverage
  # For conditions expected to be biased: check that they show bias (or poor coverage)
  results_df$param_recovery_matches <- with(results_df, {
    ifelse(expect_unbiased,
           # If expected unbiased: should be unbiased with good coverage
           b_unbiased & b_coverage_ok,
           # If expected biased: should show bias OR poor coverage
           !b_unbiased | !b_coverage_ok)
  })

  return(results_df)
}


#' Print Parameter Recovery Summary
#'
#' Pretty-prints parameter recovery results with expected vs actual patterns
#'
#' @param results_df data.frame from check_parameter_recovery()
print_parameter_recovery_summary <- function(results_df) {

  cat("=== Parameter Recovery Quality Check ===\n")
  cat("(Pooled estimates from Rubin's rules via sem.mi)\n\n")

  cat("Expected patterns:\n")
  cat("  UNBIASED: (MCAR or MAR) + (matched or saturated)\n")
  cat("  BIASED:   MNAR (any model) OR uncongenial (any mechanism)\n\n")

  # Table header
  cat(sprintf("%-6s %-11s | %7s %7s %7s | %6s %6s %6s | %8s\n",
              "Mech", "Imputation", "a_bias", "b_bias", "ab_bias",
              "cov_a", "cov_b", "cov_ab", "Expected"))
  cat(paste(rep("-", 85), collapse = ""), "\n")

  for (i in seq_len(nrow(results_df))) {
    row <- results_df[i, ]
    expected_label <- ifelse(row$expect_unbiased, "Unbiased", "Biased")

    cat(sprintf("%-6s %-11s | %+7.3f %+7.3f %+7.3f | %5.1f%% %5.1f%% %5.1f%% | %8s %s\n",
                row$mechanism,
                row$imputation_model,
                row$mean_a_bias,
                row$mean_b_bias,
                row$mean_indirect_bias,
                row$coverage_a * 100,
                row$coverage_b * 100,
                row$coverage_indirect * 100,
                expected_label,
                ifelse(row$param_recovery_matches, "", " *MISMATCH*")))
  }

  # Summary
  n_match <- sum(results_df$param_recovery_matches, na.rm = TRUE)
  n_total <- nrow(results_df)
  cat("\n")
  cat(sprintf("Parameter recovery: %d/%d conditions match expected pattern\n", n_match, n_total))

  # Highlight any mismatches
  mismatches <- results_df[!results_df$param_recovery_matches, ]
  if (nrow(mismatches) > 0) {
    cat("\nMISMATCHES (investigate these):\n")
    for (i in seq_len(nrow(mismatches))) {
      row <- mismatches[i, ]
      cat(sprintf("  - %s × %s: expected %s, got b_bias=%.3f, coverage=%.1f%%\n",
                  row$mechanism, row$imputation_model,
                  ifelse(row$expect_unbiased, "unbiased", "biased"),
                  row$mean_b_bias, row$coverage_b * 100))
    }
  }
}


#' Print Results Summary
#'
#' Pretty-prints the results table with bias formula assessment
#'
#' @param results_df data.frame from check_bias_formula()
print_results_summary <- function(results_df) {

  cat("=== Study 1 Results Summary ===\n\n")

  cat("Bias Formula Validation: Empirical Bias ≈ 0.5 × tr(RIV)\n\n")

  for (mech in unique(results_df$mechanism)) {
    cat(sprintf("Mechanism: %s\n", mech))
    cat(sprintf("%-12s %10s %10s %10s %8s %10s\n",
                "Imputation", "Emp.Bias", "Theo.Bias", "Diff", "RelDiff", "Holds?"))
    cat(paste(rep("-", 70), collapse = ""), "\n")

    subset_df <- results_df[results_df$mechanism == mech, ]
    for (i in seq_len(nrow(subset_df))) {
      row <- subset_df[i, ]
      cat(sprintf("%-12s %10.3f %10.3f %10.3f %7.1f%% %10s\n",
                  row$imputation_model,
                  row$mean_empirical_bias,
                  row$mean_theoretical_bias,
                  row$mean_bias_difference,
                  row$rel_diff * 100,
                  ifelse(row$formula_holds, "YES", "NO")))
    }
    cat("\n")
  }

  # Overall assessment
  n_correct <- sum(results_df$matches_expectation, na.rm = TRUE)
  n_total <- nrow(results_df)

  cat(sprintf("Overall: %d/%d conditions match expected pattern\n", n_correct, n_total))
  cat("\nExpected:\n")
  cat("  - matched, saturated: formula HOLDS (congenial)\n")
  cat("  - uncongenial: formula DOES NOT hold (not congenial)\n")
}


#' Print Full Results Summary
#'
#' Prints both bias formula validation AND parameter recovery quality check
#'
#' @param results_df data.frame with all checks applied
print_full_results_summary <- function(results_df) {
  print_results_summary(results_df)
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  print_parameter_recovery_summary(results_df)
}
