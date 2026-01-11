# ============================================================================
# Compute Metrics for Study 1
# ============================================================================
# Uses analytical functions from 03_analytical_functions.R

#' Compute Bias Metrics for One Replication
#'
#' @param data_complete Complete data
#' @param data_miss Data with missing values
#' @param true_params List with a, b, sigma2_M, sigma2_Y
#' @return List with metrics
compute_metrics <- function(data_complete, data_miss, true_params) {

  # Fit models using analytical functions
  theta_start <- c(0.5, 0.5, log(1), log(1))

  fit_com <- optim(
    theta_start,
    negloglik_complete,
    data = data_complete,
    method = "BFGS",
    hessian = TRUE
  )

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

  theta_obs <- fit_obs$par
  theta_obs[3:4] <- exp(theta_obs[3:4])

  # Compute tr(RIV)
  I_com <- fit_com$hessian
  I_obs <- fit_obs$hessian
  tr_RIV <- sum(diag(I_com %*% solve(I_obs))) - 4

  # Compute Q-function
  Q_at_obs <- compute_Q_analytical(data_miss, theta_obs, theta_obs)

  # Compute likelihoods
  ell_com_at_com <- -fit_com$value
  ell_com_at_obs <- -negloglik_complete(c(theta_obs[1:2], log(theta_obs[3:4])), data_complete)

  # Bias decomposition
  term1 <- Q_at_obs - ell_com_at_obs  # Imputation bias (should ≈ tr(RIV))
  term2 <- ell_com_at_obs - ell_com_at_com  # Estimation mismatch (should ≈ -0.5*tr(RIV))
  total_bias <- term1 + term2  # Total (should ≈ 0.5*tr(RIV))

  list(
    # Bias metrics
    term1 = term1,
    term2 = term2,
    total_bias = total_bias,
    tr_RIV = tr_RIV,
    theoretical_bias = 0.5 * tr_RIV,

    # Ratios (should all ≈ 1.0)
    ratio1 = term1 / tr_RIV,
    ratio2 = term2 / (-0.5 * tr_RIV),
    ratio_total = total_bias / (0.5 * tr_RIV),

    # Parameter estimates
    a_est = theta_obs[1],
    b_est = theta_obs[2],
    a_bias = theta_obs[1] - true_params$a,
    b_bias = theta_obs[2] - true_params$b
  )
}
