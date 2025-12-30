# ============================================================================
# Compute Metrics for Imputation Bias Validation
# ============================================================================
# Functions for computing:
#   - Rubin's variance components (W, B, RIV)
#   - Q-function (expected complete-data log-likelihood)
#   - Complete-data log-likelihood
#   - Imputation bias
#   - Selection entropy
# ============================================================================

#' Compute Rubin's Variance Components and RIV
#'
#' Fits model to each completed dataset and computes W, B, and RIV
#'
#' @param completed_datasets List of M completed datasets (each n x p matrix)
#' @param structure Character. "CS", "Toeplitz", or "Unstructured"
#' @param regularize Logical. Whether to regularize W^{-1} for numerical stability
#' @return List with components:
#'   \item{theta_hats}{Matrix with M rows, each row is θ̂^(m)}
#'   \item{vcov_hats}{List of M variance-covariance matrices}
#'   \item{W}{Within-imputation variance matrix}
#'   \item{B}{Between-imputation variance matrix}
#'   \item{RIV}{Relative increase in variance matrix}
#'   \item{tr_RIV}{Trace of RIV (scalar)}
#'   \item{M}{Number of imputations}
#' @details
#'   W = (1/M) Σ vcov(θ̂^(m))
#'   B = var(θ̂^(1), ..., θ̂^(M))
#'   RIV = (1 + 1/M) W^{-1} B
#'   tr(RIV) is the key quantity for imputation bias
#' @examples
#'   result <- compute_rubin_variance(completed_datasets, structure = "CS")
#'   tr_RIV <- result$tr_RIV
compute_rubin_variance <- function(completed_datasets, structure, regularize = TRUE) {

  # Source fit_models.R if not loaded
  if (!exists("fit_lavaan_model", mode = "function")) {
    source("simulations/utils/fit_models.R")
  }

  M <- length(completed_datasets)
  if (M < 2) stop("Need at least 2 imputations to compute variance components")

  # Initialize storage for parameter estimates and vcov matrices
  theta_hats <- list()
  vcov_hats <- list()
  converged_all <- TRUE

  # Fit model to each completed dataset
  for (m in 1:M) {
    fit_m <- fit_lavaan_model(
      data = completed_datasets[[m]],
      structure = structure,
      return_vcov = TRUE
    )

    if (!fit_m$converged) {
      warning(paste("Model did not converge for imputation", m))
      converged_all <- FALSE
    }

    theta_hats[[m]] <- fit_m$theta_hat
    vcov_hats[[m]] <- fit_m$vcov_theta
  }

  # Filter out failed fits (where vcov is NULL)
  vcov_is_null <- sapply(vcov_hats, is.null)
  n_null <- sum(vcov_is_null)

  if (n_null > 0) {
    warning(sprintf("%d out of %d imputations failed to produce vcov matrix; using %d successful fits",
                    n_null, M, M - n_null))

    # Keep only successful fits
    theta_hats <- theta_hats[!vcov_is_null]
    vcov_hats <- vcov_hats[!vcov_is_null]
    converged_all <- FALSE
  }

  # Check if we have at least 2 successful fits
  M_successful <- length(theta_hats)
  if (M_successful < 2) {
    stop(sprintf("Only %d out of %d imputations succeeded; need at least 2 for variance computation",
                 M_successful, M))
  }

  # Convert theta_hats to matrix (M_successful x Q)
  theta_matrix <- do.call(rbind, theta_hats)
  Q <- ncol(theta_matrix)

  # Check that all remaining vcov matrices have correct dimensions
  vcov_dims <- sapply(vcov_hats, nrow)

  if (length(unique(vcov_dims)) > 1) {
    stop(sprintf("vcov matrices have inconsistent dimensions: %s", paste(unique(vcov_dims), collapse = ", ")))
  }

  if (vcov_dims[1] != Q) {
    stop(sprintf("vcov dimension (%d) does not match theta dimension (%d)", vcov_dims[1], Q))
  }

  # Compute W: within-imputation variance (average of vcov matrices)
  # W = (1/M_successful) Σ vcov(θ̂^(m))
  W <- Reduce(`+`, vcov_hats) / M_successful

  # Verify W dimensions
  if (nrow(W) != Q || ncol(W) != Q) {
    stop(sprintf("W has wrong dimensions: %d x %d, expected %d x %d", nrow(W), ncol(W), Q, Q))
  }

  # Compute B: between-imputation variance
  # B = var(θ̂^(1), ..., θ̂^(M_successful))
  if (M_successful > 1) {
    B <- cov(theta_matrix)
  } else {
    B <- matrix(0, nrow = Q, ncol = Q)
  }

  # Verify B dimensions
  if (nrow(B) != Q || ncol(B) != Q) {
    stop(sprintf("B has wrong dimensions: %d x %d, expected %d x %d", nrow(B), ncol(B), Q, Q))
  }

  # Compute RIV = (1 + 1/M_successful) W^{-1} B
  # Need to invert W carefully
  W_inv <- tryCatch(
    {
      if (regularize) {
        # Add small ridge to diagonal for numerical stability
        ridge <- 1e-6 * mean(diag(W))
        solve(W + diag(ridge, nrow = Q))
      } else {
        solve(W)
      }
    },
    error = function(e) {
      warning("W matrix singular; using pseudoinverse")
      # Use Moore-Penrose pseudoinverse
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' required for pseudoinverse")
      }
      MASS::ginv(W)
    }
  )

  # Verify W_inv dimensions before multiplication
  if (nrow(W_inv) != Q || ncol(W_inv) != Q) {
    stop(sprintf("W_inv has wrong dimensions: %d x %d, expected %d x %d", nrow(W_inv), ncol(W_inv), Q, Q))
  }

  RIV <- (1 + 1/M_successful) * W_inv %*% B

  # Compute trace of RIV (the key quantity)
  tr_RIV <- sum(diag(RIV))

  return(list(
    theta_hats = theta_matrix,
    vcov_hats = vcov_hats,
    W = W,
    B = B,
    RIV = RIV,
    tr_RIV = tr_RIV,
    M = M,
    M_successful = M_successful,
    converged_all = converged_all
  ))
}


#' Compute Q-function (MI Estimate of Complete-Data Log-Likelihood)
#'
#' Evaluates Q̄_MI(θ̂) = (1/M) Σ ℓ_com(θ̂ | Z^obs, Z̃^mis,m)
#'
#' @param theta_hat Numeric vector or list. Parameter estimate θ̂ at which to evaluate
#'        If list, should have components mu_hat and Sigma_hat
#' @param completed_datasets List of M completed datasets
#' @return Numeric. Q̄_MI(θ̂) value
#' @details
#'   For each completed dataset, evaluates complete-data log-likelihood at θ̂,
#'   then averages across M imputations.
#' @examples
#'   Q_MI <- compute_Q_function(theta_hat, completed_datasets)
compute_Q_function <- function(theta_hat, completed_datasets) {

  # Source fit_models.R if not loaded
  if (!exists("compute_complete_loglik", mode = "function")) {
    source("simulations/utils/fit_models.R")
  }

  # Extract mu and Sigma from theta_hat
  if (is.list(theta_hat)) {
    mu_hat <- theta_hat$mu_hat
    Sigma_hat <- theta_hat$Sigma_hat
  } else {
    stop("theta_hat must be a list with components mu_hat and Sigma_hat")
  }

  M <- length(completed_datasets)

  # Compute log-likelihood on each completed dataset
  loglik_vals <- sapply(1:M, function(m) {
    compute_complete_loglik(completed_datasets[[m]], mu_hat, Sigma_hat)
  })

  # Average across imputations
  Q_MI <- mean(loglik_vals)

  return(Q_MI)
}


#' Compute Imputation Bias
#'
#' Calculates empirical imputation bias: Q̄_MI(θ̂) - ℓ_com(θ̂)
#'
#' @param Q_MI Numeric. Q-function value
#' @param ell_com Numeric. Complete-data log-likelihood value
#' @return Numeric. Imputation bias (positive indicates overestimation)
#' @details
#'   Bias = E[Q̄_MI(θ̂)] - E[ℓ_com(θ̂)]
#'   Positive bias means Q overestimates the complete-data log-likelihood
#' @examples
#'   bias <- compute_imputation_bias(Q_MI, ell_com)
compute_imputation_bias <- function(Q_MI, ell_com) {
  return(Q_MI - ell_com)
}


#' Compute Selection Entropy
#'
#' Calculates Shannon entropy of model selection frequencies
#'
#' @param selected_models Character vector. Selected model names across replications
#' @return Numeric. Entropy value (higher = more uncertainty in selection)
#' @details
#'   H = -Σ p_i log(p_i) where p_i = frequency of model i
#'   H = 0: always selects same model (no entropy)
#'   H = log(K): uniform over K models (maximum entropy)
#'   Target for calibration: H ≈ 0.9-1.1 (interesting selection problem)
#' @examples
#'   selected <- c("CS", "CS", "Toeplitz", "CS", "Unstructured", "Toeplitz")
#'   entropy <- compute_selection_entropy(selected)
compute_selection_entropy <- function(selected_models) {

  # Compute frequencies
  freq_table <- table(selected_models)
  probs <- as.numeric(freq_table) / sum(freq_table)

  # Remove zero probabilities (log(0) issue)
  probs <- probs[probs > 0]

  # Compute Shannon entropy: H = -Σ p_i log(p_i)
  entropy <- -sum(probs * log(probs))

  return(entropy)
}


#' Compute Model Selection Frequencies
#'
#' Summarizes model selection results across replications
#'
#' @param selected_models Character vector. Selected model names
#' @param true_model Character. True data-generating model (optional)
#' @return data.frame with selection frequencies and proportions
#' @details
#'   Returns frequency table plus:
#'   - P(correct): proportion selecting true model (if true_model provided)
#'   - Entropy: Shannon entropy of selection distribution
#' @examples
#'   summary <- compute_selection_frequencies(selected_models, true_model = "Toeplitz")
compute_selection_frequencies <- function(selected_models, true_model = NULL) {

  # Frequency table
  freq_table <- table(selected_models)
  n_reps <- length(selected_models)

  # Create data frame
  result_df <- data.frame(
    model = names(freq_table),
    count = as.numeric(freq_table),
    proportion = as.numeric(freq_table) / n_reps,
    stringsAsFactors = FALSE
  )

  # Add indicator for correct model
  if (!is.null(true_model)) {
    result_df$is_correct <- result_df$model == true_model
    p_correct <- result_df$proportion[result_df$is_correct]
    if (length(p_correct) == 0) p_correct <- 0
  } else {
    p_correct <- NA
  }

  # Compute entropy
  entropy <- compute_selection_entropy(selected_models)

  # Add summary attributes
  attr(result_df, "n_reps") <- n_reps
  attr(result_df, "p_correct") <- p_correct
  attr(result_df, "entropy") <- entropy
  attr(result_df, "true_model") <- true_model

  return(result_df)
}


#' Compute Summary Statistics for Bias Validation
#'
#' Aggregates results across replications for bias validation simulation
#'
#' @param bias_empirical Numeric vector. Empirical bias values across replications
#' @param tr_RIV_theoretical Numeric vector. Theoretical tr(RIV) values
#' @param Q Integer or vector. Number of parameters
#' @return data.frame with summary statistics
#' @details
#'   Computes means, SDs, and comparison between empirical and theoretical
#' @examples
#'   summary <- compute_bias_summary(bias_empirical, tr_RIV_theoretical, Q = 12)
compute_bias_summary <- function(bias_empirical, tr_RIV_theoretical, Q) {

  n_reps <- length(bias_empirical)

  if (length(tr_RIV_theoretical) != n_reps) {
    stop("bias_empirical and tr_RIV_theoretical must have same length")
  }

  # Summary statistics
  bias_mean <- mean(bias_empirical, na.rm = TRUE)
  bias_sd <- sd(bias_empirical, na.rm = TRUE)
  bias_median <- median(bias_empirical, na.rm = TRUE)

  tr_RIV_mean <- mean(tr_RIV_theoretical, na.rm = TRUE)
  tr_RIV_sd <- sd(tr_RIV_theoretical, na.rm = TRUE)

  # Comparison: empirical vs theoretical
  diff_mean <- mean(bias_empirical - tr_RIV_theoretical, na.rm = TRUE)
  diff_sd <- sd(bias_empirical - tr_RIV_theoretical, na.rm = TRUE)

  # Correlation
  cor_val <- cor(bias_empirical, tr_RIV_theoretical, use = "complete.obs")

  # RMSE
  rmse <- sqrt(mean((bias_empirical - tr_RIV_theoretical)^2, na.rm = TRUE))

  # Relative error
  rel_error_mean <- mean(abs(bias_empirical - tr_RIV_theoretical) / abs(tr_RIV_theoretical), na.rm = TRUE)

  result <- data.frame(
    n_reps = n_reps,
    Q = if (length(Q) == 1) Q else NA,
    bias_mean = bias_mean,
    bias_sd = bias_sd,
    bias_median = bias_median,
    tr_RIV_mean = tr_RIV_mean,
    tr_RIV_sd = tr_RIV_sd,
    diff_mean = diff_mean,
    diff_sd = diff_sd,
    correlation = cor_val,
    rmse = rmse,
    rel_error_mean = rel_error_mean
  )

  return(result)
}


#' Compute Information Criteria (AIC, BIC)
#'
#' Computes standard and MI-corrected information criteria
#'
#' @param logLik Numeric. Log-likelihood value
#' @param Q Integer. Number of parameters
#' @param n Integer. Sample size
#' @param tr_RIV Numeric. Trace of RIV matrix (for MI-corrected versions)
#' @return List with AIC, BIC, MI-AIC, MI-BIC values
#' @details
#'   AIC = -2 logLik + 2Q
#'   BIC = -2 logLik + Q log(n)
#'   MI-AIC = -2 Q̄_MI + 2Q + 4*tr(RIV)
#'   MI-BIC = -2 Q̄_MI + Q*log(n) + 2*tr(RIV)
#' @examples
#'   ic <- compute_information_criteria(logLik = -500, Q = 12, n = 200, tr_RIV = 3.5)
compute_information_criteria <- function(logLik, Q, n, tr_RIV = NULL) {

  # Standard criteria
  AIC_standard <- -2 * logLik + 2 * Q
  BIC_standard <- -2 * logLik + Q * log(n)

  result <- list(
    AIC = AIC_standard,
    BIC = BIC_standard
  )

  # MI-corrected criteria (if tr_RIV provided)
  if (!is.null(tr_RIV)) {
    MI_AIC <- -2 * logLik + 2 * Q + 4 * tr_RIV
    MI_BIC <- -2 * logLik + Q * log(n) + 2 * tr_RIV

    result$MI_AIC <- MI_AIC
    result$MI_BIC <- MI_BIC
    result$tr_RIV <- tr_RIV
  }

  return(result)
}
