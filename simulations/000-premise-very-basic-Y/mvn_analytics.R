# ============================================================================
# Analytic MVN Functions for Term Decomposition Verification
# ============================================================================
# All formulas are analytic (no numerical differentiation)
# Reference: simulation_plan.md for equation numbers
# ============================================================================

# ----------------------------------------------------------------------------
# Helper Functions
# ----------------------------------------------------------------------------

#' Half-vectorization of a symmetric matrix
#' @param M A symmetric matrix
#' @return Vector of lower triangular elements (including diagonal)
vech <- function(M) {
  M[lower.tri(M, diag = TRUE)]
}

#' Reconstruct symmetric matrix from half-vectorization
#' @param v Vector from vech()
#' @param p Dimension of the matrix
#' @return p x p symmetric matrix
vech_inv <- function(v, p) {
  M <- matrix(0, p, p)
  M[lower.tri(M, diag = TRUE)] <- v
  M <- M + t(M) - diag(diag(M))
  return(M)
}

#' Create the duplication matrix D_p
#' D_p maps vech(A) to vec(A) for symmetric A
#' @param p Dimension
#' @return p^2 x p(p+1)/2 duplication matrix
duplication_matrix <- function(p) {
  # Number of unique elements in symmetric p x p matrix
  n_vech <- p * (p + 1) / 2
  n_vec <- p^2

  D <- matrix(0, nrow = n_vec, ncol = n_vech)

  # Build index mapping
  vech_idx <- 1
  for (j in 1:p) {
    for (i in j:p) {
      # Position in vec (column-major)
      vec_idx_ij <- (j - 1) * p + i  # (i,j) element
      vec_idx_ji <- (i - 1) * p + j  # (j,i) element

      D[vec_idx_ij, vech_idx] <- 1
      if (i != j) {
        D[vec_idx_ji, vech_idx] <- 1
      }
      vech_idx <- vech_idx + 1
    }
  }

  return(D)
}

# ----------------------------------------------------------------------------
# Core MVN Functions
# ----------------------------------------------------------------------------

#' Complete-data log-likelihood for MVN (Eq. 7)
#' @param Y N x p data matrix
#' @param mu p-vector mean
#' @param Sigma p x p covariance matrix
#' @return Scalar log-likelihood
mvn_loglik <- function(Y, mu, Sigma) {
  Y <- as.matrix(Y)
  N <- nrow(Y)
  p <- ncol(Y)

  # Sample covariance centered at mu (Eq. 8)
  Y_centered <- sweep(Y, 2, mu)
  S_mu <- crossprod(Y_centered) / N  # (1/N) * sum of outer products

  # Log-likelihood (Eq. 7)
  log_det <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  Sigma_inv <- solve(Sigma)

  ll <- -N * p / 2 * log(2 * pi) - N / 2 * log_det - N / 2 * sum(Sigma_inv * S_mu)

  return(as.numeric(ll))
}

#' Score function for MVN (Eqs. 15, 16, 17)
#' @param Y N x p data matrix
#' @param mu p-vector mean
#' @param Sigma p x p covariance matrix
#' @return Q-vector score (Q = p + p(p+1)/2)
mvn_score <- function(Y, mu, Sigma) {
  Y <- as.matrix(Y)
  N <- nrow(Y)
  p <- ncol(Y)

  # Sample mean and covariance centered at mu
  Y_bar <- colMeans(Y)
  Y_centered <- sweep(Y, 2, mu)
  S_mu <- crossprod(Y_centered) / N

  Sigma_inv <- solve(Sigma)

  # Score w.r.t. mu (Eq. 15): N * Sigma^{-1} * (Y_bar - mu)
  score_mu <- as.vector(N * Sigma_inv %*% (Y_bar - mu))

  # Score w.r.t. Sigma (matrix form, Eq. 16)
  # dℓ/dΣ = (N/2) * Σ^{-1} * (S_μ - Σ) * Σ^{-1}
  score_Sigma_mat <- (N / 2) * Sigma_inv %*% (S_mu - Sigma) %*% Sigma_inv

  # Convert to vech form (Eq. 17): D_p' * vec(dℓ/dΣ)
  D_p <- duplication_matrix(p)
  score_vech_Sigma <- as.vector(t(D_p) %*% as.vector(score_Sigma_mat))

  # Concatenate
  return(c(score_mu, score_vech_Sigma))
}

#' Fisher information matrix for MVN (Eqs. 18, 19, 21)
#' @param N Sample size
#' @param Sigma p x p covariance matrix
#' @return Q x Q Fisher information matrix (block diagonal)
mvn_fisher_info <- function(N, Sigma) {
  p <- nrow(Sigma)
  Sigma_inv <- solve(Sigma)

  # I_μμ = N * Σ^{-1} (Eq. 18)
  I_mu <- N * Sigma_inv

  # I_vech(Σ) = (N/2) * D_p' * (Σ^{-1} ⊗ Σ^{-1}) * D_p (Eq. 19)
  D_p <- duplication_matrix(p)
  Sigma_inv_kron <- kronecker(Sigma_inv, Sigma_inv)
  I_vech_Sigma <- (N / 2) * t(D_p) %*% Sigma_inv_kron %*% D_p

  # Block diagonal (Eq. 20: cross terms are zero)
  n_mu <- p
  n_vech <- p * (p + 1) / 2
  Q <- n_mu + n_vech

  I_full <- matrix(0, Q, Q)
  I_full[1:n_mu, 1:n_mu] <- I_mu
  I_full[(n_mu + 1):Q, (n_mu + 1):Q] <- I_vech_Sigma

  return(I_full)
}

#' Complete-data MLE for MVN (Eq. 22)
#' @param Y N x p data matrix
#' @return List with mu, Sigma, and theta_vec
complete_data_mle <- function(Y) {
  Y <- as.matrix(Y)
  N <- nrow(Y)

  # MLE estimates
  mu_hat <- colMeans(Y)
  # MLE uses N divisor, not N-1
  Sigma_hat <- cov(Y) * (N - 1) / N

  # Parameter vector
  theta_vec <- c(mu_hat, vech(Sigma_hat))

  return(list(
    mu = mu_hat,
    Sigma = Sigma_hat,
    theta_vec = theta_vec
  ))
}

# ----------------------------------------------------------------------------
# Rubin's Variance Components
# ----------------------------------------------------------------------------

#' Compute within-imputation variance W_bar (Eq. 23)
#' @param completed_datasets List of M completed data matrices
#' @return Q x Q average within-imputation variance
compute_W_bar <- function(completed_datasets) {
  M <- length(completed_datasets)
  N <- nrow(completed_datasets[[1]])

  # Compute inverse Fisher info for each imputed dataset
  W_list <- lapply(completed_datasets, function(Y) {
    mle <- complete_data_mle(Y)
    I_com <- mvn_fisher_info(N, mle$Sigma)
    solve(I_com)  # Variance = I^{-1}
  })

  # Average
  W_bar <- Reduce(`+`, W_list) / M
  return(W_bar)
}

#' Compute between-imputation variance B (Eq. 24)
#' @param completed_datasets List of M completed data matrices
#' @return Q x Q between-imputation variance
compute_B <- function(completed_datasets) {
  M <- length(completed_datasets)

  # Get theta_hat for each imputed dataset
  theta_list <- lapply(completed_datasets, function(Y) {
    complete_data_mle(Y)$theta_vec
  })

  # Stack into matrix
  theta_mat <- do.call(rbind, theta_list)  # M x Q

  # Mean
  theta_bar <- colMeans(theta_mat)

  # Between-imputation variance
  theta_centered <- sweep(theta_mat, 2, theta_bar)
  B <- crossprod(theta_centered) / (M - 1)

  return(B)
}

#' Compute tr(RIV) (Eq. 27)
#' @param W_bar Within-imputation variance
#' @param B Between-imputation variance
#' @param M Number of imputations
#' @return Scalar tr(RIV)
compute_tr_RIV <- function(W_bar, B, M) {
  # RIV = (1 + 1/M) * W_bar^{-1} * B (Eq. 26)
  W_bar_inv <- solve(W_bar)
  RIV <- (1 + 1/M) * W_bar_inv %*% B

  return(sum(diag(RIV)))
}

# ----------------------------------------------------------------------------
# Covariance Structure Generators
# ----------------------------------------------------------------------------

#' Generate Toeplitz covariance matrix
#' @param p Dimension
#' @param rho Correlation decay parameter
#' @param sigma2 Variance (default 1)
#' @return p x p Toeplitz covariance matrix
generate_toeplitz <- function(p, rho, sigma2 = 1) {
  indices <- 0:(p - 1)
  first_row <- sigma2 * rho^indices
  return(toeplitz(first_row))
}

#' Generate compound symmetry (exchangeable) covariance matrix
#' @param p Dimension
#' @param rho Common correlation
#' @param sigma2 Variance (default 1)
#' @return p x p CS covariance matrix
generate_cs <- function(p, rho, sigma2 = 1) {
  Sigma <- matrix(rho * sigma2, p, p)
  diag(Sigma) <- sigma2
  return(Sigma)
}

# ----------------------------------------------------------------------------
# Lavaan Helper Functions
# ----------------------------------------------------------------------------

#' Generate saturated MVN model syntax for lavaan
#' @param var_names Character vector of variable names
#' @return Lavaan model syntax string
generate_saturated_model <- function(var_names) {
  p <- length(var_names)
  lines <- character()

  # Covariances (including variances on diagonal)
  for (i in 1:p) {
    covs <- paste(var_names[i:p], collapse = " + ")
    lines <- c(lines, paste0(var_names[i], " ~~ ", covs))
  }

  paste(lines, collapse = "\n")
}

#' Get permutation to map lavaan parameter order to our (mu, vech(Sigma)) order
#' @param fit A fitted lavaan object
#' @return Integer vector for reordering
get_lavaan_to_theta_perm <- function(fit) {
  pnames <- names(coef(fit))
  var_order <- lavaan::lavNames(fit, "ov")
  p <- length(var_order)


  # Map mean indices (intercepts: "Y1~1", "Y2~1", ...)
  mean_idx <- sapply(var_order, function(v) {
    which(pnames == paste0(v, "~1"))
  })


  # Map vech(Sigma) indices - lower triangular, column-major: (1,1), (2,1), (2,2), (3,1), ...
  cov_idx <- integer(p * (p + 1) / 2)
  k <- 1
  for (j in 1:p) {
    for (i in j:p) {
      # Try both orderings since lavaan may use either
      cov_name1 <- paste0(var_order[i], "~~", var_order[j])
      cov_name2 <- paste0(var_order[j], "~~", var_order[i])
      idx <- which(pnames == cov_name1)
      if (length(idx) == 0) idx <- which(pnames == cov_name2)
      cov_idx[k] <- idx
      k <- k + 1
    }
  }

  c(mean_idx, cov_idx)
}

#' Extract information matrix from lavaan fit in our parameter ordering
#' @param fit A fitted lavaan object
#' @return Q x Q information matrix (negative Hessian)
get_info_matrix <- function(fit) {
  perm <- get_lavaan_to_theta_perm(fit)
  H <- -lavaan::lavInspect(fit, "hessian")
  H[perm, perm]
}

#' Extract MLE parameters from lavaan fit in our ordering
#' @param fit A fitted lavaan object
#' @return List with mu, Sigma, theta_vec
get_mle_from_lavaan <- function(fit) {
  var_order <- lavaan::lavNames(fit, "ov")
  p <- length(var_order)

  # Extract estimates
  est <- lavaan::lavInspect(fit, "est")

  # For saturated models without latent variables:
  # - theta contains the covariance matrix
  # - nu contains the means (intercepts)

  # Mean vector (from nu - intercepts)
  if (!is.null(est$nu) && length(est$nu) > 0) {
    mu <- as.vector(est$nu[var_order, ])
  } else if (!is.null(est$mean)) {
    mu <- est$mean[var_order]
  } else {
    stop("Cannot find means in lavaan output")
  }
  names(mu) <- NULL

  # Covariance matrix (from theta for observed variables)
  if (!is.null(est$theta) && nrow(est$theta) > 0) {
    Sigma <- as.matrix(est$theta[var_order, var_order])
  } else if (!is.null(est$cov)) {
    Sigma <- est$cov[var_order, var_order]
  } else {
    stop("Cannot find covariance in lavaan output")
  }

  # Parameter vector in our ordering
  theta_vec <- c(mu, vech(Sigma))

  list(mu = mu, Sigma = Sigma, theta_vec = theta_vec)
}

# ----------------------------------------------------------------------------
# Analytic Q-Function
# ----------------------------------------------------------------------------

#' Compute analytic Q-function for MVN with missing data
#'
#' Q(theta|theta) = ell_obs(theta) + Q_mis
#' where Q_mis = sum over patterns of n_k * [-q_k/2*log(2pi) - 0.5*log|Sigma_m|o,k| - q_k/2]
#'
#' @param Y_miss N x p data matrix with NAs for missing values
#' @param mu p-vector mean parameter
#' @param Sigma p x p covariance matrix parameter
#' @return Scalar Q-function value
compute_Q_analytic <- function(Y_miss, mu, Sigma) {
  Y_miss <- as.matrix(Y_miss)
  N <- nrow(Y_miss)
  p <- ncol(Y_miss)

  # Identify missing patterns
  R <- !is.na(Y_miss)  # Response indicator (TRUE = observed)
  patterns <- unique(R, MARGIN = 1)

  Q_total <- 0

  for (k in seq_len(nrow(patterns))) {
    pat <- patterns[k, ]
    obs_idx <- which(pat)
    mis_idx <- which(!pat)
    q_k <- length(mis_idx)

    # Count observations with this pattern
    pattern_match <- apply(R, 1, function(r) all(r == pat))
    n_k <- sum(pattern_match)

    if (q_k == 0) {
      # Complete observation - contributes full log-likelihood
      Y_k <- Y_miss[pattern_match, , drop = FALSE]
      Q_total <- Q_total + mvn_loglik(Y_k, mu, Sigma)
    } else if (length(obs_idx) == 0) {
      # Fully missing observation - Q_mis contribution only
      # Q_mis = n_k * [-q_k/2*log(2pi) - 0.5*log|Sigma| - q_k/2]
      log_det <- determinant(Sigma, logarithm = TRUE)$modulus[1]
      Q_k_mis <- -q_k/2 * log(2*pi) - 0.5 * log_det - q_k/2
      Q_total <- Q_total + n_k * Q_k_mis
    } else {
      # Partially observed - need conditional distribution
      Y_k <- Y_miss[pattern_match, , drop = FALSE]
      Y_obs_k <- Y_k[, obs_idx, drop = FALSE]

      # Partition parameters
      mu_o <- mu[obs_idx]
      mu_m <- mu[mis_idx]
      Sigma_oo <- Sigma[obs_idx, obs_idx, drop = FALSE]
      Sigma_mm <- Sigma[mis_idx, mis_idx, drop = FALSE]
      Sigma_mo <- Sigma[mis_idx, obs_idx, drop = FALSE]
      Sigma_om <- Sigma[obs_idx, mis_idx, drop = FALSE]

      # Observed data log-likelihood contribution
      Q_obs_k <- mvn_loglik(Y_obs_k, mu_o, Sigma_oo)

      # Conditional covariance Sigma_m|o
      Sigma_oo_inv <- solve(Sigma_oo)
      Sigma_m_given_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om

      # Q_mis contribution: n_k * [-q_k/2*log(2pi) - 0.5*log|Sigma_m|o| - q_k/2]
      log_det_cond <- determinant(Sigma_m_given_o, logarithm = TRUE)$modulus[1]
      Q_k_mis <- -q_k/2 * log(2*pi) - 0.5 * log_det_cond - q_k/2

      Q_total <- Q_total + Q_obs_k + n_k * Q_k_mis
    }
  }

  return(as.numeric(Q_total))
}

#' Compute score of Q-function for MVN with missing data
#'
#' The score of Q is E[score of ℓ_com] where expectation is over Y^m|Y^o
#'
#' @param Y_miss N x p data matrix with NAs for missing values
#' @param mu p-vector mean parameter
#' @param Sigma p x p covariance matrix parameter
#' @return Q-vector score (same length as theta = (mu, vech(Sigma)))
compute_Q_score <- function(Y_miss, mu, Sigma) {
  Y_miss <- as.matrix(Y_miss)
  N <- nrow(Y_miss)
  p <- ncol(Y_miss)

  Sigma_inv <- solve(Sigma)
  D_p <- duplication_matrix(p)

  # Initialize accumulators for sufficient statistics
  # We need: sum of E[Y], sum of E[YY']
  sum_E_Y <- rep(0, p)
  sum_E_YYt <- matrix(0, p, p)

  for (i in 1:N) {
    obs_idx <- which(!is.na(Y_miss[i, ]))
    mis_idx <- which(is.na(Y_miss[i, ]))

    if (length(mis_idx) == 0) {
      # Complete observation
      Y_i <- Y_miss[i, ]
      sum_E_Y <- sum_E_Y + Y_i
      sum_E_YYt <- sum_E_YYt + outer(Y_i, Y_i)
    } else if (length(obs_idx) == 0) {
      # Fully missing - E[Y] = mu, E[YY'] = Sigma + mu mu'
      sum_E_Y <- sum_E_Y + mu
      sum_E_YYt <- sum_E_YYt + Sigma + outer(mu, mu)
    } else {
      # Partially observed
      Y_o <- Y_miss[i, obs_idx]

      # Partition
      mu_o <- mu[obs_idx]
      mu_m <- mu[mis_idx]
      Sigma_oo <- Sigma[obs_idx, obs_idx, drop = FALSE]
      Sigma_mm <- Sigma[mis_idx, mis_idx, drop = FALSE]
      Sigma_mo <- Sigma[mis_idx, obs_idx, drop = FALSE]
      Sigma_om <- Sigma[obs_idx, mis_idx, drop = FALSE]

      # Conditional moments
      Sigma_oo_inv <- solve(Sigma_oo)
      mu_m_given_o <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (Y_o - mu_o)
      Sigma_m_given_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om

      # E[Y_i]
      E_Y_i <- rep(0, p)
      E_Y_i[obs_idx] <- Y_o
      E_Y_i[mis_idx] <- as.vector(mu_m_given_o)
      sum_E_Y <- sum_E_Y + E_Y_i

      # E[Y_i Y_i']
      # For observed parts: E[Y_o Y_o'] = Y_o Y_o'
      # For missing parts: E[Y_m Y_m'] = Var(Y_m|Y_o) + E[Y_m|Y_o] E[Y_m|Y_o]'
      # For cross: E[Y_o Y_m'] = Y_o E[Y_m|Y_o]'
      E_YYt_i <- matrix(0, p, p)
      E_YYt_i[obs_idx, obs_idx] <- outer(Y_o, Y_o)
      E_YYt_i[obs_idx, mis_idx] <- outer(Y_o, as.vector(mu_m_given_o))
      E_YYt_i[mis_idx, obs_idx] <- t(E_YYt_i[obs_idx, mis_idx])  # Symmetry
      E_YYt_i[mis_idx, mis_idx] <- Sigma_m_given_o + outer(as.vector(mu_m_given_o), as.vector(mu_m_given_o))
      sum_E_YYt <- sum_E_YYt + E_YYt_i
    }
  }

  # E[Y_bar] and E[S_mu]
  E_Y_bar <- sum_E_Y / N
  # E[S_mu] = E[(1/N) sum (Y_i - mu)(Y_i - mu)']
  #         = (1/N) sum E[(Y_i - mu)(Y_i - mu)']
  #         = (1/N) sum (E[Y_i Y_i'] - mu E[Y_i]' - E[Y_i] mu' + mu mu')
  #         = (1/N) (sum_E_YYt) - mu E_Y_bar' - E_Y_bar mu' + mu mu'
  E_S_mu <- sum_E_YYt / N - outer(mu, E_Y_bar) - outer(E_Y_bar, mu) + outer(mu, mu)

  # Score w.r.t. mu: N * Sigma^{-1} * (E[Y_bar] - mu)
  score_mu <- as.vector(N * Sigma_inv %*% (E_Y_bar - mu))

  # Score w.r.t. Sigma: (N/2) * Sigma^{-1} * (E[S_mu] - Sigma) * Sigma^{-1}
  score_Sigma_mat <- (N / 2) * Sigma_inv %*% (E_S_mu - Sigma) %*% Sigma_inv

  # Convert to vech form - just extract lower triangle
  # This matches the numerical gradient behavior
  score_vech_Sigma <- score_Sigma_mat[lower.tri(score_Sigma_mat, diag = TRUE)]

  return(c(score_mu, score_vech_Sigma))
}

#' Compute tr(RIV) analytically using Missing Information Principle
#'
#' As M -> infinity: tr(RIV) = tr(I_com * I_obs^{-1}) - Q
#' where Q is the number of parameters
#'
#' @param I_com Complete-data Fisher information matrix
#' @param I_obs Observed-data Fisher information matrix
#' @return Scalar tr(RIV)
compute_tr_RIV_analytic <- function(I_com, I_obs) {
  Q <- nrow(I_com)
  I_obs_inv <- solve(I_obs)
  tr_ratio <- sum(diag(I_com %*% I_obs_inv))
  return(tr_ratio - Q)
}
