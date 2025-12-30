# ============================================================================
# Fit Models Using lavaan
# ============================================================================
# Functions for:
#   - Building lavaan model syntax for CS, Toeplitz, Unstructured
#   - Fitting models with FIML estimation
#   - Extracting parameters, covariance matrices, BIC
# ============================================================================

#' Build lavaan Model Syntax
#'
#' Creates lavaan model specification for different covariance structures
#'
#' @param p Integer. Number of variables
#' @param structure Character. One of "CS", "Toeplitz", "Unstructured"
#' @param var_names Character vector. Variable names (default: Y1, Y2, ..., Yp)
#' @return Character string with lavaan model syntax
#' @details
#'   CS: All variances equal, all covariances equal
#'   Toeplitz: Covariances depend only on distance |i-j|
#'   Unstructured: No constraints (default lavaan behavior)
#' @examples
#'   syntax_cs <- build_lavaan_model(5, "CS")
#'   syntax_toep <- build_lavaan_model(5, "Toeplitz")
build_lavaan_model <- function(p, structure, var_names = NULL) {

  if (is.null(var_names)) {
    var_names <- paste0("Y", 1:p)
  }

  if (length(var_names) != p) {
    stop("Length of var_names must equal p")
  }

  structure <- match.arg(structure, c("CS", "Toeplitz", "Unstructured"))

  # For unstructured, lavaan default (no constraints) is sufficient
  if (structure == "Unstructured") {
    # Just specify free means and covariances
    syntax <- paste0(
      "# Unstructured covariance\n",
      paste0(var_names, " ~ 1", collapse = "\n"),  # Free means
      "\n"
    )
    return(syntax)
  }

  # Compound Symmetry: equal variances, equal covariances
  if (structure == "CS") {
    syntax_parts <- c(
      "# Compound Symmetry: Equal variances + equal covariances",
      "",
      "# Free means",
      paste0(var_names, " ~ 1", collapse = "\n"),
      "",
      "# Variances (constrain equal)",
      paste0(var_names, " ~~ v*", var_names, collapse = "\n"),
      "",
      "# Covariances (constrain equal)"
    )

    # Add covariance constraints
    cov_lines <- c()
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        cov_lines <- c(cov_lines, paste0(var_names[i], " ~~ c*", var_names[j]))
      }
    }
    syntax_parts <- c(syntax_parts, cov_lines)

    syntax <- paste(syntax_parts, collapse = "\n")
    return(syntax)
  }

  # Toeplitz: covariances depend on distance
  if (structure == "Toeplitz") {
    syntax_parts <- c(
      "# Toeplitz: Covariances depend on lag |i-j|",
      "",
      "# Free means",
      paste0(var_names, " ~ 1", collapse = "\n"),
      "",
      "# Variances (constrain equal)",
      paste0(var_names, " ~~ v*", var_names, collapse = "\n"),
      "",
      "# Covariances by lag"
    )

    # Group covariances by lag (distance)
    for (lag in 1:(p-1)) {
      lag_label <- paste0("lag", lag)
      syntax_parts <- c(
        syntax_parts,
        paste0("# Lag ", lag, " covariances")
      )

      lag_lines <- c()
      for (i in 1:(p-lag)) {
        j <- i + lag
        lag_lines <- c(lag_lines, paste0(var_names[i], " ~~ ", lag_label, "*", var_names[j]))
      }
      syntax_parts <- c(syntax_parts, lag_lines, "")
    }

    syntax <- paste(syntax_parts, collapse = "\n")
    return(syntax)
  }

  stop("Should not reach here")
}


#' Fit lavaan Model to Data
#'
#' Fits a structured covariance model using lavaan with FIML
#'
#' @param data Numeric matrix or data.frame. Can contain NA values
#' @param structure Character. "CS", "Toeplitz", or "Unstructured"
#' @param return_vcov Logical. Whether to return parameter variance-covariance matrix
#' @param ... Additional arguments passed to lavaan::sem()
#' @return List with components:
#'   \item{fit}{lavaan fitted object}
#'   \item{converged}{Logical, whether optimization converged}
#'   \item{mu_hat}{Estimated mean vector}
#'   \item{Sigma_hat}{Estimated covariance matrix}
#'   \item{theta_hat}{Parameter vector (means + covariance params)}
#'   \item{vcov_theta}{Variance-covariance matrix of parameter estimates (if return_vcov = TRUE)}
#'   \item{logLik}{Log-likelihood value}
#'   \item{BIC}{BIC value}
#'   \item{nparams}{Number of free parameters}
#' @examples
#'   data_miss <- matrix(rnorm(100*5), 100, 5)
#'   data_miss[sample(500, 50)] <- NA
#'   fit_result <- fit_lavaan_model(data_miss, structure = "CS")
fit_lavaan_model <- function(data, structure, return_vcov = TRUE, ...) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' required for model fitting")
  }

  # Convert to data.frame if matrix
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  p <- ncol(data)
  var_names <- colnames(data)
  if (is.null(var_names)) {
    var_names <- paste0("Y", 1:p)
    colnames(data) <- var_names
  }

  # Build model syntax
  model_syntax <- build_lavaan_model(p, structure, var_names)

  # Fit model with FIML
  fit <- tryCatch(
    {
      lavaan::sem(
        model = model_syntax,
        data = data,
        missing = "fiml",
        meanstructure = TRUE,
        fixed.x = FALSE,
        ...
      )
    },
    error = function(e) {
      warning(paste("lavaan fitting failed:", e$message))
      return(NULL)
    }
  )

  # Handle fitting failures
  if (is.null(fit)) {
    return(list(
      converged = FALSE,
      fit = NULL,
      mu_hat = NULL,
      Sigma_hat = NULL,
      theta_hat = NULL,
      vcov_theta = NULL,
      logLik = NA,
      BIC = NA,
      nparams = NA
    ))
  }

  # Extract convergence status
  converged <- lavaan::lavInspect(fit, "converged")

  # Extract fitted means and covariance
  # lavInspect returns a list with one element per group (we have 1 group)
  mean_list <- lavaan::lavInspect(fit, "mean.ov")
  cov_list <- lavaan::lavInspect(fit, "cov.ov")

  # Handle both list and direct vector returns
  if (is.list(mean_list)) {
    mu_hat <- as.numeric(mean_list[[1]])
  } else {
    mu_hat <- as.numeric(mean_list)
  }

  if (is.list(cov_list)) {
    Sigma_hat <- as.matrix(cov_list[[1]])
  } else {
    Sigma_hat <- as.matrix(cov_list)
  }

  # Ensure names
  if (is.null(names(mu_hat))) {
    names(mu_hat) <- var_names
  }
  if (is.null(colnames(Sigma_hat))) {
    colnames(Sigma_hat) <- rownames(Sigma_hat) <- var_names
  }

  # Extract parameter estimates
  # Use lavaan::coef() which gives the free parameters that match vcov()
  # This ensures theta_hat and vcov_theta have compatible dimensions
  theta_hat <- lavaan::coef(fit)

  # Note: lavaan::coef() returns all free parameters in the model
  # For constrained models (CS, Toeplitz), this includes labeled parameters
  # The names will be lavaan's internal names, which is fine for computation

  # Get variance-covariance matrix of parameter estimates
  vcov_theta <- NULL
  if (return_vcov && converged) {
    vcov_theta <- tryCatch(
      lavaan::vcov(fit),
      error = function(e) {
        warning("Could not extract vcov matrix")
        return(NULL)
      }
    )
  }

  # Extract fit statistics
  logLik_val <- lavaan::logLik(fit)
  BIC_val <- stats::BIC(fit)  # Use stats::BIC, not lavaan::BIC
  nparams <- lavaan::fitmeasures(fit, "npar")

  return(list(
    fit = fit,
    converged = converged,
    mu_hat = mu_hat,
    Sigma_hat = Sigma_hat,
    theta_hat = theta_hat,
    vcov_theta = vcov_theta,
    logLik = as.numeric(logLik_val),
    BIC = BIC_val,
    nparams = nparams
  ))
}


#' Extract Parameter Count Q
#'
#' Returns the number of free parameters Q for a fitted model or structure
#'
#' @param p Integer. Dimension (if fit not provided)
#' @param structure Character. "CS", "Toeplitz", or "Unstructured" (if fit not provided)
#' @param fit_result List. Result from fit_lavaan_model() (optional)
#' @return Integer. Number of free parameters
#' @details
#'   CS: Q = p + 2 (means + variance + covariance)
#'   Toeplitz: Q = p + p = 2p (means + variance + p-1 lag covariances)
#'   Unstructured: Q = p + p(p+1)/2 (means + unique covariance elements)
#' @examples
#'   get_param_count(p = 10, structure = "CS")  # Returns 12
#'   get_param_count(fit_result = fit_cs)       # Extracts from fitted model
get_param_count <- function(p = NULL, structure = NULL, fit_result = NULL) {

  # If fit_result provided, extract from there
  if (!is.null(fit_result)) {
    if (!is.null(fit_result$nparams)) {
      return(fit_result$nparams)
    }
  }

  # Otherwise compute from p and structure
  if (is.null(p) || is.null(structure)) {
    stop("Must provide either fit_result or both p and structure")
  }

  structure <- match.arg(structure, c("CS", "Toeplitz", "Unstructured"))

  Q <- switch(structure,
    "CS" = p + 2,
    "Toeplitz" = 2 * p,
    "Unstructured" = p + p * (p + 1) / 2
  )

  return(as.integer(Q))
}


#' Compute Complete-Data Log-Likelihood
#'
#' Evaluates the complete-data log-likelihood at given parameters
#'
#' @param data_complete Numeric matrix. Complete data (no NAs)
#' @param mu Numeric vector. Mean parameter
#' @param Sigma Numeric matrix. Covariance parameter
#' @return Numeric. Log-likelihood value
#' @details
#'   For MVN: log L = -n/2 * (p*log(2*pi) + log|Sigma|) - 1/2 * sum((x_i - mu)' Sigma^{-1} (x_i - mu))
#' @examples
#'   data_complete <- matrix(rnorm(100*5), 100, 5)
#'   mu <- rep(0, 5)
#'   Sigma <- diag(5)
#'   ll <- compute_complete_loglik(data_complete, mu, Sigma)
compute_complete_loglik <- function(data_complete, mu, Sigma) {

  if (!is.matrix(data_complete)) {
    data_complete <- as.matrix(data_complete)
  }

  if (any(is.na(data_complete))) {
    stop("data_complete must not contain NA values")
  }

  n <- nrow(data_complete)
  p <- ncol(data_complete)

  if (length(mu) != p || nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("Dimension mismatch between data, mu, and Sigma")
  }

  # Use mvtnorm for stable computation
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' required for log-likelihood computation")
  }

  # Compute log-likelihood for each observation
  log_lik_vals <- mvtnorm::dmvnorm(data_complete, mean = mu, sigma = Sigma, log = TRUE)

  # Sum across observations
  logLik_total <- sum(log_lik_vals)

  return(logLik_total)
}
