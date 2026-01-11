# ============================================================================
# brms Imputation for Study 1
# ============================================================================
# Fits brms models with mi() syntax for Bayesian imputation
# Three imputation model types: matched, saturated, uncongenial
# ============================================================================

#' Fit brms Imputation Model
#'
#' Fits a brms multivariate model with missing data using mi() syntax
#'
#' @param data_miss data.frame with missing values (X fully observed, M/Y have NA)
#' @param model_type Character. "matched", "saturated", or "uncongenial"
#' @param backend Character. "rstan" or "cmdstanr"
#' @param chains Integer. Number of MCMC chains (default: 1 for parallel across reps)
#' @param warmup Integer. Number of warmup iterations (default: 1000)
#' @param iter Integer. Total iterations including warmup (default: 3000)
#' @param seed Integer. Random seed for Stan
#' @param refresh Integer. Progress output frequency (0 = silent)
#' @param ... Additional arguments passed to brm()
#' @return brmsfit object
#' @details
#'   Imputation model types (for DGM: X -> M -> Y full mediation):
#'
#'   - **matched**: Exactly mirrors DGM
#'     M | mi() ~ X
#'     Y | mi() ~ mi(M)
#'     Y imputed from M only (no X)
#'
#'   - **saturated**: Inclusive (adds auxiliary predictors)
#'     M | mi() ~ X
#'     Y | mi() ~ X + mi(M)
#'     Y imputed from both X and M
#'
#'   - **uncongenial**: Violates congeniality (omits M from Y model)
#'     M | mi() ~ X
#'     Y | mi() ~ X
#'     Y imputed from X only, ignoring M
#'
#' @examples
#'   fit <- fit_brms_imputation(data_miss, model_type = "saturated", seed = 123)
fit_brms_imputation <- function(data_miss,
                                 model_type = c("matched", "saturated", "uncongenial"),
                                 backend = "rstan",
                                 chains = 1,
                                 warmup = 1000,
                                 iter = 3000,
                                 seed = NULL,
                                 refresh = 0,
                                 ...) {

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' required. Install with: install.packages('brms')")
  }

  model_type <- match.arg(model_type)

  # Build brms formula based on model type
  bform <- build_imputation_formula(model_type)

  # Fit the model
  fit <- brms::brm(
    formula = bform,
    data = data_miss,
    backend = backend,
    chains = chains,
    warmup = warmup,
    iter = iter,
    seed = seed,
    refresh = refresh,
    ...
  )

  # Add metadata
  attr(fit, "imputation_model_type") <- model_type

  return(fit)
}


#' Build brms Imputation Formula
#'
#' Creates the brms formula for the specified imputation model type
#'
#' @param model_type Character. "matched", "saturated", or "uncongenial"
#' @return brmsformula object
build_imputation_formula <- function(model_type) {

  if (model_type == "matched") {
    # Matched: Y imputed from M only (mirrors full mediation DGM)
    # M | mi() ~ X
    # Y | mi() ~ mi(M)
    bform <- brms::bf(M | mi() ~ X) +
             brms::bf(Y | mi() ~ mi(M)) +
             brms::set_rescor(FALSE)

  } else if (model_type == "saturated") {
    # Saturated: Y imputed from X and M (inclusive)
    # M | mi() ~ X
    # Y | mi() ~ X + mi(M)
    bform <- brms::bf(M | mi() ~ X) +
             brms::bf(Y | mi() ~ X + mi(M)) +
             brms::set_rescor(FALSE)

  } else if (model_type == "uncongenial") {
    # Uncongenial: Y imputed from X only (omits M, violates congeniality)
    # M | mi() ~ X
    # Y | mi() ~ X
    bform <- brms::bf(M | mi() ~ X) +
             brms::bf(Y | mi() ~ X) +
             brms::set_rescor(FALSE)

  } else {
    stop(sprintf("Unknown model_type: %s", model_type))
  }

  return(bform)
}


#' Fit brms Imputation from Config
#'
#' Convenience wrapper using config parameters
#'
#' @param data_miss data.frame with missing values
#' @param model_type Character. "matched", "saturated", or "uncongenial"
#' @param config List. Configuration from config.R
#' @param seed Integer. Random seed for Stan
#' @return brmsfit object
fit_brms_imputation_from_config <- function(data_miss, model_type, config, seed = NULL) {

  fit_brms_imputation(
    data_miss = data_miss,
    model_type = model_type,
    backend = config$brms_backend,
    chains = config$brms_chains,
    warmup = config$brms_warmup,
    iter = config$brms_iter,
    seed = seed,
    refresh = config$brms_refresh
  )
}


#' Check brms Fit Convergence
#'
#' Basic convergence diagnostics for brms fit
#'
#' @param fit brmsfit object
#' @param rhat_threshold Numeric. Maximum acceptable Rhat (default: 1.05)
#' @param ess_threshold Numeric. Minimum acceptable ESS (default: 100)
#' @return List with convergence diagnostics
check_brms_convergence <- function(fit, rhat_threshold = 1.05, ess_threshold = 100) {

  # Get summary
  fit_summary <- summary(fit)

  # Extract Rhat and ESS from fixed effects
  if (!is.null(fit_summary$fixed)) {
    rhats <- fit_summary$fixed[, "Rhat"]
    ess_bulk <- fit_summary$fixed[, "Bulk_ESS"]
    ess_tail <- fit_summary$fixed[, "Tail_ESS"]
  } else {
    rhats <- NA
    ess_bulk <- NA
    ess_tail <- NA
  }

  # Check thresholds
  max_rhat <- max(rhats, na.rm = TRUE)
  min_ess_bulk <- min(ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(ess_tail, na.rm = TRUE)

  rhat_ok <- max_rhat < rhat_threshold
  ess_ok <- min_ess_bulk > ess_threshold && min_ess_tail > ess_threshold

  converged <- rhat_ok && ess_ok

  return(list(
    converged = converged,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    rhat_ok = rhat_ok,
    ess_ok = ess_ok,
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold
  ))
}


#' Print brms Imputation Summary
#'
#' Summarizes the brms imputation fit
#'
#' @param fit brmsfit object
print_brms_imputation_summary <- function(fit) {

  model_type <- attr(fit, "imputation_model_type")
  if (is.null(model_type)) model_type <- "unknown"

  convergence <- check_brms_convergence(fit)

  cat("=== brms Imputation Summary ===\n")
  cat(sprintf("Model type: %s\n", model_type))
  cat(sprintf("Chains: %d\n", fit$fit@sim$chains))
  cat(sprintf("Iterations: %d (warmup: %d)\n",
              fit$fit@sim$iter, fit$fit@sim$warmup))
  cat(sprintf("Post-warmup samples: %d\n",
              (fit$fit@sim$iter - fit$fit@sim$warmup) * fit$fit@sim$chains))

  cat("\nConvergence:\n")
  cat(sprintf("  Max Rhat: %.3f (%s)\n",
              convergence$max_rhat,
              ifelse(convergence$rhat_ok, "OK", "WARNING")))
  cat(sprintf("  Min Bulk ESS: %.0f (%s)\n",
              convergence$min_ess_bulk,
              ifelse(convergence$ess_ok, "OK", "WARNING")))
  cat(sprintf("  Overall: %s\n",
              ifelse(convergence$converged, "Converged", "Check diagnostics")))
}
