# ============================================================================
# lavaan and lavaan.mi Analysis for Study 1
# ============================================================================
# Functions for:
# - Complete-data analysis (for validation)
# - MI analysis using lavaan.mi::sem.mi()
# - W, B, tr(RIV) computation
# ============================================================================

# ----------------------------------------------------------------------------
# Model Syntax Definitions
# ----------------------------------------------------------------------------

#' Get Mediation Model Syntax
#'
#' Returns lavaan model syntax for mediation analysis
#'
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @param with_starting_values Logical. Include starting values?
#' @param pop_params List. Population parameters (for starting values)
#' @return Character string with lavaan model syntax
get_mediation_syntax <- function(model_type = c("full_mediation", "partial_mediation"),
                                  with_starting_values = FALSE,
                                  pop_params = NULL) {

  model_type <- match.arg(model_type)

  if (model_type == "full_mediation") {
    # Full mediation: No direct X -> Y path
    syntax <- '
      # Paths
      M ~ a*X
      Y ~ b*M

      # Indirect effect
      indirect := a*b
    '
  } else {
    # Partial mediation: Include direct X -> Y path
    syntax <- '
      # Paths
      M ~ a*X
      Y ~ b*M + c*X

      # Effects
      indirect := a*b
      direct := c
      total := a*b + c
    '
  }

  # Add starting values if requested
  if (with_starting_values && !is.null(pop_params)) {
    # Note: This is a simplified approach; full implementation would
    # modify the syntax to include start() modifiers
    # For now, starting values can be passed via lavaan's start argument
  }

  return(syntax)
}


# ----------------------------------------------------------------------------
# Complete Data Analysis (for validation)
# ----------------------------------------------------------------------------

#' Fit Mediation Model to Complete Data
#'
#' Fits lavaan mediation model to complete data (for validation)
#'
#' @param data data.frame with X, M, Y columns (complete, no NAs)
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @return lavaan fit object
fit_mediation_complete <- function(data, model_type = "full_mediation") {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' required. Install with: install.packages('lavaan')")
  }

  syntax <- get_mediation_syntax(model_type)

  fit <- lavaan::sem(
    model = syntax,
    data = data,
    meanstructure = TRUE
  )

  return(fit)
}


#' Extract Parameters from Lavaan Fit
#'
#' Extracts relevant parameters from complete-data lavaan fit
#'
#' @param fit lavaan fit object
#' @return data.frame with parameter estimates, SEs, CIs
extract_lavaan_params <- function(fit) {

  params <- lavaan::parameterEstimates(fit, ci = TRUE)

  # Filter to parameters of interest
  params_of_interest <- params[params$label %in% c("a", "b", "c", "indirect", "direct", "total") |
                               (params$op == "~" & params$lhs %in% c("M", "Y")), ]

  return(params_of_interest)
}


#' Check Parameter Coverage
#'
#' Checks if true parameter value falls within CI
#'
#' @param params data.frame from extract_lavaan_params()
#' @param pop_params List with true parameter values
#' @return data.frame with coverage indicators
check_parameter_coverage <- function(params, pop_params) {

  # Map lavaan labels to population parameters
  param_mapping <- list(
    a = pop_params$beta_xm,
    b = pop_params$beta_my,
    c = pop_params$beta_xy,
    indirect = pop_params$indirect_effect,
    direct = pop_params$direct_effect,
    total = pop_params$total_effect
  )

  # Add coverage indicator
  params$true_value <- sapply(params$label, function(lbl) {
    if (lbl %in% names(param_mapping)) param_mapping[[lbl]] else NA
  })

  params$covered <- with(params, !is.na(true_value) &
                                  true_value >= ci.lower &
                                  true_value <= ci.upper)

  params$bias <- params$est - params$true_value

  return(params)
}


#' Run Complete-Data Validation
#'
#' Runs complete-data analysis across replications to validate data generation
#'
#' @param n_reps Integer. Number of replications to run
#' @param config List. Configuration from config.R
#' @param seeds List. Seeds from setup_seeds()
#' @param verbose Logical. Print progress?
#' @return data.frame with parameter estimates and coverage across replications
run_complete_data_validation <- function(n_reps, config, seeds, verbose = TRUE) {

  if (verbose) cat("Running complete-data validation...\n")

  # Get population parameters
  pop_params <- get_population_params(config)

  results <- vector("list", n_reps)

  for (rep in seq_len(n_reps)) {
    if (verbose && rep %% 50 == 0) cat(sprintf("  Rep %d/%d\n", rep, n_reps))

    # Generate complete data
    seed_data <- seeds[[sprintf("rep_%03d", rep)]]$data
    data_complete <- generate_mediation_data_from_config(config, seed = seed_data)

    # Fit model
    fit <- tryCatch(
      fit_mediation_complete(data_complete, model_type = "full_mediation"),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      results[[rep]] <- NULL
      next
    }

    # Extract and check parameters
    params <- extract_lavaan_params(fit)
    params <- check_parameter_coverage(params, pop_params)
    params$rep <- rep

    results[[rep]] <- params
  }

  # Combine results
  results_df <- do.call(rbind, results[!sapply(results, is.null)])

  return(results_df)
}


#' Summarize Complete-Data Validation
#'
#' Computes summary statistics from validation results
#'
#' @param validation_results data.frame from run_complete_data_validation()
#' @param pop_params List with population parameters
#' @return data.frame with summary by parameter
summarize_validation <- function(validation_results, pop_params) {

  # Summarize by parameter
  summary_df <- do.call(rbind, lapply(c("a", "b", "indirect"), function(param) {
    param_data <- validation_results[validation_results$label == param, ]

    if (nrow(param_data) == 0) return(NULL)

    data.frame(
      parameter = param,
      true_value = param_data$true_value[1],
      mean_estimate = mean(param_data$est, na.rm = TRUE),
      sd_estimate = sd(param_data$est, na.rm = TRUE),
      mean_se = mean(param_data$se, na.rm = TRUE),
      bias = mean(param_data$bias, na.rm = TRUE),
      rel_bias_pct = 100 * mean(param_data$bias, na.rm = TRUE) / param_data$true_value[1],
      coverage = mean(param_data$covered, na.rm = TRUE),
      n_reps = sum(!is.na(param_data$est))
    )
  }))

  return(summary_df)
}


#' Print Validation Summary
#'
#' @param summary_df data.frame from summarize_validation()
print_validation_summary <- function(summary_df) {

  cat("=== Complete-Data Validation Summary ===\n\n")
  cat(sprintf("Replications: %d\n\n", summary_df$n_reps[1]))

  cat("Parameter Recovery:\n")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    cat(sprintf("  %s: true=%.3f, mean=%.3f, bias=%.4f (%.1f%%), coverage=%.1f%%\n",
                row$parameter, row$true_value, row$mean_estimate,
                row$bias, row$rel_bias_pct, row$coverage * 100))
  }

  # Check if coverage is reasonable (should be ~95% for 95% CIs)
  coverage_ok <- all(summary_df$coverage >= 0.90 & summary_df$coverage <= 0.99)
  bias_ok <- all(abs(summary_df$rel_bias_pct) < 5)

  cat(sprintf("\nValidation: %s\n",
              ifelse(coverage_ok && bias_ok, "PASSED", "CHECK RESULTS")))
}


# ----------------------------------------------------------------------------
# MI Analysis using lavaan.mi
# ----------------------------------------------------------------------------

#' Fit Mediation Model to Multiple Imputations
#'
#' Uses lavaan.mi::sem.mi() to fit model to all imputed datasets
#'
#' @param imputed_datasets List of imputed data.frames
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @return lavaan.mi fit object
fit_mediation_mi <- function(imputed_datasets, model_type = "full_mediation") {

  if (!requireNamespace("lavaan.mi", quietly = TRUE)) {
    stop("Package 'lavaan.mi' required. Install with: install.packages('lavaan.mi')")
  }

  syntax <- get_mediation_syntax(model_type)

  # Fit to all imputations
  fit_mi <- lavaan.mi::sem.mi(
    model = syntax,
    data = imputed_datasets,
    meanstructure = TRUE
  )

  return(fit_mi)
}


#' Compute RIV from lavaan.mi Fit
#'
#' Extracts W, B matrices and computes tr(RIV)
#'
#' @param fit_mi lavaan.mi fit object
#' @return List with W, B, RIV matrix, tr_RIV, and component details
compute_RIV_from_mi <- function(fit_mi) {

  # Extract variance components
  W <- vcov(fit_mi, type = "within")
  B <- vcov(fit_mi, type = "between")

  # Get M (number of imputations)
  M <- length(fit_mi@DataList)

  # Compute RIV = (1 + 1/M) * W^{-1} * B
  W_inv <- tryCatch(
    solve(W),
    error = function(e) {
      warning("W matrix singular; using regularized inverse")
      ridge <- 1e-6 * mean(diag(W))
      solve(W + diag(ridge, nrow(W)))
    }
  )

  RIV <- (1 + 1/M) * W_inv %*% B

  # Compute trace
  tr_RIV <- sum(diag(RIV))

  # Per-parameter RIV (diagonal)
  param_RIV <- diag(RIV)
  names(param_RIV) <- rownames(W)

  return(list(
    W = W,
    B = B,
    RIV = RIV,
    tr_RIV = tr_RIV,
    param_RIV = param_RIV,
    M = M
  ))
}


#' Get Pooled Parameter Estimates from lavaan.mi
#'
#' @param fit_mi lavaan.mi fit object
#' @return data.frame with pooled estimates, SEs, and FMI
get_pooled_estimates <- function(fit_mi) {

  params <- lavaan.mi::parameterEstimates.mi(fit_mi, fmi = TRUE)

  # Filter to parameters of interest
  params_of_interest <- params[params$label %in% c("a", "b", "c", "indirect", "direct", "total") |
                               (params$op == "~" & params$lhs %in% c("M", "Y")), ]

  return(params_of_interest)
}


#' Get Log-Likelihood from Imputed Datasets
#'
#' Computes Q̄ = (1/M)Σℓ(θ̄_pooled | Y^(m)) by evaluating log-likelihood
#' at the pooled Rubin's rules estimate for each imputed dataset
#'
#' @param imputed_datasets List of imputed data.frames
#' @param fit_mi lavaan.mi fit object (if NULL, uses each imputation's MLE instead)
#' @param model_type Character. "full_mediation" or "partial_mediation"
#' @return Numeric vector of log-likelihoods
#' @details
#'   The derivation assumes Q̄ is evaluated at the pooled estimate θ̄_pooled,
#'   not at each imputation's separate MLE θ̂_m.
#'
#'   Strategy: Fit the model with all parameters FIXED to pooled values,
#'   then extract log-likelihood. This is done by setting lower=upper bounds.
get_loglik_per_imputation <- function(imputed_datasets, fit_mi = NULL,
                                       model_type = "full_mediation") {

  syntax <- get_mediation_syntax(model_type)
  M <- length(imputed_datasets)

  if (!is.null(fit_mi)) {
    # CORRECT implementation: Evaluate Q̄ at pooled Rubin's rules estimate
    # Q̄ = (1/M)Σℓ(θ̄_pooled | Y^(m))

    # Get pooled parameter estimates
    pooled_coef <- coef(fit_mi)

    # Get parameter table to construct fixed-parameter syntax
    fit_template <- lavaan::sem(syntax, data = imputed_datasets[[1]],
                               do.fit = FALSE, meanstructure = TRUE)
    partable <- lavaan::parameterTable(fit_template)

    # Construct parameter names matching coef() convention
    free_idx <- partable$free > 0
    param_names <- ifelse(
      partable$label[free_idx] != "",
      partable$label[free_idx],
      paste0(partable$lhs[free_idx], partable$op[free_idx], partable$rhs[free_idx])
    )

    # Match pooled values to parameter table rows
    pooled_values <- pooled_coef[param_names]

    if (any(is.na(pooled_values))) {
      stop("Could not match all parameters to pooled estimates")
    }

    # Create modified parameter table with fixed values
    # Set ustart, lower, upper to pooled values (this fixes parameters)
    partable_fixed <- partable
    partable_fixed$ustart[free_idx] <- pooled_values
    partable_fixed$lower[free_idx] <- pooled_values
    partable_fixed$upper[free_idx] <- pooled_values

    # Evaluate log-likelihood on each imputed dataset
    logliks <- sapply(seq_len(M), function(m) {
      tryCatch({
        # Fit with fixed parameters
        fit_fixed <- lavaan::sem(model = partable_fixed,
                                data = imputed_datasets[[m]],
                                do.fit = TRUE,
                                se = "none",  # Skip SE computation
                                test = "standard",  # Need this for logLik
                                verbose = FALSE)

        # Extract log-likelihood (use logLik for reliability)
        as.numeric(logLik(fit_fixed))

      }, error = function(e) {
        warning(sprintf("Error computing Q̄ for imputation %d: %s", m, e$message))
        NA
      })
    })

  } else {
    # INCORRECT (but useful for comparison): Evaluate at each imputation's own MLE
    # This was the original implementation and does NOT match the derivation

    logliks <- sapply(seq_len(M), function(m) {
      fit <- tryCatch(
        lavaan::sem(syntax, data = imputed_datasets[[m]], meanstructure = TRUE),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA)
      lavaan::fitMeasures(fit, "logl")
    })
  }

  return(logliks)
}
