# ============================================================================
# Extract Imputations from brms Fit
# ============================================================================
# Extracts posterior draws of imputed values and creates M complete datasets
# ============================================================================

#' Extract Imputed Datasets from brms Fit
#'
#' Creates M complete datasets by extracting posterior draws of missing values
#'
#' @param brms_fit brmsfit object from fit_brms_imputation()
#' @param original_data data.frame with missing values (used as template)
#' @param M Integer. Number of imputed datasets to extract (default: all posterior samples)
#' @param thin Integer. Thinning interval - use every thin-th sample (default: 1)
#' @return List of M data.frames, each with complete data
#' @details
#'   brms stores imputed values as parameters named Ymi_<var>[i] where i is the
#'   index among missing values for that variable.
#'
#'   For example, if M has 30 missing values and Y has 25:
#'   - Ymi_M[1], Ymi_M[2], ..., Ymi_M[30]
#'   - Ymi_Y[1], Ymi_Y[2], ..., Ymi_Y[25]
#'
#'   Each posterior draw becomes one imputed dataset.
#'
#' @examples
#'   imputed_datasets <- extract_imputed_datasets(brms_fit, data_miss)
#'   length(imputed_datasets)  # Number of imputations (= posterior samples)
extract_imputed_datasets <- function(brms_fit, original_data, M = NULL, thin = 1) {

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' required. Install with: install.packages('posterior')")
  }

  # Extract posterior draws
  draws <- posterior::as_draws_df(brms_fit)

  # Get total number of posterior samples
  n_samples <- nrow(draws)

  # Determine M
  if (is.null(M)) {
    M <- floor(n_samples / thin)
  } else {
    if (M * thin > n_samples) {
      warning(sprintf(
        "Requested M=%d with thin=%d requires %d samples, but only %d available. Using M=%d.",
        M, thin, M * thin, n_samples, floor(n_samples / thin)
      ))
      M <- floor(n_samples / thin)
    }
  }

  # Sample indices (with thinning)
  sample_indices <- seq(1, n_samples, by = thin)[1:M]

  # Identify missing value locations in original data
  missing_M_idx <- which(is.na(original_data$M))
  missing_Y_idx <- which(is.na(original_data$Y))

  n_miss_M <- length(missing_M_idx)
  n_miss_Y <- length(missing_Y_idx)

  # Get column names for imputed values
  draw_names <- names(draws)
  ymi_M_cols <- grep("^Ymi_M\\[", draw_names, value = TRUE)
  ymi_Y_cols <- grep("^Ymi_Y\\[", draw_names, value = TRUE)


  # Validate column counts match missing counts
  if (length(ymi_M_cols) != n_miss_M) {
    warning(sprintf(
      "Expected %d Ymi_M columns but found %d. Check brms fit.",
      n_miss_M, length(ymi_M_cols)
    ))
  }
  if (length(ymi_Y_cols) != n_miss_Y) {
    warning(sprintf(
      "Expected %d Ymi_Y columns but found %d. Check brms fit.",
      n_miss_Y, length(ymi_Y_cols)
    ))
  }

  # Create M imputed datasets
  imputed_datasets <- vector("list", M)

  for (m in seq_len(M)) {
    # Start with copy of original data
    imp_data <- original_data

    # Get the m-th posterior sample
    draw_idx <- sample_indices[m]

    # Fill in imputed M values
    # Extract row index from column name to ensure correct matching
    for (col_name in ymi_M_cols) {
      # Extract the row index from "Ymi_M[idx]"
      row_idx <- as.integer(sub("Ymi_M\\[(\\d+)\\]", "\\1", col_name))
      imp_data$M[row_idx] <- draws[[col_name]][draw_idx]
    }

    # Fill in imputed Y values
    for (col_name in ymi_Y_cols) {
      # Extract the row index from "Ymi_Y[idx]"
      row_idx <- as.integer(sub("Ymi_Y\\[(\\d+)\\]", "\\1", col_name))
      imp_data$Y[row_idx] <- draws[[col_name]][draw_idx]
    }

    imputed_datasets[[m]] <- imp_data
  }

  # Add metadata
  attr(imputed_datasets, "M") <- M
  attr(imputed_datasets, "n_miss_M") <- n_miss_M
  attr(imputed_datasets, "n_miss_Y") <- n_miss_Y
  attr(imputed_datasets, "thin") <- thin

  return(imputed_datasets)
}


#' Validate Imputed Datasets
#'
#' Checks that imputed datasets have no missing values and reasonable distributions
#'
#' @param imputed_datasets List of imputed data.frames
#' @param original_data Original data.frame with missing values
#' @return List with validation results
validate_imputed_datasets <- function(imputed_datasets, original_data) {

  M <- length(imputed_datasets)

  # Check for remaining missing values
  n_remaining_missing <- sapply(imputed_datasets, function(d) sum(is.na(d)))

  # Check distributions of imputed values
  # Compare to observed values where available

  # Get indices of originally missing values
  missing_M_idx <- which(is.na(original_data$M))
  missing_Y_idx <- which(is.na(original_data$Y))
  observed_M_idx <- which(!is.na(original_data$M))
  observed_Y_idx <- which(!is.na(original_data$Y))

  # Compute summary stats for imputed values across datasets
  if (length(missing_M_idx) > 0) {
    imputed_M_means <- sapply(imputed_datasets, function(d) mean(d$M[missing_M_idx]))
    imputed_M_mean <- mean(imputed_M_means)
    imputed_M_sd <- sd(imputed_M_means)
  } else {
    imputed_M_mean <- NA
    imputed_M_sd <- NA
  }

  if (length(missing_Y_idx) > 0) {
    imputed_Y_means <- sapply(imputed_datasets, function(d) mean(d$Y[missing_Y_idx]))
    imputed_Y_mean <- mean(imputed_Y_means)
    imputed_Y_sd <- sd(imputed_Y_means)
  } else {
    imputed_Y_mean <- NA
    imputed_Y_sd <- NA
  }

  # Observed means for comparison
  observed_M_mean <- if (length(observed_M_idx) > 0) mean(original_data$M[observed_M_idx]) else NA
  observed_Y_mean <- if (length(observed_Y_idx) > 0) mean(original_data$Y[observed_Y_idx]) else NA

  return(list(
    M = M,
    all_complete = all(n_remaining_missing == 0),
    n_remaining_missing = n_remaining_missing,
    n_imputed_M = length(missing_M_idx),
    n_imputed_Y = length(missing_Y_idx),
    imputed_M_mean = imputed_M_mean,
    imputed_M_between_sd = imputed_M_sd,
    observed_M_mean = observed_M_mean,
    imputed_Y_mean = imputed_Y_mean,
    imputed_Y_between_sd = imputed_Y_sd,
    observed_Y_mean = observed_Y_mean
  ))
}


#' Print Imputation Extraction Summary
#'
#' @param imputed_datasets List of imputed data.frames
#' @param original_data Original data with missing values
print_imputation_summary <- function(imputed_datasets, original_data) {

  validation <- validate_imputed_datasets(imputed_datasets, original_data)

  cat("=== Imputation Extraction Summary ===\n")
  cat(sprintf("Number of imputed datasets: %d\n", validation$M))
  cat(sprintf("All datasets complete: %s\n",
              ifelse(validation$all_complete, "Yes", "NO - CHECK FOR ERRORS")))
  cat(sprintf("\nImputed values:\n"))
  cat(sprintf("  M: %d values imputed\n", validation$n_imputed_M))
  cat(sprintf("     Mean across imputations: %.3f (observed mean: %.3f)\n",
              validation$imputed_M_mean, validation$observed_M_mean))
  cat(sprintf("     Between-imputation SD: %.3f\n", validation$imputed_M_between_sd))
  cat(sprintf("  Y: %d values imputed\n", validation$n_imputed_Y))
  cat(sprintf("     Mean across imputations: %.3f (observed mean: %.3f)\n",
              validation$imputed_Y_mean, validation$observed_Y_mean))
  cat(sprintf("     Between-imputation SD: %.3f\n", validation$imputed_Y_between_sd))
}
