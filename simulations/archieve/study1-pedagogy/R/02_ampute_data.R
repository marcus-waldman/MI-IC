# ============================================================================
# Amputation Functions for Study 1
# ============================================================================
# Creates missing data patterns using mice::ampute
# Supports MCAR, MAR, and MNAR mechanisms
# X is always fully observed; M and Y have missingness
# ============================================================================

#' Ampute Mediation Data
#'
#' Creates missing data in M and Y using mice::ampute
#'
#' @param data data.frame with X, M, Y columns (complete data)
#' @param prop_missing Numeric. Proportion of observations with any missingness (0 to 1)
#' @param mechanism Character. "MCAR", "MAR", or "MNAR"
#' @param seed Integer. Random seed for reproducibility
#' @return List with:
#'   \item{data_miss}{data.frame with NA values in M and Y}
#'   \item{data_complete}{Original complete data}
#'   \item{mechanism}{Mechanism used}
#'   \item{prop_missing}{Target proportion}
#'   \item{actual_prop_missing}{Actual proportion with any missing}
#'   \item{n_miss_M}{Count of missing M values}
#'   \item{n_miss_Y}{Count of missing Y values}
#' @details
#'   Missingness patterns:
#'   - Pattern 1: M missing, Y observed (1, 0, 1)
#'   - Pattern 2: M observed, Y missing (1, 1, 0)
#'   - Pattern 3: Both M and Y missing (1, 0, 0)
#'
#'   X is always fully observed (column 1 always = 1 in patterns).
#'
#'   Mechanisms:
#'   - MCAR: Missingness independent of all data
#'   - MAR: Missingness depends on X (fully observed)
#'   - MNAR: Missingness depends on M/Y themselves
#'
#' @examples
#'   data <- generate_mediation_data(100, seed = 123)
#'   result <- ampute_mediation(data, prop_missing = 0.5, mechanism = "MAR", seed = 456)
ampute_mediation <- function(data,
                              prop_missing = 0.5,
                              mechanism = c("MCAR", "MAR", "MNAR"),
                              seed = NULL) {

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' required. Install with: install.packages('mice')")
  }

  # Input validation
  mechanism <- match.arg(mechanism)

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (!all(c("X", "M", "Y") %in% names(data))) {
    stop("data must contain columns X, M, Y")
  }

  if (prop_missing < 0 || prop_missing > 1) {
    stop("prop_missing must be between 0 and 1")
  }

  n <- nrow(data)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Handle edge case: no missingness
  if (prop_missing == 0) {
    return(list(
      data_miss = data,
      data_complete = data,
      mechanism = mechanism,
      prop_missing = 0,
      actual_prop_missing = 0,
      n_miss_M = 0,
      n_miss_Y = 0
    ))
  }

  # Define missingness patterns
  # Columns: X, M, Y (1 = observed, 0 = missing)
  # X is always observed
  patterns <- matrix(c(
    1, 0, 1,  # Pattern 1: Only M missing
    1, 1, 0,  # Pattern 2: Only Y missing
    1, 0, 0   # Pattern 3: Both M and Y missing
  ), nrow = 3, ncol = 3, byrow = TRUE)

  # Frequency of each pattern (equal by default)
  freq <- c(1/3, 1/3, 1/3)

  # Define weights based on mechanism
  # weights matrix: rows = patterns, cols = variables (X, M, Y)
  # Higher weight = more influence on missingness probability
  if (mechanism == "MCAR") {
    # MCAR: No weights (missingness independent of data)
    weights <- matrix(0, nrow = 3, ncol = 3)
    mech_type <- "MCAR"

  } else if (mechanism == "MAR") {
    # MAR: Missingness depends on X (fully observed)
    # All patterns weighted by X
    weights <- matrix(c(
      1, 0, 0,  # Pattern 1: depends on X
      1, 0, 0,  # Pattern 2: depends on X
      1, 0, 0   # Pattern 3: depends on X
    ), nrow = 3, ncol = 3, byrow = TRUE)
    mech_type <- "MAR"

  } else if (mechanism == "MNAR") {
    # MNAR: Missingness depends on the missing variable itself
    # Pattern 1 (M missing): depends on M
    # Pattern 2 (Y missing): depends on Y
    # Pattern 3 (both missing): depends on both M and Y
    weights <- matrix(c(
      0, 1, 0,  # Pattern 1: depends on M
      0, 0, 1,  # Pattern 2: depends on Y
      0, 1, 1   # Pattern 3: depends on M and Y
    ), nrow = 3, ncol = 3, byrow = TRUE)
    mech_type <- "MNAR"
  }

  # Run ampute
  ampute_result <- mice::ampute(
    data = data,
    prop = prop_missing,
    patterns = patterns,
    freq = freq,
    mech = mech_type,
    weights = weights
  )

  data_miss <- ampute_result$amp

  # Compute missingness summary
  n_miss_M <- sum(is.na(data_miss$M))
  n_miss_Y <- sum(is.na(data_miss$Y))
  n_any_missing <- sum(!complete.cases(data_miss))
  actual_prop_missing <- n_any_missing / n

  return(list(
    data_miss = data_miss,
    data_complete = data,
    mechanism = mechanism,
    prop_missing = prop_missing,
    actual_prop_missing = actual_prop_missing,
    n_miss_M = n_miss_M,
    n_miss_Y = n_miss_Y,
    prop_miss_M = n_miss_M / n,
    prop_miss_Y = n_miss_Y / n
  ))
}


#' Get Missingness Summary
#'
#' Summarizes the missingness pattern in amputed data
#'
#' @param ampute_result List returned by ampute_mediation()
#' @return data.frame with missingness statistics
get_missingness_summary <- function(ampute_result) {

  data_miss <- ampute_result$data_miss
  n <- nrow(data_miss)

  # Pattern counts
  both_obs <- sum(!is.na(data_miss$M) & !is.na(data_miss$Y))
  m_miss_only <- sum(is.na(data_miss$M) & !is.na(data_miss$Y))
  y_miss_only <- sum(!is.na(data_miss$M) & is.na(data_miss$Y))
  both_miss <- sum(is.na(data_miss$M) & is.na(data_miss$Y))

  summary_df <- data.frame(
    pattern = c("Complete (M,Y obs)", "M missing only", "Y missing only", "Both missing"),
    count = c(both_obs, m_miss_only, y_miss_only, both_miss),
    proportion = c(both_obs, m_miss_only, y_miss_only, both_miss) / n
  )

  # Add variable-level summaries
  var_summary <- data.frame(
    variable = c("X", "M", "Y"),
    n_missing = c(
      sum(is.na(data_miss$X)),
      sum(is.na(data_miss$M)),
      sum(is.na(data_miss$Y))
    ),
    prop_missing = c(
      mean(is.na(data_miss$X)),
      mean(is.na(data_miss$M)),
      mean(is.na(data_miss$Y))
    )
  )

  return(list(
    pattern_summary = summary_df,
    variable_summary = var_summary,
    mechanism = ampute_result$mechanism,
    target_prop = ampute_result$prop_missing,
    actual_prop = ampute_result$actual_prop_missing
  ))
}


#' Print Missingness Summary
#'
#' Pretty-prints the missingness summary
#'
#' @param ampute_result List returned by ampute_mediation()
print_missingness_summary <- function(ampute_result) {

  summary <- get_missingness_summary(ampute_result)

  cat("=== Missingness Summary ===\n")
  cat(sprintf("Mechanism: %s\n", summary$mechanism))
  cat(sprintf("Target proportion: %.1f%%\n", summary$target_prop * 100))
  cat(sprintf("Actual proportion with any missing: %.1f%%\n\n", summary$actual_prop * 100))

  cat("Pattern distribution:\n")
  print(summary$pattern_summary, row.names = FALSE)

  cat("\nVariable-level missingness:\n")
  print(summary$variable_summary, row.names = FALSE)
}
