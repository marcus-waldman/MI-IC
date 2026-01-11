# ============================================================================
# Data Amputation for Study 1 (Monotone Missingness Only)
# ============================================================================
# Monotone pattern: Y missing alone, or M+Y missing together (never M alone)
# This simplifies the Q-function to only 3 patterns with quadratic θ̃ terms

#' Ampute Data with Monotone MCAR
#'
#' Creates monotone missingness with pure random selection.
#' z ~ N(0,1) (no relationship to data)
#'
#' @param data data.frame with X, M, Y
#' @param cut1 Quantile for Y-only missing (default: 0.7)
#' @param cut2 Quantile for M+Y missing (default: 0.85)
#' @param seed Random seed (optional)
#' @return data.frame with missing values
ampute_mcar_monotone <- function(data, cut1 = 0.7, cut2 = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  data_miss <- data

  # z is pure noise (MCAR - no relationship to data)
  z <- rnorm(n, mean = 0, sd = 1)

  # Cut scores from standard normal quantiles
  c1 <- qnorm(cut1)
  c2 <- qnorm(cut2)

  # Apply monotone pattern
  # z > c2: both M and Y missing
  # c1 < z <= c2: Y missing only
  data_miss$Y[z > c1] <- NA
  data_miss$M[z > c2] <- NA

  data_miss
}


#' Ampute Data with Monotone MAR
#'
#' Creates monotone missingness based on X (always observed).
#' z = γ * standardize(X) + ε, where z ~ N(0,1)
#'
#' @param data data.frame with X, M, Y
#' @param gamma Strength of X's influence on missingness (default: 0.5)
#' @param cut1 Quantile for Y-only missing (default: 0.7)
#' @param cut2 Quantile for M+Y missing (default: 0.85)
#' @param seed Random seed (optional)
#' @return data.frame with missing values
ampute_mar_monotone <- function(data, gamma = 0.5, cut1 = 0.7, cut2 = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  data_miss <- data

  # Standardize X
  X_std <- (data$X - mean(data$X)) / sd(data$X)

  # z = γ * X_std + ε, where Var(ε) = 1 - γ² so z ~ N(0,1)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(1 - gamma^2))
  z <- gamma * X_std + epsilon

  # Cut scores from standard normal quantiles
  c1 <- qnorm(cut1)
  c2 <- qnorm(cut2)

  # Apply monotone pattern
  data_miss$Y[z > c1] <- NA
  data_miss$M[z > c2] <- NA

  data_miss
}


#' Ampute Data with Monotone MNAR
#'
#' Creates monotone missingness based on Y (becomes missing).
#' z = γ * standardize(Y) + ε, where z ~ N(0,1)
#'
#' @param data data.frame with X, M, Y
#' @param gamma Strength of Y's influence on missingness (default: 0.5)
#' @param cut1 Quantile for Y-only missing (default: 0.7)
#' @param cut2 Quantile for M+Y missing (default: 0.85)
#' @param seed Random seed (optional)
#' @return data.frame with missing values
ampute_mnar_monotone <- function(data, gamma = 0.5, cut1 = 0.7, cut2 = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  data_miss <- data

  # Standardize Y (the variable that will be missing)
  Y_std <- (data$Y - mean(data$Y)) / sd(data$Y)

  # z = γ * Y_std + ε, where Var(ε) = 1 - γ² so z ~ N(0,1)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(1 - gamma^2))
  z <- gamma * Y_std + epsilon

  # Cut scores from standard normal quantiles
  c1 <- qnorm(cut1)
  c2 <- qnorm(cut2)

  # Apply monotone pattern
  data_miss$Y[z > c1] <- NA
  data_miss$M[z > c2] <- NA

  data_miss
}


#' Ampute Data (Wrapper for Monotone Methods)
#'
#' @param data data.frame with X, M, Y
#' @param mechanism "MCAR", "MAR", or "MNAR" (default: "MCAR")
#' @param gamma Strength parameter for MAR/MNAR (default: 0.5)
#' @param cut1 Quantile for Y-only missing (default: 0.7)
#' @param cut2 Quantile for M+Y missing (default: 0.85)
#' @param seed Random seed (optional)
#' @return data.frame with missing values
ampute_data <- function(data, mechanism = "MCAR", gamma = 0.5,
                        cut1 = 0.7, cut2 = 0.85, seed = NULL) {
  mechanism <- toupper(mechanism)

  if (mechanism == "MCAR") {
    ampute_mcar_monotone(data, cut1 = cut1, cut2 = cut2, seed = seed)
  } else if (mechanism == "MAR") {
    ampute_mar_monotone(data, gamma = gamma, cut1 = cut1, cut2 = cut2, seed = seed)
  } else if (mechanism == "MNAR") {
    ampute_mnar_monotone(data, gamma = gamma, cut1 = cut1, cut2 = cut2, seed = seed)
  } else {
    stop("mechanism must be 'MCAR', 'MAR', or 'MNAR'")
  }
}
